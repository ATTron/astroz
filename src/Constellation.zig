//! Constellation propagation with SIMD batching and threading.
//!
//! APIs:
//! - `Constellation` struct — init from mixed TLEs, auto-classifies SGP4/SDP4, single `propagate()` call
//! - `propagateConstellation` — stateless SGP4-only batch propagation from pre-built batch elements
//! - `screenConstellation` — fused propagate+screen for conjunction assessment

const std = @import("std");
const Allocator = std.mem.Allocator;
const Sgp4 = @import("Sgp4.zig");
const Sdp4 = @import("Sdp4.zig");
pub const Sgp4Batch = @import("Sgp4Batch.zig");
pub const Sdp4Batch = @import("Sdp4Batch.zig");
const Tle = @import("Tle.zig");
const Datetime = @import("Datetime.zig");
const constants = @import("constants.zig");
pub const simdMath = @import("simdMath.zig");
const dispatch = @import("dispatch.zig");
const WCS = @import("WorldCoordinateSystem.zig");

const Error = Sgp4.Error;

pub const BatchSize: usize = Sgp4Batch.BatchSize;
const BatchResults = Sgp4.PosVelArray(BatchSize);
const Sgp4Bels = Sgp4Batch.BatchElements(BatchSize);
const Sdp4Bels = Sdp4Batch.Sdp4BatchElements(BatchSize);
const Sdp4Carry = Sdp4Batch.ResonanceCarryBatch(BatchSize);

/// Output coordinate modes for propagation
pub const OutputMode = enum(u8) {
    teme = 0,
    ecef = 1,
    geodetic = 2,
};

/// Memory layout for output arrays
pub const Layout = enum(u1) {
    /// (num_satellites, num_times, 3): threads over batch ranges
    satelliteMajor = 0,
    /// (num_times, num_satellites, 3): threads over time ranges (typically faster)
    timeMajor = 1,
};

pub const MaxThreads: usize = 128;

pub inline fn outBase(comptime layout: Layout, satIdx: usize, timeIdx: usize, numTimes: usize, numSatellites: usize) usize {
    return if (layout == .satelliteMajor)
        satIdx * numTimes * 3 + timeIdx * 3
    else
        timeIdx * numSatellites * 3 + satIdx * 3;
}

/// Scalar ECI to ECEF rotation for a single position vector
inline fn eciToEcef(pos: [3]f64, sinG: f64, cosG: f64) [3]f64 {
    return .{ pos[0] * cosG + pos[1] * sinG, pos[1] * cosG - pos[0] * sinG, pos[2] };
}

/// Cached thread count to avoid repeated syscalls
var cachedMaxThreads: usize = 0;

pub fn getMaxThreads() usize {
    if (cachedMaxThreads != 0) return cachedMaxThreads;

    var result: usize = std.Thread.getCpuCount() catch 1;
    if (@hasDecl(std, "c")) {
        const envVal = std.c.getenv("ASTROZ_THREADS");
        if (envVal) |val| {
            result = std.fmt.parseInt(usize, std.mem.sliceTo(val, 0), 10) catch result;
        }
    }
    result = @min(result, MaxThreads);
    cachedMaxThreads = result;
    return result;
}

const Constellation = @This();

// SGP4 (near earth) batches
sgp4Batches: []Sgp4Bels,
sgp4EpochOffsets: []f64,
sgp4OrigIndices: []u32,
numSgp4: usize,

// SDP4 (deep space) batches
sdp4Batches: []Sdp4Bels,
sdp4BatchEpochs: [][BatchSize]f64,
sdp4OrigIndices: []u32,
sdp4Carries: []Sdp4Carry,
numSdp4: usize,

// Reference epoch for SGP4 tsince computation
referenceEpochJd: f64,

// Totals
numSatellites: usize,
allocator: Allocator,

/// Initialize a constellation from a mixed set of TLEs.
/// Auto-classifies each TLE as SGP4 (near earth) or SDP4 (deep space),
/// groups them into SIMD batches, and prepares carry state for SDP4 resonance.
pub fn init(allocator: Allocator, tles: []const Tle, grav: constants.Sgp4GravityModel) Error!Constellation {
    const n = tles.len;

    // Temp arrays for classification
    const sgp4Idx = allocator.alloc(u32, n) catch return Error.OutOfMemory;
    defer allocator.free(sgp4Idx);
    const sdp4Idx = allocator.alloc(u32, n) catch return Error.OutOfMemory;
    defer allocator.free(sdp4Idx);
    const sdp4Els = allocator.alloc(Sdp4.Elements, n) catch return Error.OutOfMemory;
    defer allocator.free(sdp4Els);

    var nSgp4: usize = 0;
    var nSdp4: usize = 0;

    for (tles, 0..) |tle, i| {
        if (Sgp4.initElements(tle, grav)) |_| {
            sgp4Idx[nSgp4] = @intCast(i);
            nSgp4 += 1;
        } else |err| {
            if (err == Error.DeepSpaceNotSupported) {
                sdp4Els[nSdp4] = (try Sdp4.init(tle, grav)).elements;
                sdp4Idx[nSdp4] = @intCast(i);
                nSdp4 += 1;
            } else return err;
        }
    }

    // Batch SGP4
    const nSgp4B = if (nSgp4 == 0) 0 else (nSgp4 + BatchSize - 1) / BatchSize;
    const sgp4Pad = nSgp4B * BatchSize;

    const sgp4Batches = allocator.alloc(Sgp4Bels, nSgp4B) catch return Error.OutOfMemory;
    errdefer allocator.free(sgp4Batches);
    const sgp4EpochOff = allocator.alloc(f64, sgp4Pad) catch return Error.OutOfMemory;
    errdefer allocator.free(sgp4EpochOff);
    const sgp4Orig = allocator.alloc(u32, sgp4Pad) catch return Error.OutOfMemory;
    errdefer allocator.free(sgp4Orig);

    var refEpoch: f64 = 0;
    if (nSgp4 > 0) refEpoch = tleEpochJd(tles[sgp4Idx[0]]);

    for (0..nSgp4B) |bi| {
        const base = bi * BatchSize;
        var bt: [BatchSize]Tle = undefined;
        for (0..BatchSize) |lane| {
            const src = if (base + lane < nSgp4) sgp4Idx[base + lane] else sgp4Idx[nSgp4 - 1];
            bt[lane] = tles[src];
        }
        sgp4Batches[bi] = try Sgp4Batch.initBatchElements(BatchSize, bt, grav);

        const epochs: [BatchSize]f64 = simdMath.readArray(BatchSize, &sgp4Batches[bi].epochJd);
        for (0..BatchSize) |lane| {
            sgp4EpochOff[base + lane] = (refEpoch - epochs[lane]) * 1440.0;
            sgp4Orig[base + lane] = if (base + lane < nSgp4) sgp4Idx[base + lane] else sgp4Idx[nSgp4 - 1];
        }
    }

    // Batch SDP4
    const nSdp4B = if (nSdp4 == 0) 0 else (nSdp4 + BatchSize - 1) / BatchSize;
    const sdp4Pad = nSdp4B * BatchSize;

    const sdp4Batches = allocator.alloc(Sdp4Bels, nSdp4B) catch return Error.OutOfMemory;
    errdefer allocator.free(sdp4Batches);
    const sdp4Epochs = allocator.alloc([BatchSize]f64, nSdp4B) catch return Error.OutOfMemory;
    errdefer allocator.free(sdp4Epochs);
    const sdp4Orig = allocator.alloc(u32, sdp4Pad) catch return Error.OutOfMemory;
    errdefer allocator.free(sdp4Orig);
    const sdp4Carries = allocator.alloc(Sdp4Carry, nSdp4B) catch return Error.OutOfMemory;
    errdefer allocator.free(sdp4Carries);

    for (0..nSdp4B) |bi| {
        const base = bi * BatchSize;
        var be: [BatchSize]Sdp4.Elements = undefined;
        var ep: [BatchSize]f64 = undefined;
        for (0..BatchSize) |lane| {
            const src = if (base + lane < nSdp4) base + lane else nSdp4 - 1;
            be[lane] = sdp4Els[src];
            ep[lane] = sdp4Els[src].sgp4.epochJd;
            sdp4Orig[base + lane] = if (base + lane < nSdp4) sdp4Idx[src] else sdp4Idx[nSdp4 - 1];
        }
        sdp4Batches[bi] = Sdp4Batch.initFromElements(BatchSize, be, grav);
        sdp4Epochs[bi] = ep;
        sdp4Carries[bi] = Sdp4Batch.initCarry(BatchSize, &sdp4Batches[bi]);
    }

    return .{
        .sgp4Batches = sgp4Batches,
        .sgp4EpochOffsets = sgp4EpochOff,
        .sgp4OrigIndices = sgp4Orig,
        .numSgp4 = nSgp4,
        .sdp4Batches = sdp4Batches,
        .sdp4BatchEpochs = sdp4Epochs,
        .sdp4OrigIndices = sdp4Orig,
        .sdp4Carries = sdp4Carries,
        .numSdp4 = nSdp4,
        .referenceEpochJd = refEpoch,
        .numSatellites = n,
        .allocator = allocator,
    };
}

pub fn deinit(self: *Constellation) void {
    self.allocator.free(self.sgp4Batches);
    self.allocator.free(self.sgp4EpochOffsets);
    self.allocator.free(self.sgp4OrigIndices);
    self.allocator.free(self.sdp4Batches);
    self.allocator.free(self.sdp4BatchEpochs);
    self.allocator.free(self.sdp4OrigIndices);
    self.allocator.free(self.sdp4Carries);
}

/// Reset all SDP4 resonance carry states.
/// Call when switching to a non-monotonic time sequence.
pub fn resetCarry(self: *Constellation) void {
    for (self.sdp4Batches, self.sdp4Carries) |*batch, *carry| {
        carry.* = Sdp4Batch.initCarry(BatchSize, batch);
    }
}

/// Thread context for unified propagation
const PropCtx = struct {
    sgp4Batches: []const Sgp4Bels,
    sgp4EpochOffsets: []const f64,
    sgp4OrigIndices: []const u32,
    numSgp4: usize,

    sdp4Batches: []const Sdp4Bels,
    sdp4BatchEpochs: [][BatchSize]f64,
    sdp4OrigIndices: []const u32,
    sdp4Carries: []Sdp4Carry,
    numSdp4: usize,

    numSatellites: usize,
    numTimes: usize,
    tsinceBase: []const f64,
    jdFull: []const f64,
    resultsPos: []f64,
    resultsVel: ?[]f64,
    gmstSin: []const f64,
    gmstCos: []const f64,
    satelliteMask: ?[]const u8,
};

/// Propagate all satellites at the given absolute Julian Date times
pub fn propagate(
    self: *Constellation,
    jd: []const f64,
    fr: []const f64,
    resultsPos: []f64,
    resultsVel: ?[]f64,
    outputMode: OutputMode,
    layout: Layout,
) Error!void {
    const numTimes = jd.len;
    const required = self.numSatellites * numTimes * 3;
    if (resultsPos.len < required) return Error.SatelliteDecayed;
    if (resultsVel) |rv| if (rv.len < required) return Error.SatelliteDecayed;

    const pa = std.heap.page_allocator;

    const tsinceBase = pa.alloc(f64, numTimes) catch return Error.OutOfMemory;
    defer pa.free(tsinceBase);
    const jdFull = pa.alloc(f64, numTimes) catch return Error.OutOfMemory;
    defer pa.free(jdFull);

    for (0..numTimes) |t| {
        jdFull[t] = jd[t] + fr[t];
        tsinceBase[t] = (jdFull[t] - self.referenceEpochJd) * 1440.0;
    }

    var sinBuf: []f64 = &.{};
    var cosBuf: []f64 = &.{};
    defer if (sinBuf.len > 0) pa.free(sinBuf);
    defer if (cosBuf.len > 0) pa.free(cosBuf);

    if (outputMode != .teme) {
        sinBuf = pa.alloc(f64, numTimes) catch return Error.OutOfMemory;
        cosBuf = pa.alloc(f64, numTimes) catch return Error.OutOfMemory;
        for (0..numTimes) |t| {
            const gmst = WCS.julianToGmst(jdFull[t]);
            sinBuf[t] = @sin(gmst);
            cosBuf[t] = @cos(gmst);
        }
    }

    const ctx = PropCtx{
        .sgp4Batches = self.sgp4Batches,
        .sgp4EpochOffsets = self.sgp4EpochOffsets,
        .sgp4OrigIndices = self.sgp4OrigIndices,
        .numSgp4 = self.numSgp4,
        .sdp4Batches = self.sdp4Batches,
        .sdp4BatchEpochs = self.sdp4BatchEpochs,
        .sdp4OrigIndices = self.sdp4OrigIndices,
        .sdp4Carries = self.sdp4Carries,
        .numSdp4 = self.numSdp4,
        .numSatellites = self.numSatellites,
        .numTimes = numTimes,
        .tsinceBase = tsinceBase,
        .jdFull = jdFull,
        .resultsPos = resultsPos,
        .resultsVel = resultsVel,
        .gmstSin = sinBuf,
        .gmstCos = cosBuf,
        .satelliteMask = null,
    };

    dispatchImpl(layout, outputMode, resultsVel != null, ctx);
}

/// Dispatch to the correct propagateImpl variant
inline fn dispatchImpl(layout: Layout, outputMode: OutputMode, hasVel: bool, ctx: PropCtx) void {
    const key: u4 = @as(u4, @intFromEnum(layout)) << 3 |
        @as(u4, @intCast(@intFromEnum(outputMode))) << 1 |
        @intFromBool(hasVel);

    inline for (0..16) |i| {
        const mode_bits = (i >> 1) & 0b11;
        if (mode_bits < 3 and key == i) {
            const l: Layout = @enumFromInt(i >> 3);
            const m: OutputMode = @enumFromInt(mode_bits);
            const v = (i & 1) == 1;
            return propagateImpl(l, m, v, ctx);
        }
    }
}

fn propagateImpl(
    comptime layout: Layout,
    comptime mode: OutputMode,
    comptime hasVel: bool,
    ctx: PropCtx,
) void {
    const maxThreads = getMaxThreads();
    var handles: [MaxThreads]?std.Thread = .{null} ** MaxThreads;
    var idx: usize = 0;

    // SGP4 phase: thread over time ranges (timeMajor) or batch ranges (satelliteMajor)
    if (ctx.sgp4Batches.len > 0) {
        const work = if (layout == .satelliteMajor) ctx.sgp4Batches.len else ctx.numTimes;
        const nt = @min(work, maxThreads);
        const per = (work + nt - 1) / nt;
        for (0..nt) |tid| {
            const start = tid * per;
            const end = @min(start + per, work);
            if (start >= end) break;
            if (idx < MaxThreads) {
                handles[idx] = std.Thread.spawn(.{}, unifiedSgp4Range, .{
                    layout, mode, hasVel, ctx, start, end,
                }) catch {
                    unifiedSgp4Range(layout, mode, hasVel, ctx, start, end);
                    continue;
                };
                idx += 1;
            } else {
                unifiedSgp4Range(layout, mode, hasVel, ctx, start, end);
            }
        }
    }

    // SDP4 phase: thread over batch ranges (sequential time for carry)
    if (ctx.sdp4Batches.len > 0) {
        const work = ctx.sdp4Batches.len;
        const remaining = if (maxThreads > idx) maxThreads - idx else 1;
        const nt = @min(work, remaining);
        const per = (work + nt - 1) / nt;
        for (0..nt) |tid| {
            const start = tid * per;
            const end = @min(start + per, work);
            if (start >= end) break;
            if (idx < MaxThreads) {
                handles[idx] = std.Thread.spawn(.{}, unifiedSdp4Range, .{
                    layout, mode, hasVel, ctx, start, end,
                }) catch {
                    unifiedSdp4Range(layout, mode, hasVel, ctx, start, end);
                    continue;
                };
                idx += 1;
            } else {
                unifiedSdp4Range(layout, mode, hasVel, ctx, start, end);
            }
        }
    }

    for (handles) |h| if (h) |t| t.join();
}

fn unifiedSgp4Range(
    comptime layout: Layout,
    comptime mode: OutputMode,
    comptime hasVel: bool,
    ctx: PropCtx,
    rangeStart: usize,
    rangeEnd: usize,
) void {
    if (layout == .satelliteMajor) {
        // Thread over batch ranges, iterate all times per batch
        for (rangeStart..rangeEnd) |batchIdx| {
            if (!sgp4BatchActive(ctx, batchIdx)) continue;
            for (0..ctx.numTimes) |timeIdx| {
                sgp4Core(layout, mode, hasVel, ctx, batchIdx, timeIdx);
            }
        }
    } else {
        // Thread over time ranges, iterate all batches per time
        for (rangeStart..rangeEnd) |timeIdx| {
            for (ctx.sgp4Batches, 0..) |_, batchIdx| {
                sgp4Core(layout, mode, hasVel, ctx, batchIdx, timeIdx);
            }
        }
    }
}

inline fn sgp4Core(
    comptime layout: Layout,
    comptime mode: OutputMode,
    comptime hasVel: bool,
    ctx: PropCtx,
    batchIdx: usize,
    timeIdx: usize,
) void {
    const base = batchIdx * BatchSize;

    var tsince: [BatchSize]f64 = undefined;
    inline for (0..BatchSize) |lane| {
        tsince[lane] = ctx.tsinceBase[timeIdx] + ctx.sgp4EpochOffsets[base + lane];
    }

    const results = dispatch.sgp4Batch8(&ctx.sgp4Batches[batchIdx], tsince) catch {
        writeZeros(layout, hasVel, ctx, ctx.sgp4OrigIndices, base, ctx.numSgp4, timeIdx);
        return;
    };

    writeOutput(layout, mode, hasVel, ctx, ctx.sgp4OrigIndices, base, ctx.numSgp4, timeIdx, results);
}

/// Check if any lane in a batch is active (not masked out)
inline fn sgp4BatchActive(ctx: PropCtx, batchIdx: usize) bool {
    const mask = ctx.satelliteMask orelse return true;
    const base = batchIdx * BatchSize;
    inline for (0..BatchSize) |lane| {
        if (base + lane < ctx.numSgp4) {
            if (mask[ctx.sgp4OrigIndices[base + lane]] != 0) return true;
        }
    }
    return false;
}

fn unifiedSdp4Range(
    comptime layout: Layout,
    comptime mode: OutputMode,
    comptime hasVel: bool,
    ctx: PropCtx,
    batchStart: usize,
    batchEnd: usize,
) void {
    for (0..ctx.numTimes) |timeIdx| {
        const jdT = ctx.jdFull[timeIdx];

        for (batchStart..batchEnd) |batchIdx| {
            const base = batchIdx * BatchSize;

            const batchEpochs: [BatchSize]f64 = ctx.sdp4BatchEpochs[batchIdx];
            var tArr: [BatchSize]f64 = undefined;
            inline for (0..BatchSize) |lane| {
                tArr[lane] = (jdT - batchEpochs[lane]) * 1440.0;
            }

            const results = dispatch.sdp4Batch8(&ctx.sdp4Batches[batchIdx], tArr, &ctx.sdp4Carries[batchIdx]) catch {
                writeZeros(layout, hasVel, ctx, ctx.sdp4OrigIndices, base, ctx.numSdp4, timeIdx);
                continue;
            };

            writeOutput(layout, mode, hasVel, ctx, ctx.sdp4OrigIndices, base, ctx.numSdp4, timeIdx, results);
        }
    }
}

inline fn writeOutput(
    comptime layout: Layout,
    comptime mode: OutputMode,
    comptime hasVel: bool,
    ctx: PropCtx,
    origIndices: []const u32,
    base: usize,
    numReal: usize,
    timeIdx: usize,
    results: BatchResults,
) void {
    const sinG = if (mode != .teme) ctx.gmstSin[timeIdx] else 0.0;
    const cosG = if (mode != .teme) ctx.gmstCos[timeIdx] else 0.0;

    inline for (0..BatchSize) |lane| {
        if (base + lane < numReal and laneActive(ctx.satelliteMask, origIndices[base + lane])) {
            const origIdx = origIndices[base + lane];
            const ob = outBase(layout, @as(usize, origIdx), timeIdx, ctx.numTimes, ctx.numSatellites);
            ctx.resultsPos[ob..][0..3].* = switch (mode) {
                .geodetic => WCS.ecefToGeodetic(eciToEcef(results[lane][0], sinG, cosG)),
                .ecef => eciToEcef(results[lane][0], sinG, cosG),
                .teme => results[lane][0],
            };
            if (hasVel) {
                ctx.resultsVel.?[ob..][0..3].* = if (mode != .teme)
                    eciToEcef(results[lane][1], sinG, cosG)
                else
                    results[lane][1];
            }
        }
    }
}

inline fn writeZeros(
    comptime layout: Layout,
    comptime hasVel: bool,
    ctx: PropCtx,
    origIndices: []const u32,
    base: usize,
    numReal: usize,
    timeIdx: usize,
) void {
    inline for (0..BatchSize) |lane| {
        if (base + lane < numReal and laneActive(ctx.satelliteMask, origIndices[base + lane])) {
            const origIdx = origIndices[base + lane];
            const ob = outBase(layout, @as(usize, origIdx), timeIdx, ctx.numTimes, ctx.numSatellites);
            ctx.resultsPos[ob..][0..3].* = .{ 0, 0, 0 };
            if (hasVel) ctx.resultsVel.?[ob..][0..3].* = .{ 0, 0, 0 };
        }
    }
}

inline fn laneActive(satelliteMask: ?[]const u8, origIdx: u32) bool {
    if (satelliteMask) |mask| return mask[origIdx] != 0;
    return true;
}

fn tleEpochJd(tle: Tle) f64 {
    const fullYear: u16 = if (tle.firstLine.epochYear < 57)
        2000 + tle.firstLine.epochYear
    else
        1900 + tle.firstLine.epochYear;
    return Datetime.yearDoyToJulianDate(fullYear, tle.firstLine.epochDay);
}

/// Propagate SGP4 constellation with specified memory layout and threading strategy
/// Stateless: takes pre built batch elements directly
pub fn propagateConstellation(
    batches: []const Sgp4Bels,
    numSatellites: usize,
    times: []const f64,
    epochOffsets: []const f64,
    resultsPos: []f64,
    resultsVel: ?[]f64,
    outputMode: OutputMode,
    referenceJd: f64,
    satelliteMask: ?[]const u8,
    layout: Layout,
) Error!void {
    const requiredSize = numSatellites * times.len * 3;
    if (resultsPos.len < requiredSize) return Error.SatelliteDecayed;
    if (resultsVel) |rv| {
        if (rv.len < requiredSize) return Error.SatelliteDecayed;
    }

    const pa = std.heap.page_allocator;

    // Build identity origIndices: satellite i maps to output position i
    const padded = batches.len * BatchSize;
    const origIndices = pa.alloc(u32, padded) catch return Error.OutOfMemory;
    defer pa.free(origIndices);
    for (origIndices, 0..) |*v, i| v.* = @intCast(i);

    // Precompute GMST sin/cos if coordinate transform needed
    var sinBuf: []f64 = &.{};
    var cosBuf: []f64 = &.{};
    defer if (sinBuf.len > 0) pa.free(sinBuf);
    defer if (cosBuf.len > 0) pa.free(cosBuf);

    if (outputMode != .teme) {
        sinBuf = pa.alloc(f64, times.len) catch return Error.OutOfMemory;
        cosBuf = pa.alloc(f64, times.len) catch return Error.OutOfMemory;
        for (times, 0..) |time, i| {
            const gmst = WCS.julianToGmst(referenceJd + time / 1440.0);
            sinBuf[i] = @sin(gmst);
            cosBuf[i] = @cos(gmst);
        }
    }

    const ctx = PropCtx{
        .sgp4Batches = batches,
        .sgp4EpochOffsets = epochOffsets,
        .sgp4OrigIndices = origIndices,
        .numSgp4 = numSatellites,
        .sdp4Batches = &.{},
        .sdp4BatchEpochs = &.{},
        .sdp4OrigIndices = &.{},
        .sdp4Carries = &.{},
        .numSdp4 = 0,
        .numSatellites = numSatellites,
        .numTimes = times.len,
        .tsinceBase = times,
        .jdFull = &.{},
        .resultsPos = resultsPos,
        .resultsVel = resultsVel,
        .gmstSin = sinBuf,
        .gmstCos = cosBuf,
        .satelliteMask = satelliteMask,
    };

    dispatchImpl(layout, outputMode, resultsVel != null, ctx);
}

/// Stateless SDP4 only batch propagation from pre-built batch elements
/// Handles carry init, threading, and result writeback
/// `origIndices` maps batch position -> output column (length = batches.len * BatchSize)
/// `numSatellites` is the total output stride (columns in the output buffer)
pub fn propagateSdp4Constellation(
    batches: []const Sdp4Bels,
    batchEpochs: [][BatchSize]f64,
    numSdp4: usize,
    numSatellites: usize,
    origIndices: []const u32,
    jd: []const f64,
    fr: []const f64,
    resultsPos: []f64,
    resultsVel: ?[]f64,
    outputMode: OutputMode,
    layout: Layout,
) Error!void {
    const numTimes = jd.len;
    const requiredSize = numSatellites * numTimes * 3;
    if (resultsPos.len < requiredSize) return Error.SatelliteDecayed;
    if (resultsVel) |rv| if (rv.len < requiredSize) return Error.SatelliteDecayed;

    const pa = std.heap.page_allocator;

    // Compute jdFull
    const jdFull = pa.alloc(f64, numTimes) catch return Error.OutOfMemory;
    defer pa.free(jdFull);
    for (0..numTimes) |t| jdFull[t] = jd[t] + fr[t];

    // Initialize carries
    const carries = pa.alloc(Sdp4Carry, batches.len) catch return Error.OutOfMemory;
    defer pa.free(carries);
    for (batches, carries) |*b, *carry| carry.* = Sdp4Batch.initCarry(BatchSize, b);

    // Precompute GMST sin/cos if coordinate transform needed
    var sinBuf: []f64 = &.{};
    var cosBuf: []f64 = &.{};
    defer if (sinBuf.len > 0) pa.free(sinBuf);
    defer if (cosBuf.len > 0) pa.free(cosBuf);

    if (outputMode != .teme) {
        sinBuf = pa.alloc(f64, numTimes) catch return Error.OutOfMemory;
        cosBuf = pa.alloc(f64, numTimes) catch return Error.OutOfMemory;
        for (jdFull, 0..) |jdT, i| {
            const gmst = WCS.julianToGmst(jdT);
            sinBuf[i] = @sin(gmst);
            cosBuf[i] = @cos(gmst);
        }
    }

    const ctx = PropCtx{
        .sgp4Batches = &.{},
        .sgp4EpochOffsets = &.{},
        .sgp4OrigIndices = &.{},
        .numSgp4 = 0,
        .sdp4Batches = batches,
        .sdp4BatchEpochs = batchEpochs,
        .sdp4OrigIndices = origIndices,
        .sdp4Carries = carries,
        .numSdp4 = numSdp4,
        .numSatellites = numSatellites,
        .numTimes = numTimes,
        .tsinceBase = &.{},
        .jdFull = jdFull,
        .resultsPos = resultsPos,
        .resultsVel = resultsVel,
        .gmstSin = sinBuf,
        .gmstCos = cosBuf,
        .satelliteMask = null,
    };

    dispatchImpl(layout, outputMode, resultsVel != null, ctx);
}

/// Fused propagate + screen: propagate SGP4 constellation and screen against target.
/// Returns min_distances and min_t_indices per satellite (excluding target).
pub fn screenConstellation(
    batches: []const Sgp4Bels,
    numSatellites: usize,
    times: []const f64,
    epochOffsets: []const f64,
    targetIdx: usize,
    targetEpochOffset: f64,
    threshold: f64,
    referenceJd: f64,
    outMinDists: []f64,
    outMinTIndices: []u32,
) Error!void {
    if (outMinDists.len < numSatellites or outMinTIndices.len < numSatellites)
        return Error.SatelliteDecayed;

    const thresholdSq = threshold * threshold;
    const targetBatchIdx = targetIdx / BatchSize;
    const targetLane = targetIdx % BatchSize;
    const targetBatch = &batches[targetBatchIdx];

    for (0..numSatellites) |i| {
        outMinDists[i] = thresholdSq;
        outMinTIndices[i] = 0;
    }

    const page_alloc = std.heap.page_allocator;
    const sinBuf = page_alloc.alloc(f64, times.len) catch return Error.OutOfMemory;
    defer page_alloc.free(sinBuf);
    const cosBuf = page_alloc.alloc(f64, times.len) catch return Error.OutOfMemory;
    defer page_alloc.free(cosBuf);

    for (times, 0..) |time, i| {
        const gmst = WCS.julianToGmst(referenceJd + time / 1440.0);
        sinBuf[i] = @sin(gmst);
        cosBuf[i] = @cos(gmst);
    }

    for (times, 0..) |time, timeIdx| {
        var targetTimeArr: [BatchSize]f64 = undefined;
        @memset(&targetTimeArr, time + targetEpochOffset);
        const targetResults = dispatch.sgp4Batch8(targetBatch, targetTimeArr) catch continue;
        const targetPos = eciToEcef(targetResults[targetLane][0], sinBuf[timeIdx], cosBuf[timeIdx]);

        for (batches, 0..) |*batch, batchIdx| {
            const satBase = batchIdx * BatchSize;
            var batchTimes: [BatchSize]f64 = undefined;
            inline for (0..BatchSize) |i| {
                batchTimes[i] = time + epochOffsets[satBase + i];
            }

            const batchResults = dispatch.sgp4Batch8(batch, batchTimes) catch continue;

            for (0..BatchSize) |lane| {
                const satIdx = satBase + lane;
                if (satIdx >= numSatellites or satIdx == targetIdx) continue;

                const pos = eciToEcef(batchResults[lane][0], sinBuf[timeIdx], cosBuf[timeIdx]);
                const dx = targetPos[0] - pos[0];
                const dy = targetPos[1] - pos[1];
                const dz = targetPos[2] - pos[2];
                const distSq = dx * dx + dy * dy + dz * dz;

                if (distSq < outMinDists[satIdx]) {
                    outMinDists[satIdx] = distSq;
                    outMinTIndices[satIdx] = @intCast(timeIdx);
                }
            }
        }
    }

    for (0..numSatellites) |i| {
        outMinDists[i] = @sqrt(outMinDists[i]);
    }
}

const testing = std.testing;

const leoTle1 = "1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995\n2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123";
const leoTle2 = "1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998\n2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371";
const leoTle3 = "1 55910U 23035C   24187.17717543  .00022897  00000+0  15695-2 0  9993\n2 55910  43.9975 313.1680 0011485 277.1866  82.7988 15.05748091 71356";
const geoTle = "1 28626U 05004A   24186.00000000 -.00000098  00000+0  00000+0 0  9998\n2 28626   0.0163 279.8379 0003069  20.3251 343.1766  1.00270142 70992";
const gpsTle = "1 20413U 90005A   24186.00000000  .00000012  00000+0  10000-3 0  9992\n2 20413  55.4408  61.4858 0112981 129.5765 231.5553  2.00561730104446";

test "mixed init classifies correctly" {
    var tles: [5]Tle = undefined;
    defer for (&tles) |*t| t.deinit();
    tles[0] = try Tle.parse(leoTle1, testing.allocator);
    tles[1] = try Tle.parse(geoTle, testing.allocator);
    tles[2] = try Tle.parse(leoTle2, testing.allocator);
    tles[3] = try Tle.parse(gpsTle, testing.allocator);
    tles[4] = try Tle.parse(leoTle3, testing.allocator);

    var c = try Constellation.init(testing.allocator, &tles, constants.wgs72);
    defer c.deinit();

    try testing.expectEqual(@as(usize, 3), c.numSgp4);
    try testing.expectEqual(@as(usize, 2), c.numSdp4);
    try testing.expectEqual(@as(usize, 5), c.numSatellites);
}

test "propagate matches scalar - TEME timeMajor" {
    var tles: [5]Tle = undefined;
    defer for (&tles) |*t| t.deinit();
    tles[0] = try Tle.parse(leoTle1, testing.allocator);
    tles[1] = try Tle.parse(geoTle, testing.allocator);
    tles[2] = try Tle.parse(leoTle2, testing.allocator);
    tles[3] = try Tle.parse(gpsTle, testing.allocator);
    tles[4] = try Tle.parse(leoTle3, testing.allocator);

    var c = try Constellation.init(testing.allocator, &tles, constants.wgs72);
    defer c.deinit();

    const epoch = tleEpochJd(tles[0]);
    const jd1 = [_]f64{@floor(epoch)};
    const fr1 = [_]f64{epoch - @floor(epoch)};
    var pos: [5 * 3]f64 = undefined;
    var vel: [5 * 3]f64 = undefined;

    try c.propagate(&jd1, &fr1, &pos, &vel, .teme, .timeMajor);

    const leoIndices = [_]usize{ 0, 2, 4 };
    const leoTles = [_][]const u8{ leoTle1, leoTle2, leoTle3 };
    for (leoIndices, leoTles) |idx, tleStr| {
        var tle = try Tle.parse(tleStr, testing.allocator);
        defer tle.deinit();
        const sgp4 = try Sgp4.init(tle, constants.wgs72);
        const satEpoch = tleEpochJd(tle);
        const tsince = (jd1[0] + fr1[0] - satEpoch) * 1440.0;
        const scalar = try sgp4.propagate(tsince);
        const ob = idx * 3;
        try testing.expectApproxEqAbs(scalar[0][0], pos[ob + 0], 0.01);
        try testing.expectApproxEqAbs(scalar[0][1], pos[ob + 1], 0.01);
        try testing.expectApproxEqAbs(scalar[0][2], pos[ob + 2], 0.01);
        try testing.expectApproxEqAbs(scalar[1][0], vel[ob + 0], 1e-5);
        try testing.expectApproxEqAbs(scalar[1][1], vel[ob + 1], 1e-5);
        try testing.expectApproxEqAbs(scalar[1][2], vel[ob + 2], 1e-5);
    }

    const dsIndices = [_]usize{ 1, 3 };
    const dsTles = [_][]const u8{ geoTle, gpsTle };
    for (dsIndices, dsTles) |idx, tleStr| {
        var tle = try Tle.parse(tleStr, testing.allocator);
        defer tle.deinit();
        const sdp4 = try Sdp4.init(tle, constants.wgs72);
        const satEpoch = tleEpochJd(tle);
        const tsince = (jd1[0] + fr1[0] - satEpoch) * 1440.0;
        const scalar = try sdp4.propagate(tsince);
        const ob = idx * 3;
        try testing.expectApproxEqAbs(scalar[0][0], pos[ob + 0], 0.01);
        try testing.expectApproxEqAbs(scalar[0][1], pos[ob + 1], 0.01);
        try testing.expectApproxEqAbs(scalar[0][2], pos[ob + 2], 0.01);
        try testing.expectApproxEqAbs(scalar[1][0], vel[ob + 0], 1e-5);
        try testing.expectApproxEqAbs(scalar[1][1], vel[ob + 1], 1e-5);
        try testing.expectApproxEqAbs(scalar[1][2], vel[ob + 2], 1e-5);
    }
}

test "timeMajor vs satelliteMajor produce identical results" {
    var tles: [5]Tle = undefined;
    defer for (&tles) |*t| t.deinit();
    tles[0] = try Tle.parse(leoTle1, testing.allocator);
    tles[1] = try Tle.parse(geoTle, testing.allocator);
    tles[2] = try Tle.parse(leoTle2, testing.allocator);
    tles[3] = try Tle.parse(gpsTle, testing.allocator);
    tles[4] = try Tle.parse(leoTle3, testing.allocator);

    var c = try Constellation.init(testing.allocator, &tles, constants.wgs72);
    defer c.deinit();

    const epoch = tleEpochJd(tles[0]);
    const jdArr = [_]f64{ @floor(epoch), @floor(epoch) };
    const frArr = [_]f64{ epoch - @floor(epoch), epoch - @floor(epoch) + 1.0 / 1440.0 };

    var posTM: [5 * 2 * 3]f64 = undefined;
    var posSM: [5 * 2 * 3]f64 = undefined;

    c.resetCarry();
    try c.propagate(&jdArr, &frArr, &posTM, null, .teme, .timeMajor);
    c.resetCarry();
    try c.propagate(&jdArr, &frArr, &posSM, null, .teme, .satelliteMajor);

    for (0..2) |t| {
        for (0..5) |s| {
            const tmOb = t * 5 * 3 + s * 3;
            const smOb = s * 2 * 3 + t * 3;
            for (0..3) |j| {
                try testing.expectApproxEqAbs(posTM[tmOb + j], posSM[smOb + j], 1e-10);
            }
        }
    }
}

test "all-SGP4 constellation" {
    var tles: [3]Tle = undefined;
    defer for (&tles) |*t| t.deinit();
    tles[0] = try Tle.parse(leoTle1, testing.allocator);
    tles[1] = try Tle.parse(leoTle2, testing.allocator);
    tles[2] = try Tle.parse(leoTle3, testing.allocator);

    var c = try Constellation.init(testing.allocator, &tles, constants.wgs72);
    defer c.deinit();

    try testing.expectEqual(@as(usize, 3), c.numSgp4);
    try testing.expectEqual(@as(usize, 0), c.numSdp4);

    const epoch = tleEpochJd(tles[0]);
    const jdArr = [_]f64{@floor(epoch)};
    const frArr = [_]f64{epoch - @floor(epoch)};
    var pos: [3 * 3]f64 = undefined;
    try c.propagate(&jdArr, &frArr, &pos, null, .teme, .timeMajor);

    var tle = try Tle.parse(leoTle1, testing.allocator);
    defer tle.deinit();
    const sgp4 = try Sgp4.init(tle, constants.wgs72);
    const scalar = try sgp4.propagate(0.0);
    try testing.expectApproxEqAbs(scalar[0][0], pos[0], 0.01);
    try testing.expectApproxEqAbs(scalar[0][1], pos[1], 0.01);
    try testing.expectApproxEqAbs(scalar[0][2], pos[2], 0.01);
}

test "all-SDP4 constellation" {
    var tles: [2]Tle = undefined;
    defer for (&tles) |*t| t.deinit();
    tles[0] = try Tle.parse(geoTle, testing.allocator);
    tles[1] = try Tle.parse(gpsTle, testing.allocator);

    var c = try Constellation.init(testing.allocator, &tles, constants.wgs72);
    defer c.deinit();

    try testing.expectEqual(@as(usize, 0), c.numSgp4);
    try testing.expectEqual(@as(usize, 2), c.numSdp4);

    const epoch = tleEpochJd(tles[0]);
    const jdArr = [_]f64{@floor(epoch)};
    const frArr = [_]f64{epoch - @floor(epoch)};
    var pos: [2 * 3]f64 = undefined;
    try c.propagate(&jdArr, &frArr, &pos, null, .teme, .timeMajor);

    var tle = try Tle.parse(geoTle, testing.allocator);
    defer tle.deinit();
    const sdp4 = try Sdp4.init(tle, constants.wgs72);
    const scalar = try sdp4.propagate(0.0);
    try testing.expectApproxEqAbs(scalar[0][0], pos[0], 0.01);
    try testing.expectApproxEqAbs(scalar[0][1], pos[1], 0.01);
    try testing.expectApproxEqAbs(scalar[0][2], pos[2], 0.01);
}

test "ECEF output mode" {
    var tles: [3]Tle = undefined;
    defer for (&tles) |*t| t.deinit();
    tles[0] = try Tle.parse(leoTle1, testing.allocator);
    tles[1] = try Tle.parse(geoTle, testing.allocator);
    tles[2] = try Tle.parse(leoTle2, testing.allocator);

    var c = try Constellation.init(testing.allocator, &tles, constants.wgs72);
    defer c.deinit();

    const epoch = tleEpochJd(tles[0]);
    const jdArr = [_]f64{@floor(epoch)};
    const frArr = [_]f64{epoch - @floor(epoch)};

    var posTeme: [3 * 3]f64 = undefined;
    var posEcef: [3 * 3]f64 = undefined;

    c.resetCarry();
    try c.propagate(&jdArr, &frArr, &posTeme, null, .teme, .timeMajor);
    c.resetCarry();
    try c.propagate(&jdArr, &frArr, &posEcef, null, .ecef, .timeMajor);

    const gmst = WCS.julianToGmst(epoch);
    const cosG = @cos(gmst);
    const sinG = @sin(gmst);

    for (0..3) |sat| {
        const ob = sat * 3;
        const ex = posTeme[ob + 0] * cosG + posTeme[ob + 1] * sinG;
        const ey = posTeme[ob + 1] * cosG - posTeme[ob + 0] * sinG;
        try testing.expectApproxEqAbs(ex, posEcef[ob + 0], 1e-6);
        try testing.expectApproxEqAbs(ey, posEcef[ob + 1], 1e-6);
        try testing.expectApproxEqAbs(posTeme[ob + 2], posEcef[ob + 2], 1e-6);
    }
}

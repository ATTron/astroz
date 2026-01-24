//! Multi-batch constellation propagation with threading support
//! Propagates multiple batches of 4 satellites across multiple time points
//! with optional coordinate transformations (TEME, ECEF, Geodetic)

const std = @import("std");
const Sgp4Batch = @import("Sgp4Batch.zig");
const Sgp4 = @import("Sgp4.zig");
const simdMath = @import("simdMath.zig");
const WCS = @import("WorldCoordinateSystem.zig");

const Vec4 = simdMath.Vec4;
const ElementsV4 = Sgp4Batch.ElementsV4;
const Error = Sgp4.Error;

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
    /// (num_times, num_satellites, 3): threads over time ranges (this is typically faster than satelliteMajor)
    timeMajor = 1,
};

/// All the data needed for propagation (passed to threads)
const PropagateCtx = struct {
    batches: []const ElementsV4,
    numSatellites: usize,
    times: []const f64,
    epochOffsets: []const f64,
    resultsPos: []f64,
    resultsVel: ?[]f64,
    referenceJd: f64,
    satelliteMask: ?[]const u8,
};

/// Propagate constellation with specified memory layout and threading strategy
pub fn propagateConstellation(
    batches: []const ElementsV4,
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

    const ctx = PropagateCtx{
        .batches = batches,
        .numSatellites = numSatellites,
        .times = times,
        .epochOffsets = epochOffsets,
        .resultsPos = resultsPos,
        .resultsVel = resultsVel,
        .referenceJd = referenceJd,
        .satelliteMask = satelliteMask,
    };

    // Dispatch: 2 layouts × 3 modes × 2 velocity = 12 variants
    // Key encoding: [layout:1][mode:2][vel:1] = 4 bits
    const key: u4 = @as(u4, @intFromEnum(layout)) << 3 |
        @as(u4, @intCast(@intFromEnum(outputMode))) << 1 |
        @intFromBool(resultsVel != null);

    inline for (0..16) |i| {
        const mode_bits = (i >> 1) & 0b11;
        // filter out the unusable modes
        if (mode_bits < 3 and key == i) {
            const l: Layout = @enumFromInt(i >> 3);
            const m: OutputMode = @enumFromInt(mode_bits);
            const v = (i & 1) == 1;
            return propagateThreaded(l, m, v, ctx);
        }
    }
}

fn propagateThreaded(
    comptime layout: Layout,
    comptime mode: OutputMode,
    comptime hasVel: bool,
    ctx: PropagateCtx,
) Error!void {
    const numTimes = ctx.times.len;
    const maxThreads = getMaxThreads();

    // Thread over batches for satellite major, over time ranges for time major
    const workItems = if (layout == .satelliteMajor) ctx.batches.len else numTimes;

    // Memory tuning: limit threads to avoid cache thrashing
    // Each thread should have enough work to amortize thread overhead
    const minWorkPerThread: usize = 16;
    const effectiveThreads = @max(1, workItems / minWorkPerThread);
    const numThreads = @min(effectiveThreads, @min(workItems, maxThreads));
    const perThread = (workItems + numThreads - 1) / numThreads;
    const maxHandles = 64;
    const clampedThreads = @min(numThreads, maxHandles);
    var handles: [maxHandles]?std.Thread = .{null} ** maxHandles;

    for (0..clampedThreads) |tid| {
        const start = tid * perThread;
        const end = @min(start + perThread, workItems);
        if (start >= end) break;

        handles[tid] = std.Thread.spawn(.{}, propagateRange, .{
            layout,
            mode,
            hasVel,
            ctx,
            start,
            end,
        }) catch {
            propagateRange(
                layout,
                mode,
                hasVel,
                ctx,
                start,
                end,
            );
            continue;
        };
    }

    for (handles) |h| if (h) |t| t.join();
}

fn propagateRange(
    comptime layout: Layout,
    comptime mode: OutputMode,
    comptime hasVel: bool,
    ctx: PropagateCtx,
    rangeStart: usize,
    rangeEnd: usize,
) void {
    const numTimes = ctx.times.len;
    const batchStart = if (layout == .satelliteMajor) rangeStart else 0;
    const batchEnd = if (layout == .satelliteMajor) rangeEnd else ctx.batches.len;
    const timeStart = if (layout == .timeMajor) rangeStart else 0;
    const timeEnd = if (layout == .timeMajor) rangeEnd else numTimes;

    for (ctx.batches[batchStart..batchEnd], batchStart..) |*batch, batchIdx| {
        if (!batchIsActive(batchIdx, ctx.numSatellites, ctx.satelliteMask)) continue;
        const offsets = loadBatchOffsets(batchIdx, ctx.epochOffsets);
        const active = computeActive(batchIdx, ctx.numSatellites, ctx.satelliteMask);

        for (timeStart..timeEnd) |timeIdx| {
            propagateCore(
                layout,
                mode,
                hasVel,
                ctx,
                batch,
                batchIdx,
                timeIdx,
                offsets,
                active,
            );
        }
    }
}

inline fn propagateCore(
    comptime layout: Layout,
    comptime mode: OutputMode,
    comptime hasVel: bool,
    ctx: PropagateCtx,
    batch: *const ElementsV4,
    batchIdx: usize,
    timeIdx: usize,
    batchOffsets: Vec4,
    active: [4]bool,
) void {
    const time = ctx.times[timeIdx];
    const baseTime: Vec4 = @splat(time);
    const satBase = batchIdx * 4;
    const numTimes = ctx.times.len;

    const pv = Sgp4Batch.propagateBatchV4Direct(batch, baseTime + batchOffsets) catch {
        inline for (0..4) |s| {
            if (active[s]) {
                const ob = outBase(layout, satBase + s, timeIdx, numTimes, ctx.numSatellites);
                ctx.resultsPos[ob..][0..3].* = .{ 0, 0, 0 };
                if (hasVel) ctx.resultsVel.?[ob..][0..3].* = .{ 0, 0, 0 };
            }
        }
        return;
    };

    // Compute GMST once for coordinate rotation (comptime eliminated for TEME)
    const gmst = if (mode != .teme) WCS.julianToGmst(ctx.referenceJd + time / 1440.0) else undefined;

    // Transform to output coordinate system
    const pos = if (mode != .teme)
        WCS.eciToEcefV4(pv.rx, pv.ry, pv.rz, gmst)
    else
        .{ .x = pv.rx, .y = pv.ry, .z = pv.rz };

    const vel = if (hasVel) blk: {
        break :blk if (mode != .teme)
            WCS.eciToEcefV4(pv.vx, pv.vy, pv.vz, gmst)
        else
            .{ .x = pv.vx, .y = pv.vy, .z = pv.vz };
    } else undefined;

    // Write results for each active satellite in the batch
    inline for (0..4) |s| {
        if (active[s]) {
            const ob = outBase(layout, satBase + s, timeIdx, numTimes, ctx.numSatellites);
            if (mode == .geodetic) {
                ctx.resultsPos[ob..][0..3].* = WCS.ecefToGeodetic(.{ pos.x[s], pos.y[s], pos.z[s] });
            } else {
                ctx.resultsPos[ob..][0..3].* = .{ pos.x[s], pos.y[s], pos.z[s] };
            }
            if (hasVel) ctx.resultsVel.?[ob..][0..3].* = .{ vel.x[s], vel.y[s], vel.z[s] };
        }
    }
}

inline fn outBase(comptime layout: Layout, satIdx: usize, timeIdx: usize, numTimes: usize, numSatellites: usize) usize {
    return if (layout == .satelliteMajor)
        satIdx * numTimes * 3 + timeIdx * 3
    else
        timeIdx * numSatellites * 3 + satIdx * 3;
}

inline fn batchIsActive(batchIdx: usize, numSatellites: usize, satelliteMask: ?[]const u8) bool {
    if (satelliteMask) |mask| {
        const base = batchIdx * 4;
        var anyActive = false;
        inline for (0..4) |s| {
            if (base + s < numSatellites and mask[base + s] != 0) {
                anyActive = true;
            }
        }
        return anyActive;
    }
    return true;
}

inline fn loadBatchOffsets(batchIdx: usize, epochOffsets: []const f64) Vec4 {
    const base = batchIdx * 4;
    return Vec4{ epochOffsets[base], epochOffsets[base + 1], epochOffsets[base + 2], epochOffsets[base + 3] };
}

inline fn computeActive(batchIdx: usize, numSatellites: usize, satelliteMask: ?[]const u8) [4]bool {
    var active: [4]bool = undefined;
    inline for (0..4) |s| {
        const satIdx = batchIdx * 4 + s;
        active[s] = if (satelliteMask) |mask|
            (satIdx < numSatellites and mask[satIdx] != 0)
        else
            (satIdx < numSatellites);
    }
    return active;
}

fn getMaxThreads() usize {
    if (@hasDecl(std, "c")) {
        const envVal = std.c.getenv("ASTROZ_THREADS");
        if (envVal) |val| {
            return std.fmt.parseInt(usize, std.mem.sliceTo(val, 0), 10) catch
                std.Thread.getCpuCount() catch 1;
        }
    }
    return std.Thread.getCpuCount() catch 1;
}

//! Multi-batch constellation propagation with threading support
//! Propagates multiple batches of N satellites (N=4 or N=8 based on AVX512) across multiple time points
//! with optional coordinate transformations (TEME, ECEF, Geodetic)

const std = @import("std");
const Sgp4Batch = @import("Sgp4Batch.zig");
const Sgp4 = @import("Sgp4.zig");
const simdMath = @import("simdMath.zig");
const WCS = @import("WorldCoordinateSystem.zig");

const Error = Sgp4.Error;

/// Compile-time batch size: 8 for AVX512, 4 for AVX2/SSE
pub const BatchSize: usize = Sgp4Batch.BatchSize;

/// Generic vector type for current batch size
const Vec = simdMath.VecN(BatchSize);

/// Batch elements type for current batch size
const BatchElements = Sgp4Batch.BatchElements(BatchSize);

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
    batches: []const BatchElements,
    numSatellites: usize,
    times: []const f64,
    epochOffsets: []const f64,
    resultsPos: []f64,
    resultsVel: ?[]f64,
    referenceJd: f64,
    satelliteMask: ?[]const u8,
    gmstSin: []const f64, // Precomputed sin(GMST) values
    gmstCos: []const f64, // Precomputed cos(GMST) values
};

/// Propagate constellation with specified memory layout and threading strategy
pub fn propagateConstellation(
    batches: []const BatchElements,
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

    // Precompute GMST sin/cos for non-TEME modes (GMST depends only on time, not satellite)
    const page_alloc = std.heap.page_allocator;
    var sinBuf: []f64 = &.{};
    var cosBuf: []f64 = &.{};
    defer if (sinBuf.len > 0) page_alloc.free(sinBuf);
    defer if (cosBuf.len > 0) page_alloc.free(cosBuf);

    if (outputMode != .teme) {
        sinBuf = page_alloc.alloc(f64, times.len) catch return Error.OutOfMemory;
        cosBuf = page_alloc.alloc(f64, times.len) catch return Error.OutOfMemory;
        for (times, 0..) |time, i| {
            const gmst = WCS.julianToGmst(referenceJd + time / 1440.0);
            sinBuf[i] = @sin(gmst);
            cosBuf[i] = @cos(gmst);
        }
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
        .gmstSin = sinBuf,
        .gmstCos = cosBuf,
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

    const numThreads = @min(workItems, maxThreads);
    const perThread = (workItems + numThreads - 1) / numThreads;
    var handles: [MaxThreads]?std.Thread = .{null} ** MaxThreads;

    for (0..numThreads) |tid| {
        const start = tid * perThread;
        const end = @min(start + perThread, workItems);
        if (start >= end) break;
        handles[tid] = std.Thread.spawn(.{}, propagateRange, .{
            layout, mode, hasVel, ctx, start, end,
        }) catch {
            propagateRange(layout, mode, hasVel, ctx, start, end);
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
    batch: *const BatchElements,
    batchIdx: usize,
    timeIdx: usize,
    batchOffsets: Vec,
    active: [BatchSize]bool,
) void {
    const time = ctx.times[timeIdx];
    const baseTime: Vec = @splat(time);
    const satBase = batchIdx * BatchSize;
    const numTimes = ctx.times.len;

    const pv = Sgp4Batch.propagateBatchDirect(BatchSize, batch, baseTime + batchOffsets) catch {
        inline for (0..BatchSize) |s| {
            if (active[s]) {
                const ob = outBase(layout, satBase + s, timeIdx, numTimes, ctx.numSatellites);
                ctx.resultsPos[ob..][0..3].* = .{ 0, 0, 0 };
                if (hasVel) ctx.resultsVel.?[ob..][0..3].* = .{ 0, 0, 0 };
            }
        }
        return;
    };

    // Transform to output coordinate system using precomputed sin/cos
    const pos = if (mode != .teme)
        eciToEcefDirect(BatchSize, pv.rx, pv.ry, pv.rz, ctx.gmstSin[timeIdx], ctx.gmstCos[timeIdx])
    else
        .{ .x = pv.rx, .y = pv.ry, .z = pv.rz };

    const vel = if (hasVel)
        (if (mode != .teme)
            eciToEcefDirect(BatchSize, pv.vx, pv.vy, pv.vz, ctx.gmstSin[timeIdx], ctx.gmstCos[timeIdx])
        else
            .{ .x = pv.vx, .y = pv.vy, .z = pv.vz })
    else
        undefined;

    // Write results for each active satellite in the batch
    inline for (0..BatchSize) |s| {
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

/// Rotate N ECI position vectors to ECEF using precomputed sin/cos values
inline fn eciToEcefDirect(comptime N: usize, x: simdMath.VecN(N), y: simdMath.VecN(N), z: simdMath.VecN(N), sinG: f64, cosG: f64) WCS.Vec3N(N) {
    const VecT = simdMath.VecN(N);
    const cosVec: VecT = @splat(cosG);
    const sinVec: VecT = @splat(sinG);
    return .{
        .x = x * cosVec + y * sinVec,
        .y = y * cosVec - x * sinVec,
        .z = z,
    };
}

inline fn batchIsActive(batchIdx: usize, numSatellites: usize, satelliteMask: ?[]const u8) bool {
    if (satelliteMask) |mask| {
        const base = batchIdx * BatchSize;
        var anyActive = false;
        inline for (0..BatchSize) |s| {
            if (base + s < numSatellites and mask[base + s] != 0) {
                anyActive = true;
            }
        }
        return anyActive;
    }
    return true;
}

inline fn loadBatchOffsets(batchIdx: usize, epochOffsets: []const f64) Vec {
    const base = batchIdx * BatchSize;
    var arr: [BatchSize]f64 = undefined;
    inline for (0..BatchSize) |i| {
        arr[i] = epochOffsets[base + i];
    }
    return arr;
}

inline fn computeActive(batchIdx: usize, numSatellites: usize, satelliteMask: ?[]const u8) [BatchSize]bool {
    var active: [BatchSize]bool = undefined;
    inline for (0..BatchSize) |s| {
        const satIdx = batchIdx * BatchSize + s;
        active[s] = if (satelliteMask) |mask|
            (satIdx < numSatellites and mask[satIdx] != 0)
        else
            (satIdx < numSatellites);
    }
    return active;
}

/// Cached thread count to avoid repeated syscalls
var cachedMaxThreads: usize = 0;

const MaxThreads: usize = 128;

fn getMaxThreads() usize {
    // Fast path: return cached value
    if (cachedMaxThreads != 0) return cachedMaxThreads;

    // Slow path: compute and cache
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

/// Fused propagate+screen: propagate constellation and screen against target on-the-fly
/// Returns min_distances and min_t_indices per satellite (excluding target)
pub fn screenConstellation(
    batches: []const BatchElements,
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

    // Store squared distances, sqrt at the end
    for (0..numSatellites) |i| {
        outMinDists[i] = thresholdSq; // Only track if closer than threshold
        outMinTIndices[i] = 0;
    }

    // Precompute GMST sin/cos for ECEF output
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
        const targetTime: Vec = @splat(time + targetEpochOffset);
        const targetPv = Sgp4Batch.propagateBatchDirect(BatchSize, targetBatch, targetTime) catch continue;
        const targetPos = eciToEcefDirect(BatchSize, targetPv.rx, targetPv.ry, targetPv.rz, sinBuf[timeIdx], cosBuf[timeIdx]);
        const tArr: [BatchSize]f64 = targetPos.x;
        const tx = tArr[targetLane];
        const tyArr: [BatchSize]f64 = targetPos.y;
        const ty = tyArr[targetLane];
        const tzArr: [BatchSize]f64 = targetPos.z;
        const tz = tzArr[targetLane];

        for (batches, 0..) |*batch, batchIdx| {
            const satBase = batchIdx * BatchSize;
            var batchTimes: Vec = undefined;
            inline for (0..BatchSize) |i| {
                batchTimes[i] = time + epochOffsets[satBase + i];
            }

            const pv = Sgp4Batch.propagateBatchDirect(BatchSize, batch, batchTimes) catch continue;
            const pos = eciToEcefDirect(BatchSize, pv.rx, pv.ry, pv.rz, sinBuf[timeIdx], cosBuf[timeIdx]);
            const xArr: [BatchSize]f64 = pos.x;
            const yArr: [BatchSize]f64 = pos.y;
            const zArr: [BatchSize]f64 = pos.z;

            for (0..BatchSize) |lane| {
                const satIdx = satBase + lane;
                if (satIdx >= numSatellites or satIdx == targetIdx) continue;

                const dx = tx - xArr[lane];
                const dy = ty - yArr[lane];
                const dz = tz - zArr[lane];
                const distSq = dx * dx + dy * dy + dz * dz;

                if (distSq < outMinDists[satIdx]) {
                    outMinDists[satIdx] = distSq;
                    outMinTIndices[satIdx] = @intCast(timeIdx);
                }
            }
        }
    }

    // Convert squared distances to distances
    for (0..numSatellites) |i| {
        outMinDists[i] = @sqrt(outMinDists[i]);
    }
}

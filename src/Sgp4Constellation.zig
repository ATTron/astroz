//! Multi-batch constellation propagation with threading support.
//! Propagates multiple batches of 4 satellites across multiple time points
//! with optional coordinate transformations (TEME, ECEF, Geodetic)

const std = @import("std");
const Sgp4Batch = @import("Sgp4Batch.zig");
const Sgp4 = @import("Sgp4.zig");
const simdMath = @import("simdMath.zig");

const Vec4 = simdMath.Vec4;
const ElementsV4 = Sgp4Batch.ElementsV4;
const Error = Sgp4.Error;

/// Output coordinate modes for propagation
pub const OutputMode = enum(u8) {
    teme = 0,
    ecef = 1,
    geodetic = 2,
};

/// Propagate constellation with per-satellite epoch offsets and coordinate output
pub fn propagateConstellationWithOffsets(
    batches: []const ElementsV4,
    num_satellites: usize,
    times: []const f64,
    epoch_offsets: []const f64,
    results_pos: []f64,
    results_vel: ?[]f64,
    output_mode: OutputMode,
    reference_jd: f64,
    satellite_mask: ?[]const u8,
) Error!void {
    const numTimes = times.len;

    const requiredPosSize = num_satellites * numTimes * 3;
    if (results_pos.len < requiredPosSize) {
        return Error.SatelliteDecayed;
    }
    if (results_vel) |rv| {
        if (rv.len < requiredPosSize) {
            return Error.SatelliteDecayed;
        }
    }

    const has_vel = results_vel != null;
    inline for (.{ .teme, .ecef, .geodetic }) |mode| {
        if (output_mode == mode) {
            if (has_vel) {
                try propagateInner(mode, true, batches, num_satellites, times, epoch_offsets, results_pos, results_vel, reference_jd, satellite_mask);
            } else {
                try propagateInner(mode, false, batches, num_satellites, times, epoch_offsets, results_pos, null, reference_jd, satellite_mask);
            }
        }
    }
}

fn propagateInner(
    comptime mode: OutputMode,
    comptime has_vel: bool,
    batches: []const ElementsV4,
    num_satellites: usize,
    times: []const f64,
    epoch_offsets: []const f64,
    results_pos: []f64,
    results_vel: ?[]f64,
    reference_jd: f64,
    satellite_mask: ?[]const u8,
) Error!void {
    const max_threads = getMaxThreads();
    const num_threads = @min(batches.len, max_threads);

    if (num_threads <= 1) {
        for (batches, 0..) |*batch, batchIdx| {
            propagateBatch(mode, has_vel, batch, batchIdx, num_satellites, times, epoch_offsets, results_pos, results_vel, reference_jd, satellite_mask);
        }
        return;
    }

    const per_thread = (batches.len + num_threads - 1) / num_threads;
    const max_handles = 64;
    const clamped_threads = @min(num_threads, max_handles);
    var handles: [max_handles]?std.Thread = .{null} ** max_handles;
    var spawned: usize = 0;

    for (0..clamped_threads) |tid| {
        const start = tid * per_thread;
        const end = @min(start + per_thread, batches.len);
        if (start >= end) break;
        handles[tid] = std.Thread.spawn(.{}, batchRangeWorker(mode, has_vel), .{
            batches[start..end], start, num_satellites, times, epoch_offsets, results_pos, results_vel, reference_jd, satellite_mask,
        }) catch null;
        if (handles[tid] != null) spawned += 1;
    }

    // If no threads spawned successfully, fall back to sequential
    if (spawned == 0) {
        for (batches, 0..) |*batch, batchIdx| {
            propagateBatch(mode, has_vel, batch, batchIdx, num_satellites, times, epoch_offsets, results_pos, results_vel, reference_jd, satellite_mask);
        }
        return;
    }

    for (&handles) |*h| {
        if (h.*) |t| t.join();
    }
}

fn getMaxThreads() usize {
    if (@hasDecl(std, "c")) {
        const env_val = std.c.getenv("ASTROZ_THREADS");
        if (env_val) |val| {
            return std.fmt.parseInt(usize, std.mem.sliceTo(val, 0), 10) catch
                std.Thread.getCpuCount() catch 1;
        }
    }
    return std.Thread.getCpuCount() catch 1;
}

fn batchRangeWorker(comptime mode: OutputMode, comptime has_vel: bool) fn ([]const ElementsV4, usize, usize, []const f64, []const f64, []f64, ?[]f64, f64, ?[]const u8) void {
    return struct {
        fn worker(batch_slice: []const ElementsV4, start_idx: usize, num_satellites: usize, times: []const f64, epoch_offsets: []const f64, results_pos: []f64, results_vel: ?[]f64, reference_jd: f64, satellite_mask: ?[]const u8) void {
            for (batch_slice, start_idx..) |*batch, batchIdx| {
                propagateBatch(mode, has_vel, batch, batchIdx, num_satellites, times, epoch_offsets, results_pos, results_vel, reference_jd, satellite_mask);
            }
        }
    }.worker;
}

inline fn propagateBatch(
    comptime mode: OutputMode,
    comptime has_vel: bool,
    batch: *const ElementsV4,
    batchIdx: usize,
    num_satellites: usize,
    times: []const f64,
    epoch_offsets: []const f64,
    results_pos: []f64,
    results_vel: ?[]f64,
    reference_jd: f64,
    satellite_mask: ?[]const u8,
) void {
    const WCS = @import("WorldCoordinateSystem.zig");
    const numTimes = times.len;

    if (satellite_mask) |mask| {
        const base = batchIdx * 4;
        var any_active = false;
        inline for (0..4) |s| {
            if (base + s < num_satellites and mask[base + s] != 0) {
                any_active = true;
            }
        }
        if (!any_active) return;
    }

    const offsetBase = batchIdx * 4;
    const batchOffsets = Vec4{
        epoch_offsets[offsetBase],
        epoch_offsets[offsetBase + 1],
        epoch_offsets[offsetBase + 2],
        epoch_offsets[offsetBase + 3],
    };

    var active: [4]bool = undefined;
    inline for (0..4) |s| {
        const satIdx = batchIdx * 4 + s;
        active[s] = if (satellite_mask) |mask|
            (satIdx < num_satellites and mask[satIdx] != 0)
        else
            (satIdx < num_satellites);
    }

    for (0..numTimes) |timeIdx| {
        const baseTime: Vec4 = @splat(times[timeIdx]);
        const satResults = Sgp4Batch.propagateSatellitesV4Vec(batch, baseTime + batchOffsets) catch {
            inline for (0..4) |s| {
                if (active[s]) {
                    const outBase = (batchIdx * 4 + s) * numTimes * 3 + timeIdx * 3;
                    results_pos[outBase] = 0.0;
                    results_pos[outBase + 1] = 0.0;
                    results_pos[outBase + 2] = 0.0;
                    if (has_vel) {
                        results_vel.?[outBase] = 0.0;
                        results_vel.?[outBase + 1] = 0.0;
                        results_vel.?[outBase + 2] = 0.0;
                    }
                }
            }
            continue;
        };

        var pos: [4][3]f64 = undefined;
        var vel: [4][3]f64 = undefined;

        if (mode == .ecef or mode == .geodetic) {
            const gmst = WCS.julianToGmst(reference_jd + times[timeIdx] / 1440.0);
            const ecef = WCS.eciToEcefV4(
                .{ satResults[0][0][0], satResults[1][0][0], satResults[2][0][0], satResults[3][0][0] },
                .{ satResults[0][0][1], satResults[1][0][1], satResults[2][0][1], satResults[3][0][1] },
                .{ satResults[0][0][2], satResults[1][0][2], satResults[2][0][2], satResults[3][0][2] },
                gmst,
            );
            inline for (0..4) |s| {
                if (mode == .geodetic) {
                    pos[s] = WCS.ecefToGeodetic(.{ ecef.x[s], ecef.y[s], ecef.z[s] });
                } else {
                    pos[s] = .{ ecef.x[s], ecef.y[s], ecef.z[s] };
                }
            }
            if (has_vel) {
                const ve = WCS.eciToEcefV4(
                    .{ satResults[0][1][0], satResults[1][1][0], satResults[2][1][0], satResults[3][1][0] },
                    .{ satResults[0][1][1], satResults[1][1][1], satResults[2][1][1], satResults[3][1][1] },
                    .{ satResults[0][1][2], satResults[1][1][2], satResults[2][1][2], satResults[3][1][2] },
                    gmst,
                );
                inline for (0..4) |s| {
                    vel[s] = .{ ve.x[s], ve.y[s], ve.z[s] };
                }
            }
        } else {
            inline for (0..4) |s| {
                pos[s] = satResults[s][0];
                if (has_vel) vel[s] = satResults[s][1];
            }
        }

        inline for (0..4) |s| {
            if (active[s]) {
                const outBase = (batchIdx * 4 + s) * numTimes * 3 + timeIdx * 3;
                results_pos[outBase] = pos[s][0];
                results_pos[outBase + 1] = pos[s][1];
                results_pos[outBase + 2] = pos[s][2];
                if (has_vel) {
                    results_vel.?[outBase] = vel[s][0];
                    results_vel.?[outBase + 1] = vel[s][1];
                    results_vel.?[outBase + 2] = vel[s][2];
                }
            }
        }
    }
}

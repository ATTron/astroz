//! SGP4 C API exports

const std = @import("std");
const astroz = @import("astroz");
const Sgp4 = astroz.Sgp4;
const Tle = astroz.Tle;
const constants = astroz.constants;
const dispatch = astroz.dispatch;

/// Matches the lane count of `dispatch.sgp4Times8`.
const BatchSize = 8;

const err = @import("error.zig");
const tleApi = @import("tle.zig");

pub const Handle = *anyopaque;

const gpa = std.heap.page_allocator;

fn propagateErr(e: Sgp4.Error) err.Code {
    return switch (e) {
        error.DeepSpaceNotSupported => .deepSpaceNotSupported,
        error.InvalidEccentricity => .invalidEccentricity,
        error.NegativeMeanMotion => .negativeMeanMotion,
        error.PerturbedEccentricity => .perturbedEccentricity,
        error.NegativeSemilatusRectum => .negativeSemilatusRectum,
        error.SatelliteDecayed => .satelliteDecayed,
        error.OutOfMemory => .allocFailed,
    };
}

/// gravity model: 0 = WGS84 (default), 1 = WGS72
pub fn init(tleHandle: tleApi.Handle, gravModel: i32, out: *Handle) err.Code {
    const tlePtr: *Tle = @ptrCast(@alignCast(tleHandle));
    const grav = if (gravModel == 1) constants.wgs72 else constants.wgs84;

    const sgp4 = Sgp4.init(tlePtr.*, grav) catch |e| return propagateErr(e);

    const ptr = gpa.create(Sgp4) catch return .allocFailed;
    ptr.* = sgp4;
    out.* = @ptrCast(ptr);
    return .ok;
}

pub fn free(handle: Handle) void {
    const ptr: *Sgp4 = @ptrCast(@alignCast(handle));
    gpa.destroy(ptr);
}

/// returns position (km) and velocity (km/s)
pub fn propagate(handle: Handle, tsince: f64, pos: *[3]f64, vel: *[3]f64) err.Code {
    const ptr: *Sgp4 = @ptrCast(@alignCast(handle));
    const result = ptr.propagate(tsince) catch |e| return propagateErr(e);
    pos.* = result[0];
    vel.* = result[1];
    return .ok;
}

/// Propagates `count` times (minutes since epoch) in SIMD batches of `BatchSize`.
/// `results` receives `[x, y, z, vx, vy, vz]` per time.
pub fn propagateBatch(handle: Handle, times: [*]const f64, results: [*]f64, count: u32) err.Code {
    const ptr: *Sgp4 = @ptrCast(@alignCast(handle));

    var i: usize = 0;
    while (i < count) : (i += BatchSize) {
        const remaining = count - i;
        var batch: [BatchSize]f64 = undefined;
        inline for (0..BatchSize) |j| batch[j] = times[i + @min(j, remaining - 1)];

        const batchResult = dispatch.sgp4Times8(ptr, batch) catch |e| return propagateErr(e);

        const writeCount = @min(remaining, BatchSize);
        for (0..writeCount) |j| {
            const base = (i + j) * 6;
            results[base + 0] = batchResult[j][0][0];
            results[base + 1] = batchResult[j][0][1];
            results[base + 2] = batchResult[j][0][2];
            results[base + 3] = batchResult[j][1][0];
            results[base + 4] = batchResult[j][1][1];
            results[base + 5] = batchResult[j][1][2];
        }
    }

    return .ok;
}

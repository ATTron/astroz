//! SGP4 C API exports

const astroz = @import("astroz");
const Sgp4 = astroz.Sgp4;
const Tle = astroz.Tle;
const constants = astroz.constants;

const allocator = @import("allocator.zig");
const err = @import("error.zig");
const tle_api = @import("tle.zig");

pub const Handle = *anyopaque;

/// gravity model: 0 = WGS84 (default), 1 = WGS72
pub fn init(tle_handle: tle_api.Handle, grav_model: i32, out: *Handle) err.Code {
    const tle_ptr: *Tle = @ptrCast(@alignCast(tle_handle));
    const grav = if (grav_model == 1) constants.wgs72 else constants.wgs84;

    const sgp4 = Sgp4.init(tle_ptr.*, grav) catch |e| {
        return switch (e) {
            Sgp4.Error.DeepSpaceNotSupported => .deep_space_not_supported,
            Sgp4.Error.InvalidEccentricity => .invalid_eccentricity,
            Sgp4.Error.SatelliteDecayed => .satellite_decayed,
        };
    };

    const ptr = allocator.get().create(Sgp4) catch return .alloc_failed;
    ptr.* = sgp4;
    out.* = @ptrCast(ptr);
    return .ok;
}

pub fn free(handle: Handle) void {
    const ptr: *Sgp4 = @ptrCast(@alignCast(handle));
    allocator.get().destroy(ptr);
}

/// returns position (km) and velocity (km/s)
pub fn propagate(handle: Handle, tsince: f64, pos: *[3]f64, vel: *[3]f64) err.Code {
    const ptr: *Sgp4 = @ptrCast(@alignCast(handle));
    const result = ptr.propagate(tsince) catch |e| {
        return switch (e) {
            Sgp4.Error.SatelliteDecayed => .satellite_decayed,
            else => err.fromError(e),
        };
    };
    pos.* = result[0];
    vel.* = result[1];
    return .ok;
}

/// Batch propagation - single FFI call for multiple time steps
/// Uses SIMD which will process 4 time steps at a time
/// times: array of tsince values (minutes since epoch)
/// results: output array of [pos_x, pos_y, pos_z, vel_x, vel_y, vel_z] for each time
/// count: number of time steps
pub fn propagateBatch(handle: Handle, times: [*]const f64, results: [*]f64, count: u32) err.Code {
    const ptr: *Sgp4 = @ptrCast(@alignCast(handle));

    // 4 times per cycle using SIMD
    const batches = count / 4;
    for (0..batches) |batch| {
        const baseTime = batch * 4;
        const batchTimes: [4]f64 = times[baseTime..][0..4].*;

        const batchResult = ptr.propagateV4(batchTimes) catch |e| {
            return switch (e) {
                Sgp4.Error.SatelliteDecayed => .satellite_decayed,
                else => err.fromError(e),
            };
        };

        for (0..4) |i| {
            const base = (baseTime + i) * 6;
            results[base + 0] = batchResult[i][0][0];
            results[base + 1] = batchResult[i][0][1];
            results[base + 2] = batchResult[i][0][2];
            results[base + 3] = batchResult[i][1][0];
            results[base + 4] = batchResult[i][1][1];
            results[base + 5] = batchResult[i][1][2];
        }
    }

    // handle leftovers with scalar functions
    const remainderStart = batches * 4;
    for (remainderStart..count) |i| {
        const result = ptr.propagate(times[i]) catch |e| {
            return switch (e) {
                Sgp4.Error.SatelliteDecayed => .satellite_decayed,
                else => err.fromError(e),
            };
        };

        const base = i * 6;
        results[base + 0] = result[0][0];
        results[base + 1] = result[0][1];
        results[base + 2] = result[0][2];
        results[base + 3] = result[1][0];
        results[base + 4] = result[1][1];
        results[base + 5] = result[1][2];
    }

    return .ok;
}

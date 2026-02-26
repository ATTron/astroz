//! SGP4 C API exports

const astroz = @import("astroz");
const Sgp4 = astroz.Sgp4;
const Tle = astroz.Tle;
const constants = astroz.constants;
const dispatch = astroz.dispatch;

const BatchSize = astroz.Constellation.BatchSize;

const allocator = @import("allocator.zig");
const err = @import("error.zig");
const tleApi = @import("tle.zig");

pub const Handle = *anyopaque;

/// gravity model: 0 = WGS84 (default), 1 = WGS72
pub fn init(tleHandle: tleApi.Handle, gravModel: i32, out: *Handle) err.Code {
    const tlePtr: *Tle = @ptrCast(@alignCast(tleHandle));
    const grav = if (gravModel == 1) constants.wgs72 else constants.wgs84;

    const sgp4 = Sgp4.init(tlePtr.*, grav) catch |e| {
        return switch (e) {
            Sgp4.Error.DeepSpaceNotSupported => .deepSpaceNotSupported,
            Sgp4.Error.InvalidEccentricity => .invalidEccentricity,
            Sgp4.Error.SatelliteDecayed => .satelliteDecayed,
        };
    };

    const ptr = allocator.get().create(Sgp4) catch return .allocFailed;
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
            Sgp4.Error.SatelliteDecayed => .satelliteDecayed,
            else => err.fromError(e),
        };
    };
    pos.* = result[0];
    vel.* = result[1];
    return .ok;
}

/// Batch propagation - single FFI call for multiple time steps
/// Uses SIMD to process BatchSize time steps at a time
/// times: array of tsince values (minutes since epoch)
/// results: output array of [posX, posY, posZ, velX, velY, velZ] for each time
/// count: number of time steps
pub fn propagateBatch(handle: Handle, times: [*]const f64, results: [*]f64, count: u32) err.Code {
    const ptr: *Sgp4 = @ptrCast(@alignCast(handle));

    var i: usize = 0;
    while (i < count) : (i += BatchSize) {
        const remaining = count - i;
        var batch: [BatchSize]f64 = undefined;
        inline for (0..BatchSize) |j| batch[j] = times[i + @min(j, remaining - 1)];

        const batchResult = dispatch.sgp4Times8(ptr, batch) catch |e| {
            return switch (e) {
                Sgp4.Error.SatelliteDecayed => .satelliteDecayed,
                else => err.fromError(e),
            };
        };

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

// Multi-Satellite Batch API
pub const BatchHandle = *anyopaque;

const Elements8 = astroz.Constellation.Sgp4Batch.BatchElements(8);
const PVA = Sgp4.PosVelArray(8);

/// Initialize batch propagator for up to 4 satellites (padded to 8 for SIMD)
/// All satellites must use the same gravity model
pub fn initBatch(tleHandles: [*]const tleApi.Handle, gravModel: i32, out: *BatchHandle) err.Code {
    const grav = if (gravModel == 1) constants.wgs72 else constants.wgs84;

    var tles: [8]Tle = undefined;
    for (0..4) |i| {
        const tlePtr: *Tle = @ptrCast(@alignCast(tleHandles[i]));
        tles[i] = tlePtr.*;
    }
    // Pad unused lanes with last real TLE
    for (4..8) |i| tles[i] = tles[3];

    const elements = astroz.Constellation.Sgp4Batch.initBatchElements(8, tles, grav) catch |e| {
        return switch (e) {
            Sgp4.Error.DeepSpaceNotSupported => .deepSpaceNotSupported,
            Sgp4.Error.InvalidEccentricity => .invalidEccentricity,
            Sgp4.Error.SatelliteDecayed => .satelliteDecayed,
        };
    };

    const ptr = allocator.get().create(Elements8) catch return .allocFailed;
    ptr.* = elements;
    out.* = @ptrCast(ptr);
    return .ok;
}

/// Free batch propagator handle
pub fn freeBatch(handle: BatchHandle) void {
    const ptr: *Elements8 = @ptrCast(@alignCast(handle));
    allocator.get().destroy(ptr);
}

/// Propagate 4 satellites at a single time point (via dispatch)
pub fn propagateSatellites(handle: BatchHandle, tsince: f64, results: [*]f64) err.Code {
    const ptr: *Elements8 = @ptrCast(@alignCast(handle));

    var times: [8]f64 = undefined;
    @memset(&times, tsince);
    const result: PVA = dispatch.sgp4Batch8(ptr, times) catch |e| {
        return switch (e) {
            Sgp4.Error.SatelliteDecayed => .satelliteDecayed,
            else => err.fromError(e),
        };
    };

    for (0..4) |i| {
        results[i * 6 ..][0..3].* = result[i][0];
        results[i * 6 + 3 ..][0..3].* = result[i][1];
    }

    return .ok;
}

/// Propagate 4 satellites across multiple time points (via dispatch)
pub fn propagateSatellitesBatch(handle: BatchHandle, times: [*]const f64, results: [*]f64, count: u32) err.Code {
    const ptr: *Elements8 = @ptrCast(@alignCast(handle));

    for (0..count) |t| {
        var tArr: [8]f64 = undefined;
        @memset(&tArr, times[t]);
        const result: PVA = dispatch.sgp4Batch8(ptr, tArr) catch |e| {
            return switch (e) {
                Sgp4.Error.SatelliteDecayed => .satelliteDecayed,
                else => err.fromError(e),
            };
        };

        const timeBase = t * 24;
        for (0..4) |i| {
            results[timeBase + i * 6 ..][0..3].* = result[i][0];
            results[timeBase + i * 6 + 3 ..][0..3].* = result[i][1];
        }
    }

    return .ok;
}

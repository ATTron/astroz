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
            else => err.fromError(e),
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

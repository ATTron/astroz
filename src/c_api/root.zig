//! C API for astroz - Python/FFI bindings
//! Build with: zig build c-api

const allocator = @import("allocator.zig");
pub const err = @import("error.zig");
const tle = @import("tle.zig");
const sgp4 = @import("sgp4.zig");

const version = .{ .major = 0, .minor = 3, .patch = 0 };

pub export fn astroz_version() callconv(.c) u32 {
    return (@as(u32, version.major) << 16) | (@as(u32, version.minor) << 8) | version.patch;
}

pub export fn astroz_init() callconv(.c) void {
    allocator.init();
}

pub export fn astroz_deinit() callconv(.c) void {
    allocator.deinit();
}

pub export fn tle_parse(str: [*:0]const u8, out: *tle.Handle) callconv(.c) i32 {
    return @intFromEnum(tle.parse(str, out));
}
pub export fn tle_free(handle: tle.Handle) callconv(.c) void {
    tle.free(handle);
}
pub export fn tle_get_satellite_number(h: tle.Handle) callconv(.c) u32 {
    return tle.getSatelliteNumber(h);
}
pub export fn tle_get_epoch(h: tle.Handle) callconv(.c) f64 {
    return tle.getEpoch(h);
}
pub export fn tle_get_inclination(h: tle.Handle) callconv(.c) f32 {
    return tle.getInclination(h);
}
pub export fn tle_get_eccentricity(h: tle.Handle) callconv(.c) f32 {
    return tle.getEccentricity(h);
}
pub export fn tle_get_mean_motion(h: tle.Handle) callconv(.c) f64 {
    return tle.getMeanMotion(h);
}

pub export fn sgp4_init(tle_h: tle.Handle, grav: i32, out: *sgp4.Handle) callconv(.c) i32 {
    return @intFromEnum(sgp4.init(tle_h, grav, out));
}
pub export fn sgp4_free(handle: sgp4.Handle) callconv(.c) void {
    sgp4.free(handle);
}
pub export fn sgp4_propagate(h: sgp4.Handle, tsince: f64, pos: *[3]f64, vel: *[3]f64) callconv(.c) i32 {
    return @intFromEnum(sgp4.propagate(h, tsince, pos, vel));
}

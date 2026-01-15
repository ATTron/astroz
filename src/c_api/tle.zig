//! TLE C API exports

const std = @import("std");
const astroz = @import("astroz");
const Tle = astroz.Tle;

const allocator = @import("allocator.zig");
const err = @import("error.zig");

pub const Handle = *anyopaque;

pub fn parse(tle_str: [*:0]const u8, out: *Handle) err.Code {
    const tle = Tle.parse(std.mem.span(tle_str), allocator.get()) catch |e| {
        return switch (e) {
            Tle.Error.BadTleLength => .bad_tle_length,
            else => err.fromError(e),
        };
    };
    const ptr = allocator.get().create(Tle) catch return .alloc_failed;
    ptr.* = tle;
    out.* = @ptrCast(ptr);
    return .ok;
}

pub fn free(handle: Handle) void {
    const ptr: *Tle = @ptrCast(@alignCast(handle));
    ptr.deinit();
    allocator.get().destroy(ptr);
}

pub fn getSatelliteNumber(handle: Handle) u32 {
    const ptr: *Tle = @ptrCast(@alignCast(handle));
    return ptr.firstLine.satelliteNumber;
}

pub fn getEpoch(handle: Handle) f64 {
    const ptr: *Tle = @ptrCast(@alignCast(handle));
    return ptr.firstLine.epoch;
}

pub fn getInclination(handle: Handle) f32 {
    const ptr: *Tle = @ptrCast(@alignCast(handle));
    return ptr.secondLine.inclination;
}

pub fn getRightAscension(handle: Handle) f32 {
    const ptr: *Tle = @ptrCast(@alignCast(handle));
    return ptr.secondLine.rightAscension;
}

pub fn getEccentricity(handle: Handle) f32 {
    const ptr: *Tle = @ptrCast(@alignCast(handle));
    return ptr.secondLine.eccentricity;
}

pub fn getMeanMotion(handle: Handle) f64 {
    const ptr: *Tle = @ptrCast(@alignCast(handle));
    return ptr.secondLine.mMotion;
}

pub fn getBstar(handle: Handle) f32 {
    const ptr: *Tle = @ptrCast(@alignCast(handle));
    return ptr.firstLine.bstarDrag;
}

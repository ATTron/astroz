//! TLE C API exports

const std = @import("std");
const astroz = @import("astroz");
const Tle = astroz.Tle;

const gpa = std.heap.page_allocator;
const err = @import("error.zig");

pub const Handle = *anyopaque;

pub fn parse(tleStr: [*:0]const u8, out: *Handle) err.Code {
    const tle = Tle.parse(std.mem.span(tleStr), gpa) catch |e| {
        return switch (e) {
            Tle.Error.BadTleLength => .badTleLength,
            else => err.fromError(e),
        };
    };
    const ptr = gpa.create(Tle) catch return .allocFailed;
    ptr.* = tle;
    out.* = @ptrCast(ptr);
    return .ok;
}

pub fn free(handle: Handle) void {
    const ptr: *Tle = @ptrCast(@alignCast(handle));
    ptr.deinit();
    gpa.destroy(ptr);
}

pub fn getSatelliteNumber(handle: Handle) u32 {
    const ptr: *Tle = @ptrCast(@alignCast(handle));
    return ptr.satelliteNumber;
}

pub fn getEpoch(handle: Handle) f64 {
    const ptr: *Tle = @ptrCast(@alignCast(handle));
    return ptr.epoch;
}

pub fn getInclination(handle: Handle) f64 {
    const ptr: *Tle = @ptrCast(@alignCast(handle));
    return ptr.inclination;
}

pub fn getEccentricity(handle: Handle) f64 {
    const ptr: *Tle = @ptrCast(@alignCast(handle));
    return ptr.eccentricity;
}

pub fn getMeanMotion(handle: Handle) f64 {
    const ptr: *Tle = @ptrCast(@alignCast(handle));
    return ptr.mMotion;
}

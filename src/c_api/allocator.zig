//! Global allocator for C API callers
//! Python/FFI users don't deal with allocators ;)

const std = @import("std");

var gpa: ?std.heap.GeneralPurposeAllocator(.{}) = null;

pub fn init() void {
    if (gpa == null) {
        gpa = std.heap.GeneralPurposeAllocator(.{}){};
    }
}

pub fn deinit() void {
    if (gpa) |*g| {
        _ = g.deinit();
        gpa = null;
    }
}

pub fn get() std.mem.Allocator {
    return if (gpa) |*g| g.allocator() else @panic("astroz_init() not called");
}

//! C API for astroz - Python/FFI bindings
//! Build with: zig build c-api

const allocator = @import("allocator.zig");
pub const err = @import("error.zig");

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

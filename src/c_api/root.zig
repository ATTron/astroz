//! C API for astroz - Python/FFI bindings
//! Build with: zig build c-api

const version = .{ .major = 0, .minor = 3, .patch = 0 };

pub export fn astroz_version() callconv(.c) u32 {
    return (@as(u32, version.major) << 16) | (@as(u32, version.minor) << 8) | version.patch;
}

//! Global allocator for C API callers
//! Python/FFI users don't deal with allocators ;)

const std = @import("std");

// the OS reclaims all memory on process exit anyway

pub fn init() void {
    // no-op: page_allocator doesn't need initialization
}

pub fn deinit() void {
    // no-op: let OS reclaim on process exit
}

pub fn get() std.mem.Allocator {
    return std.heap.page_allocator;
}

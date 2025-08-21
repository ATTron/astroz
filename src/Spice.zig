//! SPICE file structure
const std = @import("std");
const log = std.log;

pub const Spice = @This();

pub const LSK = struct {
    data: []u64,
};

lsk: ?LSK,

pub fn new(lsk: LSK) Spice {
    return .{
        .lsk = lsk,
    };
}

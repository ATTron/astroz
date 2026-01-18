//! SIMD Math for SGP4

const std = @import("std");
const constants = @import("../constants.zig");

const simdMath = @This();

// 4 computes per CPU tick
pub const Vec4 = @Vector(4, f64);

pub const twoPiVec: Vec4 = @splat(constants.twoPi);
pub const piVec: Vec4 = @splat(std.math.pi);
pub const halfPiVec: Vec4 = @splat(std.math.pi / 2.0);
pub const invTwoPiVec: Vec4 = @splat(1.0 / constants.twoPi);

pub fn sinSIMD(x: Vec4) Vec4 {
    return Vec4{
        @sin(x[0]),
        @sin(x[1]),
        @sin(x[2]),
        @sin(x[3]),
    };
}

pub fn cosSIMD(x: Vec4) Vec4 {
    return Vec4{
        @cos(x[0]),
        @cos(x[1]),
        @cos(x[2]),
        @cos(x[3]),
    };
}

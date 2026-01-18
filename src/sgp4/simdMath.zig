//! SIMD Math for SGP4

const std = @import("std");
const constants = @import("../constants.zig");

const simdMath = @This();

// compute 4 time segements per CPU tick
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

pub fn atan2SIMD(y: Vec4, x: Vec4) Vec4 {
    return Vec4{
        std.math.atan2(y[0], x[0]),
        std.math.atan2(y[1], x[1]),
        std.math.atan2(y[2], x[2]),
        std.math.atan2(y[3], x[3]),
    };
}

pub fn pow15V4(x: Vec4) Vec4 {
    return x * @sqrt(x);
}

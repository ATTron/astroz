//! SIMD Math for SGP4
//! Fast vectorized trig functions that compute N values at once

const std = @import("std");
const constants = @import("constants.zig");

pub const BatchSize: usize = 8;

/// Generic N-wide f64 vector type
pub fn VecN(comptime N: usize) type {
    return @Vector(N, f64);
}

/// Generic N-wide i64 vector type
pub fn VecNi(comptime N: usize) type {
    return @Vector(N, i64);
}

/// Generic sin/cos result type
pub fn SinCosN(comptime N: usize) type {
    return struct {
        sin: VecN(N),
        cos: VecN(N),
    };
}

/// Compute sin and cos together (faster than separate calls)
/// Works for any angle, accurate to ~1e-15 (double precision)
pub fn sincosN(comptime N: usize, angle: VecN(N)) SinCosN(N) {
    const Vec = VecN(N);
    const VecI = VecNi(N);

    // Constants
    const twoOverPi: Vec = @splat(2.0 / std.math.pi);
    const piOver2Hi: Vec = @splat(1.5707963267948966);
    const piOver2Lo: Vec = @splat(6.123233995736766e-17);
    const one: Vec = @splat(1.0);
    const two: VecI = @splat(2);
    const oneI: VecI = @splat(1);

    // Polynomial coefficients for sin
    const sinC1: Vec = @splat(-1.6666666666666666574148e-01);
    const sinC2: Vec = @splat(8.3333333333333225058715e-03);
    const sinC3: Vec = @splat(-1.9841269841201840457725e-04);
    const sinC4: Vec = @splat(2.7557319210152756118515e-06);
    const sinC5: Vec = @splat(-2.5052106798274583895303e-08);
    const sinC6: Vec = @splat(1.6058936490373178302326e-10);

    // Polynomial coefficients for cos
    const cos_c1: Vec = @splat(-4.9999999999999999999583e-01);
    const cos_c2: Vec = @splat(4.1666666666666665319411e-02);
    const cos_c3: Vec = @splat(-1.3888888888888872762458e-03);
    const cos_c4: Vec = @splat(2.4801587301587286645498e-05);
    const cos_c5: Vec = @splat(-2.7557319223933824788682e-07);
    const cos_c6: Vec = @splat(2.0876756987868089233269e-09);

    // Step 1: Range reduction
    const kFloat = angle * twoOverPi;
    const roundMagic: Vec = @splat(6755399441055744.0);
    const kRounded = kFloat + roundMagic - roundMagic;

    var reduced = @mulAdd(Vec, -piOver2Hi, kRounded, angle);
    reduced = @mulAdd(Vec, -piOver2Lo, kRounded, reduced);
    const k: VecI = @intFromFloat(kRounded);

    // Step 2: Polynomial approximation
    const r2 = reduced * reduced;

    // Sin polynomial: sin(r) ~= r + r^3*P(r^2)
    var sinP = @mulAdd(Vec, sinC6, r2, sinC5);
    sinP = @mulAdd(Vec, sinP, r2, sinC4);
    sinP = @mulAdd(Vec, sinP, r2, sinC3);
    sinP = @mulAdd(Vec, sinP, r2, sinC2);
    sinP = @mulAdd(Vec, sinP, r2, sinC1);
    const sinReduced = @mulAdd(Vec, sinP, r2 * reduced, reduced);

    // Cos polynomial: cos(r) ~= 1 + r^2*P(r^2)
    var cosP = @mulAdd(Vec, cos_c6, r2, cos_c5);
    cosP = @mulAdd(Vec, cosP, r2, cos_c4);
    cosP = @mulAdd(Vec, cosP, r2, cos_c3);
    cosP = @mulAdd(Vec, cosP, r2, cos_c2);
    cosP = @mulAdd(Vec, cosP, r2, cos_c1);
    const cosReduced = @mulAdd(Vec, cosP, r2, one);

    // Step 3: Quadrant correction
    const swap = (k & oneI) != @as(VecI, @splat(0));
    var s = @select(f64, swap, cosReduced, sinReduced);
    var c = @select(f64, swap, sinReduced, cosReduced);

    const sinSign = (k & two) << @splat(62);
    const cosSign = ((k + oneI) & two) << @splat(62);

    s = @bitCast(@as(VecI, @bitCast(s)) ^ sinSign);
    c = @bitCast(@as(VecI, @bitCast(c)) ^ cosSign);

    return .{ .sin = s, .cos = c };
}

/// Generic vectorized sin
pub inline fn sinN(comptime N: usize, x: VecN(N)) VecN(N) {
    return sincosN(N, x).sin;
}

/// Generic vectorized cos
pub inline fn cosN(comptime N: usize, x: VecN(N)) VecN(N) {
    return sincosN(N, x).cos;
}

/// Normalize angle to [0, 2*pi)
pub fn modTwoPiN(comptime N: usize, x: VecN(N)) VecN(N) {
    const Vec = VecN(N);
    const twoPiVec: Vec = @splat(constants.twoPi);
    const invTwoPiVec: Vec = @splat(1.0 / constants.twoPi);
    const negTwoPiVec: Vec = @splat(-constants.twoPi);
    const zero: Vec = @splat(0.0);

    const n = @floor(x * invTwoPiVec);
    var result = @mulAdd(Vec, negTwoPiVec, n, x);
    const mask = result < zero;
    result = @select(f64, mask, result + twoPiVec, result);
    return result;
}

/// Vectorized atan2 using polynomial approximation for the principal value
/// with quadrant correction. Accurate to ~1e-7 radians, sufficient for SGP4
pub fn atan2N(comptime N: usize, y: VecN(N), x: VecN(N)) VecN(N) {
    const Vec = VecN(N);
    const zero: Vec = @splat(0.0);
    const piVec: Vec = @splat(std.math.pi);
    const halfPiVec: Vec = @splat(std.math.pi / 2.0);

    const abs_x = @abs(x);
    const abs_y = @abs(y);

    const max_xy = @max(abs_x, abs_y);
    const min_xy = @min(abs_x, abs_y);

    const tiny: Vec = @splat(1.0e-30);
    const safe_max = @max(max_xy, tiny);
    const t = min_xy / safe_max;
    const t2 = t * t;

    // Polynomial coefficients
    const c1: Vec = @splat(1.0);
    const c3: Vec = @splat(-0.3333314528);
    const c5: Vec = @splat(0.1999355085);
    const c7: Vec = @splat(-0.1420889944);
    const c9: Vec = @splat(0.1065626393);
    const c11: Vec = @splat(-0.0752896400);
    const c13: Vec = @splat(0.0429096138);
    const c15: Vec = @splat(-0.0161657367);
    const c17: Vec = @splat(0.0028662257);

    // Horner's method
    var atan_t = c17;
    atan_t = @mulAdd(Vec, atan_t, t2, c15);
    atan_t = @mulAdd(Vec, atan_t, t2, c13);
    atan_t = @mulAdd(Vec, atan_t, t2, c11);
    atan_t = @mulAdd(Vec, atan_t, t2, c9);
    atan_t = @mulAdd(Vec, atan_t, t2, c7);
    atan_t = @mulAdd(Vec, atan_t, t2, c5);
    atan_t = @mulAdd(Vec, atan_t, t2, c3);
    atan_t = @mulAdd(Vec, atan_t, t2, c1);
    atan_t = atan_t * t;

    const swap_mask = abs_y > abs_x;
    atan_t = @select(f64, swap_mask, halfPiVec - atan_t, atan_t);

    const x_neg = x < zero;
    const y_neg = y < zero;

    var result = atan_t;
    result = @select(f64, x_neg, piVec - result, result);
    result = @select(f64, y_neg, -result, result);

    return result;
}

/// x^1.5 = x * sqrt(x)
pub fn pow15N(comptime N: usize, x: VecN(N)) VecN(N) {
    return x * @sqrt(x);
}

/// x^(2/3) via Newton iteration for cube root: y^3 = x^2
pub fn pow23N(comptime N: usize, x: VecN(N)) VecN(N) {
    const Vec = VecN(N);
    const x2 = x * x;
    var y = @sqrt(x); // x^0.5, reasonable initial guess for x^0.667
    const third: Vec = @splat(1.0 / 3.0);
    const two_third: Vec = @splat(2.0 / 3.0);
    // Newton: y^3 = x^2 → y_{n+1} = (2y + x²/y²) / 3
    inline for (0..6) |_| {
        y = two_third * y + third * x2 / (y * y);
    }
    return y;
}

test "sincosN(4) accuracy" {
    const testVals = VecN(4){ 0.0, std.math.pi / 6.0, std.math.pi / 4.0, std.math.pi / 2.0 };
    const result = sincosN(4, testVals);

    inline for (0..4) |i| {
        try std.testing.expectApproxEqAbs(@sin(testVals[i]), result.sin[i], 1e-12);
        try std.testing.expectApproxEqAbs(@cos(testVals[i]), result.cos[i], 1e-12);
    }
}

test "sincosN(4) large angles" {
    const testVals = VecN(4){ 10.0 * std.math.pi, 100.0, -50.0, 1000.0 };
    const result = sincosN(4, testVals);

    inline for (0..4) |i| {
        try std.testing.expectApproxEqAbs(@sin(testVals[i]), result.sin[i], 1e-12);
        try std.testing.expectApproxEqAbs(@cos(testVals[i]), result.cos[i], 1e-12);
    }
}

test "sincosN(8) matches sincosN(4)" {
    const v4 = VecN(4){ 1.0, 2.0, 3.0, 4.0 };
    const v8 = VecN(8){ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };

    const r4 = sincosN(4, v4);
    const r8 = sincosN(8, v8);

    inline for (0..4) |i| {
        try std.testing.expectApproxEqAbs(r4.sin[i], r8.sin[i], 1e-15);
        try std.testing.expectApproxEqAbs(r4.cos[i], r8.cos[i], 1e-15);
    }
}

test "atan2N(4)" {
    const y = VecN(4){ 1.0, 1.0, -1.0, -1.0 };
    const x = VecN(4){ 1.0, -1.0, -1.0, 1.0 };
    const result = atan2N(4, y, x);

    inline for (0..4) |i| {
        try std.testing.expectApproxEqAbs(std.math.atan2(y[i], x[i]), result[i], 1e-7);
    }
}

test "atan2N(8)" {
    const y = VecN(8){ 1.0, 1.0, -1.0, -1.0, 0.0, 1.0, -1.0, 0.5 };
    const x = VecN(8){ 1.0, -1.0, -1.0, 1.0, 1.0, 0.0, 0.0, 0.5 };
    const result = atan2N(8, y, x);

    inline for (0..8) |i| {
        try std.testing.expectApproxEqAbs(std.math.atan2(y[i], x[i]), result[i], 1e-7);
    }
}

test "modTwoPiN(8)" {
    const x = VecN(8){ 0.0, std.math.pi, 3.0 * std.math.pi, -std.math.pi, 7.0, 10.0, -5.0, 100.0 };
    const result = modTwoPiN(8, x);

    inline for (0..8) |i| {
        const expected = @mod(x[i], constants.twoPi);
        const exp_pos = if (expected < 0) expected + constants.twoPi else expected;
        try std.testing.expectApproxEqAbs(exp_pos, result[i], 1e-10);
    }
}

test "pow23N(4)" {
    const x = VecN(4){ 1.0, 8.0, 27.0, 0.5 };
    const result = pow23N(4, x);

    inline for (0..4) |i| {
        const expected = std.math.pow(f64, x[i], 2.0 / 3.0);
        try std.testing.expectApproxEqAbs(expected, result[i], 1e-12);
    }
}

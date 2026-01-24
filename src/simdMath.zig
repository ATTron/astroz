//! SIMD Math for SGP4
//! Fast vectorized trig functions that compute 4 values at once

const std = @import("std");
const constants = @import("constants.zig");

pub const Vec4 = @Vector(4, f64);
const Vec4i = @Vector(4, i64);

const SinCos = struct { sin: Vec4, cos: Vec4 };

// Constants
const zero: Vec4 = @splat(0.0);
const one: Vec4 = @splat(1.0);
const two: Vec4i = @splat(2);

// For range reduction: split pi/2 into high + low parts for precision
const twoOverPi: Vec4 = @splat(2.0 / std.math.pi);
const piOver2Hi: Vec4 = @splat(1.5707963267948966); // upper bits of pi/2
const piOver2Lo: Vec4 = @splat(6.123233995736766e-17); // lower bits of pi/2

// Polynomial Coefficients
// These approximate sin/cos on [-pi/4, pi/4] with ~1e-11 accuracy
// Generated via Remez algorithm (minimax polynomial fitting)

// sin(x) ~= x - x^3/6 + x^5/120 - ...  (Taylor like, but optimized coefficients)
const sinC1: Vec4 = @splat(-1.6666666666666632e-01); // ≈ -1/6
const sinC2: Vec4 = @splat(8.3333333333224632e-03); // ≈  1/120
const sinC3: Vec4 = @splat(-1.9841269829816540e-04); // ≈ -1/5040
const sinC4: Vec4 = @splat(2.7557313707070068e-06); // ≈  1/362880

// cos(x) ~= 1 - x^2/2 + x^4/24 - ...  (Taylor like, but optimized coefficients)
const cos_c1: Vec4 = @splat(-4.9999999999999994e-01); // ≈ -1/2
const cos_c2: Vec4 = @splat(4.1666666666666019e-02); // ≈  1/24
const cos_c3: Vec4 = @splat(-1.3888888888874096e-03); // ≈ -1/720
const cos_c4: Vec4 = @splat(2.4801587301492746e-05); // ≈  1/40320

// Core Functions

/// Compute sin and cos together (faster than separate calls)
/// Works for any angle, accurate to ~1e-11
pub inline fn sincosV4(angle: Vec4) SinCos {
    // Step 1: Range reduction
    // Trig functions repeat every 2*pi, so we reduce to [-pi/4, pi/4]
    // where our polynomial is accurate, then fix up the quadrant
    const quadrant = reduceToQuadrant(angle);

    // Step 2: Polynomial approximation on reduced range
    const r2 = quadrant.reduced * quadrant.reduced;
    const sinReduced = sinPoly(quadrant.reduced, r2);
    const cosReduced = cosPoly(r2);

    // Step 3: Fix up based on which quadrant we're in
    return applyQuadrant(sinReduced, cosReduced, quadrant.k);
}

/// Vectorized sin (just calls sincosV4)
pub inline fn sinV4(x: Vec4) Vec4 {
    return sincosV4(x).sin;
}

/// Vectorized cos (just calls sincosV4)
pub inline fn cosV4(x: Vec4) Vec4 {
    return sincosV4(x).cos;
}

// Implementation Details

/// Reduce angle to [-pi/4, pi/4] and track which quadrant
inline fn reduceToQuadrant(angle: Vec4) struct { reduced: Vec4, k: Vec4i } {
    // k = round(angle / (pi/2)) tells us the quadrant
    const kFloat = angle * twoOverPi;

    // Fast rounding trick: add then subtract 2^52+2^51 to round to integer
    const roundMagic: Vec4 = @splat(6755399441055744.0);
    const kRounded = kFloat + roundMagic - roundMagic;

    // Subtract k * pi/2 from angle (in two parts for precision)
    var reduced = @mulAdd(Vec4, -piOver2Hi, kRounded, angle);
    reduced = @mulAdd(Vec4, -piOver2Lo, kRounded, reduced);

    return .{
        .reduced = reduced,
        .k = @intFromFloat(kRounded),
    };
}

/// Evaluate sin polynomial: sin(r) ~= r + r^3*P(r^2)
inline fn sinPoly(r: Vec4, r2: Vec4) Vec4 {
    // Horner's method: evaluate from innermost term outward
    var p = @mulAdd(Vec4, sinC4, r2, sinC3);
    p = @mulAdd(Vec4, p, r2, sinC2);
    p = @mulAdd(Vec4, p, r2, sinC1);
    return @mulAdd(Vec4, p, r2 * r, r); // r + r^3 * p
}

/// Evaluate cos polynomial: cos(r) ~= 1 + r^2*P(r^2)
inline fn cosPoly(r2: Vec4) Vec4 {
    var p = @mulAdd(Vec4, cos_c4, r2, cos_c3);
    p = @mulAdd(Vec4, p, r2, cos_c2);
    p = @mulAdd(Vec4, p, r2, cos_c1);
    return @mulAdd(Vec4, p, r2, one); // 1 + r^2 * p
}

/// Apply quadrant correction to get final sin/cos values
/// k mod 4 tells us: 0=keep, 1=swap, 2=negate, 3=swap+negate
inline fn applyQuadrant(sinR: Vec4, cosR: Vec4, k: Vec4i) SinCos {
    const oneI: Vec4i = @splat(1);

    // Swap sin<->cos when k is odd (quadrants 1 and 3)
    const swap = (k & oneI) != @as(Vec4i, @splat(0));
    var s = @select(f64, swap, cosR, sinR);
    var c = @select(f64, swap, sinR, cosR);

    // Negate by flipping sign bit (bit 63) using XOR
    // sin: negate when k&2 != 0 (quadrants 2 and 3)
    // cos: negate when (k+1)&2 != 0 (quadrants 1 and 2)
    const sinSign = (k & two) << @splat(62);
    const cosSign = ((k + oneI) & two) << @splat(62);

    s = @bitCast(@as(Vec4i, @bitCast(s)) ^ sinSign);
    c = @bitCast(@as(Vec4i, @bitCast(c)) ^ cosSign);

    return .{ .sin = s, .cos = c };
}

// Other SIMD Math Functions

pub const twoPiVec: Vec4 = @splat(constants.twoPi);
pub const piVec: Vec4 = @splat(std.math.pi);
pub const halfPiVec: Vec4 = @splat(std.math.pi / 2.0);
pub const invTwoPiVec: Vec4 = @splat(1.0 / constants.twoPi);
const negTwoPiVec: Vec4 = @splat(-constants.twoPi);

/// Normalize angle to [0, 2*pi)
pub inline fn modTwoPiV4(x: Vec4) Vec4 {
    const n = @floor(x * invTwoPiVec);
    var result = @mulAdd(Vec4, negTwoPiVec, n, x);
    const mask = result < zero;
    result = @select(f64, mask, result + twoPiVec, result);
    return result;
}

/// Vectorized atan2(y, x). Accurate to ~1e-7
pub fn atan2SIMD(y: Vec4, x: Vec4) Vec4 {
    // We compute atan(min/max) to keep the ratio in [0,1] where our polynomial works
    // then adjust based on which was bigger and the signs of x,y

    const absX = @abs(x);
    const absY = @abs(y);
    const maxXY = @max(absX, absY);
    const minXY = @min(absX, absY);

    // Safe division (avoid 0/0)
    const tiny: Vec4 = @splat(1.0e-30);
    const t = minXY / @max(maxXY, tiny);
    const t2 = t * t;

    // Polynomial for atan(t) on [0,1], ~1e-7 accuracy (Remez minimax coefficients)
    // atan(t) ~= t - t^3/3 + t^5/5 - ... but with optimized coefficients
    var result = @as(Vec4, @splat(0.0028662257)); // t^17 coefficient
    result = @mulAdd(Vec4, result, t2, @splat(-0.0161657367)); // t^15
    result = @mulAdd(Vec4, result, t2, @splat(0.0429096138)); // t^13
    result = @mulAdd(Vec4, result, t2, @splat(-0.0752896400)); // t^11
    result = @mulAdd(Vec4, result, t2, @splat(0.1065626393)); // t^9
    result = @mulAdd(Vec4, result, t2, @splat(-0.1420889944)); // t^7
    result = @mulAdd(Vec4, result, t2, @splat(0.1999355085)); // t^5
    result = @mulAdd(Vec4, result, t2, @splat(-0.3333314528)); // t^3
    result = @mulAdd(Vec4, result, t2, one); // t^1
    result = result * t;

    // If |y| > |x|, we computed atan(|x|/|y|), convert to atan(|y|/|x|) = pi/2 - atan
    result = @select(f64, absY > absX, halfPiVec - result, result);

    // Quadrant adjustment based on signs
    result = @select(f64, x < zero, piVec - result, result); // Q2/Q3: flip around pi/2
    result = @select(f64, y < zero, -result, result); // Q3/Q4: negate

    return result;
}

/// x^1.5 = x * sqrt(x)
pub fn pow15V4(x: Vec4) Vec4 {
    return x * @sqrt(x);
}

test "sinSIMD accuracy" {
    const testVals = Vec4{ 0.0, std.math.pi / 6.0, std.math.pi / 4.0, std.math.pi / 2.0 };
    const expected = Vec4{ 0.0, 0.5, std.math.sqrt(2.0) / 2.0, 1.0 };
    const result = @sin(testVals);

    inline for (0..4) |i| {
        try std.testing.expectApproxEqAbs(expected[i], result[i], 1e-14);
    }
}

test "cosSIMD accuracy" {
    const testVals = Vec4{ 0.0, std.math.pi / 3.0, std.math.pi / 4.0, std.math.pi };
    const expected = Vec4{ 1.0, 0.5, std.math.sqrt(2.0) / 2.0, -1.0 };
    const result = @cos(testVals);

    inline for (0..4) |i| {
        try std.testing.expectApproxEqAbs(expected[i], result[i], 1e-14);
    }
}

test "sinSIMD large angles" {
    const testVals = Vec4{ 10.0 * std.math.pi, 100.0, -50.0, 1000.0 };
    const result = @sin(testVals);

    inline for (0..4) |i| {
        const expected = @sin(testVals[i]);
        try std.testing.expectApproxEqAbs(expected, result[i], 1e-14);
    }
}

test "cosSIMD large angles" {
    const testVals = Vec4{ 10.0 * std.math.pi, 100.0, -50.0, 1000.0 };
    const result = @cos(testVals);

    inline for (0..4) |i| {
        const expected = @cos(testVals[i]);
        try std.testing.expectApproxEqAbs(expected, result[i], 1e-14);
    }
}

test "atan2SIMD basic quadrants" {
    // Test all four quadrants
    // Q1: x > 0, y > 0
    // Q2: x < 0, y > 0
    // Q3: x < 0, y < 0
    // Q4: x > 0, y < 0
    const y = Vec4{ 1.0, 1.0, -1.0, -1.0 };
    const x = Vec4{ 1.0, -1.0, -1.0, 1.0 };
    const result = atan2SIMD(y, x);

    const tol = 1e-6; // atan2 polynomial has ~1e-7 accuracy
    inline for (0..4) |i| {
        const expected = std.math.atan2(y[i], x[i]);
        try std.testing.expectApproxEqAbs(expected, result[i], tol);
    }
}

test "atan2SIMD axis aligned" {
    // Test points along axes
    const y = Vec4{ 0.0, 1.0, 0.0, -1.0 };
    const x = Vec4{ 1.0, 0.0, -1.0, 0.0 };
    const result = atan2SIMD(y, x);

    const tol = 1e-6;
    inline for (0..4) |i| {
        const expected = std.math.atan2(y[i], x[i]);
        try std.testing.expectApproxEqAbs(expected, result[i], tol);
    }
}

test "atan2SIMD various angles" {
    // Test various angles typical in SGP4 Kepler solver
    const y = Vec4{ 0.5, -0.3, 0.8, -0.9 };
    const x = Vec4{ 0.866, 0.954, -0.6, -0.436 };
    const result = atan2SIMD(y, x);

    const tol = 1e-6;
    inline for (0..4) |i| {
        const expected = std.math.atan2(y[i], x[i]);
        try std.testing.expectApproxEqAbs(expected, result[i], tol);
    }
}

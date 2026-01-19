//! SIMD Math for SGP4
//!
//! Uses Zig's built-in @sin/@cos which auto-vectorize on modern CPUs.
//! The builtin functions map directly to LLVM intrinsics which use
//! platform-optimal implementations (e.g., libmvec on Linux x86-64).

const std = @import("std");
const constants = @import("constants.zig");

const simdMath = @This();

// compute 4 time segments per CPU tick
pub const Vec4 = @Vector(4, f64);

pub const twoPiVec: Vec4 = @splat(constants.twoPi);
pub const piVec: Vec4 = @splat(std.math.pi);
pub const halfPiVec: Vec4 = @splat(std.math.pi / 2.0);
pub const invTwoPiVec: Vec4 = @splat(1.0 / constants.twoPi);

/// Vectorized sine - uses Zig's @sin builtin which auto-vectorizes
pub fn sinSIMD(x: Vec4) Vec4 {
    return @sin(x);
}

/// Vectorized cosine - uses Zig's @cos builtin which auto-vectorizes
pub fn cosSIMD(x: Vec4) Vec4 {
    return @cos(x);
}

/// Vectorized atan2 using polynomial approximation for the principal value,
/// with quadrant correction. Accurate to ~1e-7 radians, sufficient for SGP4.
///
/// This avoids scalar fallback and maintains SIMD flow through the Kepler solver.
pub fn atan2SIMD(y: Vec4, x: Vec4) Vec4 {
    const abs_x = @abs(x);
    const abs_y = @abs(y);

    // Compute atan(min/max) to keep argument in [0, 1] for better polynomial accuracy
    const max_xy = @max(abs_x, abs_y);
    const min_xy = @min(abs_x, abs_y);

    // Avoid division by zero - if both are zero, result is undefined anyway
    const tiny: Vec4 = @splat(1.0e-30);
    const safe_max = @max(max_xy, tiny);
    const t = min_xy / safe_max;
    const t2 = t * t;

    // Polynomial approximation for atan(t) where t in [0, 1]
    // atan(t) ≈ t - t³/3 + t⁵/5 - t⁷/7 + ...
    // Using optimized coefficients for better accuracy in this range
    const c1: Vec4 = @splat(1.0);
    const c3: Vec4 = @splat(-0.3333314528);
    const c5: Vec4 = @splat(0.1999355085);
    const c7: Vec4 = @splat(-0.1420889944);
    const c9: Vec4 = @splat(0.1065626393);
    const c11: Vec4 = @splat(-0.0752896400);
    const c13: Vec4 = @splat(0.0429096138);
    const c15: Vec4 = @splat(-0.0161657367);
    const c17: Vec4 = @splat(0.0028662257);

    // Horner's method for polynomial evaluation
    var atan_t = c17;
    atan_t = atan_t * t2 + c15;
    atan_t = atan_t * t2 + c13;
    atan_t = atan_t * t2 + c11;
    atan_t = atan_t * t2 + c9;
    atan_t = atan_t * t2 + c7;
    atan_t = atan_t * t2 + c5;
    atan_t = atan_t * t2 + c3;
    atan_t = atan_t * t2 + c1;
    atan_t = atan_t * t;

    // If |y| > |x|, we computed atan(|x|/|y|), need atan(|y|/|x|) = π/2 - atan(|x|/|y|)
    const swap_mask = abs_y > abs_x;
    atan_t = @select(f64, swap_mask, halfPiVec - atan_t, atan_t);

    // Quadrant correction based on signs of x and y
    // Q1: x > 0, y >= 0: result = atan_t
    // Q2: x < 0, y >= 0: result = π - atan_t
    // Q3: x < 0, y < 0: result = -π + atan_t
    // Q4: x > 0, y < 0: result = -atan_t

    const zero: Vec4 = @splat(0.0);
    const x_neg = x < zero;
    const y_neg = y < zero;

    // Start with the absolute angle
    var result = atan_t;

    // If x < 0, flip angle around π/2 (i.e., result = π - result)
    result = @select(f64, x_neg, piVec - result, result);

    // If y < 0, negate the result
    result = @select(f64, y_neg, -result, result);

    return result;
}

pub fn pow15V4(x: Vec4) Vec4 {
    return x * @sqrt(x);
}

// =============================================================================
// Tests
// =============================================================================

test "sinSIMD accuracy" {
    const testVals = Vec4{ 0.0, std.math.pi / 6.0, std.math.pi / 4.0, std.math.pi / 2.0 };
    const expected = Vec4{ 0.0, 0.5, std.math.sqrt(2.0) / 2.0, 1.0 };
    const result = sinSIMD(testVals);

    inline for (0..4) |i| {
        try std.testing.expectApproxEqAbs(expected[i], result[i], 1e-14);
    }
}

test "cosSIMD accuracy" {
    const testVals = Vec4{ 0.0, std.math.pi / 3.0, std.math.pi / 4.0, std.math.pi };
    const expected = Vec4{ 1.0, 0.5, std.math.sqrt(2.0) / 2.0, -1.0 };
    const result = cosSIMD(testVals);

    inline for (0..4) |i| {
        try std.testing.expectApproxEqAbs(expected[i], result[i], 1e-14);
    }
}

test "sinSIMD large angles" {
    const testVals = Vec4{ 10.0 * std.math.pi, 100.0, -50.0, 1000.0 };
    const result = sinSIMD(testVals);

    inline for (0..4) |i| {
        const expected = @sin(testVals[i]);
        try std.testing.expectApproxEqAbs(expected, result[i], 1e-14);
    }
}

test "cosSIMD large angles" {
    const testVals = Vec4{ 10.0 * std.math.pi, 100.0, -50.0, 1000.0 };
    const result = cosSIMD(testVals);

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

test "atan2SIMD axis-aligned" {
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

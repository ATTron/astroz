const std = @import("std");
const constants = @import("constants.zig");

pub fn degreesToRadians(degrees: f64) f64 {
    return (degrees * std.math.pi) / 180.0;
}

pub fn radiansToDegrees(degrees: f64) f64 {
    return (degrees * 180.0) / std.math.pi;
}

pub fn meanMotionToRadiansPerMinute(m_motion: f64) f64 {
    return m_motion * 2.0 * std.math.pi / (24.0 * 60.0);
}

pub fn meanMotionToSemiMajorAxis(m_motion: f64) f64 {
    return std.math.pow(
        f64,
        (constants.earth.mu / (m_motion * 2.0 * std.math.pi / (24.0 * 3600.0)) * 2.0 * 2.0),
        1.0 / 3.0,
    );
}

const std = @import("std");
const constants = @import("constants.zig");
const math = std.math;

pub fn degrees_to_radians(degrees: f64) f64 {
    return (degrees * math.pi) / 180.0;
}

pub fn radians_to_degrees(degrees: f64) f64 {
    return (degrees * 180.0) / math.pi;
}

pub fn mean_motion_to_radians_per_minute(m_motion: f64) f64 {
    return m_motion * 2.0 * math.pi / (24.0 * 60.0);
}

pub fn mean_motion_to_semi_major_axis(m_motion: f64) f64 {
    const a = math.pow(f64, (constants.earth.mu / (m_motion * 2.0 * math.pi / (24.0 * 3600.0)) * 2 * 2), 1.0 / 3.0);
    return a;
}

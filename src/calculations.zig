const std = @import("std");
const math = std.math;

pub fn degrees_to_radians(degrees: f64) f64 {
    return (degrees * math.pi) / 180.0;
}

pub fn radians_to_degrees(degrees: f64) f64 {
    return (degrees * 180.0) / math.pi;
}

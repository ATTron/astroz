//! Orbital mechanics C API exports

const std = @import("std");
const astroz = @import("astroz");
const constants = astroz.constants;

const err = @import("error.zig");

/// hohmann transfer result (C-compatible)
pub const HohmannResult = extern struct {
    semi_major_axis: f64,
    delta_v1: f64,
    delta_v2: f64,
    total_delta_v: f64,
    transfer_time: f64,
    transfer_time_days: f64,
};

/// mu: gravitational parameter (km³/s²)
/// r1: initial orbit radius (km)
/// r2: final orbit radius (km)
pub fn hohmann(mu: f64, r1: f64, r2: f64, out: *HohmannResult) err.Code {
    if (r1 <= 0 or r2 <= 0) return .value_error;
    if (@abs(r1 - r2) < 1000) return .value_error;

    const sma = (r1 + r2) / 2.0;
    const v1_circ = @sqrt(mu / r1);
    const v2_circ = @sqrt(mu / r2);
    const v1_trans = @sqrt(mu / r1) * @sqrt(2.0 * r2 / (r1 + r2));
    const v2_trans = @sqrt(mu / r2) * @sqrt(2.0 * r1 / (r1 + r2));
    const dv1 = v1_trans - v1_circ;
    const dv2 = v2_circ - v2_trans;
    const transfer_time = std.math.pi * @sqrt(sma * sma * sma / mu);

    out.* = .{
        .semi_major_axis = sma,
        .delta_v1 = dv1,
        .delta_v2 = dv2,
        .total_delta_v = @abs(dv1) + @abs(dv2),
        .transfer_time = transfer_time,
        .transfer_time_days = transfer_time / constants.seconds_per_day,
    };
    return .ok;
}

/// mu: gravitational parameter (km³/s²)
/// radius: distance from central body (km)
/// sma: semi-major axis (km), or 0 for circular orbit
pub fn orbitalVelocity(mu: f64, radius: f64, sma: f64) f64 {
    if (radius <= 0) return -1.0;
    if (sma == 0) return @sqrt(mu / radius);
    if (sma <= 0) return -1.0;
    return @sqrt(mu * (2.0 / radius - 1.0 / sma));
}

/// mu: gravitational parameter (km³/s²)
/// sma: semi-major axis (km)
pub fn orbitalPeriod(mu: f64, sma: f64) f64 {
    if (sma <= 0) return -1.0;
    return 2.0 * std.math.pi * @sqrt(sma * sma * sma / mu);
}

/// mu: gravitational parameter (km³/s²)
/// radius: distance from central body (km)
pub fn escapeVelocity(mu: f64, radius: f64) f64 {
    if (radius <= 0) return -1.0;
    return @sqrt(2.0 * mu / radius);
}

//! Orbital mechanics C API exports

const astroz = @import("astroz");
const calc = astroz.calculations;
const constants = astroz.constants;

const err = @import("error.zig");

pub const HohmannResult = extern struct {
    semi_major_axis: f64,
    delta_v1: f64,
    delta_v2: f64,
    total_delta_v: f64,
    transfer_time: f64,
    transfer_time_days: f64,
};

pub fn hohmann(mu: f64, r1: f64, r2: f64, out: *HohmannResult) err.Code {
    if (r1 <= 0 or r2 <= 0) return .valueError;
    if (@abs(r1 - r2) < 1000) return .valueError;

    const result = calc.hohmannTransfer(mu, r1, r2);
    out.* = .{
        .semi_major_axis = result.semiMajorAxis,
        .delta_v1 = result.deltaV1,
        .delta_v2 = result.deltaV2,
        .total_delta_v = result.totalDeltaV,
        .transfer_time = result.transferTime,
        .transfer_time_days = result.transferTime / constants.secondsPerDay,
    };
    return .ok;
}

pub fn orbitalVelocity(mu: f64, radius: f64, sma: f64) f64 {
    if (radius <= 0) return -1.0;
    if (sma == 0) return calc.orbitalVelocity(mu, radius, null);
    if (sma <= 0) return -1.0;
    return calc.orbitalVelocity(mu, radius, sma);
}

pub fn orbitalPeriod(mu: f64, sma: f64) f64 {
    if (sma <= 0) return -1.0;
    return calc.orbitalPeriod(mu, sma);
}

pub fn escapeVelocity(mu: f64, radius: f64) f64 {
    if (radius <= 0) return -1.0;
    return calc.escapeVelocity(mu, radius);
}

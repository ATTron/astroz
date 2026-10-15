//! Orbital mechanics helpers bound to a gravitational parameter, with input validation.

const std = @import("std");

const calculations = @import("calculations.zig");
const Vector3D = calculations.Vector3D;
const constants = @import("constants.zig");

pub const ValidationError = error{ InvalidDeltaV, InvalidDistance, InvalidSMA, ValueError };

const OrbitalMechanics = @This();

pub const LambertResult = struct {
    departureVelocity: Vector3D,
    arrivalVelocity: Vector3D,
    transferAngle: f64,
    semiMajorAxis: f64,
    timeOfFlight: f64,
};

pub const TransferResult = struct {
    semiMajorAxis: f64,
    deltaV1: f64,
    deltaV2: f64,
    totalDeltaV: f64,
    transferTime: f64,
    transferTimeDays: f64,
};

pub const BiEllipticTransferResult = struct {
    semiMajorAxis: f64,
    deltaV1: f64,
    deltaV2: f64,
    deltaV3: f64,
    totalDeltaV: f64,
    totalTime: f64,
    totalTimeDays: f64,
};

mu: f64,

pub fn init(mu: f64) OrbitalMechanics {
    return .{ .mu = mu };
}

pub fn orbitalVelocity(self: *OrbitalMechanics, radius: f64, semiMajorAxis: ?f64) ValidationError!f64 {
    if (radius <= 0) return ValidationError.InvalidDistance;
    if (semiMajorAxis) |sma| if (sma <= 0) return ValidationError.InvalidSMA;
    return calculations.orbitalVelocity(self.mu, radius, semiMajorAxis);
}

pub fn orbitalPeriod(self: *OrbitalMechanics, semiMajorAxis: f64) ValidationError!f64 {
    if (semiMajorAxis <= 0) return ValidationError.InvalidSMA;
    return calculations.orbitalPeriod(self.mu, semiMajorAxis);
}

pub fn escapeVelocity(self: *OrbitalMechanics, radius: f64) ValidationError!f64 {
    if (radius <= 0) return ValidationError.InvalidDistance;
    return calculations.escapeVelocity(self.mu, radius);
}

/// Hohmann transfer between two circular orbits.
pub fn hohmannTransfer(self: *OrbitalMechanics, initialRadius: f64, finalRadius: f64) ValidationError!TransferResult {
    if (initialRadius <= 0 or finalRadius <= 0) return ValidationError.ValueError;
    if (@abs(initialRadius - finalRadius) < 1000) return ValidationError.ValueError;

    const r = calculations.hohmannTransfer(self.mu, initialRadius, finalRadius);
    return .{
        .semiMajorAxis = r.semiMajorAxis,
        .deltaV1 = r.deltaV1,
        .deltaV2 = r.deltaV2,
        .totalDeltaV = r.totalDeltaV,
        .transferTime = r.transferTime,
        .transferTimeDays = r.transferTime / constants.secondsPerDay,
    };
}

pub fn biEllipicTransfer(self: *OrbitalMechanics, initialRadius: f64, finalRadius: f64, aphelionRadius: f64) !BiEllipticTransferResult {
    if (initialRadius <= 0 or finalRadius <= 0 or aphelionRadius <= 0) {
        return ValidationError.ValueError;
    }

    if (aphelionRadius <= @max(initialRadius, finalRadius)) {
        return ValidationError.ValueError;
    }

    const firstTransferSemiMajorAxis = (initialRadius + aphelionRadius) / 2;
    const velocityInitialCircular = try self.orbitalVelocity(initialRadius, null);
    const velocityInitialTransfer = try self.orbitalVelocity(initialRadius, firstTransferSemiMajorAxis);
    const deltaVFirstBurn = velocityInitialTransfer - velocityInitialCircular;

    const secondTransferSemiMajorAxis = (finalRadius + aphelionRadius) / 2;
    const velocityAphelionFirstTransfer = try self.orbitalVelocity(aphelionRadius, firstTransferSemiMajorAxis);
    const velocityAphelionSecondTransfer = try self.orbitalVelocity(aphelionRadius, secondTransferSemiMajorAxis);
    const deltaVSecondBurn = velocityAphelionSecondTransfer - velocityAphelionFirstTransfer;

    const velocityFinalCircular = try self.orbitalVelocity(finalRadius, null);
    const velocityFinalTransfer = try self.orbitalVelocity(finalRadius, secondTransferSemiMajorAxis);
    const deltaVThirdBurn = velocityFinalCircular - velocityFinalTransfer;

    const totalDeltaV = @abs(deltaVFirstBurn) + @abs(deltaVSecondBurn) + @abs(deltaVThirdBurn);

    const firstTransferTime = std.math.pi * @sqrt(std.math.pow(f64, firstTransferSemiMajorAxis, 3) / self.mu);
    const secondTransferTime = std.math.pi * @sqrt(std.math.pow(f64, secondTransferSemiMajorAxis, 3) / self.mu);
    const totalTime = firstTransferTime + secondTransferTime;

    return .{
        .semiMajorAxis = firstTransferSemiMajorAxis,
        .deltaV1 = deltaVFirstBurn,
        .deltaV2 = deltaVSecondBurn,
        .deltaV3 = deltaVThirdBurn,
        .totalDeltaV = totalDeltaV,
        .totalTime = totalTime,
        .totalTimeDays = totalTime / constants.secondsPerDay,
    };
}

/// Solves Lambert's problem: the orbit from `r1` to `r2` in `timeOfFlight`
/// seconds, taking the prograde (short or long way, by the sense of r1 x r2
/// about +z) zero-revolution transfer. Universal variables with Stumpff
/// functions (Curtis, Orbital Mechanics for Engineering Students, Alg. 5.2).
pub fn lambert(self: OrbitalMechanics, r1Vec: Vector3D, r2Vec: Vector3D, timeOfFlight: f64) !LambertResult {
    const r1 = r1Vec.magnitude();
    const r2 = r2Vec.magnitude();
    if (timeOfFlight <= 0 or r1 <= 0 or r2 <= 0) return ValidationError.ValueError;

    const cosTheta = std.math.clamp(r1Vec.dot(r2Vec) / (r1 * r2), -1.0, 1.0);
    var theta = std.math.acos(cosTheta);
    if (r1Vec.x() * r2Vec.y() - r1Vec.y() * r2Vec.x() < 0) theta = 2.0 * std.math.pi - theta;
    // A is zero for a 180-degree transfer, whose plane r1 and r2 don't define
    const a = @sin(theta) * @sqrt(r1 * r2 / (1.0 - @cos(theta)));
    if (!std.math.isFinite(a) or @abs(a) < 1e-9) return ValidationError.ValueError;

    const geo: LambertGeometry = .{ .r1 = r1, .r2 = r2, .a = a, .sqrtMu = @sqrt(self.mu) };

    // Time of flight grows monotonically with z, so bracket the root and bisect.
    // z < 0 is hyperbolic, 0 parabolic, and z < (2 pi)^2 keeps to zero revolutions;
    // stop just short of that, where C(z) runs out of precision.
    var lo: f64 = -4.0 * std.math.pi * std.math.pi;
    var hi: f64 = 4.0 * std.math.pi * std.math.pi - 1e-4;
    // Lower lo until the time there is short enough; null means y < 0, already too low
    while (geo.time(lo)) |t| : (lo *= 2.0) {
        if (t < timeOfFlight) break;
        if (lo < -1e6) return ValidationError.ValueError;
    }
    const tHi = geo.time(hi) orelse return ValidationError.ValueError;
    if (tHi < timeOfFlight) return ValidationError.ValueError; // would need more than one revolution
    var z: f64 = 0.5 * (lo + hi);
    for (0..200) |_| {
        z = 0.5 * (lo + hi);
        const t = geo.time(z) orelse {
            lo = z;
            continue;
        };
        if (t < timeOfFlight) lo = z else hi = z;
        if (hi - lo < 1e-12 * @max(1.0, @abs(z))) break;
    }

    const yz = geo.y(z);
    const f = 1.0 - yz / r1;
    const g = a * @sqrt(yz / self.mu);
    const gDot = 1.0 - yz / r2;
    const v1 = r2Vec.sub(r1Vec.mul(f)).mul(1.0 / g);
    const v2 = r2Vec.mul(gDot).sub(r1Vec).mul(1.0 / g);
    const v1Sq = v1.dot(v1);

    return .{
        .departureVelocity = v1,
        .arrivalVelocity = v2,
        .transferAngle = theta,
        .semiMajorAxis = 1.0 / (2.0 / r1 - v1Sq / self.mu),
        .timeOfFlight = timeOfFlight,
    };
}

const LambertGeometry = struct {
    r1: f64,
    r2: f64,
    a: f64,
    sqrtMu: f64,

    fn y(g: LambertGeometry, z: f64) f64 {
        return g.r1 + g.r2 + g.a * (z * stumpffS(z) - 1.0) / @sqrt(stumpffC(z));
    }

    /// Time of flight for universal-anomaly parameter z, or null where y < 0.
    fn time(g: LambertGeometry, z: f64) ?f64 {
        const yz = g.y(z);
        if (yz < 0) return null;
        return (std.math.pow(f64, yz / stumpffC(z), 1.5) * stumpffS(z) + g.a * @sqrt(yz)) / g.sqrtMu;
    }
};

fn stumpffS(z: f64) f64 {
    if (z > 1e-6) {
        const s = @sqrt(z);
        return (s - @sin(s)) / (s * s * s);
    }
    if (z < -1e-6) {
        const s = @sqrt(-z);
        return (std.math.sinh(s) - s) / (s * s * s);
    }
    return 1.0 / 6.0 - z / 120.0;
}

fn stumpffC(z: f64) f64 {
    if (z > 1e-6) return (1.0 - @cos(@sqrt(z))) / z;
    if (z < -1e-6) return (std.math.cosh(@sqrt(-z)) - 1.0) / -z;
    return 0.5 - z / 24.0;
}

test "lambert matches Curtis example 5.2" {
    const om = OrbitalMechanics.init(398600.0);
    const res = try om.lambert(Vector3D.new(5000, 10000, 2100), Vector3D.new(-14600, 2500, 7000), 3600);
    const v1 = res.departureVelocity;
    const v2 = res.arrivalVelocity;
    try std.testing.expectApproxEqAbs(-5.9925, v1.x(), 1e-3);
    try std.testing.expectApproxEqAbs(1.9254, v1.y(), 1e-3);
    try std.testing.expectApproxEqAbs(3.2456, v1.z(), 1e-3);
    try std.testing.expectApproxEqAbs(-3.3125, v2.x(), 1e-3);
    try std.testing.expectApproxEqAbs(-4.1966, v2.y(), 1e-3);
    try std.testing.expectApproxEqAbs(-0.38529, v2.z(), 1e-3);
    try std.testing.expectApproxEqRel(20003.0, res.semiMajorAxis, 1e-3);
}

test "lambert departure state reaches r2 after the time of flight" {
    const mu = constants.sun.mu;
    const om = OrbitalMechanics.init(mu);
    const au = constants.auKm;
    for ([_][2]f64{ .{ 135.0, 250.0 }, .{ 60.0, 120.0 }, .{ 250.0, 300.0 } }) |case| {
        const th = case[0] * constants.deg2rad;
        const r1 = Vector3D.new(au, 0, 0);
        const r2 = Vector3D.new(1.524 * au * @cos(th), 1.524 * au * @sin(th), 0.01 * au);
        const tof = case[1] * constants.secondsPerDay;
        const res = try om.lambert(r1, r2, tof);
        // RK4 on two-body gravity, 1000 steps
        var state = [6]f64{ r1.x(), r1.y(), r1.z(), res.departureVelocity.x(), res.departureVelocity.y(), res.departureVelocity.z() };
        const n = 20000;
        const h = tof / n;
        for (0..n) |_| state = rk4TwoBody(state, h, mu);
        const miss = @sqrt((state[0] - r2.x()) * (state[0] - r2.x()) + (state[1] - r2.y()) * (state[1] - r2.y()) + (state[2] - r2.z()) * (state[2] - r2.z()));
        try std.testing.expect(miss < 10.0); // km after hundreds of days
    }
}

fn rk4TwoBody(s: [6]f64, h: f64, mu: f64) [6]f64 {
    const deriv = struct {
        fn d(x: [6]f64, m: f64) [6]f64 {
            const r3 = std.math.pow(f64, x[0] * x[0] + x[1] * x[1] + x[2] * x[2], 1.5);
            return .{ x[3], x[4], x[5], -m * x[0] / r3, -m * x[1] / r3, -m * x[2] / r3 };
        }
    }.d;
    const k1 = deriv(s, mu);
    var t: [6]f64 = undefined;
    for (&t, s, k1) |*o, a, k| o.* = a + 0.5 * h * k;
    const k2 = deriv(t, mu);
    for (&t, s, k2) |*o, a, k| o.* = a + 0.5 * h * k;
    const k3 = deriv(t, mu);
    for (&t, s, k3) |*o, a, k| o.* = a + h * k;
    const k4 = deriv(t, mu);
    var out: [6]f64 = undefined;
    for (&out, s, k1, k2, k3, k4) |*o, a, b, c, d, e| o.* = a + h / 6.0 * (b + 2 * c + 2 * d + e);
    return out;
}

//! Base struct that provides various orbital mechanic functions

const std = @import("std");
const log = std.log;

const calculations = @import("calculations.zig");
const Vector3D = calculations.Vector3D;
const constants = @import("constants.zig");
const CelestialBody = constants.CelestialBody;

pub const ValidationError = error{ InvalidDeltaV, ValueError };

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

    pub fn init(semiMajorAxis: f64, deltaV1: f64, deltaV2: f64, totalDeltaV: f64, transferTime: f64, transferTimeDays: f64) !TransferResult {
        const expectedTotal = @abs(deltaV1) + @abs(deltaV2);
        if (@abs(totalDeltaV - expectedTotal) > 1e-6) {
            return ValidationError.InvalidDeltaV;
        }

        return .{
            .semiMajorAxis = semiMajorAxis,
            .deltaV1 = deltaV1,
            .deltaV2 = deltaV2,
            .totalDeltaV = totalDeltaV,
            .transferTime = transferTime,
            .transferTimeDays = transferTimeDays,
        };
    }
};

pub const BiEllipticTransferResult = struct {
    semiMajorAxis: f64,
    deltaV1: f64,
    deltaV2: f64,
    deltaV3: f64,
    totalDeltaV: f64,
    totalTime: f64,
    totalTimeDays: f64,

    pub fn init(semiMajorAxis: f64, deltaV1: f64, deltaV2: f64, deltaV3: f64, totalDeltaV: f64, totalTime: f64, totalTimeDays: f64) !BiEllipticTransferResult {
        return .{
            .semiMajorAxis = semiMajorAxis,
            .deltaV1 = deltaV1,
            .deltaV2 = deltaV2,
            .deltaV3 = deltaV3,
            .totalDeltaV = totalDeltaV,
            .totalTime = totalTime,
            .totalTimeDays = totalTimeDays,
        };
    }
};

centralBody: CelestialBody,
mu: f64,

pub fn init(mu: f64, centralBody: CelestialBody) OrbitalMechanics {
    return .{
        .centralBody = centralBody,
        .mu = mu,
    };
}

/// calculate orbital velocity
pub fn orbitalVelocity(self: *OrbitalMechanics, r: f64, a: ?f64) f64 {
    if (r <= 0) {
        log.warn("Distance cannot be negative: {d}", .{r});
        return 0;
    }

    if (a == null) {
        return @sqrt(self.mu / r);
    } else {
        const sma = a.?;
        if (sma <= 0) {
            log.warn("Semi major axis cannot be negative: {d}", .{sma});
            return 0;
        }
        return @sqrt(self.mu * (2.0 / r - 1.0 / sma));
    }
}

/// calculate orbital period
pub fn orbitalPeriod(self: *OrbitalMechanics, a: f64) f64 {
    if (a <= 0) {
        log.warn("Semi major axis cannot be negative: {d}", .{a});
        return 0;
    }

    return 2.0 * std.math.pi * @sqrt((a * a * a) / self.mu);
}

/// calculate escape velocity from a given radius
pub fn escapeVelocity(self: *OrbitalMechanics, r: f64) f64 {
    if (r <= 0) {
        log.warn("Distance cannot be negative: {d}", .{r});
    }
    return @sqrt(2 * self.mu / r);
}

/// calculate hohmann transfer parameters between two circular orbits
/// note: this calculates heliocentric transfer delta-v
/// ex. earth - mars
/// heliocentric delta-v: ~5.6 km/s (this calculation)
/// LEO to mars transfer: ~3.6 km/s (includes earth's orbital velocity advantage)
pub fn hohmannTransfer(self: *OrbitalMechanics, r1: f64, r2: f64) !TransferResult {
    if (r1 <= 0 or r2 <= 0) {
        return ValidationError.ValueError;
    }
    if (@abs(r1 - r2) < 1000) {
        return ValidationError.ValueError;
    }

    const aTransfer = (r1 + r2) / 2.0;

    const v1Circular = @sqrt(self.mu / r1);
    const v2Circular = @sqrt(self.mu / r2);

    const v1Transfer = @sqrt(self.mu / r1) * @sqrt(2.0 * r2 / (r1 + r2));
    const v2Transfer = @sqrt(self.mu / r2) * @sqrt(2.0 * r1 / (r1 + r2));

    const deltaV1 = v1Transfer - v1Circular;
    const deltaV2 = v2Circular - v2Transfer;

    const totalDeltaV = @abs(deltaV1) + @abs(deltaV2);

    const transferTime = std.math.pi * @sqrt((aTransfer * aTransfer * aTransfer) / self.mu);

    return TransferResult.init(
        aTransfer,
        deltaV1,
        deltaV2,
        totalDeltaV,
        transferTime,
        transferTime / (24.0 * 3600.0),
    );
}

pub fn biEllipicTransfer(self: *OrbitalMechanics, r1: f64, r2: f64, rAphelion: f64) !BiEllipticTransferResult {
    if (r1 <= 0 or r2 <= 0 or rAphelion <= 0) {
        return ValidationError.ValueError;
    }

    if (rAphelion <= @max(r1, r2)) {
        return ValidationError.ValueError;
    }

    const a1 = (r1 + rAphelion) / 2;
    const v1Circular = self.orbitalVelocity(r1, null);
    const v1Transfer = self.orbitalVelocity(r1, a1);
    const deltaV1 = v1Transfer - v1Circular;

    const a2 = (r2 + rAphelion) / 2;
    const vAphelion1 = self.orbitalVelocity(rAphelion, a1);
    const vAphelion2 = self.orbitalVelocity(rAphelion, a2);
    const deltaV2 = vAphelion2 - vAphelion1;

    const v2Circular = self.orbitalVelocity(r2, null);
    const v2Transfer = self.orbitalVelocity(r2, a2);
    const deltaV3 = v2Circular - v2Transfer;

    const totalDeltaV = @abs(deltaV1) + @abs(deltaV2) + @abs(deltaV3);

    const t1 = std.math.pi * @sqrt(std.math.pow(f64, a1, 3) / self.mu);
    const t2 = std.math.pi * @sqrt(std.math.pow(f64, a2, 3) / self.mu);
    const totalTime = t1 + t2;

    return .{
        .semiMajorAxis = a1,
        .deltaV1 = deltaV1,
        .deltaV2 = deltaV2,
        .deltaV3 = deltaV3,
        .totalDeltaV = totalDeltaV,
        .totalTime = totalTime,
        .totalTimeDays = totalTime / (24 * 3600),
    };
}

/// simplified lambert solver using universal variable approach
pub fn lambertSolverSimple(self: OrbitalMechanics, r1Vec: Vector3D, r2Vec: Vector3D, tof: f64) !LambertResult {
    if (tof <= 0) {
        return ValidationError.ValueError;
    }

    const r1 = r1Vec.magnitude();
    const r2 = r2Vec.magnitude();

    if (r1 <= 0 or r2 <= 0) {
        return ValidationError.ValueError;
    }

    var cosDnu = r1Vec.dot(r2Vec) / (r1 * r2);
    cosDnu = @max(-1.0, @min(1.0, cosDnu));
    const transferAngle = std.math.acos(cosDnu);

    const c = @sqrt(r1 * r1 + r2 * r2 - 2.0 * r1 * r2 * cosDnu);
    const s = (r1 + r2 + c) / 2.0;

    if (s <= 0) {
        return ValidationError.ValueError;
    }

    // initial guess for semi-major axis (minimum energy case)
    const aMin = s / 2.0;
    var a = aMin;

    // for a more robust solver, we would iterate to match flight time
    const expectedTof = std.math.pi * @sqrt((a * a * a) / self.mu);

    const tofRatio = tof / expectedTof;
    if (tofRatio > 0.1 and tofRatio < 10.0) {
        a = aMin * std.math.pow(f64, tofRatio, 2.0 / 3.0);
    }

    const f = 1.0 - a / r1 * (1.0 - cosDnu);
    const sinTransferAngle = @sin(transferAngle);

    if (@abs(sinTransferAngle) < 1e-12) {
        return ValidationError.ValueError;
    }

    const g = r1 * r2 * sinTransferAngle / @sqrt(self.mu * a);

    // in case we have a very small g value
    if (@abs(g) < 1e-15) {
        return ValidationError.ValueError;
    }

    const deltaR = r2Vec.sub(r1Vec);
    const v1Vec = deltaR.sub(r1Vec.mul(f)).mul(1.0 / g);

    const fdot = @sqrt(self.mu * a) / (r1 * r2) * @sin(transferAngle) * (1.0 - cosDnu) / a - 1.0 / r1;
    const gdot = 1.0 - a / r2 * (1.0 - cosDnu);

    const v2Vec = r1Vec.mul(fdot).add(v1Vec.mul(gdot));

    return .{
        .departureVelocity = v1Vec,
        .arrivalVelocity = v2Vec,
        .transferAngle = transferAngle,
        .semiMajorAxis = a,
        .timeOfFlight = tof,
    };
}

//! Base struct that provides various orbital mechanic functions

const std = @import("std");
const log = std.log;
const calculations = @import("calculations.zig");
const constants = @import("constants.zig");
const Vector3D = calculations.Vector3D;
const CelestialBody = constants.CelestialBody;
const ValidationError = error{ InvalidDeltaV, ValueError };

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

pub const LambertResult = struct {
    departureVelocity: Vector3D,
    arrivalVelocity: Vector3D,
    transferAngle: f64,
    semiMajorAxis: f64,
    timeOfFlight: f64,
};

centralBody: CelestialBody,
mu: f64,

pub fn init(mu: f64, centralBody: CelestialBody) OrbitalMechanics {
    return .{
        .centralBody = centralBody,
        .mu = mu,
    };
}

/// Calulcate orbital velocity using vis-viva equation
pub fn orbitalVelocity(self: *OrbitalMechanics, r: f64, a: ?f64) f64 {
    if (r <= 0) {
        log.warn("Distance cannot be negative.");
    }

    if (a != null) {
        return @sqrt(self.mu / r);
    } else {
        if (a <= 0) {
            log.warn("Semi major axis cannot be negative.");
        }
        return @sqrt(self.mu * (2 / r - 1 / a));
    }
}

/// Calculate orbital period using Kepler's third law
pub fn orbitalPeriod(self: *OrbitalMechanics, a: f64) f64 {
    if (a <= 0) {
        log.warn("Semi major axis cannot be negative.");
    }

    return 2 * std.math.pi * @sqrt(std.math.pow(f64, a, 3) / self.mu);
}

/// Calculate escape velocity from a given radius
pub fn escapeVelocity(self: *OrbitalMechanics, r: f64) f64 {
    if (r <= 0) {
        log.warn("Distance cannot be negative.");
    }
    return @sqrt(2 * self.mu / r);
}

pub fn hohmannTransfer(self: *OrbitalMechanics, r1: f64, r2: f64) !TransferResult {
    if (r1 <= 0 or r2 <= 0) {
        return ValidationError.ValueError;
    }
    if (@abs(r1 - r2) < 1000) {
        return ValidationError.ValueError;
    }

    const aTransfer = (r1 + r2) / 2;

    const v1Circular = self.orbitalVelocity(r1, null);
    const v2Circular = self.orbitalVelocity(r2, null);

    const v1Transfer = self.orbitalVelocity(r1, aTransfer);
    const v2Transfer = self.orbitalVelocity(r2, aTransfer);

    const deltaV1 = v1Transfer - v1Circular;
    const deltaV2 = v2Transfer - v2Circular;

    const totalDeltaV = @abs(deltaV1) + @abs(deltaV2);

    const transferTime = std.math.pi * @sqrt(std.math.pow(f64, aTransfer, 3) / self.mu);

    return TransferResult.init(
        aTransfer,
        deltaV1,
        deltaV2,
        totalDeltaV,
        transferTime,
        transferTime / (24 * 3600),
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
        .deltaV1 = deltaV1,
        .deltaV2 = deltaV2,
        .deltaV3 = deltaV3,
        .totalDeltaV = totalDeltaV,
        .totalTime = totalTime,
        .totalTimeDays = totalTime / (24 * 3600),
    };
}

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
    cosDnu = @max(-1, @min(1, cosDnu));
    const sinDnu = @sqrt(1 - std.math.pow(f64, cosDnu, 2));
    const transferAngle = std.math.acos(cosDnu);

    const c = @sqrt(std.math.pow(f64, r1, 2) + std.math.pow(f64, r2, 2) - 2 * r1 * r2 * cosDnu);

    const s = (r1 + r2 + c) / 2;

    const aMin = s / 2;

    const a = aMin;

    const f = 1 - r2 / a * (1 - cosDnu);
    const g = r1 * r2 * sinDnu / @sqrt(self.mu * a);

    const r2MinusFr1 = r2Vec - (r1Vec * f);
    const v1Array = r2MinusFr1 / g;
    const v1Vec: Vector3D = .{ .x = v1Array[0], .y = v1Array[1], .z = v1Array[2] };

    const fdot = @sqrt(self.mu / a) * @tan(transferAngle / 2) * ((s - c) / r1 - (s - c) / r2);
    const gdot = 1 - r1 / a * (1 - cosDnu);

    const v2Array = (r1Vec * fdot + v1Vec * gdot);
    const v2Vec: Vector3D = .{ .x = v2Array[0], .y = v2Array[1], .z = v2Array[2] };

    return .{
        .departureVelocity = v1Vec,
        .arrivalVelocity = v2Vec,
        .transferAngle = transferAngle,
        .semiMajorAxis = a,
        .timeOfFlight = tof,
    };
}

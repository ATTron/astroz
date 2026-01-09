//! Base struct that provides various orbital mechanic functions

const std = @import("std");
const log = std.log;

const calculations = @import("calculations.zig");
const Vector3D = calculations.Vector3D;
const constants = @import("constants.zig");
const CelestialBody = constants.CelestialBody;

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

    pub fn init(
        semiMajorAxis: f64,
        deltaV1: f64,
        deltaV2: f64,
        deltaV3: f64,
        totalDeltaV: f64,
        totalTime: f64,
        totalTimeDays: f64,
    ) !BiEllipticTransferResult {
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
pub fn orbitalVelocity(self: *OrbitalMechanics, radius: f64, semiMajorAxis: ?f64) ValidationError!f64 {
    if (radius <= 0) {
        log.warn("Distance cannot be negative: {d}", .{radius});
        return ValidationError.InvalidDistance;
    }

    if (semiMajorAxis == null) {
        return @sqrt(self.mu / radius);
    } else {
        const sma = semiMajorAxis.?;
        if (sma <= 0) {
            log.warn("Semi major axis cannot be negative: {d}", .{sma});
            return ValidationError.InvalidSMA;
        }
        return @sqrt(self.mu * (2.0 / radius - 1.0 / sma));
    }
}

/// calculate orbital period
pub fn orbitalPeriod(self: *OrbitalMechanics, semiMajorAxis: f64) ValidationError!f64 {
    if (semiMajorAxis <= 0) {
        log.warn("Semi major axis cannot be negative: {d}", .{semiMajorAxis});
        return ValidationError.InvalidSMA;
    }

    return 2.0 * std.math.pi * @sqrt((semiMajorAxis * semiMajorAxis * semiMajorAxis) / self.mu);
}

/// calculate escape velocity from a given radius
pub fn escapeVelocity(self: *OrbitalMechanics, radius: f64) ValidationError!f64 {
    if (radius <= 0) {
        log.warn("Distance cannot be negative: {d}", .{radius});
        return ValidationError.InvalidDistance;
    }
    return @sqrt(2 * self.mu / radius);
}

/// calculate hohmann transfer parameters between two circular orbits
/// note: this calculates heliocentric transfer delta-v
/// ex. earth - mars
/// heliocentric delta-v: ~5.6 km/s (this calculation)
/// LEO to mars transfer: ~3.6 km/s (includes earth's orbital velocity advantage)
pub fn hohmannTransfer(self: *OrbitalMechanics, initialRadius: f64, finalRadius: f64) !TransferResult {
    if (initialRadius <= 0 or finalRadius <= 0) {
        return ValidationError.ValueError;
    }
    if (@abs(initialRadius - finalRadius) < 1000) {
        return ValidationError.ValueError;
    }

    const transferSemiMajorAxis = (initialRadius + finalRadius) / 2.0;

    const velocityInitialCircular = @sqrt(self.mu / initialRadius);
    const velocityFinalCircular = @sqrt(self.mu / finalRadius);

    const velocityInitialTransfer = @sqrt(self.mu / initialRadius) * @sqrt(2.0 * finalRadius / (initialRadius + finalRadius));
    const velocityFinalTransfer = @sqrt(self.mu / finalRadius) * @sqrt(2.0 * initialRadius / (initialRadius + finalRadius));

    const deltaVInitialBurn = velocityInitialTransfer - velocityInitialCircular;
    const deltaVFinalBurn = velocityFinalCircular - velocityFinalTransfer;

    const totalDeltaV = @abs(deltaVInitialBurn) + @abs(deltaVFinalBurn);

    const transferTime = std.math.pi * @sqrt((transferSemiMajorAxis * transferSemiMajorAxis * transferSemiMajorAxis) / self.mu);

    return TransferResult.init(
        transferSemiMajorAxis,
        deltaVInitialBurn,
        deltaVFinalBurn,
        totalDeltaV,
        transferTime,
        transferTime / (24.0 * 3600.0),
    );
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
        .totalTimeDays = totalTime / (24 * 3600),
    };
}

/// simplified lambert solver using universal variable approach
pub fn lambertSolverSimple(self: OrbitalMechanics, initialPositionVec: Vector3D, finalPositionVec: Vector3D, timeOfFlight: f64) !LambertResult {
    if (timeOfFlight <= 0) {
        return ValidationError.ValueError;
    }

    const initialRadius = initialPositionVec.magnitude();
    const finalRadius = finalPositionVec.magnitude();

    if (initialRadius <= 0 or finalRadius <= 0) {
        return ValidationError.ValueError;
    }

    var cosTransferAngle = initialPositionVec.dot(finalPositionVec) / (initialRadius * finalRadius);
    cosTransferAngle = @max(-1.0, @min(1.0, cosTransferAngle));
    const transferAngle = std.math.acos(cosTransferAngle);

    const chordLength = @sqrt(initialRadius * initialRadius + finalRadius * finalRadius - 2.0 * initialRadius * finalRadius * cosTransferAngle);
    const semiperimeter = (initialRadius + finalRadius + chordLength) / 2.0;

    if (semiperimeter <= 0) {
        return ValidationError.ValueError;
    }

    const minSemiMajorAxis = semiperimeter / 2.0;
    var semiMajorAxis = minSemiMajorAxis;

    const expectedTimeOfFlight = std.math.pi * @sqrt((semiMajorAxis * semiMajorAxis * semiMajorAxis) / self.mu);

    const timeOfFlightRatio = timeOfFlight / expectedTimeOfFlight;
    if (timeOfFlightRatio > 0.1 and timeOfFlightRatio < 10.0) {
        semiMajorAxis = minSemiMajorAxis * std.math.pow(f64, timeOfFlightRatio, 2.0 / 3.0);
    }

    const lagrangeF = 1.0 - semiMajorAxis / initialRadius * (1.0 - cosTransferAngle);
    const sinTransferAngle = @sin(transferAngle);

    if (@abs(sinTransferAngle) < 1e-12) {
        return ValidationError.ValueError;
    }

    const lagrangeG = initialRadius * finalRadius * sinTransferAngle / @sqrt(self.mu * semiMajorAxis);

    if (@abs(lagrangeG) < 1e-15) {
        return ValidationError.ValueError;
    }

    const deltaPosition = finalPositionVec.sub(initialPositionVec);
    const initialVelocityVec = deltaPosition.sub(initialPositionVec.mul(lagrangeF)).mul(1.0 / lagrangeG);

    const lagrangeFDot = @sqrt(self.mu * semiMajorAxis) / (initialRadius * finalRadius) * @sin(transferAngle) * (1.0 - cosTransferAngle) / semiMajorAxis - 1.0 / initialRadius;
    const lagrangeGDot = 1.0 - semiMajorAxis / finalRadius * (1.0 - cosTransferAngle);

    const finalVelocityVec = initialPositionVec.mul(lagrangeFDot).add(initialVelocityVec.mul(lagrangeGDot));

    return .{
        .departureVelocity = initialVelocityVec,
        .arrivalVelocity = finalVelocityVec,
        .transferAngle = transferAngle,
        .semiMajorAxis = semiMajorAxis,
        .timeOfFlight = timeOfFlight,
    };
}

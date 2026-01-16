//! Base struct that provides various orbital mechanic functions

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

/// calculate orbital velocity
pub fn orbitalVelocity(self: *OrbitalMechanics, radius: f64, semiMajorAxis: ?f64) ValidationError!f64 {
    if (radius <= 0) return ValidationError.InvalidDistance;
    if (semiMajorAxis) |sma| if (sma <= 0) return ValidationError.InvalidSMA;
    return calculations.orbitalVelocity(self.mu, radius, semiMajorAxis);
}

/// calculate orbital period
pub fn orbitalPeriod(self: *OrbitalMechanics, semiMajorAxis: f64) ValidationError!f64 {
    if (semiMajorAxis <= 0) return ValidationError.InvalidSMA;
    return calculations.orbitalPeriod(self.mu, semiMajorAxis);
}

/// calculate escape velocity from a given radius
pub fn escapeVelocity(self: *OrbitalMechanics, radius: f64) ValidationError!f64 {
    if (radius <= 0) return ValidationError.InvalidDistance;
    return calculations.escapeVelocity(self.mu, radius);
}

/// calculate hohmann transfer parameters between two circular orbits
pub fn hohmannTransfer(self: *OrbitalMechanics, initialRadius: f64, finalRadius: f64) ValidationError!TransferResult {
    if (initialRadius <= 0 or finalRadius <= 0) return ValidationError.ValueError;
    if (@abs(initialRadius - finalRadius) < 1000) return ValidationError.ValueError;

    const r = calculations.hohmannTransfer(self.mu, initialRadius, finalRadius);
    return .{
        .semiMajorAxis = r.semi_major_axis,
        .deltaV1 = r.delta_v1,
        .deltaV2 = r.delta_v2,
        .totalDeltaV = r.total_delta_v,
        .transferTime = r.transfer_time,
        .transferTimeDays = r.transfer_time / constants.seconds_per_day,
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
        .totalTimeDays = totalTime / constants.seconds_per_day,
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

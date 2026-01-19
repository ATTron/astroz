//! Base struct for mission planning
const std = @import("std");
const log = std.log;

const calculations = @import("calculations.zig");
const constants = @import("constants.zig");
const OrbitalMechanics = @import("OrbitalMechanics.zig");
const ValidationError = OrbitalMechanics.ValidationError;

const Mission = @This();

pub const TrajectoryPoint = struct {
    time: f64,
    body: []const u8,
    position: calculations.Vector3D,
    label: []const u8,
};

allocator: std.mem.Allocator,
parameters: MissionParameters,
orbitalMechanics: OrbitalMechanics,
trajectoryPredictions: std.ArrayList(TrajectoryPoint),

pub const PlanetaryPositions = struct {
    planet: constants.CelestialBody,
    position: calculations.Vector3D,
    velocity: ?calculations.Vector3D,
    radius: f64,
    theta: f64,
    time: f64,
};

pub const MissionParameters = struct {
    departureBody: constants.CelestialBody,
    arrivalBody: constants.CelestialBody,
    departureTime: f64,
    arrivalTime: ?f64,
    transferType: []const u8,

    pub fn init(
        departureBody: constants.CelestialBody,
        arrivalBody: constants.CelestialBody,
        departureTime: f64,
        arrivalTime: ?f64,
        transferType: []const u8,
    ) !MissionParameters {
        if (departureTime < 0) {
            return ValidationError.ValueError;
        }

        return .{
            .departureBody = departureBody,
            .arrivalBody = arrivalBody,
            .departureTime = departureTime,
            .arrivalTime = arrivalTime,
            .transferType = transferType,
        };
    }
};

pub const TransferType = union(enum) {
    hohmann: OrbitalMechanics.TransferResult,
    biElliptic: OrbitalMechanics.BiEllipticTransferResult,
};

pub const MissionPlan = struct {
    transfer: TransferType,
    synodicPeriodYears: f64,
    synodicPeriodDays: f64,
    departureBody: constants.CelestialBody,
    arrivalBody: constants.CelestialBody,
    departureOrbitalPeriodDays: f64,
    arrivalOrbitalPeriodDays: f64,
    departureEccentricity: f64,
    arrivalEccentricity: f64,
};

pub fn init(allocator: std.mem.Allocator, parameters: MissionParameters, orbitalMechanics: OrbitalMechanics) Mission {
    return .{
        .allocator = allocator,
        .parameters = parameters,
        .orbitalMechanics = orbitalMechanics,
        .trajectoryPredictions = std.ArrayList(TrajectoryPoint){},
    };
}

pub fn deinit(self: *Mission) void {
    self.trajectoryPredictions.deinit(self.allocator);
}

pub fn printWaypoints(self: *Mission, limit: ?u8) void {
    var count: usize = 0;
    for (self.trajectoryPredictions.items) |point| {
        if (std.mem.eql(u8, point.label, "waypoint")) {
            std.debug.print("Day {d:>6.1}: {s} at ({d:>8.0}, {d:>8.0}) km\n", .{
                point.time,
                point.body,
                point.position.x(),
                point.position.y(),
            });
        }
        if (limit != null and count < limit.?) {
            count += 1;
        }
    }
}

// TODO: probably want a sane limit of like 100 or something but for now this will just dump everything
pub fn printTrajectories(self: *Mission, limit: ?u8) void {
    var count: usize = 0;
    for (self.trajectoryPredictions.items) |point| {
        std.debug.print("Day {d:>6.1}: Transfer at ({d:>8.0}, {d:>8.0}) km\n", .{
            point.time,
            point.position.x(),
            point.position.y(),
        });
        if (limit != null and count < limit.?) {
            count += 1;
        }
    }
}

fn logTransferInfo(
    departureName: []const u8,
    arrivalName: []const u8,
    transfer: OrbitalMechanics.TransferResult,
    arrivalLeadAngle: f64,
) void {
    log.info("{s}-{s} Hohmann Transfer:", .{ departureName, arrivalName });
    log.info("Total Delta-V: {d:.2} km/s", .{transfer.totalDeltaV});
    log.info("Transfer Time: {d:.1} days", .{transfer.transferTimeDays});
    log.info("{s} lead angle required: {d:.1} degrees", .{
        arrivalName,
        arrivalLeadAngle * constants.rad2deg,
    });
}

pub fn propagateTransfer(self: *Mission, totalDays: f64, timeStepDays: f64) !void {
    const departureBody = self.parameters.departureBody;
    const arrivalBody = self.parameters.arrivalBody;

    const departureRadius = departureBody.semiMajorAxis;
    const arrivalRadius = arrivalBody.semiMajorAxis;

    const transfer = try self.orbitalMechanics.hohmannTransfer(departureRadius, arrivalRadius);

    const transferEccentricity = (arrivalRadius - departureRadius) / (arrivalRadius + departureRadius);
    const arrivalAngularVelocity = 2.0 * std.math.pi / arrivalBody.period;
    const arrivalAngleAtTransfer = arrivalAngularVelocity * transfer.transferTimeDays;
    const arrivalLeadAngle = std.math.pi - arrivalAngleAtTransfer;

    logTransferInfo(
        departureBody.name,
        arrivalBody.name,
        transfer,
        arrivalLeadAngle,
    );

    self.trajectoryPredictions.clearRetainingCapacity();

    var day: f64 = 0.0;
    while (day <= totalDays) : (day += timeStepDays) {
        const departureAngle = (day / departureBody.period) * 2.0 * std.math.pi;
        const departureX = departureRadius * @cos(departureAngle);
        const departureY = departureRadius * @sin(departureAngle);

        try self.trajectoryPredictions.append(self.allocator, TrajectoryPoint{
            .time = day,
            .body = departureBody.name,
            .position = calculations.Vector3D.new(departureX, departureY, 0.0),
            .label = "planet",
        });

        const arrivalAngle = arrivalLeadAngle + (day / arrivalBody.period) * 2.0 * std.math.pi;
        const arrivalX = arrivalRadius * @cos(arrivalAngle);
        const arrivalY = arrivalRadius * @sin(arrivalAngle);

        try self.trajectoryPredictions.append(self.allocator, TrajectoryPoint{
            .time = day,
            .body = arrivalBody.name,
            .position = calculations.Vector3D.new(arrivalX, arrivalY, 0.0),
            .label = "planet",
        });

        if (day <= transfer.transferTimeDays) {
            const transferAngle = (day / transfer.transferTimeDays) * std.math.pi;
            const transferRadius = transfer.semiMajorAxis * (1.0 - transferEccentricity * transferEccentricity) /
                (1.0 + transferEccentricity * @cos(transferAngle));

            const transferX = transferRadius * @cos(transferAngle);
            const transferY = transferRadius * @sin(transferAngle);

            try self.trajectoryPredictions.append(self.allocator, TrajectoryPoint{
                .time = day,
                .body = "Transfer",
                .position = calculations.Vector3D.new(transferX, transferY, 0.0),
                .label = "trajectory",
            });
        }
    }

    try self.trajectoryPredictions.append(self.allocator, TrajectoryPoint{
        .time = 0.0,
        .body = "Departure",
        .position = calculations.Vector3D.new(departureRadius, 0.0, 0.0),
        .label = "waypoint",
    });

    const finalArrivalAngle = arrivalLeadAngle +
        (transfer.transferTimeDays / arrivalBody.period) * 2.0 * std.math.pi;

    const arrivalWaypointX = arrivalRadius * @cos(finalArrivalAngle);
    const arrivalWaypointY = arrivalRadius * @sin(finalArrivalAngle);

    try self.trajectoryPredictions.append(self.allocator, TrajectoryPoint{
        .time = transfer.transferTimeDays,
        .body = "Arrival",
        .position = calculations.Vector3D.new(arrivalWaypointX, arrivalWaypointY, 0.0),
        .label = "waypoint",
    });

    log.info("Generated {d} trajectory points", .{self.trajectoryPredictions.items.len});
}

pub fn planetaryPositions(self: *Mission, timeYears: f64) std.ArrayList(PlanetaryPositions) {
    var positions = std.ArrayList(PlanetaryPositions){};

    for (constants.allBodies) |planet| {
        if (std.mem.eql(u8, planet.name, "sun")) {
            continue;
        }

        const orbitSemiMajorAxis = planet.semiMajorAxis * 1000;
        const orbitEccentricity = planet.eccentricity;

        const orbitalPeriodYears = planet.period / 365.25;

        const meanMotion = 2 * std.math.pi / orbitalPeriodYears;
        const meanAnomaly = meanMotion * timeYears;

        const eccentricAnomaly = calculations.solveKeplerEquation(meanAnomaly, orbitEccentricity);

        const trueAnomaly = 2 * std.math.atan2(@sqrt(1 + orbitEccentricity) * @sin(eccentricAnomaly / 2), @sqrt(1 - orbitEccentricity) * @cos(eccentricAnomaly / 2));
        const radius = orbitSemiMajorAxis * (1 - orbitEccentricity * @cos(eccentricAnomaly));

        const xOrbit = radius * @cos(trueAnomaly);
        const yOrbit = radius * @sin(trueAnomaly);

        const inclinationRad = planet.inclination * constants.deg2rad;
        const x = xOrbit;
        const y = yOrbit * @cos(inclinationRad);
        const z = yOrbit * @sin(inclinationRad);

        const position: calculations.Vector3D = calculations.Vector3D.new(x, y, z);
        const planetaryPos = PlanetaryPositions{
            .planet = planet,
            .position = position,
            .velocity = null,
            .radius = radius,
            .theta = trueAnomaly,
            .time = timeYears,
        };

        positions.append(self.allocator, planetaryPos) catch {};
    }

    return positions;
}

/// Calculate a complete mission plan including transfer and synodic period
/// The synodic period is the time between optimal launch windows
pub fn planMission(self: *Mission, params: MissionParameters) !MissionPlan {
    const departureRadius = params.departureBody.semiMajorAxis;
    const arrivalRadius = params.arrivalBody.semiMajorAxis;

    const transfer = if (std.mem.eql(u8, params.transferType, "hohmann"))
        TransferType{ .hohmann = try self.orbitalMechanics.hohmannTransfer(departureRadius, arrivalRadius) }
    else if (std.mem.eql(u8, params.transferType, "biElliptic")) blk: {
        const aphelionRadius = 3 * @max(departureRadius, arrivalRadius);
        break :blk TransferType{ .biElliptic = try self.orbitalMechanics.biEllipicTransfer(departureRadius, arrivalRadius, aphelionRadius) };
    } else return ValidationError.ValueError;

    const departurePeriodYears = params.departureBody.period / 365.25;
    const arrivalPeriodYears = params.arrivalBody.period / 365.25;

    const synodicPeriod = if (@abs(departurePeriodYears - arrivalPeriodYears) < 1e-6)
        std.math.inf(f64)
    else
        @abs(1 / (1 / departurePeriodYears - 1 / arrivalPeriodYears));

    return .{
        .transfer = transfer,
        .synodicPeriodYears = synodicPeriod,
        .synodicPeriodDays = if (std.math.isInf(synodicPeriod)) synodicPeriod else synodicPeriod * 365.25,
        .departureBody = params.departureBody,
        .arrivalBody = params.arrivalBody,
        .departureOrbitalPeriodDays = params.departureBody.period,
        .arrivalOrbitalPeriodDays = params.arrivalBody.period,
        .departureEccentricity = params.departureBody.eccentricity,
        .arrivalEccentricity = params.arrivalBody.eccentricity,
    };
}

test "mission planning Earth-Mars transfer" {
    const testing = std.testing;
    const ta = std.testing.allocator;

    const orbitalMechanics = OrbitalMechanics.init(constants.sun.mu);

    const params = try MissionParameters.init(
        constants.earth,
        constants.mars,
        365.0,
        null,
        "hohmann",
    );

    var mission = Mission.init(ta, params, orbitalMechanics);
    defer mission.deinit();

    try testing.expectEqual(constants.earth, params.departureBody);
    try testing.expectEqual(constants.mars, params.arrivalBody);
    try testing.expectEqual(@as(f64, 365.0), params.departureTime);

    const tYears = 1.0;
    var positions = mission.planetaryPositions(tYears);
    defer positions.deinit(ta);

    try testing.expect(positions.items.len > 0);

    const earthRadius = constants.earth.semiMajorAxis;
    const marsRadius = constants.mars.semiMajorAxis;

    const transfer = try mission.orbitalMechanics.hohmannTransfer(earthRadius, marsRadius);

    try testing.expect(transfer.semiMajorAxis > earthRadius);
    try testing.expect(transfer.semiMajorAxis < marsRadius);

    try testing.expectApproxEqRel(5.6, transfer.totalDeltaV, 0.01);
    try testing.expectApproxEqRel(259, transfer.transferTimeDays, 0.01);
}

test "mission planning parameter validation" {
    const testing = std.testing;

    try testing.expectError(ValidationError.ValueError, MissionParameters.init(constants.earth, constants.mars, -1.0, // Invalid negative time
        null, "hohmann"));

    // Zero time should be valid for examples/demonstrations
    const validParams = try MissionParameters.init(
        constants.earth,
        constants.mars,
        0.0,
        null,
        "hohmann",
    );
    try testing.expectEqual(@as(f64, 0.0), validParams.departureTime);
}

test "bi-elliptic vs hohmann transfer comparison" {
    const testing = std.testing;

    var orbitalMechanics = OrbitalMechanics.init(constants.earth.mu);

    // test transfers from LEO to GEO
    const rLEO = 6700.0; // ~300km altitude
    const rGEO = 42164.0; // geostationary orbit

    const hohmann = try orbitalMechanics.hohmannTransfer(rLEO, rGEO);

    const rAphelion = 2.0 * rGEO;
    const biElliptic = try orbitalMechanics.biEllipicTransfer(rLEO, rGEO, rAphelion);

    try testing.expect(hohmann.totalDeltaV > 0);
    try testing.expect(biElliptic.totalDeltaV > 0);

    try testing.expect(hohmann.totalDeltaV < biElliptic.totalDeltaV);
    // bi-elliptic should take longer
    try testing.expect(biElliptic.totalTimeDays > hohmann.transferTimeDays);
}

test "propagate transfer generates correct trajectory types" {
    const testing = std.testing;
    const ta = std.testing.allocator;

    const orbitalMechanics = OrbitalMechanics.init(constants.sun.mu);
    const params = try MissionParameters.init(constants.earth, constants.mars, 0.0, null, "hohmann");
    var mission = Mission.init(ta, params, orbitalMechanics);
    defer mission.deinit();

    try mission.propagateTransfer(520.0, 5.0);

    try testing.expect(mission.trajectoryPredictions.items.len > 0);

    var waypointCount: usize = 0;
    var planetCount: usize = 0;
    var transferCount: usize = 0;

    for (mission.trajectoryPredictions.items) |point| {
        if (std.mem.eql(u8, point.label, "waypoint")) waypointCount += 1;
        if (std.mem.eql(u8, point.label, "planet")) planetCount += 1;
        if (std.mem.eql(u8, point.label, "trajectory")) transferCount += 1;
    }

    try testing.expectEqual(@as(usize, 2), waypointCount);
    try testing.expect(planetCount > 0);
    try testing.expect(transferCount > 0);
}

test "propagate transfer trajectory validation" {
    const testing = std.testing;
    const ta = std.testing.allocator;

    const orbitalMechanics = OrbitalMechanics.init(constants.sun.mu);
    const params = try MissionParameters.init(
        constants.earth,
        constants.mars,
        0.0,
        null,
        "hohmann",
    );
    var mission = Mission.init(ta, params, orbitalMechanics);
    defer mission.deinit();

    try mission.propagateTransfer(520.0, 10.0);

    // Validate trajectory structure
    var departureWaypoint: ?TrajectoryPoint = null;
    var arrivalWaypoint: ?TrajectoryPoint = null;
    var firstTransferPoint: ?TrajectoryPoint = null;
    var lastTransferPoint: ?TrajectoryPoint = null;

    for (mission.trajectoryPredictions.items) |point| {
        if (std.mem.eql(u8, point.label, "waypoint")) {
            if (point.time == 0.0) {
                departureWaypoint = point;
            } else {
                arrivalWaypoint = point;
            }
        }
        if (std.mem.eql(u8, point.body, "Transfer")) {
            if (firstTransferPoint == null) {
                firstTransferPoint = point;
            }
            lastTransferPoint = point;
        }
    }

    try testing.expect(departureWaypoint != null);
    try testing.expect(arrivalWaypoint != null);
    try testing.expect(firstTransferPoint != null);
    try testing.expect(lastTransferPoint != null);

    try testing.expectEqual(@as(f64, 0.0), departureWaypoint.?.time);
    try testing.expect(arrivalWaypoint.?.time > departureWaypoint.?.time);

    try testing.expect(firstTransferPoint.?.time <= 10.0);
    try testing.expect(lastTransferPoint.?.time <= arrivalWaypoint.?.time);
    try testing.expect(!std.math.isNan(departureWaypoint.?.position.x()));
    try testing.expect(!std.math.isNan(departureWaypoint.?.position.y()));
    try testing.expect(!std.math.isNan(arrivalWaypoint.?.position.x()));
    try testing.expect(!std.math.isNan(arrivalWaypoint.?.position.y()));

    const departureRadius = departureWaypoint.?.position.magnitude();
    try testing.expectApproxEqRel(149598000.0, departureRadius, 0.01);

    const arrivalRadius = arrivalWaypoint.?.position.magnitude();
    try testing.expectApproxEqRel(227939000.0, arrivalRadius, 0.01);
}

test "Lambert solver integration" {
    const testing = std.testing;

    const orbitalMechanics = OrbitalMechanics.init(constants.sun.mu);

    const earthPos = calculations.Vector3D.new(constants.earth.semiMajorAxis, 0.0, 0.0);

    // place mars at 45 degrees ahead in its orbit for a realistic transfer scenario
    const marsAngle = std.math.pi / 4.0; // 45 degrees
    const marsPos = calculations.Vector3D.new(
        constants.mars.semiMajorAxis * @cos(marsAngle),
        constants.mars.semiMajorAxis * @sin(marsAngle),
        0.0,
    );

    // approximate transfer time for Hohmann transfer to Mars (~259 days)
    const transferTime = 259.0 * constants.secondsPerDay;

    const lambertResult = try orbitalMechanics.lambertSolverSimple(earthPos, marsPos, transferTime);

    const depVelMag = lambertResult.departureVelocity.magnitude();
    const arrVelMag = lambertResult.arrivalVelocity.magnitude();

    try testing.expect(!std.math.isNan(depVelMag));
    try testing.expect(!std.math.isNan(arrVelMag));
    try testing.expect(depVelMag > 0 and depVelMag < 100.0);
    try testing.expect(arrVelMag > 0 and arrVelMag < 100.0);
    try testing.expect(lambertResult.transferAngle > 0);
    try testing.expect(lambertResult.semiMajorAxis > constants.earth.semiMajorAxis);
}

test "planMission with Hohmann transfer" {
    const testing = std.testing;
    const ta = std.testing.allocator;

    const orbitalMechanics = OrbitalMechanics.init(constants.sun.mu);
    const params = try MissionParameters.init(
        constants.earth,
        constants.mars,
        0.0,
        null,
        "hohmann",
    );

    var mission = Mission.init(ta, params, orbitalMechanics);
    defer mission.deinit();

    const plan = try mission.planMission(params);

    // Verify it's a Hohmann transfer
    try testing.expect(plan.transfer == .hohmann);

    // Check transfer parameters
    const hohmannTransfer = plan.transfer.hohmann;
    try testing.expect(hohmannTransfer.totalDeltaV > 0);
    try testing.expectApproxEqRel(5.6, hohmannTransfer.totalDeltaV, 0.1);
    try testing.expectApproxEqRel(259.0, hohmannTransfer.transferTimeDays, 1.0);

    // Check synodic period (Earth-Mars synodic period is ~780 days)
    try testing.expect(plan.synodicPeriodDays > 700 and plan.synodicPeriodDays < 800);

    // Verify body information is preserved
    try testing.expectEqual(constants.earth, plan.departureBody);
    try testing.expectEqual(constants.mars, plan.arrivalBody);
}

test "planMission with bi-elliptic transfer" {
    const testing = std.testing;
    const ta = std.testing.allocator;

    const orbitalMechanics = OrbitalMechanics.init(constants.sun.mu);
    const params = try MissionParameters.init(
        constants.earth,
        constants.mars,
        0.0,
        null,
        "biElliptic",
    );

    var mission = Mission.init(ta, params, orbitalMechanics);
    defer mission.deinit();

    const plan = try mission.planMission(params);

    // Verify it's a bi-elliptic transfer
    try testing.expect(plan.transfer == .biElliptic);

    // Check transfer parameters
    const biEllipticTransfer = plan.transfer.biElliptic;
    try testing.expect(biEllipticTransfer.totalDeltaV > 0);
    try testing.expect(biEllipticTransfer.totalTimeDays > 0);

    // Bi-elliptic should take longer than Hohmann
    try testing.expect(biEllipticTransfer.totalTimeDays > 259.0);
}

test "planMission with invalid transfer type" {
    const testing = std.testing;
    const ta = std.testing.allocator;

    const orbitalMechanics = OrbitalMechanics.init(constants.sun.mu);
    const params = try MissionParameters.init(
        constants.earth,
        constants.mars,
        0.0,
        null,
        "invalidType",
    );

    var mission = Mission.init(ta, params, orbitalMechanics);
    defer mission.deinit();

    try testing.expectError(ValidationError.ValueError, mission.planMission(params));
}

test "planetaryPositions calculates valid positions" {
    const testing = std.testing;
    const ta = std.testing.allocator;

    const orbitalMechanics = OrbitalMechanics.init(constants.sun.mu);
    const params = try MissionParameters.init(constants.earth, constants.mars, 0.0, null, "hohmann");
    var mission = Mission.init(ta, params, orbitalMechanics);
    defer mission.deinit();

    var positions = mission.planetaryPositions(1.0);
    defer positions.deinit(ta);

    try testing.expect(positions.items.len > 0);

    for (positions.items) |pos| {
        try testing.expect(pos.radius > 0);
        try testing.expect(!std.math.isNan(pos.position.x()));
        try testing.expect(!std.math.isNan(pos.position.y()));
        try testing.expect(!std.math.isNan(pos.position.z()));
        try testing.expect(!std.math.isNan(pos.theta));
    }
}

test "propagateTransfer edge cases" {
    const testing = std.testing;
    const ta = std.testing.allocator;

    const orbitalMechanics = OrbitalMechanics.init(constants.sun.mu);
    const params = try MissionParameters.init(constants.earth, constants.mars, 0.0, null, "hohmann");
    var mission = Mission.init(ta, params, orbitalMechanics);
    defer mission.deinit();

    try mission.propagateTransfer(0.0, 1.0);
    try testing.expect(mission.trajectoryPredictions.items.len >= 2);

    mission.trajectoryPredictions.clearRetainingCapacity();

    try mission.propagateTransfer(10.0, 20.0);
    try testing.expect(mission.trajectoryPredictions.items.len >= 2);
}

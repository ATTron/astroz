//! Base struct for mission planning

const std = @import("std");
const log = std.log;

const OrbitalMechanics = @import("OrbitalMechanics.zig");
const ValidationError = @import("OrbitalMechanics.zig").ValidationError;
const constants = @import("constants.zig");
const calculations = @import("calculations.zig");

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

    pub fn init(departureBody: constants.CelestialBody, arrivalBody: constants.CelestialBody, departureTime: f64, arrivalTime: ?f64, transferType: []const u8) !MissionParameters {
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

pub fn propagateTransfer(self: *Mission, totalDays: f64, timeStepDays: f64) !void {
    self.trajectoryPredictions.clearRetainingCapacity();

    const departureBody = self.parameters.departureBody;
    const arrivalBody = self.parameters.arrivalBody;

    const departureRadius = departureBody.semiMajorAxis;
    const arrivalRadius = arrivalBody.semiMajorAxis;

    const departurePeriod = departureBody.period;
    const arrivalPeriod = arrivalBody.period;
    const centralMu = self.orbitalMechanics.centralBody.mu;
    const math = std.math;

    const aTransfer = (departureRadius + arrivalRadius) / 2.0;
    const v1Circular = @sqrt(centralMu / departureRadius);
    const v2Circular = @sqrt(centralMu / arrivalRadius);
    const v1Transfer = @sqrt(centralMu / departureRadius) * @sqrt(2.0 * arrivalRadius / (departureRadius + arrivalRadius));
    const v2Transfer = @sqrt(centralMu / arrivalRadius) * @sqrt(2.0 * departureRadius / (departureRadius + arrivalRadius));

    const deltaV1 = v1Transfer - v1Circular;
    const deltaV2 = v2Circular - v2Transfer;
    const totalDeltaV = @abs(deltaV1) + @abs(deltaV2);
    const transferTime = math.pi * @sqrt((aTransfer * aTransfer * aTransfer) / centralMu);
    const transferTimeDays = transferTime / (24.0 * 3600.0);

    log.info("{s}-{s} Hohmann Transfer:", .{ departureBody.name, arrivalBody.name });
    log.info("Total Delta-V: {d:.2} km/s", .{totalDeltaV});
    log.info("Transfer Time: {d:.1} days", .{transferTimeDays});

    const arrivalAngularVelocity = 2.0 * math.pi / arrivalPeriod;
    const arrivalAngleAtTransfer = arrivalAngularVelocity * transferTimeDays;
    const arrivalLeadAngle = math.pi - arrivalAngleAtTransfer;

    log.info("{s} lead angle required: {d:.1} degrees", .{ arrivalBody.name, arrivalLeadAngle * 180.0 / math.pi });

    var day: f64 = 0.0;
    while (day <= totalDays) : (day += timeStepDays) {
        const departureAngle = (day / departurePeriod) * 2.0 * math.pi;
        const departureX = departureRadius * @cos(departureAngle);
        const departureY = departureRadius * @sin(departureAngle);
        const departurePos = calculations.Vector3D.new(departureX, departureY, 0.0);

        try self.trajectoryPredictions.append(self.allocator, TrajectoryPoint{
            .time = day,
            .body = departureBody.name,
            .position = departurePos,
            .label = "planet",
        });

        const arrivalAngle = arrivalLeadAngle + (day / arrivalPeriod) * 2.0 * math.pi;
        const arrivalX = arrivalRadius * @cos(arrivalAngle);
        const arrivalY = arrivalRadius * @sin(arrivalAngle);
        const arrivalPos = calculations.Vector3D.new(arrivalX, arrivalY, 0.0);

        try self.trajectoryPredictions.append(self.allocator, TrajectoryPoint{
            .time = day,
            .body = arrivalBody.name,
            .position = arrivalPos,
            .label = "planet",
        });

        if (day <= transferTimeDays) {
            const transferAngle = (day / transferTimeDays) * math.pi; // 0 to π
            const a = aTransfer;
            const e = (arrivalRadius - departureRadius) / (arrivalRadius + departureRadius);
            const r = a * (1.0 - e * e) / (1.0 + e * @cos(transferAngle));

            const transferX = r * @cos(transferAngle);
            const transferY = r * @sin(transferAngle);
            const transferPos = calculations.Vector3D.new(transferX, transferY, 0.0);

            try self.trajectoryPredictions.append(self.allocator, TrajectoryPoint{
                .time = day,
                .body = "Transfer",
                .position = transferPos,
                .label = "trajectory",
            });
        }
    }

    const departureWaypointPos = calculations.Vector3D.new(departureRadius, 0.0, 0.0);
    try self.trajectoryPredictions.append(self.allocator, TrajectoryPoint{
        .time = 0.0,
        .body = "Departure",
        .position = departureWaypointPos,
        .label = "waypoint",
    });

    const finalArrivalAngle = arrivalLeadAngle + (transferTimeDays / arrivalPeriod) * 2.0 * math.pi;
    const arrivalWaypointX = arrivalRadius * @cos(finalArrivalAngle);
    const arrivalWaypointY = arrivalRadius * @sin(finalArrivalAngle);
    const arrivalWaypointPos = calculations.Vector3D.new(arrivalWaypointX, arrivalWaypointY, 0.0);

    try self.trajectoryPredictions.append(self.allocator, TrajectoryPoint{
        .time = transferTimeDays,
        .body = "Arrival",
        .position = arrivalWaypointPos,
        .label = "waypoint",
    });

    log.info("Generated {d} trajectory points", .{self.trajectoryPredictions.items.len});
}

pub fn planetaryPositions(self: *Mission, tYears: f64) std.ArrayList(PlanetaryPositions) {
    var positions = std.ArrayList(PlanetaryPositions){};

    for (constants.allBodies) |planet| {
        if (std.mem.eql(u8, planet.name, "sun")) {
            continue;
        }

        const a = planet.semiMajorAxis * 1000;
        const e = planet.eccentricity;
        const periodYears = planet.period / 365.25;
        const n = 2 * std.math.pi / periodYears;
        const M = n * tYears;
        var E = M;

        for (0..10) |_| {
            E = E - (E - e * @sin(E) - M) / (1 - e * @cos(E));
        }

        const nu = 2 * std.math.atan2(@sqrt(1 + e) * @sin(E / 2), @sqrt(1 - e) * @cos(E / 2));
        const r = a * (1 - e * @cos(E));

        const xOrbit = r * @cos(nu);
        const yOrbit = r * @sin(nu);

        const inclinationRad = std.math.degreesToRadians(planet.inclination);
        const x = xOrbit;
        const y = yOrbit * @cos(inclinationRad);
        const z = yOrbit * @sin(inclinationRad);

        const position: calculations.Vector3D = calculations.Vector3D.new(x, y, z);
        const planetaryPos = PlanetaryPositions{
            .planet = planet,
            .position = position,
            .velocity = null,
            .radius = r,
            .theta = nu,
            .time = tYears,
        };

        positions.append(self.allocator, planetaryPos) catch {};
    }

    return positions;
}

pub fn planMission(self: *Mission, params: MissionParameters) void {
    const rDeparture = params.departureBody.semiMajorAxis * 1000;
    const rArrival = params.arrivalBody.semiMajorAxis * 1000;

    var transfer = null;
    if (std.mem.eql(u8, params.transferType, "hohmann")) {
        transfer = self.orbitalMechanics.hohmannTransfer(rDeparture, rArrival);
    } else if (std.mem.eql(u8, params.transferType, "bi_elliptic")) {
        const rAphelion = 3 * @max(rDeparture, rArrival);
        transfer = self.orbitalMechanics.biEllipicTransfer(rDeparture, rArrival, rAphelion);
    } else {
        return ValidationError.ValueError;
    }

    const tDeparture = params.departureBody.period / 365.25;
    const tArrival = params.arrivalBody / 365.25;

    var synodicPeriod = null;
    if (@abs(tDeparture - tArrival) < 1e-6) {
        synodicPeriod = std.math.inf(f32);
    } else {
        synodicPeriod = @abs(1 / (1 / tDeparture - 1 / tArrival));
    }

    return .{
        .transfer = transfer,
        .synodicPeriodYears = synodicPeriod,
        .synodicPeriodDays = if (synodicPeriod != std.math.inf(f32)) synodicPeriod * 365.25 else std.math.inf(f32),
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

    const orbitalMechanics = OrbitalMechanics.init(constants.sun.mu, constants.sun);

    const params = try MissionParameters.init(constants.earth, constants.mars, 365.0, null, "hohmann");

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
    const validParams = try MissionParameters.init(constants.earth, constants.mars, 0.0, null, "hohmann");
    try testing.expectEqual(@as(f64, 0.0), validParams.departureTime);
}

test "bi-elliptic vs hohmann transfer comparison" {
    const testing = std.testing;

    var orbitalMechanics = OrbitalMechanics.init(constants.earth.mu, constants.earth);

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

test "planetary orbital mechanics integration" {
    const testing = std.testing;

    const earthPeriod = constants.earth.period;
    const earthSMA = constants.earth.semiMajorAxis;

    var earthOM = OrbitalMechanics.init(constants.sun.mu, constants.sun);
    const earthVelocity = earthOM.orbitalVelocity(earthSMA, null);

    // average orbital velocity should be ~29.78 km/s
    try testing.expectApproxEqRel(29.78, earthVelocity, 0.01);

    const marsSMA = constants.mars.semiMajorAxis; // km
    const marsVelocity = earthOM.orbitalVelocity(marsSMA, null);

    try testing.expect(marsVelocity < earthVelocity);
    try testing.expect(marsVelocity > 20.0 and marsVelocity < 30.0);

    const calculatedEarthPeriod = earthOM.orbitalPeriod(earthSMA) / (24.0 * 3600.0);
    try testing.expectApproxEqRel(earthPeriod, calculatedEarthPeriod, 0.01);
}

test "propagate transfer for different planetary pairs" {
    const testing = std.testing;
    const ta = std.testing.allocator;

    const orbitalMechanics = OrbitalMechanics.init(constants.sun.mu, constants.sun);

    const earthMarsParams = try MissionParameters.init(constants.earth, constants.mars, 0.0, null, "hohmann");
    var earthMarsMission = Mission.init(ta, earthMarsParams, orbitalMechanics);
    defer earthMarsMission.deinit();

    try earthMarsMission.propagateTransfer(520.0, 5.0);

    try testing.expect(earthMarsMission.trajectoryPredictions.items.len > 0);

    var waypointCount: usize = 0;
    var earthCount: usize = 0;
    var marsCount: usize = 0;
    var transferCount: usize = 0;

    for (earthMarsMission.trajectoryPredictions.items) |point| {
        if (std.mem.eql(u8, point.label, "waypoint")) waypointCount += 1;
        if (std.mem.eql(u8, point.body, "earth")) earthCount += 1;
        if (std.mem.eql(u8, point.body, "mars")) marsCount += 1;
        if (std.mem.eql(u8, point.body, "Transfer")) transferCount += 1;
    }

    try testing.expectEqual(@as(usize, 2), waypointCount); // Departure and arrival
    try testing.expect(earthCount > 0);
    try testing.expect(marsCount > 0);
    try testing.expect(transferCount > 0);

    const marsJupiterParams = try MissionParameters.init(constants.mars, constants.jupiter, 0.0, null, "hohmann");
    var marsJupiterMission = Mission.init(ta, marsJupiterParams, orbitalMechanics);
    defer marsJupiterMission.deinit();

    try marsJupiterMission.propagateTransfer(1000.0, 10.0);

    try testing.expect(marsJupiterMission.trajectoryPredictions.items.len > 0);

    var marsJupiterCount: usize = 0;
    var jupiterCount: usize = 0;

    for (marsJupiterMission.trajectoryPredictions.items) |point| {
        if (std.mem.eql(u8, point.body, "mars")) marsJupiterCount += 1;
        if (std.mem.eql(u8, point.body, "jupiter")) jupiterCount += 1;
    }

    try testing.expect(marsJupiterCount > 0);
    try testing.expect(jupiterCount > 0);

    const venusEarthParams = try MissionParameters.init(constants.venus, constants.earth, 0.0, null, "hohmann");
    var venusEarthMission = Mission.init(ta, venusEarthParams, orbitalMechanics);
    defer venusEarthMission.deinit();

    try venusEarthMission.propagateTransfer(300.0, 2.0);

    try testing.expect(venusEarthMission.trajectoryPredictions.items.len > 0);

    var venusCount: usize = 0;
    var venusEarthCount: usize = 0;

    for (venusEarthMission.trajectoryPredictions.items) |point| {
        if (std.mem.eql(u8, point.body, "venus")) venusCount += 1;
        if (std.mem.eql(u8, point.body, "earth")) venusEarthCount += 1;
    }

    try testing.expect(venusCount > 0);
    try testing.expect(venusEarthCount > 0);
}

test "propagate transfer trajectory validation" {
    const testing = std.testing;
    const ta = std.testing.allocator;

    const orbitalMechanics = OrbitalMechanics.init(constants.sun.mu, constants.sun);
    const params = try MissionParameters.init(constants.earth, constants.mars, 0.0, null, "hohmann");
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

    const orbitalMechanics = OrbitalMechanics.init(constants.sun.mu, constants.sun);

    const earthPos = calculations.Vector3D.new(constants.earth.semiMajorAxis, 0.0, 0.0);

    // place mars at 45 degrees ahead in its orbit for a realistic transfer scenario
    const marsAngle = std.math.pi / 4.0; // 45 degrees
    const marsPos = calculations.Vector3D.new(constants.mars.semiMajorAxis * @cos(marsAngle), constants.mars.semiMajorAxis * @sin(marsAngle), 0.0);

    // approximate transfer time for Hohmann transfer to Mars (~259 days)
    const transferTime = 259.0 * 24.0 * 3600.0;

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

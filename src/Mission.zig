//! Base struct for mission planning

const std = @import("std");
const log = std.log;

const OrbitalMechanics = @import("OrbitalMechanics.zig");
const ValidationError = @import("OrbitalMechanics.zig").ValidationError;
const constants = @import("constants.zig");
const calculations = @import("calculations.zig");

const Mission = @This();

allocator: std.mem.Allocator,
parameters: MissionParameters,
orbitalMechanics: OrbitalMechanics,

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
        if (departureTime <= 0) {
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

pub fn planetaryPositions(self: *Mission, tYears: f64) std.ArrayList(PlanetaryPositions) {
    var positions = std.ArrayList(PlanetaryPositions).init(self.allocator);

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

        positions.append(planetaryPos) catch {};
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

    var mission = Mission{
        .allocator = ta,
        .parameters = undefined,
        .orbitalMechanics = orbitalMechanics,
    };

    const params = try MissionParameters.init(constants.earth, constants.mars, 365.0, // 1 year departure time
        null, "hohmann");

    try testing.expectEqual(constants.earth, params.departureBody);
    try testing.expectEqual(constants.mars, params.arrivalBody);
    try testing.expectEqual(@as(f64, 365.0), params.departureTime);

    const tYears = 1.0; // 1 year from J2000
    var positions = mission.planetaryPositions(tYears);
    defer positions.deinit();

    try testing.expect(positions.items.len > 0);

    const earthRadius = constants.earth.semiMajorAxis; // Already in km
    const marsRadius = constants.mars.semiMajorAxis; // Already in km

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

    try testing.expectError(ValidationError.ValueError, MissionParameters.init(constants.earth, constants.mars, 0.0, // Invalid zero time
        null, "hohmann"));
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

    // mars orbital characteristics
    const marsSMA = constants.mars.semiMajorAxis; // km
    const marsVelocity = earthOM.orbitalVelocity(marsSMA, null);

    // mars should be slower than Earth
    try testing.expect(marsVelocity < earthVelocity);
    try testing.expect(marsVelocity > 20.0 and marsVelocity < 30.0);

    // test orbital periods using Kepler's third law
    const calculatedEarthPeriod = earthOM.orbitalPeriod(earthSMA) / (24.0 * 3600.0);
    try testing.expectApproxEqRel(earthPeriod, calculatedEarthPeriod, 0.01);
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

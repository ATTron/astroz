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

pub fn planetaryPositions(self: *Mission, tYears: f64) std.AutoHashMap(constants.CelestialBody, PlanetaryPositions) {
    const positions = std.AutoHashMap(
        constants.CelestialBody,
        PlanetaryPositions,
    ).init(self.allocator);

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

        positions.put(planet.name, position);
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

test "mission planning" {
    const ta = std.testing.allocator;
    std.debug.print("This is the testing allocator: {any}", .{ta});
}

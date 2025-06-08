//! Base struct for mission planning

const std = @import("std");
const log = std.log;

const OrbitalMechanics = @import("OrbitalMechanics.zig");
const ValidationError = @import("OrbitalMechanics.zig").ValidationError;
const constants = @import("constants.zig");
const calculations = @import("calculations.zig");

const Mission = @This();

allocator: std.mem.Allocator,

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
    transferTime: []const u8,

    pub fn init(departureBody: constants.CelestialBody, arrivalBody: constants.CelestialBody, departureTime: f64, arrivalTime: ?f64, transferTime: []const u8) !MissionParameters {
        if (departureTime <= 0) {
            return ValidationError.ValueError;
        }

        return .{
            .departureBody = departureBody,
            .arrivalBody = arrivalBody,
            .departureTime = departureTime,
            .arrivalTime = arrivalTime,
            .transferTime = transferTime,
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

        const position: calculations.Vector3D = .{ .x = x, .y = y, .z = z };

        positions.put(planet.name, position);
    }

    return positions;
}

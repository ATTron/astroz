//! World Coordinate System is commonly used
const std = @import("std");
const constants = @import("constants.zig");
const calculations = @import("calculations.zig");
const Tle = @import("Tle.zig");
const Spacecraft = @import("Spacecraft.zig");

const Matrix3x3 = [3][3]f64;

const Vector3 = [3]f64;

const WorldCoordinateSystem = @This();

x: f64,
y: f64,
z: f64,

/// Currently this function assumes a fully parsed TLE already
pub fn fromTle(tle: Tle, t0: f64, celestial_object: constants.CelestialBody) WorldCoordinateSystem {
    std.log.info("TLE PARSING, {}", .{tle});
    const orbital_elements = calculations.tleToOrbitalElements(tle);

    const eci = orbitalElementsToECI(orbital_elements);
    const ecef = eciToECEF(eci, t0, celestial_object);

    return .{ .x = ecef[0], .y = ecef[1], .z = ecef[2] };
}

fn orbitalElementsToECI(elements: Spacecraft.OrbitalElements) Vector3 {
    const r = elements.a * (1 - elements.e * elements.e) / (1 + elements.e * @cos(elements.true_anomaly));
    const x_orbit = r * @cos(elements.true_anomaly);
    const y_orbit = r * @sin(elements.true_anomaly);

    const R_w = Matrix3x3{
        .{ @cos(elements.arg_periapsis), -@sin(elements.arg_periapsis), 0 },
        .{ @sin(elements.arg_periapsis), @cos(elements.arg_periapsis), 0 },
        .{ 0, 0, 1 },
    };

    const R_i = Matrix3x3{
        .{ 1, 0, 0 },
        .{ 0, @cos(elements.i), -@sin(elements.i) },
        .{ 0, @sin(elements.i), @cos(elements.i) },
    };

    const R_o = Matrix3x3{
        .{ @cos(elements.raan), -@sin(elements.raan), 0 },
        .{ @sin(elements.raan), @cos(elements.raan), 0 },
        .{ 0, 0, 1 },
    };

    const R = calculations.multiplyMatrices(R_o, calculations.multiplyMatrices(R_i, R_w));

    return .{
        R[0][0] * x_orbit + R[0][1] * y_orbit,
        R[1][0] * x_orbit + R[1][1] * y_orbit,
        R[2][0] * x_orbit + R[2][1] * y_orbit,
    };
}

fn eciToECEF(eci: Vector3, time_since_epoch: f64, celestial_object: constants.CelestialBody) Vector3 {
    const m = celestial_object.rotation_rate * time_since_epoch;
    return .{
        eci[0] * @cos(m) + eci[1] * @sin(m),
        -eci[0] * @sin(m) + eci[1] * @cos(m),
        eci[2],
    };
}

test WorldCoordinateSystem {
    const raw_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;
    const expected_ecs: WorldCoordinateSystem = .{ .x = 4.628063569540487e3, .y = -5.164768842168279e3, .z = 7.220776206732921e0 };

    var test_tle = try Tle.parse(raw_tle, std.testing.allocator);
    defer test_tle.deinit();
    const wcs = WorldCoordinateSystem.fromTle(test_tle, 0.0);

    try std.testing.expectEqualDeep(expected_ecs, wcs);
}

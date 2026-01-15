//! World Coordinate System is commonly used
const std = @import("std");

const calculations = @import("calculations.zig");
const constants = @import("constants.zig");
const Spacecraft = @import("Spacecraft.zig");
const Tle = @import("Tle.zig");

const Matrix3x3 = [3][3]f64;

pub const Vector3 = [3]f64;

pub const wgs84_a = constants.wgs84.radiusEarthKm;
pub const wgs84_f = constants.wgs84_flattening;
pub const wgs84_e2 = constants.wgs84_eccentricity_sq;

const WorldCoordinateSystem = @This();

x: f64,
y: f64,
z: f64,

/// Currently this function assumes a fully parsed TLE already
pub fn fromTle(tle: Tle, t0: f64, celestialObject: constants.CelestialBody) WorldCoordinateSystem {
    std.log.info("TLE PARSING, {}", .{tle});
    const orbitalElements = calculations.tleToOrbitalElements(tle);

    const eci = orbitalElementsToECI(orbitalElements);
    const ecef = eciToECEF(eci, t0, celestialObject);

    return .{ .x = ecef[0], .y = ecef[1], .z = ecef[2] };
}

fn orbitalElementsToECI(elements: calculations.OrbitalElements) Vector3 {
    const r = elements.a * (1 - elements.e * elements.e) / (1 + elements.e * @cos(elements.trueAnomaly));
    const xOrbit = r * @cos(elements.trueAnomaly);
    const yOrbit = r * @sin(elements.trueAnomaly);

    const Rw = Matrix3x3{
        .{ @cos(elements.argPeriapsis), -@sin(elements.argPeriapsis), 0 },
        .{ @sin(elements.argPeriapsis), @cos(elements.argPeriapsis), 0 },
        .{ 0, 0, 1 },
    };

    const Ri = Matrix3x3{
        .{ 1, 0, 0 },
        .{ 0, @cos(elements.i), -@sin(elements.i) },
        .{ 0, @sin(elements.i), @cos(elements.i) },
    };

    const Ro = Matrix3x3{
        .{ @cos(elements.raan), -@sin(elements.raan), 0 },
        .{ @sin(elements.raan), @cos(elements.raan), 0 },
        .{ 0, 0, 1 },
    };

    const R = calculations.multiplyMatrices(Ro, calculations.multiplyMatrices(Ri, Rw));

    return .{
        R[0][0] * xOrbit + R[0][1] * yOrbit,
        R[1][0] * xOrbit + R[1][1] * yOrbit,
        R[2][0] * xOrbit + R[2][1] * yOrbit,
    };
}

fn eciToECEF(eci: Vector3, timeSinceEpoch: f64, celestialObject: constants.CelestialBody) Vector3 {
    const m = celestialObject.rotationRate * timeSinceEpoch;
    return .{
        eci[0] * @cos(m) + eci[1] * @sin(m),
        -eci[0] * @sin(m) + eci[1] * @cos(m),
        eci[2],
    };
}

/// convert ECI to ECEF using GMST (Greenwich Mean Sidereal Time).
pub fn eciToEcefGmst(eci: Vector3, gmst: f64) Vector3 {
    const cos_gmst = @cos(gmst);
    const sin_gmst = @sin(gmst);
    return .{
        eci[0] * cos_gmst + eci[1] * sin_gmst,
        -eci[0] * sin_gmst + eci[1] * cos_gmst,
        eci[2],
    };
}

/// convert ECEF to geodetic coordinates (WGS84).
pub fn ecefToGeodetic(ecef: Vector3) Vector3 {
    const x = ecef[0];
    const y = ecef[1];
    const z = ecef[2];

    const lon = std.math.atan2(y, x);
    const p = @sqrt(x * x + y * y);

    var lat = std.math.atan2(z, p * (1.0 - wgs84_e2));
    for (0..10) |_| {
        const lat_prev = lat;
        const sin_lat = @sin(lat);
        const N = wgs84_a / @sqrt(1.0 - wgs84_e2 * sin_lat * sin_lat);
        lat = std.math.atan2(z + wgs84_e2 * N * sin_lat, p);
        if (@abs(lat - lat_prev) < 1e-12) break;
    }

    const sin_lat = @sin(lat);
    const cos_lat = @cos(lat);
    const N = wgs84_a / @sqrt(1.0 - wgs84_e2 * sin_lat * sin_lat);
    const alt = p / cos_lat - N;

    return .{ lat, lon, alt };
}

/// convert ECEF to geodetic coordinates (WGS84) in degrees.
pub fn ecefToGeodeticDeg(ecef: Vector3) Vector3 {
    const lla = ecefToGeodetic(ecef);
    return .{
        lla[0] * constants.rad2deg,
        lla[1] * constants.rad2deg,
        lla[2],
    };
}

/// Compute GMST (Greenwich Mean Sidereal Time) from Julian date.
pub fn julianToGmst(jd: f64) f64 {
    const jd_j2000 = jd - constants.j2000_jd;
    const t = jd_j2000 / constants.julian_days_per_century;
    var gmst = 280.46061837 + 360.98564736629 * jd_j2000 +
        0.000387933 * t * t - t * t * t / 38710000.0;
    gmst = @mod(gmst, 360.0);
    if (gmst < 0) gmst += 360.0;
    return gmst * constants.deg2rad;
}

test WorldCoordinateSystem {
    const raw_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;

    var test_tle = try Tle.parse(raw_tle, std.testing.allocator);
    defer test_tle.deinit();
    const wcs = WorldCoordinateSystem.fromTle(test_tle, 0.0, constants.earth);

    try std.testing.expectApproxEqAbs(wcs.x, 4628.0, 1.0);
    try std.testing.expectApproxEqAbs(wcs.y, -5164.8, 1.0);
    try std.testing.expectApproxEqAbs(wcs.z, 7.2, 1.0);
}

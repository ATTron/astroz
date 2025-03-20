const std = @import("std");
const Tle = @import("Tle.zig");
const OrbitalElements = @import("Spacecraft.zig").OrbitalElements;
const constants = @import("constants.zig");

pub fn degreesToRadians(degrees: f64) f64 {
    return (degrees * std.math.pi) / 180.0;
}

pub fn radiansToDegrees(degrees: f64) f64 {
    return (degrees * 180.0) / std.math.pi;
}

pub fn meanMotionToRadiansPerMinute(mMotion: f64) f64 {
    return mMotion * 2.0 * std.math.pi / (24.0 * 60.0);
}

pub fn meanMotionToSemiMajorAxis(mMotion: f64) f64 {
    return std.math.pow(
        f64,
        (constants.earth.mu / (mMotion * 2.0 * std.math.pi / (24.0 * 3600.0)) * 2.0 * 2.0),
        1.0 / 3.0,
    );
}

pub fn multiplyMatrices(a: [3][3]f64, b: [3][3]f64) [3][3]f64 {
    var result: [3][3]f64 = undefined;
    for (0..3) |i| {
        for (0..3) |j| {
            result[i][j] = 0;
            for (0..3) |k| {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return result;
}

/// convert tle to orbital elements
pub fn tleToOrbitalElements(tle: Tle) OrbitalElements {
    const inclination = degreesToRadians(tle.second_line.inclination);
    const raan = degreesToRadians(tle.second_line.right_ascension);
    const eccentricity = tle.second_line.eccentricity;
    const ArgumentOfPerigee = degreesToRadians(tle.second_line.perigee);
    const mAnomaly = degreesToRadians(tle.second_line.m_anomaly);
    const mMotion = meanMotionToRadiansPerMinute(tle.second_line.m_motion);

    const a = meanMotionToSemiMajorAxis(mMotion);

    const E = solveKeplerEquation(mAnomaly, eccentricity);
    const trueAnomaly = 2.0 * std.math.atan(@sqrt((1 + eccentricity) / (1 - eccentricity)) * @tan(E / 2.0));

    return .{
        .a = a,
        .e = eccentricity,
        .i = inclination,
        .raan = raan,
        .argPeriapsis = ArgumentOfPerigee,
        .trueAnomaly = trueAnomaly,
    };
}

/// convert orbital elements to a usable state vector
pub fn orbitalElementsToStateVector(elements: OrbitalElements, mu: f64) [6]f64 {
    const p = elements.a * (1 - elements.e * elements.e);
    const r = p / (1 + elements.e * @cos(elements.trueAnomaly));

    const rOrbital = [3]f64{
        r * @cos(elements.trueAnomaly),
        r * @sin(elements.trueAnomaly),
        0,
    };
    const vOrbital = [3]f64{
        -@sqrt(mu / p) * @sin(elements.trueAnomaly),
        @sqrt(mu / p) * (elements.e + @cos(elements.trueAnomaly)),
        0,
    };

    const cosRaan = @cos(elements.raan);
    const sinRaan = @sin(elements.raan);
    const cosI = @cos(elements.i);
    const sinI = @sin(elements.i);
    const cosArg = @cos(elements.argPeriapsis);
    const sinArg = @sin(elements.argPeriapsis);

    const rot = [3][3]f64{
        .{ cosRaan * cosArg - sinRaan * sinArg * cosI, -cosRaan * sinArg - sinRaan * cosArg * cosI, sinRaan * sinI },
        .{ sinRaan * cosArg + cosRaan * sinArg * cosI, -sinRaan * sinArg + cosRaan * cosArg * cosI, -cosRaan * sinI },
        .{ sinArg * sinI, cosArg * sinI, cosI },
    };

    var rInertial: [3]f64 = undefined;
    var vInertial: [3]f64 = undefined;

    rInertial[0] = rot[0][0] * rOrbital[0] + rot[0][1] * rOrbital[1] + rot[0][2] * rOrbital[2];
    rInertial[1] = rot[1][0] * rOrbital[0] + rot[1][1] * rOrbital[1] + rot[1][2] * rOrbital[2];
    rInertial[2] = rot[2][0] * rOrbital[0] + rot[2][1] * rOrbital[1] + rot[2][2] * rOrbital[2];

    vInertial[0] = rot[0][0] * vOrbital[0] + rot[0][1] * vOrbital[1] + rot[0][2] * vOrbital[2];
    vInertial[1] = rot[1][0] * vOrbital[0] + rot[1][1] * vOrbital[1] + rot[1][2] * vOrbital[2];
    vInertial[2] = rot[2][0] * vOrbital[0] + rot[2][1] * vOrbital[1] + rot[2][2] * vOrbital[2];

    return .{ rInertial[0], rInertial[1], rInertial[2], vInertial[0], vInertial[1], vInertial[2] };
}

pub fn stateVectorToOrbitalElements(r: [3]f64, v: [3]f64, mu: f64) OrbitalElements {
    const rMag = magnitude(r);
    const vMag = magnitude(v);
    const h = crossProduct(r, v);
    const hMag = magnitude(h);
    const n = crossProduct(.{ 0, 0, 1 }, h);
    const nMag = magnitude(n);

    const eVec = .{
        (vMag * vMag - mu / rMag) * r[0] / mu - dotProduct(r, v) * v[0] / mu,
        (vMag * vMag - mu / rMag) * r[1] / mu - dotProduct(r, v) * v[1] / mu,
        (vMag * vMag - mu / rMag) * r[2] / mu - dotProduct(r, v) * v[2] / mu,
    };
    const e = magnitude(eVec);

    const eps = vMag * vMag / 2.0 - mu / rMag;
    const a = if (@abs(e - 1.0) > 1e-10) -mu / (2.0 * eps) else std.math.inf(f64);
    const i = std.math.acos(h[2] / hMag);
    const raan = std.math.atan2(n[1], n[0]);
    const argPeriapsis = std.math.acos(dotProduct(n, eVec) / (nMag * e));
    const trueAnomaly = std.math.acos(dotProduct(eVec, r) / (e * rMag));

    return .{
        .a = a,
        .e = e,
        .i = i,
        .raan = raan,
        .argPeriapsis = if (r[2] >= 0) argPeriapsis else 2 * std.math.pi - argPeriapsis,
        .trueAnomaly = if (dotProduct(r, v) >= 0) trueAnomaly else 2 * std.math.pi - trueAnomaly,
    };
}

fn dotProduct(a: [3]f64, b: [3]f64) f64 {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

fn crossProduct(a: [3]f64, b: [3]f64) [3]f64 {
    return .{
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    };
}

fn magnitude(v: [3]f64) f64 {
    return std.math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

fn solveKeplerEquation(mAnomaly: f64, eccentricity: f64) f64 {
    var E = mAnomaly;
    const tolerance: f64 = 1e-6;
    var d: f64 = 1.0;

    while (@abs(d) > tolerance) {
        d = E - eccentricity * @sin(E) - mAnomaly;
        E -= d / (1 - eccentricity * @cos(E));
    }

    return E;
}

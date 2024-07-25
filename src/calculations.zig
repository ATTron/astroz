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

pub fn meanMotionToRadiansPerMinute(m_motion: f64) f64 {
    return m_motion * 2.0 * std.math.pi / (24.0 * 60.0);
}

pub fn meanMotionToSemiMajorAxis(m_motion: f64) f64 {
    return std.math.pow(
        f64,
        (constants.earth.mu / (m_motion * 2.0 * std.math.pi / (24.0 * 3600.0)) * 2.0 * 2.0),
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
    const argument_of_perigee = degreesToRadians(tle.second_line.perigee);
    const m_anomaly = degreesToRadians(tle.second_line.m_anomaly);
    const m_motion = meanMotionToRadiansPerMinute(tle.second_line.m_motion);

    const a = meanMotionToSemiMajorAxis(m_motion);

    const E = solveKeplerEquation(m_anomaly, eccentricity);
    const true_anomaly = 2.0 * std.math.atan(@sqrt((1 + eccentricity) / (1 - eccentricity)) * @tan(E / 2.0));

    return .{
        .a = a,
        .e = eccentricity,
        .i = inclination,
        .raan = raan,
        .arg_periapsis = argument_of_perigee,
        .true_anomaly = true_anomaly,
    };
}

/// convert orbital elements to a usable state vector
pub fn orbitalElementsToStateVector(elements: OrbitalElements, mu: f64) [6]f64 {
    const p = elements.a * (1 - elements.e * elements.e);
    const r = p / (1 + elements.e * @cos(elements.true_anomaly));

    const r_orbital = [3]f64{
        r * @cos(elements.true_anomaly),
        r * @sin(elements.true_anomaly),
        0,
    };
    const v_orbital = [3]f64{
        -@sqrt(mu / p) * @sin(elements.true_anomaly),
        @sqrt(mu / p) * (elements.e + @cos(elements.true_anomaly)),
        0,
    };

    const cos_raan = @cos(elements.raan);
    const sin_raan = @sin(elements.raan);
    const cos_i = @cos(elements.i);
    const sin_i = @sin(elements.i);
    const cos_arg = @cos(elements.arg_periapsis);
    const sin_arg = @sin(elements.arg_periapsis);

    const rot = [3][3]f64{
        .{ cos_raan * cos_arg - sin_raan * sin_arg * cos_i, -cos_raan * sin_arg - sin_raan * cos_arg * cos_i, sin_raan * sin_i },
        .{ sin_raan * cos_arg + cos_raan * sin_arg * cos_i, -sin_raan * sin_arg + cos_raan * cos_arg * cos_i, -cos_raan * sin_i },
        .{ sin_arg * sin_i, cos_arg * sin_i, cos_i },
    };

    var r_inertial: [3]f64 = undefined;
    var v_inertial: [3]f64 = undefined;

    r_inertial[0] = rot[0][0] * r_orbital[0] + rot[0][1] * r_orbital[1] + rot[0][2] * r_orbital[2];
    r_inertial[1] = rot[1][0] * r_orbital[0] + rot[1][1] * r_orbital[1] + rot[1][2] * r_orbital[2];
    r_inertial[2] = rot[2][0] * r_orbital[0] + rot[2][1] * r_orbital[1] + rot[2][2] * r_orbital[2];

    v_inertial[0] = rot[0][0] * v_orbital[0] + rot[0][1] * v_orbital[1] + rot[0][2] * v_orbital[2];
    v_inertial[1] = rot[1][0] * v_orbital[0] + rot[1][1] * v_orbital[1] + rot[1][2] * v_orbital[2];
    v_inertial[2] = rot[2][0] * v_orbital[0] + rot[2][1] * v_orbital[1] + rot[2][2] * v_orbital[2];

    return .{ r_inertial[0], r_inertial[1], r_inertial[2], v_inertial[0], v_inertial[1], v_inertial[2] };
}

pub fn stateVectorToOrbitalElements(r: [3]f64, v: [3]f64, mu: f64) OrbitalElements {
    const r_mag = magnitude(r);
    const v_mag = magnitude(v);
    const h = crossProduct(r, v);
    const h_mag = magnitude(h);
    const n = crossProduct(.{ 0, 0, 1 }, h);
    const n_mag = magnitude(n);

    const e_vec = .{
        (v_mag * v_mag - mu / r_mag) * r[0] / mu - dotProduct(r, v) * v[0] / mu,
        (v_mag * v_mag - mu / r_mag) * r[1] / mu - dotProduct(r, v) * v[1] / mu,
        (v_mag * v_mag - mu / r_mag) * r[2] / mu - dotProduct(r, v) * v[2] / mu,
    };
    const e = magnitude(e_vec);

    const eps = v_mag * v_mag / 2.0 - mu / r_mag;
    const a = if (@abs(e - 1.0) > 1e-10) -mu / (2.0 * eps) else std.math.inf(f64);
    const i = std.math.acos(h[2] / h_mag);
    const raan = std.math.atan2(n[1], n[0]);
    const arg_periapsis = std.math.acos(dotProduct(n, e_vec) / (n_mag * e));
    const true_anomaly = std.math.acos(dotProduct(e_vec, r) / (e * r_mag));

    return .{
        .a = a,
        .e = e,
        .i = i,
        .raan = raan,
        .arg_periapsis = if (r[2] >= 0) arg_periapsis else 2 * std.math.pi - arg_periapsis,
        .true_anomaly = if (dotProduct(r, v) >= 0) true_anomaly else 2 * std.math.pi - true_anomaly,
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

fn solveKeplerEquation(m_anomaly: f64, eccentricity: f64) f64 {
    var E = m_anomaly;
    const tolerance: f64 = 1e-6;
    var d: f64 = 1.0;

    while (@abs(d) > tolerance) {
        d = E - eccentricity * @sin(E) - m_anomaly;
        E -= d / (1 - eccentricity * @cos(E));
    }

    return E;
}

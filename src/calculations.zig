//! Various functions to perform the astrodynamics calculations

const std = @import("std");
const Tle = @import("Tle.zig");
const constants = @import("constants.zig");

/// State Vector - Used for position and velocity knowledge
pub const StateV = [6]f64;

/// Contains time and state vector to be used during propagation
pub const StateTime = struct {
    time: f64,
    state: StateV,
};

// an easy win would be to combine this with StateTime to form a *spacecraftState* struct to deal with all of this
/// Responsible for spacecraft orientation
pub const AttitudeState = struct {
    quaternion: [4]f64,
    angularVelocity: [3]f64,
};

/// Needed for propagation.
pub const OrbitalElements = struct {
    a: f64,
    e: f64,
    i: f64,
    raan: f64,
    argPeriapsis: f64,
    trueAnomaly: f64,
};

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
    const inclination = degreesToRadians(tle.secondLine.inclination);
    const raan = degreesToRadians(tle.secondLine.rightAscension);
    const eccentricity = tle.secondLine.eccentricity;
    const ArgumentOfPerigee = degreesToRadians(tle.secondLine.perigee);
    const mAnomaly = degreesToRadians(tle.secondLine.mAnomaly);
    const mMotion = meanMotionToRadiansPerMinute(tle.secondLine.mMotion);

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

pub fn triad(v1_body: [3]f64, v2_body: [3]f64, v1_ref: [3]f64, v2_ref: [3]f64) [3][3]f64 {
    const t1Body = normalize(v1_body);
    const t2Body = normalize(cross(v1_body, v2_body));
    const t3Body = cross(t1Body, t2Body);

    const t1Ref = normalize(v1_ref);
    const t2Ref = normalize(cross(v1_ref, v2_ref));
    const t3Ref = cross(t1Ref, t2Ref);

    const bodyMatrix = [3][3]f64{
        .{ t1Body[0], t2Body[0], t3Body[0] },
        .{ t1Body[1], t2Body[1], t3Body[1] },
        .{ t1Body[2], t2Body[2], t3Body[2] },
    };

    const refMatrix = [3][3]f64{
        .{ t1Ref[0], t2Ref[0], t3Ref[0] },
        .{ t1Ref[1], t2Ref[1], t3Ref[1] },
        .{ t1Ref[2], t2Ref[2], t3Ref[2] },
    };

    return multiplyMatrices(bodyMatrix, transposeMatrix(refMatrix));
}

pub fn normalize(v: [3]f64) [3]f64 {
    const magni = @sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    return .{ v[0] / magni, v[1] / magni, v[2] / magni };
}

pub fn cross(a: [3]f64, b: [3]f64) [3]f64 {
    return .{
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    };
}

pub fn transposeMatrix(m: [3][3]f64) [3][3]f64 {
    return .{
        .{ m[0][0], m[1][0], m[2][0] },
        .{ m[0][1], m[1][1], m[2][1] },
        .{ m[0][2], m[1][2], m[2][2] },
    };
}

pub fn matrixToQuaternion(m: [3][3]f64) [4]f64 {
    var q: [4]f64 = undefined;
    const trace = m[0][0] + m[1][1] + m[2][2];

    if (trace > 0) {
        const s = 0.5 / @sqrt(trace + 1.0);
        q[0] = 0.25 / s;
        q[1] = (m[2][1] - m[1][2]) * s;
        q[2] = (m[0][2] - m[2][0]) * s;
        q[3] = (m[1][0] - m[0][1]) * s;
    } else {
        if (m[0][0] > m[1][1] and m[0][0] > m[2][2]) {
            const s = 2.0 * @sqrt(1.0 + m[0][0] - m[1][1] - m[2][2]);
            q[0] = (m[2][1] - m[1][2]) / s;
            q[1] = 0.25 * s;
            q[2] = (m[0][1] + m[1][0]) / s;
            q[3] = (m[0][2] + m[2][0]) / s;
        } else if (m[1][1] > m[2][2]) {
            const s = 2.0 * @sqrt(1.0 + m[1][1] - m[0][0] - m[2][2]);
            q[0] = (m[0][2] - m[2][0]) / s;
            q[1] = (m[0][1] + m[1][0]) / s;
            q[2] = 0.25 * s;
            q[3] = (m[1][2] + m[2][1]) / s;
        } else {
            const s = 2.0 * @sqrt(1.0 + m[2][2] - m[0][0] - m[1][1]);
            q[0] = (m[1][0] - m[0][1]) / s;
            q[1] = (m[0][2] + m[2][0]) / s;
            q[2] = (m[1][2] + m[2][1]) / s;
            q[3] = 0.25 * s;
        }
    }

    return q;
}

pub fn addAttitudeStates(a: AttitudeState, b: AttitudeState) AttitudeState {
    return .{
        .quaternion = .{
            a.quaternion[0] + b.quaternion[0],
            a.quaternion[1] + b.quaternion[1],
            a.quaternion[2] + b.quaternion[2],
            a.quaternion[3] + b.quaternion[3],
        },
        .angularVelocity = .{
            a.angularVelocity[0] + b.angularVelocity[0],
            a.angularVelocity[1] + b.angularVelocity[1],
            a.angularVelocity[2] + b.angularVelocity[2],
        },
    };
}

pub fn scaleAttitudeState(state: AttitudeState, scalar: f64) AttitudeState {
    return .{
        .quaternion = .{
            state.quaternion[0] * scalar,
            state.quaternion[1] * scalar,
            state.quaternion[2] * scalar,
            state.quaternion[3] * scalar,
        },
        .angularVelocity = .{
            state.angularVelocity[0] * scalar,
            state.angularVelocity[1] * scalar,
            state.angularVelocity[2] * scalar,
        },
    };
}

pub fn addScaledAttitudeState(state: AttitudeState, delta: AttitudeState, scalar: f64) AttitudeState {
    return addAttitudeStates(state, scaleAttitudeState(delta, scalar));
}

pub fn vectorAdd(a: StateV, b: StateV) StateV {
    var result: StateV = undefined;
    for (0..6) |i| {
        result[i] = a[i] + b[i];
    }
    return result;
}

pub fn scalarMultiply(scalar: f64, vector: StateV) StateV {
    var result: StateV = undefined;
    for (0..6) |i| {
        result[i] = scalar * vector[i];
    }
    return result;
}

pub fn impulse(state: StateV, delta_v: [3]f64) StateV {
    return .{
        state[0],              state[1],              state[2],
        state[3] + delta_v[0], state[4] + delta_v[1], state[5] + delta_v[2],
    };
}

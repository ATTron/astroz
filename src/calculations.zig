//! Various functions to perform the astrodynamics calculations

const std = @import("std");

const constants = @import("constants.zig");
const Tle = @import("Tle.zig");

pub const Vector3D = struct {
    data: @Vector(3, f64),

    pub fn new(vx: f64, vy: f64, vz: f64) Vector3D {
        return .{ .data = @Vector(3, f64){ vx, vy, vz } };
    }

    pub fn x(self: Vector3D) f64 {
        return self.data[0];
    }

    pub fn y(self: Vector3D) f64 {
        return self.data[1];
    }

    pub fn z(self: Vector3D) f64 {
        return self.data[2];
    }

    pub fn magnitude(self: Vector3D) f64 {
        return @sqrt(std.math.pow(f64, self.x(), 2) + std.math.pow(f64, self.y(), 2) + std.math.pow(f64, self.z(), 2));
    }

    pub fn array(self: Vector3D) [3]f64 {
        return .{ self.x(), self.y(), self.z() };
    }

    pub fn dot(self: Vector3D, other: Vector3D) f64 {
        return self.x() * other.x() + self.y() * other.y() + self.z() * other.z();
    }

    pub fn add(self: Vector3D, other: Vector3D) Vector3D {
        return .{ .data = @Vector(3, f64){ self.x() + other.x(), self.y() + other.y(), self.z() + other.z() } };
    }

    pub fn sub(self: Vector3D, other: Vector3D) Vector3D {
        return .{ .data = @Vector(3, f64){ self.x() - other.x(), self.y() - other.y(), self.z() - other.z() } };
    }

    pub fn mul(self: Vector3D, scalar: f64) Vector3D {
        return .{ .data = @Vector(3, f64){ self.x() * scalar, self.y() * scalar, self.z() * scalar } };
    }
};

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

pub fn meanMotionToRadiansPerMinute(mMotion: f64) f64 {
    return mMotion * constants.twoPi / constants.minutesPerDay;
}

/// Orbital velocity at given radius. Use sma=null for circular orbit.
pub fn orbitalVelocity(mu: f64, radius: f64, sma: ?f64) f64 {
    if (sma) |a| {
        return @sqrt(mu * (2.0 / radius - 1.0 / a));
    }
    return @sqrt(mu / radius);
}

/// Orbital period from semi-major axis.
pub fn orbitalPeriod(mu: f64, sma: f64) f64 {
    return 2.0 * std.math.pi * @sqrt(sma * sma * sma / mu);
}

/// Escape velocity at given radius.
pub fn escapeVelocity(mu: f64, radius: f64) f64 {
    return @sqrt(2.0 * mu / radius);
}

/// Hohmann transfer result.
pub const HohmannTransfer = struct {
    semiMajorAxis: f64,
    deltaV1: f64,
    deltaV2: f64,
    totalDeltaV: f64,
    transferTime: f64,
};

/// Calculate Hohmann transfer between two circular orbits.
pub fn hohmannTransfer(mu: f64, r1: f64, r2: f64) HohmannTransfer {
    const sma = (r1 + r2) / 2.0;
    const v1Circ = @sqrt(mu / r1);
    const v2Circ = @sqrt(mu / r2);
    const v1Trans = @sqrt(mu / r1) * @sqrt(2.0 * r2 / (r1 + r2));
    const v2Trans = @sqrt(mu / r2) * @sqrt(2.0 * r1 / (r1 + r2));
    const dv1 = v1Trans - v1Circ;
    const dv2 = v2Circ - v2Trans;
    return .{
        .semiMajorAxis = sma,
        .deltaV1 = dv1,
        .deltaV2 = dv2,
        .totalDeltaV = @abs(dv1) + @abs(dv2),
        .transferTime = std.math.pi * @sqrt(sma * sma * sma / mu),
    };
}

pub fn meanMotionToSemiMajorAxis(mMotion: f64) f64 {
    const secondsPerDay = constants.minutesPerDay * 60.0;
    return std.math.pow(
        f64,
        (constants.earth.mu / (mMotion * constants.twoPi / secondsPerDay) * 2.0 * 2.0),
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
    const inclination = tle.secondLine.inclination * constants.deg2rad;
    const raan = tle.secondLine.rightAscension * constants.deg2rad;
    const eccentricity = tle.secondLine.eccentricity;
    const argPeriapsis = tle.secondLine.perigee * constants.deg2rad;
    const mAnomaly = tle.secondLine.mAnomaly * constants.deg2rad;
    const mMotion = meanMotionToRadiansPerMinute(tle.secondLine.mMotion);

    const a = meanMotionToSemiMajorAxis(mMotion);

    const E = solveKeplerEquation(mAnomaly, eccentricity);
    const trueAnomaly = 2.0 * std.math.atan(@sqrt((1 + eccentricity) / (1 - eccentricity)) * @tan(E / 2.0));

    return .{
        .a = a,
        .e = eccentricity,
        .i = inclination,
        .raan = raan,
        .argPeriapsis = argPeriapsis,
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
    const rMag = mag(r);
    const vMag = mag(v);
    const h = cross(r, v);
    const hMag = mag(h);
    const n = cross(.{ 0, 0, 1 }, h);
    const nMag = mag(n);

    const rv = dot(r, v);
    const eVec = .{
        (vMag * vMag - mu / rMag) * r[0] / mu - rv * v[0] / mu,
        (vMag * vMag - mu / rMag) * r[1] / mu - rv * v[1] / mu,
        (vMag * vMag - mu / rMag) * r[2] / mu - rv * v[2] / mu,
    };
    const e = mag(eVec);

    const eps = vMag * vMag / 2.0 - mu / rMag;
    const a = if (@abs(e - 1.0) > 1e-10) -mu / (2.0 * eps) else std.math.inf(f64);
    const i = std.math.acos(h[2] / hMag);
    const raan = std.math.atan2(n[1], n[0]);
    const argPeriapsis = std.math.acos(dot(n, eVec) / (nMag * e));
    const trueAnomaly = std.math.acos(dot(eVec, r) / (e * rMag));

    return .{
        .a = a,
        .e = e,
        .i = i,
        .raan = raan,
        .argPeriapsis = if (r[2] >= 0) argPeriapsis else constants.twoPi - argPeriapsis,
        .trueAnomaly = if (rv >= 0) trueAnomaly else constants.twoPi - trueAnomaly,
    };
}

/// Solve Kepler's equation: E - e*sin(E) = M
/// Returns eccentric anomaly E given mean anomaly M and eccentricity e
/// Uses Newton-Raphson iteration with optional damping for high eccentricity
pub fn solveKeplerEquation(M: f64, e: f64) f64 {
    return solveKeplerNewton(M, e, 1e-6, 50, null);
}

/// Newton-Raphson solver for Kepler's equation with configurable tolerance and damping.
/// - M: mean anomaly (radians)
/// - e: eccentricity
/// - tol: convergence tolerance
/// - maxIter: maximum iterations
/// - damp: optional damping factor (clamps correction magnitude), null for no damping
pub fn solveKeplerNewton(M: f64, e: f64, tol: f64, maxIter: u32, damp: ?f64) f64 {
    var E = M;
    var iter: u32 = 0;

    while (iter < maxIter) {
        const sinE = @sin(E);
        const cosE = @cos(E);
        const f = E - e * sinE - M;
        if (@abs(f) < tol) break;

        var dE = f / (1.0 - e * cosE);
        if (damp) |d| {
            if (@abs(dE) >= d) dE = if (dE > 0) d else -d;
        }
        E -= dE;
        iter += 1;
    }

    return E;
}

pub fn triad(v1Body: [3]f64, v2Body: [3]f64, v1Ref: [3]f64, v2Ref: [3]f64) [3][3]f64 {
    const t1Body = normalize(v1Body);
    const t2Body = normalize(cross(v1Body, v2Body));
    const t3Body = cross(t1Body, t2Body);

    const t1Ref = normalize(v1Ref);
    const t2Ref = normalize(cross(v1Ref, v2Ref));
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
    const m = mag(v);
    return .{ v[0] / m, v[1] / m, v[2] / m };
}

pub fn cross(a: [3]f64, b: [3]f64) [3]f64 {
    return .{
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    };
}

pub fn dot(a: [3]f64, b: [3]f64) f64 {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

pub fn mag(v: [3]f64) f64 {
    return @sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/// Position magnitude from state vector [x, y, z, vx, vy, vz]
pub fn posMag(state: StateV) f64 {
    return @sqrt(state[0] * state[0] + state[1] * state[1] + state[2] * state[2]);
}

/// Velocity magnitude from state vector [x, y, z, vx, vy, vz]
pub fn velMag(state: StateV) f64 {
    return @sqrt(state[3] * state[3] + state[4] * state[4] + state[5] * state[5]);
}

/// Extract position vector from state
pub fn posVec(state: StateV) [3]f64 {
    return .{ state[0], state[1], state[2] };
}

/// Extract velocity vector from state
pub fn velVec(state: StateV) [3]f64 {
    return .{ state[3], state[4], state[5] };
}

pub fn transposeMatrix(m: [3][3]f64) [3][3]f64 {
    return .{
        .{ m[0][0], m[1][0], m[2][0] },
        .{ m[0][1], m[1][1], m[2][1] },
        .{ m[0][2], m[1][2], m[2][2] },
    };
}

// TODO: this can probably be improved somehow
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

/// propagate attitude state using RK4 integration
pub fn propagateAttitude(state: AttitudeState, inertiaTensor: [3][3]f64, dt: f64) AttitudeState {
    const k1 = attitudeDerivative(state, inertiaTensor);
    const k2 = attitudeDerivative(addScaledAttitudeState(state, k1, 0.5 * dt), inertiaTensor);
    const k3 = attitudeDerivative(addScaledAttitudeState(state, k2, 0.5 * dt), inertiaTensor);
    const k4 = attitudeDerivative(addScaledAttitudeState(state, k3, dt), inertiaTensor);

    return addScaledAttitudeState(
        state,
        addAttitudeStates(
            addAttitudeStates(k1, scaleAttitudeState(k2, 2)),
            addAttitudeStates(scaleAttitudeState(k3, 2), k4),
        ),
        dt / 6.0,
    );
}

fn attitudeDerivative(state: AttitudeState, I: [3][3]f64) AttitudeState {
    const q = state.quaternion;
    const w = state.angularVelocity;

    return .{
        .quaternion = .{
            0.5 * (-q[1] * w[0] - q[2] * w[1] - q[3] * w[2]),
            0.5 * (q[0] * w[0] + q[2] * w[2] - q[3] * w[1]),
            0.5 * (q[0] * w[1] - q[1] * w[2] + q[3] * w[0]),
            0.5 * (q[0] * w[2] + q[1] * w[1] - q[2] * w[0]),
        },
        .angularVelocity = .{
            (I[1][1] - I[2][2]) / I[0][0] * w[1] * w[2],
            (I[2][2] - I[0][0]) / I[1][1] * w[2] * w[0],
            (I[0][0] - I[1][1]) / I[2][2] * w[0] * w[1],
        },
    };
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

pub fn impulse(state: StateV, deltaV: [3]f64) StateV {
    return .{
        state[0],             state[1],             state[2],
        state[3] + deltaV[0], state[4] + deltaV[1], state[5] + deltaV[2],
    };
}

test "Vector3D integration tests" {
    const testing = std.testing;

    const pos = Vector3D.new(7000.0, 0.0, 0.0);
    const vel = Vector3D.new(0.0, 7.5, 0.0);

    const r = pos.magnitude();
    const v = vel.magnitude();

    try testing.expectApproxEqAbs(@as(f64, 7000.0), r, 1e-10);
    try testing.expectApproxEqAbs(@as(f64, 7.5), v, 1e-10);

    const specificEnergy = (v * v / 2.0) - (constants.earth.mu / r);
    try testing.expect(specificEnergy < 0);

    const sunDir = Vector3D.new(1.0, 0.0, 0.0);
    const earthDir = Vector3D.new(0.0, 1.0, 0.0);
    const vectorDotProduct = sunDir.dot(earthDir);
    try testing.expectApproxEqAbs(@as(f64, 0.0), vectorDotProduct, 1e-10);

    const rVec = Vector3D.new(1000.0, 2000.0, 3000.0);
    const vVec = Vector3D.new(5.0, -2.0, 1.0);

    const hVec = cross(.{ rVec.x(), rVec.y(), rVec.z() }, .{ vVec.x(), vVec.y(), vVec.z() });
    const hMagnitude = mag(hVec);
    try testing.expect(hMagnitude > 0);

    const testElements = OrbitalElements{
        .a = 6700.0,
        .e = 0.001,
        .i = 51.6 * constants.deg2rad,
        .raan = 339.7 * constants.deg2rad,
        .argPeriapsis = 92.8 * constants.deg2rad,
        .trueAnomaly = 267.1 * constants.deg2rad,
    };

    const stateVector = orbitalElementsToStateVector(testElements, constants.earth.mu);

    const position = [3]f64{ stateVector[0], stateVector[1], stateVector[2] };
    const velocity = [3]f64{ stateVector[3], stateVector[4], stateVector[5] };

    const positionMag = mag(position);
    const velocityMag = mag(velocity);

    try testing.expect(positionMag > 6000.0 and positionMag < 8000.0);
    try testing.expect(velocityMag > 6.0 and velocityMag < 8.0);

    const backToElements = stateVectorToOrbitalElements(position, velocity, constants.earth.mu);
    try testing.expectApproxEqAbs(testElements.i, backToElements.i, 1e-2);
    try testing.expectApproxEqAbs(testElements.e, backToElements.e, 1e-2);
}

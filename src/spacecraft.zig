const std = @import("std");
const constants = @import("constants.zig");
const calculations = @import("calculations.zig");
const TLE = @import("tle.zig").TLE;
const CelestialBody = constants.CelestialBody;
const Datetime = @import("time.zig").Datetime;

/// State Vector - Used for position and velocity knowledge
const StateV = [6]f64;

/// Contains time and state vector to be used during propagation
const StateTime = struct {
    time: f64,
    state: StateV,
};

// an easy win would be to combine this with StateTime to form a *spacecraftState* struct to deal with all of this
/// Responsible for spacecraft orientation
const AttitudeState = struct {
    quaternion: [4]f64,
    angular_velocity: [3]f64,
};

/// Satellite details used in calculations
pub const SatelliteParameters = struct {
    drag: f64,
    cross_section: f64,
    width: f64,
    height: f64,
    depth: f64,
};

/// The impulse maneuver type
pub const Impulse = struct {
    time: f64,
    delta_v: [3]f64,
    mode: enum { Absolute, Prograde, Phase, Plane_Change },
    phase_change: ?f64 = null,
    plane_change: ?struct {
        delta_inclination: f64,
        delta_raan: f64,
    } = null,
};

/// Needed for propagation.
pub const OrbitalElements = struct {
    a: f64,
    e: f64,
    i: f64,
    raan: f64,
    arg_periapsis: f64,
    true_anomaly: f64,
};

/// Enum that helps determine the values in the SatelliteParameters struct
pub const SatelliteSize = enum {
    Cube,
    Medium,
    Large,

    pub fn getDragAndCrossSectional(self: SatelliteSize) SatelliteParameters {
        return switch (self) {
            .Cube => .{
                .drag = 2.2,
                .cross_section = 0.05,
                .width = 0.1,
                .height = 0.1,
                .depth = 0.3,
            },
            .Medium => .{
                .drag = 2.2,
                .cross_section = 5.0,
                .width = 1.4,
                .height = 1.4,
                .depth = 1.6,
            },
            .Large => .{
                .drag = 2.2,
                .cross_section = 50.0,
                .width = 3.2,
                .height = 3.2,
                .depth = 4.0,
            },
        };
    }
};

/// This spacecraft type is the base class that takes the inputs needed to determine future orbit paths
pub const Spacecraft = struct {
    name: []const u8,
    tle: TLE,
    mass: f64,
    size: SatelliteParameters,
    quaternion: [4]f64,
    angular_velocity: [3]f64,
    inertia_tensor: [3][3]f64,
    body_vectors: [2][3]f64,
    reference_vectors: [2][3]f64,
    orbiting_object: CelestialBody = constants.earth,
    orbit_predictions: std.ArrayList(StateTime),
    allocator: std.mem.Allocator,

    pub fn init(name: []const u8, tle: TLE, mass: f64, size: SatelliteSize, orbiting_object: ?CelestialBody, allocator: std.mem.Allocator) Spacecraft {
        return .{
            .name = name,
            .tle = tle,
            .mass = mass,
            .size = size.getDragAndCrossSectional(),
            .quaternion = .{ 1.0, 0.0, 0.0, 0.0 },
            .angular_velocity = .{ 0.0, 0.0, 0.0 },
            .inertia_tensor = .{
                .{ 1.0, 0.0, 0.0 },
                .{ 0.0, 1.0, 0.0 },
                .{ 0.0, 0.0, 1.0 },
            },
            .body_vectors = .{
                .{ 1.0, 0.0, 0.0 },
                .{ 0.0, 1.0, 0.0 },
            },
            .reference_vectors = .{
                .{ 1.0, 0.0, 0.0 },
                .{ 0.0, 1.0, 0.0 },
            },
            .orbiting_object = orbiting_object.?,
            .orbit_predictions = std.ArrayList(StateTime).init(allocator),
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Spacecraft) void {
        self.orbit_predictions.deinit();
    }

    pub fn updateAttitude(self: *Spacecraft) void {
        const attitude_matrix = triad(self.body_vectors[0], self.body_vectors[1], self.reference_vectors[0], self.reference_vectors[1]);
        self.quaternion = matrixToQuaternion(attitude_matrix);
    }

    pub fn propagateAttitude(self: *Spacecraft, dt: f64) void {
        const state = AttitudeState{ .quaternion = self.quaternion, .angular_velocity = self.angular_velocity };
        const new_state = rk4Attitude(self, state, dt);
        self.quaternion = new_state.quaternion;
        self.angular_velocity = new_state.angular_velocity;
    }

    /// This will call the proper propagation methods based on a TLE epoch and recalculation time
    /// Most of the functions this calls are private and will need code inspection to see
    pub fn propagate(self: *Spacecraft, t0: f64, days: f64, h: f64, impulse_list: ?[]const Impulse) !void {
        const y0_oe = self.tleToOrbitalElements();
        const y0 = orbitalElementsToStateVector(y0_oe, self.orbiting_object.mu);
        var t = t0;
        const tf = self.tle.first_line.epoch + days * 86400.0;
        var y = y0;
        const initial_energy = self.calculateEnergy(y);
        std.log.info("Initial energy established: {d}\n", .{initial_energy});
        try self.orbit_predictions.append(StateTime{ .time = t, .state = y });
        var step: usize = 0;
        var impulse_index: usize = 0;

        while (t < tf) : (step += 1) {
            if (impulse_list) |impulses| {
                while (impulse_index < impulses.len and impulses[impulse_index].time <= t + h) {
                    const dt = impulses[impulse_index].time - t;
                    if (dt > 0) {
                        y = self.rk4(y, dt);
                        t += dt;
                        try self.orbit_predictions.append(StateTime{ .time = t, .state = y });
                    }

                    switch (impulses[impulse_index].mode) {
                        .Absolute => {
                            y = impulse(y, impulses[impulse_index].delta_v);
                            std.log.info("Impulse completed at t={d:.2}s", .{t});
                        },
                        .Prograde => {
                            const velocity_magnitude = @sqrt(y[3] * y[3] + y[4] * y[4] + y[5] * y[5]);
                            const delta_v_magnitude = impulses[impulse_index].delta_v[0]; // Use the x-component as magnitude
                            const delta_v = [3]f64{
                                y[3] / velocity_magnitude * delta_v_magnitude,
                                y[4] / velocity_magnitude * delta_v_magnitude,
                                y[5] / velocity_magnitude * delta_v_magnitude,
                            };
                            y = impulse(y, delta_v);
                            std.log.info("Prograde impulse completed at t={d:.2}s", .{t});
                        },
                        .Phase => {
                            const r = @sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
                            const phase_change = impulses[impulse_index].phase_change orelse 0;
                            const transfer_orbits = impulses[impulse_index].delta_v[0];
                            const delta_v_magnitude = self.calculatePhaseChange(r, phase_change, transfer_orbits);
                            const velocity_magnitude = @sqrt(y[3] * y[3] + y[4] * y[4] + y[5] * y[5]);
                            const delta_v = [3]f64{
                                y[3] / velocity_magnitude * delta_v_magnitude,
                                y[4] / velocity_magnitude * delta_v_magnitude,
                                y[5] / velocity_magnitude * delta_v_magnitude,
                            };
                            y = impulse(y, delta_v);

                            const period = 2 * std.math.pi * @sqrt(std.math.pow(f64, r, 3) / self.orbiting_object.mu);
                            const second_impulse_time = t + period * transfer_orbits;
                            while (t < second_impulse_time) : (t += h) {
                                y = self.rk4(y, h);
                                try self.orbit_predictions.append(StateTime{ .time = t + h, .state = y });
                            }
                            y = impulse(y, .{ -delta_v[0], -delta_v[1], -delta_v[2] });
                            try self.orbit_predictions.append(StateTime{ .time = t, .state = y });
                            std.log.info("Phase change completed at t={d:.2}s", .{t});
                        },
                        .Plane_Change => {
                            const plane_change = impulses[impulse_index];
                            y = self.applyPlaneChange(y, plane_change);
                            try self.orbit_predictions.append(StateTime{ .time = t, .state = y });
                            std.log.info("Plane change completed at t={d:.2}s", .{t});
                        },
                    }

                    try self.orbit_predictions.append(StateTime{ .time = t, .state = y });

                    impulse_index += 1;
                }
            }

            // Regular propagation step
            const step_size = @min(h, tf - t);
            y = self.rk4(y, step_size);
            t += step_size;
            const r = @sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
            const energy = self.calculateEnergy(y);

            try self.orbit_predictions.append(StateTime{ .time = t, .state = y });

            if (energy > 0 or std.math.isNan(energy) or r > 100_000) {
                std.log.warn("Abnormal orbit detected at step {d}", .{step});
                std.log.warn("Time: {d} seconds", .{t - t0});
                std.log.warn("State vector: {any}", .{y});
                std.log.warn("Radius: {d} km", .{r});
                std.log.warn("Energy: {d}", .{energy});
                break;
            }
        }
    }

    fn rk4(self: Spacecraft, y: StateV, h: f64) StateV {
        const k1 = scalarMultiply(h, self.derivative(y));
        const k2 = scalarMultiply(h, self.derivative(vectorAdd(y, scalarMultiply(0.5, k1))));
        const k3 = scalarMultiply(h, self.derivative(vectorAdd(y, scalarMultiply(0.5, k2))));
        const k4 = scalarMultiply(h, self.derivative(vectorAdd(y, k3)));
        return vectorAdd(
            y,
            scalarMultiply(
                1.0 / 6.0,
                vectorAdd(
                    vectorAdd(k1, scalarMultiply(2, k2)),
                    vectorAdd(scalarMultiply(2, k3), k4),
                ),
            ),
        );
    }

    fn triad(v1_body: [3]f64, v2_body: [3]f64, v1_ref: [3]f64, v2_ref: [3]f64) [3][3]f64 {
        const t1_body = normalize(v1_body);
        const t2_body = normalize(cross(v1_body, v2_body));
        const t3_body = cross(t1_body, t2_body);

        const t1_ref = normalize(v1_ref);
        const t2_ref = normalize(cross(v1_ref, v2_ref));
        const t3_ref = cross(t1_ref, t2_ref);

        const body_matrix = [3][3]f64{
            .{ t1_body[0], t2_body[0], t3_body[0] },
            .{ t1_body[1], t2_body[1], t3_body[1] },
            .{ t1_body[2], t2_body[2], t3_body[2] },
        };

        const ref_matrix = [3][3]f64{
            .{ t1_ref[0], t2_ref[0], t3_ref[0] },
            .{ t1_ref[1], t2_ref[1], t3_ref[1] },
            .{ t1_ref[2], t2_ref[2], t3_ref[2] },
        };

        return multiplyMatrices(body_matrix, transposeMatrix(ref_matrix));
    }

    fn normalize(v: [3]f64) [3]f64 {
        const magni = @sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        return .{ v[0] / magni, v[1] / magni, v[2] / magni };
    }

    fn cross(a: [3]f64, b: [3]f64) [3]f64 {
        return .{
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        };
    }

    fn multiplyMatrices(a: [3][3]f64, b: [3][3]f64) [3][3]f64 {
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

    fn transposeMatrix(m: [3][3]f64) [3][3]f64 {
        return .{
            .{ m[0][0], m[1][0], m[2][0] },
            .{ m[0][1], m[1][1], m[2][1] },
            .{ m[0][2], m[1][2], m[2][2] },
        };
    }

    fn matrixToQuaternion(m: [3][3]f64) [4]f64 {
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

    fn rk4Attitude(spacecraft: *Spacecraft, state: AttitudeState, dt: f64) AttitudeState {
        const k1 = attitudeDerivative(spacecraft, state);
        const k2 = attitudeDerivative(spacecraft, addScaledAttitudeState(state, k1, 0.5 * dt));
        const k3 = attitudeDerivative(spacecraft, addScaledAttitudeState(state, k2, 0.5 * dt));
        const k4 = attitudeDerivative(spacecraft, addScaledAttitudeState(state, k3, dt));

        return addScaledAttitudeState(state, addAttitudeStates(addAttitudeStates(k1, scaleAttitudeState(k2, 2)), addAttitudeStates(scaleAttitudeState(k3, 2), k4)), dt / 6.0);
    }

    fn attitudeDerivative(spacecraft: *Spacecraft, state: AttitudeState) AttitudeState {
        const q = state.quaternion;
        const w = state.angular_velocity;

        const q_dot = [4]f64{
            0.5 * (-q[1] * w[0] - q[2] * w[1] - q[3] * w[2]),
            0.5 * (q[0] * w[0] + q[2] * w[2] - q[3] * w[1]),
            0.5 * (q[0] * w[1] - q[1] * w[2] + q[3] * w[0]),
            0.5 * (q[0] * w[2] + q[1] * w[1] - q[2] * w[0]),
        };

        const I = spacecraft.inertia_tensor;
        const w_dot = [3]f64{
            (I[1][1] - I[2][2]) / I[0][0] * w[1] * w[2],
            (I[2][2] - I[0][0]) / I[1][1] * w[2] * w[0],
            (I[0][0] - I[1][1]) / I[2][2] * w[0] * w[1],
        };

        return AttitudeState{
            .quaternion = q_dot,
            .angular_velocity = w_dot,
        };
    }

    fn addAttitudeStates(a: AttitudeState, b: AttitudeState) AttitudeState {
        return AttitudeState{
            .quaternion = .{
                a.quaternion[0] + b.quaternion[0],
                a.quaternion[1] + b.quaternion[1],
                a.quaternion[2] + b.quaternion[2],
                a.quaternion[3] + b.quaternion[3],
            },
            .angular_velocity = .{
                a.angular_velocity[0] + b.angular_velocity[0],
                a.angular_velocity[1] + b.angular_velocity[1],
                a.angular_velocity[2] + b.angular_velocity[2],
            },
        };
    }

    fn scaleAttitudeState(state: AttitudeState, scalar: f64) AttitudeState {
        return AttitudeState{
            .quaternion = .{
                state.quaternion[0] * scalar,
                state.quaternion[1] * scalar,
                state.quaternion[2] * scalar,
                state.quaternion[3] * scalar,
            },
            .angular_velocity = .{
                state.angular_velocity[0] * scalar,
                state.angular_velocity[1] * scalar,
                state.angular_velocity[2] * scalar,
            },
        };
    }

    fn addScaledAttitudeState(state: AttitudeState, delta: AttitudeState, scalar: f64) AttitudeState {
        return addAttitudeStates(state, scaleAttitudeState(delta, scalar));
    }

    fn vectorAdd(a: StateV, b: StateV) StateV {
        var result: StateV = undefined;
        for (0..6) |i| {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    fn crossProduct(a: [3]f64, b: [3]f64) [3]f64 {
        return .{
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        };
    }

    fn dotProduct(a: [3]f64, b: [3]f64) f64 {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    fn magnitude(v: [3]f64) f64 {
        return std.math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }

    fn scalarMultiply(scalar: f64, vector: StateV) StateV {
        var result: StateV = undefined;
        for (0..6) |i| {
            result[i] = scalar * vector[i];
        }
        return result;
    }

    fn derivative(self: Spacecraft, state: StateV) StateV {
        const x = state[0];
        const y = state[1];
        const z = state[2];
        const vx = state[3];
        const vy = state[4];
        const vz = state[5];
        const r = @sqrt(x * x + y * y + z * z);
        const r2 = r * r;
        const r3 = r2 * r;
        const r5 = r3 * r2;

        var ax = -self.orbiting_object.mu * x / r3;
        var ay = -self.orbiting_object.mu * y / r3;
        var az = -self.orbiting_object.mu * z / r3;

        const j2_factor = -1.5 * self.orbiting_object.j2_perturbation * self.orbiting_object.mu * self.orbiting_object.eq_radius.? * self.orbiting_object.eq_radius.? / r5;
        ax += j2_factor * x * (5 * z * z / r2 - 1);
        ay += j2_factor * y * (5 * z * z / r2 - 1);
        az += j2_factor * z * (5 * z * z / r2 - 3);

        const altitude = r - self.orbiting_object.eq_radius.?;
        if (altitude < 1000) { // Only apply drag below 1000 km
            const v = @sqrt(vx * vx + vy * vy + vz * vz);
            const rho = self.atmosphericDensity(altitude);
            const drag_factor = -0.5 * self.size.drag *
                self.size.cross_section *
                rho * v / self.mass;
            ax += drag_factor * vx;
            ay += drag_factor * vy;
            az += drag_factor * vz;
        }

        return .{ vx, vy, vz, ax, ay, az };
    }

    fn atmosphericDensity(self: Spacecraft, altitude: f64) f64 {
        return self.orbiting_object.sea_level_density * std.math.exp(-altitude / self.orbiting_object.scale_height);
    }

    fn calculateEnergy(self: Spacecraft, state: StateV) f64 {
        const r = @sqrt(state[0] * state[0] + state[1] * state[1] + state[2] * state[2]);
        const v_squared = state[3] * state[3] + state[4] * state[4] + state[5] * state[5];
        return 0.5 * v_squared - self.orbiting_object.mu / r;
    }

    fn applyPlaneChange(self: *Spacecraft, y: [6]f64, plane_change: Impulse) [6]f64 {
        const mu = self.orbiting_object.mu;
        const r = [3]f64{ y[0], y[1], y[2] };
        const v = [3]f64{ y[3], y[4], y[5] };

        var elements = stateVectorToOrbitalElements(r, v, mu);
        elements.i += plane_change.plane_change.?.delta_inclination;
        elements.raan += plane_change.plane_change.?.delta_raan;

        return orbitalElementsToStateVector(elements, mu);
    }

    fn tleToOrbitalElements(self: Spacecraft) OrbitalElements {
        const inclination = calculations.degreesToRadians(self.tle.second_line.inclination);
        const raan = calculations.degreesToRadians(self.tle.second_line.right_ascension);
        const eccentricity = self.tle.second_line.eccentricity;
        const argument_of_perigee = calculations.degreesToRadians(self.tle.second_line.perigee);
        const m_anomaly = calculations.degreesToRadians(self.tle.second_line.m_anomaly);
        const m_motion = calculations.meanMotionToRadiansPerMinute(self.tle.second_line.m_motion);

        const a = calculations.meanMotionToSemiMajorAxis(m_motion);

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

        return OrbitalElements{
            .a = a,
            .e = e,
            .i = i,
            .raan = raan,
            .arg_periapsis = if (r[2] >= 0) arg_periapsis else 2 * std.math.pi - arg_periapsis,
            .true_anomaly = if (dotProduct(r, v) >= 0) true_anomaly else 2 * std.math.pi - true_anomaly,
        };
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

    fn impulse(state: StateV, delta_v: [3]f64) StateV {
        return .{
            state[0],              state[1],              state[2],
            state[3] + delta_v[0], state[4] + delta_v[1], state[5] + delta_v[2],
        };
    }

    fn calculatePhaseChange(self: Spacecraft, radius: f64, phase_change: f64, transfer_orbits: f64) f64 {
        const v_circular = @sqrt(self.orbiting_object.mu / radius);
        const period = 2 * std.math.pi * @sqrt(std.math.pow(f64, radius, 3) / self.orbiting_object.mu);
        const delta_t = phase_change * period / (2 * std.math.pi);
        const a_transfer = std.math.pow(f64, (period + delta_t) / (2 * std.math.pi) * @sqrt(self.orbiting_object.mu), 2.0 / 3.0);
        const v_transfer = @sqrt(self.orbiting_object.mu * (2 / radius - 1 / a_transfer));
        return (v_transfer - v_circular) * 2 / transfer_orbits;
    }
};

test "init spacecraft" {
    const raw_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;
    var test_tle = try TLE.parse(raw_tle, std.testing.allocator);
    defer test_tle.deinit();

    var test_sc = Spacecraft.init("dummy_sc", test_tle, 300.000, SatelliteSize.Cube, constants.earth, std.testing.allocator);
    defer test_sc.deinit();

    try test_sc.propagate(
        test_sc.tle.first_line.epoch,
        3, // days to predict
        1, // steps, i.e. predict every simulated second
        null,
    );

    for (test_sc.orbit_predictions.items) |iter| {
        const r = std.math.sqrt(iter.state[0] * iter.state[0] + iter.state[1] * iter.state[1] + iter.state[2] * iter.state[2]);

        try std.testing.expect(r > test_sc.orbiting_object.eq_radius.?);
    }

    // write out to file for creating test mapping data
    // const file = try std.fs.cwd().createFile("./test/orbit_data.csv", .{});
    // defer file.close();
    // const writer = file.writer();
    //
    // try writer.writeAll("time,x,y,z\n");
    //
    // for (test_sc.orbit_predictions.items) |item| {
    //     try writer.print("{d},{d},{d},{d}\n", .{ item.time, item.state[0], item.state[1], item.state[2] });
    // }
    //
    // std.debug.print("Orbit data written to orbit_data.csv\n", .{});
}

test "prop spacecraft w/ impulse" {
    const raw_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;
    var test_tle = try TLE.parse(raw_tle, std.testing.allocator);
    defer test_tle.deinit();

    var test_sc = Spacecraft.init("dummy_sc", test_tle, 300.000, SatelliteSize.Cube, constants.earth, std.testing.allocator);
    defer test_sc.deinit();

    const impulses = [_]Impulse{
        .{ .time = 2635014.50, .delta_v = .{ 0.2, 0.0, 0.0 }, .mode = .Prograde },
        .{ .time = 2638026.50, .delta_v = .{ 0.2, 0.0, 0.0 }, .mode = .Prograde },
        .{ .time = 2638103.50, .delta_v = .{ 0.2, 0.0, 0.0 }, .mode = .Prograde },
    };

    try test_sc.propagate(
        test_sc.tle.first_line.epoch,
        3, // days to predict
        1, // steps, i.e. predict every simulated second
        &impulses,
    );

    for (test_sc.orbit_predictions.items) |iter| {
        const r = std.math.sqrt(iter.state[0] * iter.state[0] + iter.state[1] * iter.state[1] + iter.state[2] * iter.state[2]);

        try std.testing.expect(r > test_sc.orbiting_object.eq_radius.?);
    }
}

test "prop spacecraft w/ phase" {
    const raw_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;
    var test_tle = try TLE.parse(raw_tle, std.testing.allocator);
    defer test_tle.deinit();

    var test_sc = Spacecraft.init("dummy_sc", test_tle, 300.000, SatelliteSize.Cube, constants.earth, std.testing.allocator);
    defer test_sc.deinit();

    const phase_maneuver = Impulse{
        .time = 2500000.0,
        .delta_v = .{ 1.0, 0.0, 0.0 },
        .mode = .Phase,
        .phase_change = std.math.pi / 2.0,
    };

    const impulses = [_]Impulse{phase_maneuver};

    try test_sc.propagate(
        test_sc.tle.first_line.epoch,
        3, // days to predict
        1, // steps, i.e. predict every simulated second
        &impulses,
    );

    for (test_sc.orbit_predictions.items) |iter| {
        const r = std.math.sqrt(iter.state[0] * iter.state[0] + iter.state[1] * iter.state[1] + iter.state[2] * iter.state[2]);

        try std.testing.expect(r > test_sc.orbiting_object.eq_radius.?);
    }

    // const file = try std.fs.cwd().createFile("./test/orbit_data_with_phase.csv", .{});
    // defer file.close();
    // const writer = file.writer();
    //
    // try writer.writeAll("time,x,y,z\n");
    //
    // for (test_sc.orbit_predictions.items) |item| {
    //     try writer.print("{d},{d},{d},{d}\n", .{ item.time, item.state[0], item.state[1], item.state[2] });
    // }
    //
    // std.debug.print("Orbit data written to orbit_data.csv\n", .{});
}

test "prop spacecraft w/ plane change" {
    const raw_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;
    var test_tle = try TLE.parse(raw_tle, std.testing.allocator);
    defer test_tle.deinit();

    var test_sc = Spacecraft.init("dummy_sc", test_tle, 300.000, SatelliteSize.Cube, constants.earth, std.testing.allocator);
    defer test_sc.deinit();

    const plane_change_maneuver = Impulse{
        .time = 2500000.0,
        .delta_v = .{ 0.0, 0.0, 0.0 }, // Not used for plane changes
        .mode = .Plane_Change,
        .plane_change = .{
            .delta_inclination = std.math.pi / 18.0, // 10-degree inclination change
            .delta_raan = std.math.pi / 36.0, // 5-degree RAAN change
        },
    };

    const impulses = [_]Impulse{plane_change_maneuver};

    try test_sc.propagate(
        test_sc.tle.first_line.epoch,
        3, // days to predict
        1, // steps, i.e. predict every simulated second
        &impulses,
    );

    for (test_sc.orbit_predictions.items) |iter| {
        const r = std.math.sqrt(iter.state[0] * iter.state[0] + iter.state[1] * iter.state[1] + iter.state[2] * iter.state[2]);

        try std.testing.expect(r > test_sc.orbiting_object.eq_radius.?);
    }

    // const file = try std.fs.cwd().createFile("./test/orbit_data_with_plane_change.csv", .{});
    // defer file.close();
    // const writer = file.writer();
    //
    // try writer.writeAll("time,x,y,z\n");
    //
    // for (test_sc.orbit_predictions.items) |item| {
    //     try writer.print("{d},{d},{d},{d}\n", .{ item.time, item.state[0], item.state[1], item.state[2] });
    // }
    //
    // std.debug.print("Orbit data written to orbit_data.csv\n", .{});
}

test "orientation determination testing" {
    const raw_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;
    var test_tle = try TLE.parse(raw_tle, std.testing.allocator);
    defer test_tle.deinit();
    var spacecraft = Spacecraft.init("dummy_sc", test_tle, 300.000, SatelliteSize.Cube, constants.earth, std.testing.allocator);

    spacecraft.angular_velocity = .{ 0.1, 0.05, 0.02 };

    // const file = try std.fs.cwd().createFile("test/attitude_data.csv", .{});
    // defer file.close();
    // const writer = file.writer();
    //
    // try writer.writeAll("Time,qw,qx,qy,qz,wx,wy,wz\n");

    const dt = 60.0; // 1 minute time step
    const simulation_time = 3 * 24 * 60 * 60.0; // 3 days in seconds
    const orbital_period = 90 * 60.0; // 90 minutes orbital period
    var t: f64 = 0;
    while (t < simulation_time) : (t += dt) {
        const angle = 0.5 * @sin(2 * std.math.pi * t / orbital_period);

        spacecraft.body_vectors[0] = .{ @cos(angle), 0, @sin(angle) };
        spacecraft.body_vectors[1] = .{ 0, 1, 0 };

        spacecraft.updateAttitude();
        spacecraft.propagateAttitude(dt);

        // Simulate simple circular orbit
        const orbit_radius = 7000.0;
        // const x = orbit_radius * @cos(2 * std.math.pi * t / orbital_period);
        // const y = orbit_radius * @sin(2 * std.math.pi * t / orbital_period);
        // const z = 0.0;

        try std.testing.expect(orbit_radius > spacecraft.orbiting_object.m_radius.?);

        // try writer.print("{d},{d},{d},{d},{d},{d},{d},{d},{d},{d},{d}\n", .{
        //     t,
        //     spacecraft.quaternion[0],
        //     spacecraft.quaternion[1],
        //     spacecraft.quaternion[2],
        //     spacecraft.quaternion[3],
        //     spacecraft.angular_velocity[0],
        //     spacecraft.angular_velocity[1],
        //     spacecraft.angular_velocity[2],
        //     x,
        //     y,
        //     z,
        // });
    }
}

test "orientation determination with dramatic changes" {
    const raw_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;
    var test_tle = try TLE.parse(raw_tle, std.testing.allocator);
    defer test_tle.deinit();
    var spacecraft = Spacecraft.init("dummy_sc", test_tle, 300.000, SatelliteSize.Cube, constants.earth, std.testing.allocator);
    defer spacecraft.deinit();

    // Initial angular velocity
    spacecraft.angular_velocity = .{ 0.0, 0.0, 0.0 };

    // const file = try std.fs.cwd().createFile("test/dramatic_attitude_data.csv", .{});
    // defer file.close();
    // const writer = file.writer();

    // try writer.writeAll("Time,qw,qx,qy,qz,wx,wy,wz,x,y,z\n");

    const dt = 120.0; // 2 mins time step
    const simulation_time = 3 * 24 * 60 * 60.0; // 3 days in seconds
    const orbital_period = 90 * 60.0; // 90 minutes orbital period
    var t: f64 = 0;

    while (t < simulation_time) : (t += dt) {
        // Simulate a dramatic torque effect
        const torque_x = 0.001 * @sin(2 * std.math.pi * t / (orbital_period * 2));
        const torque_y = 0.0005 * @cos(2 * std.math.pi * t / (orbital_period * 3));
        const torque_z = 0.0002 * @sin(2 * std.math.pi * t / orbital_period);

        // Update angular velocity based on torque (simplified)
        spacecraft.angular_velocity[0] += torque_x * dt;
        spacecraft.angular_velocity[1] += torque_y * dt;
        spacecraft.angular_velocity[2] += torque_z * dt;

        // Update attitude
        spacecraft.updateAttitude();
        spacecraft.propagateAttitude(dt);

        // Simulate simple circular orbit
        const orbit_radius = 7000.0;
        // const x = orbit_radius * @cos(2 * std.math.pi * t / orbital_period);
        // const y = orbit_radius * @sin(2 * std.math.pi * t / orbital_period);
        // const z = 0.0;

        try std.testing.expect(orbit_radius > spacecraft.orbiting_object.m_radius.?);

        // try writer.print("{d},{d},{d},{d},{d},{d},{d},{d},{d},{d},{d}\n", .{
        //     t,
        //     spacecraft.quaternion[0],
        //     spacecraft.quaternion[1],
        //     spacecraft.quaternion[2],
        //     spacecraft.quaternion[3],
        //     spacecraft.angular_velocity[0],
        //     spacecraft.angular_velocity[1],
        //     spacecraft.angular_velocity[2],
        //     x,
        //     y,
        //     z,
        // });
    }
}

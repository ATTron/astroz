const std = @import("std");
const math = std.math;
const constants = @import("constants.zig");
const calculations = @import("calculations.zig");
const TLE = @import("tle.zig").TLE;
const Celistal_Body = constants.Celestial_Body;
const Datetime = @import("time.zig").Datetime;

/// State Vector - Used for position and velcity knowledge
const StateV = [6]f64;
/// Contains time and state vector to be used during propagation
const StateTime = struct {
    time: f64,
    state: StateV,
};

/// Satellite details used in calculations
pub const Satellite_Parameters = struct {
    drag: f64,
    cross_section: f64,
};

/// The impulse maneuver type
pub const Impulse = struct {
    time: f64,
    delta_v: [3]f64,
    mode: enum { Absolute, Prograde, Phase },
    phase_change: ?f64 = null,
};

/// Enum that helps determine the values in the Satellite_Parameters struct
pub const Satellite_Size = enum {
    Cube,
    Medium,
    Large,

    pub fn get_drag_and_cross_sectional(self: @This()) Satellite_Parameters {
        return switch (self) {
            .Cube => .{
                .drag = 2.2,
                .cross_section = 0.05,
            },
            .Medium => .{
                .drag = 2.2,
                .cross_section = 5.0,
            },
            .Large => .{
                .drag = 2.2,
                .cross_section = 50.0,
            },
        };
    }
};

/// This spacecraft type is the base class that takes the inputs needed to determine future orbit paths
pub const Spacecraft = struct {
    name: []const u8,
    tle: TLE,
    mass: f64,
    size: Satellite_Parameters,
    orbiting_object: Celistal_Body = constants.earth,
    orbit_predictions: std.ArrayList(StateTime),
    allocator: std.mem.Allocator,

    const Self = @This();

    pub fn create(name: []const u8, tle: TLE, mass: f64, size: Satellite_Size, orbiting_object: ?Celistal_Body, allocator: std.mem.Allocator) Self {
        return .{
            .name = name,
            .tle = tle,
            .mass = mass,
            .size = size.get_drag_and_cross_sectional(),
            .orbiting_object = orbiting_object.?,
            .orbit_predictions = std.ArrayList(StateTime).init(allocator),
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Self) void {
        self.orbit_predictions.deinit();
    }

    /// This will call the proper propagation methods based on a TLE epoch and recalculation time
    /// Most of the functions this calls are private and will need code inspection to see
    pub fn propagate(self: *Self, t0: f64, days: f64, h: f64, impulse_list: ?[]const Impulse) !void {
        const y0 = self.tle_to_state_vector();
        var t = t0;
        const tf = self.tle.first_line.epoch + days * 86400.0;
        var y = y0;
        const initial_energy = self.calculate_energy(y);
        std.log.debug("Initial energy established: {d}\n", .{initial_energy});
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
                        },
                        .Phase => {
                            const r = @sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
                            const phase_change = impulses[impulse_index].phase_change orelse 0;
                            const transfer_orbits = impulses[impulse_index].delta_v[0];
                            const delta_v_magnitude = self.calculate_phase_change(r, phase_change, transfer_orbits);
                            const velocity_magnitude = @sqrt(y[3] * y[3] + y[4] * y[4] + y[5] * y[5]);
                            const delta_v = [3]f64{
                                y[3] / velocity_magnitude * delta_v_magnitude,
                                y[4] / velocity_magnitude * delta_v_magnitude,
                                y[5] / velocity_magnitude * delta_v_magnitude,
                            };
                            y = impulse(y, delta_v);

                            const period = 2 * math.pi * @sqrt(math.pow(f64, r, 3) / self.orbiting_object.mu);
                            const second_impulse_time = t + period * transfer_orbits;
                            while (t < second_impulse_time) : (t += h) {
                                y = self.rk4(y, h);
                                try self.orbit_predictions.append(StateTime{ .time = t + h, .state = y });
                            }
                            y = impulse(y, .{ -delta_v[0], -delta_v[1], -delta_v[2] });
                            try self.orbit_predictions.append(StateTime{ .time = t, .state = y });
                            std.log.info("Phase change completed at t={d:.2}s", .{t});
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

            const energy = self.calculate_energy(y);
            const r = @sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
            try self.orbit_predictions.append(StateTime{ .time = t, .state = y });

            if (energy > 0 or std.math.isNan(energy) or r > 100000) {
                std.log.warn("Abnormal orbit detected at step {d}", .{step});
                std.log.warn("Time: {d} seconds", .{t - t0});
                std.log.warn("State vector: {any}", .{y});
                std.log.warn("Radius: {d} km", .{r});
                std.log.warn("Energy: {d}", .{energy});
                break;
            }
        }
    }

    fn rk4(self: Self, y: StateV, h: f64) StateV {
        const k1 = scalar_multiply(h, self.derivative(y));
        const k2 = scalar_multiply(h, self.derivative(vector_add(y, scalar_multiply(0.5, k1))));
        const k3 = scalar_multiply(h, self.derivative(vector_add(y, scalar_multiply(0.5, k2))));
        const k4 = scalar_multiply(h, self.derivative(vector_add(y, k3)));
        return vector_add(
            y,
            scalar_multiply(
                1.0 / 6.0,
                vector_add(
                    vector_add(k1, scalar_multiply(2, k2)),
                    vector_add(scalar_multiply(2, k3), k4),
                ),
            ),
        );
    }

    fn vector_add(a: StateV, b: StateV) StateV {
        var result: StateV = undefined;
        for (0..6) |i| {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    fn scalar_multiply(scalar: f64, vector: StateV) StateV {
        var result: StateV = undefined;
        for (0..6) |i| {
            result[i] = scalar * vector[i];
        }
        return result;
    }

    fn derivative(self: Self, state: StateV) StateV {
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

        const j2_factor = -1.5 * self.orbiting_object.j2_pertrubation * self.orbiting_object.mu * self.orbiting_object.eq_radius.? * self.orbiting_object.eq_radius.? / r5;
        ax += j2_factor * x * (5 * z * z / r2 - 1);
        ay += j2_factor * y * (5 * z * z / r2 - 1);
        az += j2_factor * z * (5 * z * z / r2 - 3);

        const altitude = r - self.orbiting_object.eq_radius.?;
        if (altitude < 1000) { // Only apply drag below 1000 km
            const v = @sqrt(vx * vx + vy * vy + vz * vz);
            const rho = self.atmospheric_density(altitude);
            const drag_factor = -0.5 * self.size.drag *
                self.size.cross_section *
                rho * v / self.mass;
            ax += drag_factor * vx;
            ay += drag_factor * vy;
            az += drag_factor * vz;
        }

        return .{ vx, vy, vz, ax, ay, az };
    }

    fn atmospheric_density(self: Self, altitude: f64) f64 {
        return self.orbiting_object.sea_level_density * math.exp(-altitude / self.orbiting_object.scale_height);
    }

    fn calculate_energy(self: Self, state: StateV) f64 {
        const r = @sqrt(state[0] * state[0] + state[1] * state[1] + state[2] * state[2]);
        const v_squared = state[3] * state[3] + state[4] * state[4] + state[5] * state[5];
        return 0.5 * v_squared - self.orbiting_object.mu / r;
    }
    fn correct_energy(self: Self, state: StateV, target_energy: f64) StateV {
        const current_energy = self.calculate_energy(state);
        std.log.debug("Current energy: {d}, Target energy: {d}", .{ current_energy, target_energy });

        if (std.math.isNan(current_energy) or std.math.isNan(target_energy)) {
            std.log.err("NaN energy detected in correct_energy", .{});
            return state;
        }

        if (current_energy >= 0) {
            std.log.err("Positive energy detected: {d}", .{current_energy});
            return state;
        }

        const energy_ratio = @sqrt(target_energy / current_energy);

        if (std.math.isNan(energy_ratio) or energy_ratio <= 0) {
            std.log.err("Invalid energy ratio: {d}", .{energy_ratio});
            return state;
        }

        var corrected_state = state;
        corrected_state[3] *= energy_ratio;
        corrected_state[4] *= energy_ratio;
        corrected_state[5] *= energy_ratio;
        return corrected_state;
    }

    /// This will calculate and return the position and velocity state vector from the TLE of the spacecraft
    pub fn tle_to_state_vector(self: Self) [6]f64 {
        const inclination = calculations.degrees_to_radians(self.tle.second_line.inclination);
        const raan = calculations.degrees_to_radians(self.tle.second_line.right_ascension);
        const eccentricity = self.tle.second_line.eccentricity;
        const argument_of_perigee = calculations.degrees_to_radians(self.tle.second_line.perigee);
        const m_anomaly = calculations.degrees_to_radians(self.tle.second_line.m_anomaly);
        const m_motion = calculations.mean_motion_to_radians_per_minute(self.tle.second_line.m_motion);

        const a = calculations.mean_motion_to_semi_major_axis(m_motion);

        const E = solve_kepler_equation(m_anomaly, eccentricity);
        const true_anomaly = 2.0 * math.atan(@sqrt((1 + eccentricity) / (1 - eccentricity)) * @tan(E / 2.0));

        const r = a * (1 - eccentricity * @cos(E));

        const x_orbital = r * @cos(true_anomaly);
        const y_orbital = r * @sin(true_anomaly);

        const p = a * (1 - eccentricity * eccentricity);
        const vx_orbital = -@sqrt(self.orbiting_object.mu / p) * @sin(E);
        const vy_orbital = @sqrt(self.orbiting_object.mu / p) * (@sqrt(1 - eccentricity * eccentricity) * @cos(E));

        var position = [3]f64{ 0, 0, 0 };
        var velocity = [3]f64{ 0, 0, 0 };

        const cos_raan = @cos(raan);
        const sin_raan = @sin(raan);
        const cos_argp = @cos(argument_of_perigee);
        const sin_argp = @sin(argument_of_perigee);
        const cos_inc = @cos(inclination);
        const sin_inc = @sin(inclination);

        position[0] = (cos_raan * cos_argp - sin_raan * sin_argp * cos_inc) * x_orbital + (-cos_raan * sin_argp - sin_raan * cos_argp * cos_inc) * y_orbital;
        position[1] = (sin_raan * cos_argp + cos_raan * sin_argp * cos_inc) * x_orbital + (-sin_raan * sin_argp + cos_raan * cos_argp * cos_inc) * y_orbital;
        position[2] = (sin_argp * sin_inc) * x_orbital + (cos_argp * sin_inc) * y_orbital;

        velocity[0] = (cos_raan * cos_argp - sin_raan * sin_argp * cos_inc) * vx_orbital + (-cos_raan * sin_argp - sin_raan * cos_argp * cos_inc) * vy_orbital;
        velocity[1] = (sin_raan * cos_argp + cos_raan * sin_argp * cos_inc) * vx_orbital + (-sin_raan * sin_argp + cos_raan * cos_argp * cos_inc) * vy_orbital;
        velocity[2] = (sin_argp * sin_inc) * vx_orbital + (cos_argp * sin_inc) * vy_orbital;

        return .{ position[0], position[1], position[2], velocity[0], velocity[1], velocity[2] };
    }

    fn solve_kepler_equation(m_anomaly: f64, eccentricity: f64) f64 {
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

    fn calculate_phase_change(self: Self, radius: f64, phase_change: f64, transfer_orbits: f64) f64 {
        const v_circular = @sqrt(self.orbiting_object.mu / radius);
        const period = 2 * math.pi * @sqrt(math.pow(f64, radius, 3) / self.orbiting_object.mu);
        const delta_t = phase_change * period / (2 * math.pi);
        const a_transfer = math.pow(f64, (period + delta_t) / (2 * math.pi) * @sqrt(self.orbiting_object.mu), 2.0 / 3.0);
        const v_transfer = @sqrt(self.orbiting_object.mu * (2 / radius - 1 / a_transfer));
        return (v_transfer - v_circular) * 2 / transfer_orbits;
    }
};

test "create spacecraft" {
    const raw_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;
    var test_tle = try TLE.parse(raw_tle, std.testing.allocator);
    defer test_tle.deinit();

    var test_sc = Spacecraft.create("dummy_sc", test_tle, 300.000, Satellite_Size.Cube, constants.earth, std.testing.allocator);
    defer test_sc.deinit();

    try test_sc.propagate(
        test_sc.tle.first_line.epoch,
        3, // days to predict
        1, // steps, i.e. predict every simulated second
        null,
    );

    for (test_sc.orbit_predictions.items) |iter| {
        const r = math.sqrt(iter.state[0] * iter.state[0] + iter.state[1] * iter.state[1] + iter.state[2] * iter.state[2]);

        try std.testing.expect(r > test_sc.orbiting_object.eq_radius.?);
    }

    // write out to file for creating test mapping data
    // const file = try std.fs.cwd().createFile("./test/files/orbit_data.csv", .{});
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

    var test_sc = Spacecraft.create("dummy_sc", test_tle, 300.000, Satellite_Size.Cube, constants.earth, std.testing.allocator);
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
        const r = math.sqrt(iter.state[0] * iter.state[0] + iter.state[1] * iter.state[1] + iter.state[2] * iter.state[2]);

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

    var test_sc = Spacecraft.create("dummy_sc", test_tle, 300.000, Satellite_Size.Cube, constants.earth, std.testing.allocator);
    defer test_sc.deinit();

    const phase_maneuver = Impulse{
        .time = 2500000.0,
        .delta_v = .{ 1.0, 0.0, 0.0 },
        .mode = .Phase,
        .phase_change = math.pi / 2.0,
    };

    const impulses = [_]Impulse{phase_maneuver};

    try test_sc.propagate(
        test_sc.tle.first_line.epoch,
        3, // days to predict
        1, // steps, i.e. predict every simulated second
        &impulses,
    );

    for (test_sc.orbit_predictions.items) |iter| {
        const r = math.sqrt(iter.state[0] * iter.state[0] + iter.state[1] * iter.state[1] + iter.state[2] * iter.state[2]);

        try std.testing.expect(r > test_sc.orbiting_object.eq_radius.?);
    }
}

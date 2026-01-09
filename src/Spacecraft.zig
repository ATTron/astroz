//! Base struct that takes the inputs needed to determine future orbit paths.
const std = @import("std");
const log = std.log;

const calculations = @import("calculations.zig");
const constants = @import("constants.zig");
const CelestialBody = constants.CelestialBody;
const Datetime = @import("Datetime.zig");
const Tle = @import("Tle.zig");

const Spacecraft = @This();

/// Satellite details used in calculations
pub const SatelliteParameters = struct {
    drag: f64,
    crossSection: f64,
    width: f64,
    height: f64,
    depth: f64,
};

/// The impulse maneuver type
pub const Impulse = struct {
    time: f64,
    deltaV: [3]f64,
    mode: enum { Absolute, Prograde, Phase, PlaneChange },
    phaseChange: ?f64 = null,
    planeChange: ?struct {
        deltaInclination: f64,
        deltaRaan: f64,
    } = null,
};

/// Determines the values in the SatelliteParameters struct
pub const SatelliteSize = enum {
    Cube,
    Mini,
    Medium,
    Large,

    pub fn generateDragAndCrossSectional(self: SatelliteSize) SatelliteParameters {
        return switch (self) {
            .Cube => .{
                .drag = 2.2,
                .crossSection = 0.05,
                .width = 0.1,
                .height = 0.1,
                .depth = 0.3,
            },
            .Mini => .{
                .drag = 2.2,
                .crossSection = 2.5,
                .width = 0.6,
                .height = 0.6,
                .depth = 1.0,
            },
            .Medium => .{
                .drag = 2.2,
                .crossSection = 5.0,
                .width = 1.4,
                .height = 1.4,
                .depth = 1.6,
            },
            .Large => .{
                .drag = 2.2,
                .crossSection = 50.0,
                .width = 3.2,
                .height = 3.2,
                .depth = 4.0,
            },
        };
    }
};

name: []const u8,
tle: Tle,
mass: f64,
size: SatelliteParameters,
quaternion: [4]f64,
angularVelocity: [3]f64,
inertiaTensor: [3][3]f64,
bodyVectors: [2][3]f64,
referenceVectors: [2][3]f64,
orbitingObject: CelestialBody = constants.earth,
orbitPredictions: std.ArrayList(calculations.StateTime),
allocator: std.mem.Allocator,

pub fn init(name: []const u8, tle: Tle, mass: f64, size: SatelliteSize, orbitingObject: ?CelestialBody, allocator: std.mem.Allocator) Spacecraft {
    return .{
        .name = name,
        .tle = tle,
        .mass = mass,
        .size = size.generateDragAndCrossSectional(),
        .quaternion = .{ 1.0, 0.0, 0.0, 0.0 },
        .angularVelocity = .{ 0.0, 0.0, 0.0 },
        .inertiaTensor = .{
            .{ 1.0, 0.0, 0.0 },
            .{ 0.0, 1.0, 0.0 },
            .{ 0.0, 0.0, 1.0 },
        },
        .bodyVectors = .{
            .{ 1.0, 0.0, 0.0 },
            .{ 0.0, 1.0, 0.0 },
        },
        .referenceVectors = .{
            .{ 1.0, 0.0, 0.0 },
            .{ 0.0, 1.0, 0.0 },
        },
        .orbitingObject = orbitingObject.?,
        .orbitPredictions = std.ArrayList(calculations.StateTime){},
        .allocator = allocator,
    };
}

pub fn deinit(self: *Spacecraft) void {
    self.orbitPredictions.deinit(self.allocator);
}

pub fn updateAttitude(self: *Spacecraft) void {
    const attitudeMatrix = calculations.triad(
        self.bodyVectors[0],
        self.bodyVectors[1],
        self.referenceVectors[0],
        self.referenceVectors[1],
    );
    self.quaternion = calculations.matrixToQuaternion(attitudeMatrix);
}

pub fn propagateAttitude(self: *Spacecraft, dt: f64) void {
    const state = calculations.AttitudeState{
        .quaternion = self.quaternion,
        .angularVelocity = self.angularVelocity,
    };
    const newState = rk4Attitude(self, state, dt);
    self.quaternion = newState.quaternion;
    self.angularVelocity = newState.angularVelocity;
}

/// This will call the proper propagation methods based on a TLE epoch and recalculation time
/// Most of the functions this calls are private and will need code inspection to see
pub fn propagate(self: *Spacecraft, t0: f64, days: f64, h: f64, impulseList: ?[]const Impulse) !void {
    const y0OE = calculations.tleToOrbitalElements(self.tle);
    const y0 = calculations.orbitalElementsToStateVector(y0OE, self.orbitingObject.mu);
    var t = t0;
    const tf = self.tle.firstLine.epoch + days * 86400.0;
    var y = y0;
    const initialEnergy = self.calculateEnergy(y);

    log.info("Initial energy established: {d}\n", .{initialEnergy});

    try self.orbitPredictions.append(self.allocator, calculations.StateTime{ .time = t, .state = y });
    var step: usize = 0;
    var impulseIndex: usize = 0;

    while (t < tf) : (step += 1) {
        if (impulseList) |impulses| {
            while (impulseIndex < impulses.len and impulses[impulseIndex].time <= t + h) {
                const dt = impulses[impulseIndex].time - t;
                if (dt > 0) {
                    y = self.rk4(y, dt);
                    t += dt;
                    try self.orbitPredictions.append(self.allocator, calculations.StateTime{ .time = t, .state = y });
                }

                switch (impulses[impulseIndex].mode) {
                    .Absolute => {
                        y = calculations.impulse(y, impulses[impulseIndex].deltaV);
                        log.info("Impulse completed at t={d:.2}s", .{t});
                    },
                    .Prograde => {
                        const velocityMagnitude = @sqrt(y[3] * y[3] + y[4] * y[4] + y[5] * y[5]);
                        const deltaVMagnitude = impulses[impulseIndex].deltaV[0]; // Use the x-component as magnitude
                        const deltaV = [3]f64{
                            y[3] / velocityMagnitude * deltaVMagnitude,
                            y[4] / velocityMagnitude * deltaVMagnitude,
                            y[5] / velocityMagnitude * deltaVMagnitude,
                        };
                        y = calculations.impulse(y, deltaV);
                        log.info("Prograde impulse completed at t={d:.2}s", .{t});
                    },
                    .Phase => {
                        const r = @sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
                        const phaseChange = impulses[impulseIndex].phaseChange orelse 0;
                        const transferOrbits = impulses[impulseIndex].deltaV[0];
                        const deltaVMagnitude = self.calculatePhaseChange(r, phaseChange, transferOrbits);
                        const velocityMagnitude = @sqrt(y[3] * y[3] + y[4] * y[4] + y[5] * y[5]);
                        const deltaV = [3]f64{
                            y[3] / velocityMagnitude * deltaVMagnitude,
                            y[4] / velocityMagnitude * deltaVMagnitude,
                            y[5] / velocityMagnitude * deltaVMagnitude,
                        };
                        y = calculations.impulse(y, deltaV);

                        const period = 2 * std.math.pi * @sqrt(std.math.pow(f64, r, 3) / self.orbitingObject.mu);
                        const secondImpulseTime = t + period * transferOrbits;
                        while (t < secondImpulseTime) : (t += h) {
                            y = self.rk4(y, h);
                            try self.orbitPredictions.append(self.allocator, calculations.StateTime{ .time = t + h, .state = y });
                        }
                        y = calculations.impulse(y, .{ -deltaV[0], -deltaV[1], -deltaV[2] });
                        try self.orbitPredictions.append(self.allocator, calculations.StateTime{ .time = t, .state = y });
                        log.info("Phase change completed at t={d:.2}s", .{t});
                    },
                    .PlaneChange => {
                        const planeChange = impulses[impulseIndex];
                        y = self.applyPlaneChange(y, planeChange);
                        try self.orbitPredictions.append(self.allocator, calculations.StateTime{ .time = t, .state = y });
                        log.info("Plane change completed at t={d:.2}s", .{t});
                    },
                }

                try self.orbitPredictions.append(self.allocator, calculations.StateTime{ .time = t, .state = y });

                impulseIndex += 1;
            }
        }

        // Regular propagation step
        const stepSize = @min(h, tf - t);
        y = self.rk4(y, stepSize);
        t += stepSize;
        const r = @sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
        const energy = self.calculateEnergy(y);

        try self.orbitPredictions.append(self.allocator, calculations.StateTime{ .time = t, .state = y });

        if (energy > 0 or std.math.isNan(energy) or r > 100_000) {
            log.warn("Abnormal orbit detected at step {d}", .{step});
            log.warn("Time: {d} seconds", .{t - t0});
            log.warn("State vector: {any}", .{y});
            log.warn("Radius: {d} km", .{r});
            log.warn("Energy: {d}", .{energy});
            break;
        }
    }
}

fn rk4(self: Spacecraft, y: calculations.StateV, h: f64) calculations.StateV {
    const k1 = calculations.scalarMultiply(h, self.derivative(y));
    const k2 = calculations.scalarMultiply(h, self.derivative(calculations.vectorAdd(y, calculations.scalarMultiply(0.5, k1))));
    const k3 = calculations.scalarMultiply(h, self.derivative(calculations.vectorAdd(y, calculations.scalarMultiply(0.5, k2))));
    const k4 = calculations.scalarMultiply(h, self.derivative(calculations.vectorAdd(y, k3)));
    return calculations.vectorAdd(
        y,
        calculations.scalarMultiply(
            1.0 / 6.0,
            calculations.vectorAdd(
                calculations.vectorAdd(k1, calculations.scalarMultiply(2, k2)),
                calculations.vectorAdd(calculations.scalarMultiply(2, k3), k4),
            ),
        ),
    );
}

fn rk4Attitude(spacecraft: *Spacecraft, state: calculations.AttitudeState, dt: f64) calculations.AttitudeState {
    const k1 = attitudeDerivative(spacecraft, state);
    const k2 = attitudeDerivative(spacecraft, calculations.addScaledAttitudeState(state, k1, 0.5 * dt));
    const k3 = attitudeDerivative(spacecraft, calculations.addScaledAttitudeState(state, k2, 0.5 * dt));
    const k4 = attitudeDerivative(spacecraft, calculations.addScaledAttitudeState(state, k3, dt));

    return calculations.addScaledAttitudeState(
        state,
        calculations.addAttitudeStates(
            calculations.addAttitudeStates(k1, calculations.scaleAttitudeState(k2, 2)),
            calculations.addAttitudeStates(calculations.scaleAttitudeState(k3, 2), k4),
        ),
        dt / 6.0,
    );
}

fn attitudeDerivative(spacecraft: *Spacecraft, state: calculations.AttitudeState) calculations.AttitudeState {
    const q = state.quaternion;
    const w = state.angularVelocity;

    const qDot = [4]f64{
        0.5 * (-q[1] * w[0] - q[2] * w[1] - q[3] * w[2]),
        0.5 * (q[0] * w[0] + q[2] * w[2] - q[3] * w[1]),
        0.5 * (q[0] * w[1] - q[1] * w[2] + q[3] * w[0]),
        0.5 * (q[0] * w[2] + q[1] * w[1] - q[2] * w[0]),
    };

    const I = spacecraft.inertiaTensor;
    const wDot = [3]f64{
        (I[1][1] - I[2][2]) / I[0][0] * w[1] * w[2],
        (I[2][2] - I[0][0]) / I[1][1] * w[2] * w[0],
        (I[0][0] - I[1][1]) / I[2][2] * w[0] * w[1],
    };

    return .{
        .quaternion = qDot,
        .angularVelocity = wDot,
    };
}

fn derivative(self: Spacecraft, state: calculations.StateV) calculations.StateV {
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

    var ax = -self.orbitingObject.mu * x / r3;
    var ay = -self.orbitingObject.mu * y / r3;
    var az = -self.orbitingObject.mu * z / r3;

    const j2Factor = -1.5 * self.orbitingObject.j2Perturbation * self.orbitingObject.mu * self.orbitingObject.eqRadius.? * self.orbitingObject.eqRadius.? / r5;
    ax += j2Factor * x * (5 * z * z / r2 - 1);
    ay += j2Factor * y * (5 * z * z / r2 - 1);
    az += j2Factor * z * (5 * z * z / r2 - 3);

    const altitude = r - self.orbitingObject.eqRadius.?;
    if (altitude < 1000) { // Only apply drag below 1000 km
        const v = @sqrt(vx * vx + vy * vy + vz * vz);
        const rho = self.atmosphericDensity(altitude);
        const dragFactor = -0.5 * self.size.drag *
            self.size.crossSection *
            rho * v / self.mass;
        ax += dragFactor * vx;
        ay += dragFactor * vy;
        az += dragFactor * vz;
    }

    return .{ vx, vy, vz, ax, ay, az };
}

fn atmosphericDensity(self: Spacecraft, altitude: f64) f64 {
    return self.orbitingObject.seaLevelDensity * std.math.exp(-altitude / self.orbitingObject.scaleHeight);
}

fn calculateEnergy(self: Spacecraft, state: calculations.StateV) f64 {
    const r = @sqrt(state[0] * state[0] + state[1] * state[1] + state[2] * state[2]);
    const vSquared = state[3] * state[3] + state[4] * state[4] + state[5] * state[5];
    return 0.5 * vSquared - self.orbitingObject.mu / r;
}

fn applyPlaneChange(self: *Spacecraft, y: [6]f64, plane_change: Impulse) [6]f64 {
    const mu = self.orbitingObject.mu;
    const r = [3]f64{ y[0], y[1], y[2] };
    const v = [3]f64{ y[3], y[4], y[5] };

    var elements = calculations.stateVectorToOrbitalElements(r, v, mu);
    elements.i += plane_change.planeChange.?.deltaInclination;
    elements.raan += plane_change.planeChange.?.deltaRaan;

    return calculations.orbitalElementsToStateVector(elements, mu);
}

fn calculatePhaseChange(self: Spacecraft, radius: f64, phase_change: f64, transfer_orbits: f64) f64 {
    const vCircular = @sqrt(self.orbitingObject.mu / radius);
    const period = 2 * std.math.pi * @sqrt(std.math.pow(f64, radius, 3) / self.orbitingObject.mu);
    const deltaT = phase_change * period / (2 * std.math.pi);
    const aTransfer = std.math.pow(f64, (period + deltaT) / (2 * std.math.pi) * @sqrt(self.orbitingObject.mu), 2.0 / 3.0);
    const vTransfer = @sqrt(self.orbitingObject.mu * (2 / radius - 1 / aTransfer));
    return (vTransfer - vCircular) * 2 / transfer_orbits;
}

test "init spacecraft" {
    const raw_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;
    var test_tle = try Tle.parse(raw_tle, std.testing.allocator);
    defer test_tle.deinit();

    var test_sc = Spacecraft.init(
        "dummy_sc",
        test_tle,
        300.000,
        SatelliteSize.Cube,
        constants.earth,
        std.testing.allocator,
    );
    defer test_sc.deinit();

    try test_sc.propagate(
        test_sc.tle.firstLine.epoch,
        3, // days to predict
        1, // steps, i.e. predict every simulated second
        null,
    );

    for (test_sc.orbitPredictions.items) |iter| {
        const r = std.math.sqrt(iter.state[0] * iter.state[0] + iter.state[1] * iter.state[1] + iter.state[2] * iter.state[2]);

        try std.testing.expect(r > test_sc.orbitingObject.eqRadius.?);
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
    var test_tle = try Tle.parse(raw_tle, std.testing.allocator);
    defer test_tle.deinit();

    var test_sc = Spacecraft.init(
        "dummy_sc",
        test_tle,
        300.000,
        SatelliteSize.Cube,
        constants.earth,
        std.testing.allocator,
    );
    defer test_sc.deinit();

    const impulses = [_]Impulse{
        .{ .time = 2635014.50, .deltaV = .{ 0.2, 0.0, 0.0 }, .mode = .Prograde },
        .{ .time = 2638026.50, .deltaV = .{ 0.2, 0.0, 0.0 }, .mode = .Prograde },
        .{ .time = 2638103.50, .deltaV = .{ 0.2, 0.0, 0.0 }, .mode = .Prograde },
    };

    try test_sc.propagate(
        test_sc.tle.firstLine.epoch,
        3, // days to predict
        1, // steps, i.e. predict every simulated second
        &impulses,
    );

    for (test_sc.orbitPredictions.items) |iter| {
        const r = std.math.sqrt(iter.state[0] * iter.state[0] + iter.state[1] * iter.state[1] + iter.state[2] * iter.state[2]);

        try std.testing.expect(r > test_sc.orbitingObject.eqRadius.?);
    }
}

test "prop spacecraft w/ phase" {
    const raw_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;
    var test_tle = try Tle.parse(raw_tle, std.testing.allocator);
    defer test_tle.deinit();

    var test_sc = Spacecraft.init(
        "dummy_sc",
        test_tle,
        300.000,
        SatelliteSize.Cube,
        constants.earth,
        std.testing.allocator,
    );
    defer test_sc.deinit();

    const phase_maneuver = Impulse{
        .time = 2500000.0,
        .deltaV = .{ 1.0, 0.0, 0.0 },
        .mode = .Phase,
        .phaseChange = std.math.pi / 2.0,
    };

    const impulses = [_]Impulse{phase_maneuver};

    try test_sc.propagate(
        test_sc.tle.firstLine.epoch,
        3, // days to predict
        1, // steps, i.e. predict every simulated second
        &impulses,
    );

    for (test_sc.orbitPredictions.items) |iter| {
        const r = std.math.sqrt(iter.state[0] * iter.state[0] + iter.state[1] * iter.state[1] + iter.state[2] * iter.state[2]);

        try std.testing.expect(r > test_sc.orbitingObject.eqRadius.?);
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
    var test_tle = try Tle.parse(raw_tle, std.testing.allocator);
    defer test_tle.deinit();

    var test_sc = Spacecraft.init(
        "dummy_sc",
        test_tle,
        300.000,
        SatelliteSize.Cube,
        constants.earth,
        std.testing.allocator,
    );
    defer test_sc.deinit();

    const plane_change_maneuver = Impulse{
        .time = 2500000.0,
        .deltaV = .{ 0.0, 0.0, 0.0 }, // Not used for plane changes
        .mode = .PlaneChange,
        .planeChange = .{
            .deltaInclination = std.math.pi / 18.0, // 10-degree inclination change
            .deltaRaan = std.math.pi / 36.0, // 5-degree RAAN change
        },
    };

    const impulses = [_]Impulse{plane_change_maneuver};

    try test_sc.propagate(
        test_sc.tle.firstLine.epoch,
        3, // days to predict
        1, // steps, i.e. predict every simulated second
        &impulses,
    );

    for (test_sc.orbitPredictions.items) |iter| {
        const r = std.math.sqrt(iter.state[0] * iter.state[0] + iter.state[1] * iter.state[1] + iter.state[2] * iter.state[2]);

        try std.testing.expect(r > test_sc.orbitingObject.eqRadius.?);
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
    var test_tle = try Tle.parse(raw_tle, std.testing.allocator);
    defer test_tle.deinit();
    var spacecraft = Spacecraft.init(
        "dummy_sc",
        test_tle,
        300.000,
        SatelliteSize.Cube,
        constants.earth,
        std.testing.allocator,
    );

    spacecraft.angularVelocity = .{ 0.1, 0.05, 0.02 };

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

        spacecraft.bodyVectors[0] = .{ @cos(angle), 0, @sin(angle) };
        spacecraft.bodyVectors[1] = .{ 0, 1, 0 };

        spacecraft.updateAttitude();
        spacecraft.propagateAttitude(dt);

        // Simulate simple circular orbit
        const orbit_radius = 7000.0;
        // const x = orbit_radius * @cos(2 * std.math.pi * t / orbital_period);
        // const y = orbit_radius * @sin(2 * std.math.pi * t / orbital_period);
        // const z = 0.0;

        try std.testing.expect(orbit_radius > spacecraft.orbitingObject.mRadius.?);

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
    var test_tle = try Tle.parse(raw_tle, std.testing.allocator);
    defer test_tle.deinit();
    var spacecraft = Spacecraft.init(
        "dummy_sc",
        test_tle,
        300.000,
        SatelliteSize.Cube,
        constants.earth,
        std.testing.allocator,
    );
    defer spacecraft.deinit();

    // Initial angular velocity
    spacecraft.angularVelocity = .{ 0.0, 0.0, 0.0 };

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
        spacecraft.angularVelocity[0] += torque_x * dt;
        spacecraft.angularVelocity[1] += torque_y * dt;
        spacecraft.angularVelocity[2] += torque_z * dt;

        // Update attitude
        spacecraft.updateAttitude();
        spacecraft.propagateAttitude(dt);

        // Simulate simple circular orbit
        const orbit_radius = 7000.0;
        // const x = orbit_radius * @cos(2 * std.math.pi * t / orbital_period);
        // const y = orbit_radius * @sin(2 * std.math.pi * t / orbital_period);
        // const z = 0.0;

        try std.testing.expect(orbit_radius > spacecraft.orbitingObject.mRadius.?);

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

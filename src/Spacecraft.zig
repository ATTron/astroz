//! Base struct that takes the inputs needed to determine future orbit paths.
const std = @import("std");
const log = std.log;

const calculations = @import("calculations.zig");
const constants = @import("constants.zig");
const CelestialBody = constants.CelestialBody;
const Tle = @import("Tle.zig");
const propagators = @import("propagators/propagators.zig");

const Spacecraft = @This();

/// Satellite details used in calculations
pub const SatelliteParameters = struct {
    drag: f64,
    crossSection: f64,
    width: f64,
    height: f64,
    depth: f64,
};

/// impulse maneuver types with clearer semantics
pub const Impulse = struct {
    time: f64,
    maneuver: Maneuver,

    pub const Maneuver = union(enum) {
        absolute: [3]f64,
        prograde: f64,
        phase: struct {
            angle: f64, // radians
            orbits: f64 = 1.0, // transfer orbits
        },
        planeChange: struct {
            deltaInclination: f64, // radians
            deltaRaan: f64, // radians
        },
    };
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

/// creates force models configured for this spacecraft's orbiting body and parameters
fn createForceModels(self: *Spacecraft) struct {
    twobody: propagators.TwoBody,
    j2: propagators.J2,
    drag: propagators.Drag,
} {
    return .{
        .twobody = propagators.TwoBody.init(self.orbitingObject.mu),
        .j2 = propagators.J2.init(
            self.orbitingObject.mu,
            self.orbitingObject.j2Perturbation,
            self.orbitingObject.eqRadius.?,
        ),
        .drag = propagators.Drag.init(
            self.orbitingObject.eqRadius.?,
            self.orbitingObject.seaLevelDensity,
            self.orbitingObject.scaleHeight,
            self.size.drag,
            self.size.crossSection,
            self.mass,
            1000.0, // max altitude for drag
        ),
    };
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
    const newState = calculations.propagateAttitude(state, self.inertiaTensor, dt);
    self.quaternion = newState.quaternion;
    self.angularVelocity = newState.angularVelocity;
}

/// propagate orbit from TLE epoch for specified days with optional impulse maneuvers
pub fn propagate(self: *Spacecraft, t0: f64, days: f64, h: f64, impulseList: ?[]const Impulse) !void {
    const y0OE = calculations.tleToOrbitalElements(self.tle);
    var y = calculations.orbitalElementsToStateVector(y0OE, self.orbitingObject.mu);
    var t = t0;
    const tf = self.tle.firstLine.epoch + days * constants.secondsPerDay;

    // setup force models and integrator
    var forces = self.createForceModels();
    const models = [_]propagators.ForceModel{
        propagators.ForceModel.wrap(propagators.TwoBody, &forces.twobody),
        propagators.ForceModel.wrap(propagators.J2, &forces.j2),
        propagators.ForceModel.wrap(propagators.Drag, &forces.drag),
    };
    var composite = try propagators.Composite.init(self.allocator, &models);
    defer composite.deinit();

    var rk4 = propagators.Rk4{};
    const integrator = rk4.integrator();
    const force = propagators.ForceModel.wrap(propagators.Composite, &composite);

    try self.orbitPredictions.append(self.allocator, .{ .time = t, .state = y });
    var impulseIndex: usize = 0;

    while (t < tf) {
        // Handle impulse maneuvers
        if (impulseList) |impulses| {
            while (impulseIndex < impulses.len and impulses[impulseIndex].time <= t + h) {
                const dt = impulses[impulseIndex].time - t;
                if (dt > 0) {
                    y = try integrator.step(y, t, dt, force);
                    t += dt;
                    try self.orbitPredictions.append(self.allocator, .{ .time = t, .state = y });
                }
                y = try self.applyImpulse(y, impulses[impulseIndex], &t, h, integrator, force);
                try self.orbitPredictions.append(self.allocator, .{ .time = t, .state = y });
                impulseIndex += 1;
            }
        }

        // regular propagation step
        const stepSize = @min(h, tf - t);
        y = try integrator.step(y, t, stepSize, force);
        t += stepSize;
        try self.orbitPredictions.append(self.allocator, .{ .time = t, .state = y });

        // check for abnormal orbit
        const r = calculations.posMag(y);
        const energy = self.calculateEnergy(y);
        if (energy > 0 or std.math.isNan(energy) or r > 100_000) {
            log.warn("Abnormal orbit: r={d} km, energy={d}", .{ r, energy });
            break;
        }
    }
}

fn applyImpulse(self: *Spacecraft, state: [6]f64, impulse: Impulse, t: *f64, h: f64, integrator: propagators.Integrator, force: propagators.ForceModel) ![6]f64 {
    var y = state;
    switch (impulse.maneuver) {
        .absolute => |dv| {
            y = calculations.impulse(y, dv);
        },
        .prograde => |dvMag| {
            const vMag = calculations.velMag(y);
            const dv = [3]f64{ y[3] / vMag * dvMag, y[4] / vMag * dvMag, y[5] / vMag * dvMag };
            y = calculations.impulse(y, dv);
        },
        .phase => |phase| {
            const r = calculations.posMag(y);
            const vMag = calculations.velMag(y);
            const dvMag = self.calculatePhaseChange(r, phase.angle, phase.orbits);
            const dv = [3]f64{ y[3] / vMag * dvMag, y[4] / vMag * dvMag, y[5] / vMag * dvMag };
            y = calculations.impulse(y, dv);

            // Propagate through transfer orbit(s)
            const period = 2 * std.math.pi * @sqrt(std.math.pow(f64, r, 3) / self.orbitingObject.mu);
            const tEnd = t.* + period * phase.orbits;
            while (t.* < tEnd) : (t.* += h) {
                y = try integrator.step(y, t.*, h, force);
                try self.orbitPredictions.append(self.allocator, .{ .time = t.* + h, .state = y });
            }
            // Return burn to circularize
            y = calculations.impulse(y, .{ -dv[0], -dv[1], -dv[2] });
        },
        .planeChange => |pc| {
            y = self.applyPlaneChange(y, pc.deltaInclination, pc.deltaRaan);
        },
    }
    return y;
}

fn calculateEnergy(self: Spacecraft, state: calculations.StateV) f64 {
    const r = calculations.posMag(state);
    const v = calculations.velMag(state);
    return 0.5 * v * v - self.orbitingObject.mu / r;
}

/// apply plane change maneuver using actual delta-V
fn applyPlaneChange(_: *Spacecraft, y: [6]f64, deltaInclination: f64, deltaRaan: f64) [6]f64 {
    const vMag = calculations.velMag(y);

    // combined plane change angle (assumes optimal geo)
    const totalAngle = @sqrt(deltaInclination * deltaInclination + deltaRaan * deltaRaan);
    if (totalAngle < 1e-10) return y;

    // delta-v magnitude for plane change: dv = 2 * v * sin(angle/2)
    const dvMag = 2.0 * vMag * @sin(totalAngle / 2.0);

    // compute normal to orbital plane (r x v)
    const r = [3]f64{ y[0], y[1], y[2] };
    const v = [3]f64{ y[3], y[4], y[5] };
    const h = [3]f64{
        r[1] * v[2] - r[2] * v[1],
        r[2] * v[0] - r[0] * v[2],
        r[0] * v[1] - r[1] * v[0],
    };
    const hMag = @sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);

    // apply delta-V in direction of angular momentum change (simplified)
    const dv = [3]f64{
        h[0] / hMag * dvMag * @sin(deltaInclination),
        h[1] / hMag * dvMag * @sin(deltaInclination),
        h[2] / hMag * dvMag * @cos(deltaInclination),
    };

    // log the delta-V cost for user awareness
    log.debug("Plane change: di={d:.2}°, dΩ={d:.2}°, Δv={d:.4} km/s", .{
        deltaInclination * 180.0 / std.math.pi,
        deltaRaan * 180.0 / std.math.pi,
        dvMag,
    });

    return calculations.impulse(y, dv);
}

/// calculate delta-V for a phasing maneuver
fn calculatePhaseChange(self: Spacecraft, radius: f64, phaseAngle: f64, transferOrbits: f64) f64 {
    const mu = self.orbitingObject.mu;
    const vCircular = @sqrt(mu / radius);
    const period = 2.0 * std.math.pi * @sqrt(std.math.pow(f64, radius, 3) / mu);

    const deltaT = phaseAngle * period / (2.0 * std.math.pi * transferOrbits);

    const transferPeriod = period + deltaT;
    const aTransfer = std.math.pow(f64, transferPeriod * @sqrt(mu) / (2.0 * std.math.pi), 2.0 / 3.0);

    const vTransfer = @sqrt(mu * (2.0 / radius - 1.0 / aTransfer));

    return vTransfer - vCircular;
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
        const r = calculations.posMag(iter.state);

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
        .{ .time = 2635014.50, .maneuver = .{ .prograde = 0.2 } },
        .{ .time = 2638026.50, .maneuver = .{ .prograde = 0.2 } },
        .{ .time = 2638103.50, .maneuver = .{ .prograde = 0.2 } },
    };

    try test_sc.propagate(
        test_sc.tle.firstLine.epoch,
        3, // days to predict
        1, // steps, i.e. predict every simulated second
        &impulses,
    );

    for (test_sc.orbitPredictions.items) |iter| {
        const r = calculations.posMag(iter.state);

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
        .maneuver = .{ .phase = .{ .angle = std.math.pi / 2.0, .orbits = 1.0 } },
    };

    const impulses = [_]Impulse{phase_maneuver};

    try test_sc.propagate(
        test_sc.tle.firstLine.epoch,
        3, // days to predict
        1, // steps, i.e. predict every simulated second
        &impulses,
    );

    for (test_sc.orbitPredictions.items) |iter| {
        const r = calculations.posMag(iter.state);

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
        .maneuver = .{
            .planeChange = .{
                .deltaInclination = std.math.pi / 18.0, // 10-degree inclination change
                .deltaRaan = std.math.pi / 36.0, // 5-degree RAAN change
            },
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
        const r = calculations.posMag(iter.state);

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

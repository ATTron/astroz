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
    time: f64, // seconds after the start of propagation
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
        .orbitPredictions = std.ArrayList(calculations.StateTime).empty,
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

/// propagate from the TLE state starting at time t0 (J2000 seconds) for the given days,
/// replacing any previous predictions. Impulse times are seconds after t0, in ascending order.
pub fn propagate(self: *Spacecraft, t0: f64, days: f64, h: f64, impulseList: ?[]const Impulse) !void {
    const impulses = impulseList orelse &.{};
    const duration = days * constants.secondsPerDay;
    var prevTime: f64 = 0;
    for (impulses) |impulse| {
        if (std.math.isNan(impulse.time) or impulse.time < 0 or impulse.time > duration) return error.ImpulseOutOfRange;
        if (impulse.time < prevTime) return error.ImpulsesNotSorted;
        prevTime = impulse.time;
    }

    const y0OE = calculations.tleToOrbitalElements(self.tle);
    var y = calculations.orbitalElementsToStateVector(y0OE, self.orbitingObject.mu);
    var t = t0;
    const tf = t0 + duration;

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

    self.orbitPredictions.clearRetainingCapacity();
    try self.orbitPredictions.append(self.allocator, .{ .time = t, .state = y });
    var next: usize = 0;

    while (t < tf) {
        // coast to any burn due within this step, then apply it
        while (next < impulses.len and t0 + impulses[next].time <= t + h) : (next += 1) {
            const burnAt = t0 + impulses[next].time;
            // only happens when a phase transfer carried us past this burn
            if (burnAt < t) return error.ImpulseDuringTransfer;
            if (burnAt > t) {
                y = try integrator.step(y, t, burnAt - t, force);
                t = burnAt;
                try self.orbitPredictions.append(self.allocator, .{ .time = t, .state = y });
            }
            y = try self.applyImpulse(y, impulses[next], &t, h, tf, integrator, force);
            try self.orbitPredictions.append(self.allocator, .{ .time = t, .state = y });
        }
        if (t >= tf) break;

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

fn applyImpulse(self: *Spacecraft, state: [6]f64, impulse: Impulse, t: *f64, h: f64, tf: f64, integrator: propagators.Integrator, force: propagators.ForceModel) ![6]f64 {
    var y = state;
    switch (impulse.maneuver) {
        .absolute => |dv| {
            y = calculations.impulse(y, dv);
        },
        .prograde => |dvMag| {
            const dv = progradeVec(y, dvMag);
            y = calculations.impulse(y, dv);
        },
        .phase => |phase| {
            const r = calculations.posMag(y);
            const dvMag = self.calculatePhaseChange(r, phase.angle, phase.orbits);
            const dv = progradeVec(y, dvMag);
            y = calculations.impulse(y, dv);

            // Propagate through transfer orbit(s)
            const period = 2 * std.math.pi * @sqrt(std.math.pow(f64, r, 3) / self.orbitingObject.mu);
            const tEnd = t.* + period * phase.orbits;
            if (tEnd > tf) return error.ImpulseOutOfRange;
            while (t.* < tEnd) {
                const step = @min(h, tEnd - t.*);
                y = try integrator.step(y, t.*, step, force);
                t.* += step;
                try self.orbitPredictions.append(self.allocator, .{ .time = t.*, .state = y });
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

fn progradeVec(y: [6]f64, dvMag: f64) [3]f64 {
    const vMag = calculations.velMag(y);
    return .{ y[3] / vMag * dvMag, y[4] / vMag * dvMag, y[5] / vMag * dvMag };
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

const testTle =
    \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
    \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
;

fn testSpacecraft(tle: Tle) Spacecraft {
    return Spacecraft.init("test_sc", tle, 300.0, SatelliteSize.Cube, constants.earth, std.testing.allocator);
}

test "propagate" {
    var tle = try Tle.parse(testTle, std.testing.allocator);
    defer tle.deinit();
    const epoch = tle.epoch;
    const hour = 3600.0;
    const days = 3;

    var baseline = testSpacecraft(tle);
    defer baseline.deinit();
    try baseline.propagate(epoch, days, 1, null);
    const base = baseline.orbitPredictions.items;

    const prograde = [_]Impulse{
        .{ .time = 48 * hour, .maneuver = .{ .prograde = 0.2 } },
        .{ .time = 49 * hour, .maneuver = .{ .prograde = 0.2 } },
        .{ .time = 50 * hour, .maneuver = .{ .prograde = 0.2 } },
    };
    const phase = [_]Impulse{
        .{ .time = 11 * hour, .maneuver = .{ .phase = .{ .angle = std.math.pi / 2.0, .orbits = 1.0 } } },
    };
    const planeChange = [_]Impulse{
        // 10 deg inclination, 5 deg RAAN
        .{ .time = 11 * hour, .maneuver = .{ .planeChange = .{ .deltaInclination = std.math.pi / 18.0, .deltaRaan = std.math.pi / 36.0 } } },
    };

    for ([_][]const Impulse{ &.{}, &prograde, &phase, &planeChange }) |impulses| {
        var sc = testSpacecraft(tle);
        defer sc.deinit();
        try sc.propagate(epoch, days, 1, impulses);
        const points = sc.orbitPredictions.items;

        // stays above the surface, time never runs backwards, ends exactly at the end of the run
        for (points[1..], points[0 .. points.len - 1]) |p, prev| {
            try std.testing.expect(calculations.posMag(p.state) > sc.orbitingObject.eqRadius.?);
            try std.testing.expect(p.time >= prev.time);
        }
        try std.testing.expectEqual(epoch + days * constants.secondsPerDay, points[points.len - 1].time);

        if (impulses.len == 0) continue;

        // matches the unperturbed orbit up to the first burn, then diverges
        var i: usize = 1;
        while (points[i].time < epoch + impulses[0].time) : (i += 1) {
            try std.testing.expectEqual(base[i].state, points[i].state);
        }
        try std.testing.expect(i > 1);
        try std.testing.expect(!std.meta.eql(base[base.len - 1].state, points[points.len - 1].state));
    }

    // a new run replaces the previous predictions instead of appending to them
    try baseline.propagate(epoch, 1.0 / 24.0, 1, null);
    try std.testing.expect(baseline.orbitPredictions.items.len < base.len);

    // rejected schedules
    const burn = Impulse.Maneuver{ .prograde = 0.1 };
    const quarterTurn = Impulse.Maneuver{ .phase = .{ .angle = std.math.pi / 2.0 } };
    const cases = [_]struct { anyerror, []const Impulse }{
        .{ error.ImpulseOutOfRange, &.{.{ .time = -1, .maneuver = burn }} },
        .{ error.ImpulseOutOfRange, &.{.{ .time = std.math.nan(f64), .maneuver = burn }} },
        .{ error.ImpulseOutOfRange, &.{.{ .time = 7 * hour, .maneuver = burn }} },
        .{ error.ImpulsesNotSorted, &.{ .{ .time = 2 * hour, .maneuver = burn }, .{ .time = hour, .maneuver = burn } } },
        // transfer orbit (~1.5 h) would run past the end of the 6 h run
        .{ error.ImpulseOutOfRange, &.{.{ .time = 5.5 * hour, .maneuver = quarterTurn }} },
        // second burn lands inside the first burn's transfer orbit
        .{ error.ImpulseDuringTransfer, &.{ .{ .time = hour, .maneuver = quarterTurn }, .{ .time = 1.5 * hour, .maneuver = burn } } },
    };
    for (cases) |c| {
        try std.testing.expectError(c[0], baseline.propagate(epoch, 0.25, 1, c[1]));
    }
}

test "attitude update and propagation" {
    var tle = try Tle.parse(testTle, std.testing.allocator);
    defer tle.deinit();
    var sc = testSpacecraft(tle);
    defer sc.deinit();

    // one 90 min orbit of spin-up under a varying torque; numerical accuracy is
    // covered by the propagateAttitude test, this only checks the wiring
    const orbitalPeriod = 90 * 60.0;
    const dt = 60.0;
    sc.angularVelocity = .{ 0, 0, 0 };
    var t: f64 = 0;
    while (t < orbitalPeriod) : (t += dt) {
        sc.angularVelocity[0] += 0.001 * @sin(2 * std.math.pi * t / orbitalPeriod) * dt;
        sc.angularVelocity[2] += 0.0002 * @cos(2 * std.math.pi * t / orbitalPeriod) * dt;
        sc.updateAttitude();
        sc.propagateAttitude(dt);
    }
    const q = sc.quaternion;
    try std.testing.expectApproxEqAbs(1.0, @sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]), 1e-9);
}

//! Propagator orchestrates integration over a time span
const std = @import("std");
const Integrator = @import("Integrator.zig").Integrator;
const ForceModel = @import("ForceModel.zig").ForceModel;

pub const StateTime = struct {
    time: f64, // seconds from epoch
    state: [6]f64, // [x, y, z, vx, vy, vz] in km and km/s
};

pub const Propagator = struct {
    integrator: Integrator,
    force: ?ForceModel,
    allocator: std.mem.Allocator,

    pub fn init(allocator: std.mem.Allocator, integrator: Integrator, force: ?ForceModel) Propagator {
        return .{
            .integrator = integrator,
            .force = force,
            .allocator = allocator,
        };
    }

    pub fn propagate(
        self: *Propagator,
        initial: [6]f64,
        t0: f64,
        duration: f64,
        dt: f64,
    ) !std.ArrayList(StateTime) {
        var trajectory = std.ArrayList(StateTime){};
        errdefer trajectory.deinit(self.allocator);

        var state = initial;
        var t = t0;
        const t_end = t0 + duration;

        // Store initial state
        try trajectory.append(self.allocator, .{ .time = t, .state = state });

        // Step through time
        while (t < t_end) {
            const step_dt = @min(dt, t_end - t);
            state = try self.integrator.step(state, t, step_dt, self.force);
            t += step_dt;
            try trajectory.append(self.allocator, .{ .time = t, .state = state });
        }

        return trajectory;
    }

    pub fn propagateTo(
        self: *Propagator,
        initial: [6]f64,
        t0: f64,
        t_final: f64,
        dt: f64,
    ) ![6]f64 {
        var state = initial;
        var t = t0;

        while (t < t_final) {
            const step_dt = @min(dt, t_final - t);
            state = try self.integrator.step(state, t, step_dt, self.force);
            t += step_dt;
        }

        return state;
    }
};

test "propagator compiles" {
    // Just verify the types compile correctly
    _ = Propagator;
    _ = StateTime;
}

test "rk4 twobody orbit" {
    const Rk4 = @import("Integrator.zig").Rk4;
    const TwoBody = @import("ForceModel.zig").TwoBody;

    const mu = 398600.4418; // Earth
    var twobody = TwoBody.init(mu);
    var rk4 = Rk4{};

    var propagator = Propagator.init(
        std.testing.allocator,
        rk4.integrator(),
        ForceModel.wrap(TwoBody, &twobody),
    );

    // Circular orbit at ~7000 km radius
    // v_circular = sqrt(mu/r) ~= 7.55 km/s
    const r0 = 7000.0;
    const v0 = @sqrt(mu / r0);
    const initial = [6]f64{ r0, 0.0, 0.0, 0.0, v0, 0.0 };

    // Propagate for one orbit period: T = 2*pi*sqrt(a^3/mu)
    const period = 2.0 * std.math.pi * @sqrt(r0 * r0 * r0 / mu);

    var trajectory = try propagator.propagate(initial, 0, period, 10.0);
    defer trajectory.deinit(std.testing.allocator);

    // Should return to approximately the same position after one orbit
    const final = trajectory.items[trajectory.items.len - 1].state;
    const final_r = @sqrt(final[0] * final[0] + final[1] * final[1] + final[2] * final[2]);

    // Radius should be conserved (within numerical tolerance)
    try std.testing.expectApproxEqRel(r0, final_r, 1e-4);

    // Should be back near starting position
    try std.testing.expectApproxEqAbs(r0, final[0], 50.0); // within 50 km
    try std.testing.expectApproxEqAbs(0.0, final[1], 50.0);
}

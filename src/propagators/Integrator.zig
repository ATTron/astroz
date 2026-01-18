//! Integrators for orbit propagation
const std = @import("std");
const Sgp4 = @import("../Sgp4.zig");
const constants = @import("../constants.zig");
const Tle = @import("../Tle.zig");
const ForceModel = @import("ForceModel.zig").ForceModel;

pub const Integrator = struct {
    ptr: *anyopaque,
    vtable: *const VTable,

    pub const VTable = struct {
        step: *const fn (ptr: *anyopaque, state: [6]f64, t: f64, dt: f64, force: ?ForceModel) anyerror![6]f64,
    };

    pub fn step(self: Integrator, state: [6]f64, t: f64, dt: f64, force: ?ForceModel) ![6]f64 {
        return self.vtable.step(self.ptr, state, t, dt, force);
    }
};

pub const Rk4 = struct {
    pub fn integrator(self: *Rk4) Integrator {
        return .{ .ptr = self, .vtable = &vtable };
    }

    const vtable = Integrator.VTable{ .step = step };

    fn step(ptr: *anyopaque, state: [6]f64, t: f64, dt: f64, force: ?ForceModel) anyerror![6]f64 {
        _ = ptr;
        const fm = force orelse return error.ForceModelRequired;

        const k1 = derivative(state, t, fm);
        const k2 = derivative(addScaled(state, k1, 0.5 * dt), t + 0.5 * dt, fm);
        const k3 = derivative(addScaled(state, k2, 0.5 * dt), t + 0.5 * dt, fm);
        const k4 = derivative(addScaled(state, k3, dt), t + dt, fm);

        var result: [6]f64 = undefined;
        for (0..6) |i| {
            result[i] = state[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        return result;
    }

    fn derivative(state: [6]f64, t: f64, force: ForceModel) [6]f64 {
        const accel = force.acceleration(state, t);
        return .{ state[3], state[4], state[5], accel[0], accel[1], accel[2] };
    }

    fn addScaled(state: [6]f64, delta: [6]f64, scale: f64) [6]f64 {
        var result: [6]f64 = undefined;
        for (0..6) |i| {
            result[i] = state[i] + scale * delta[i];
        }
        return result;
    }
};

pub const Sgp4Integrator = struct {
    sgp4: Sgp4,

    pub const Error = Sgp4.Error;

    pub fn init(tle: Tle, grav: constants.Sgp4GravityModel) Error!Sgp4Integrator {
        return .{ .sgp4 = try Sgp4.init(tle, grav) };
    }

    pub fn integrator(self: *Sgp4Integrator) Integrator {
        return .{ .ptr = self, .vtable = &vtable };
    }

    pub fn propagate(self: *const Sgp4Integrator, tsince_minutes: f64) Error![2][3]f64 {
        return self.sgp4.propagate(tsince_minutes);
    }

    const vtable = Integrator.VTable{ .step = step };

    fn step(ptr: *anyopaque, state: [6]f64, t: f64, dt: f64, force: ?ForceModel) anyerror![6]f64 {
        _ = state;
        _ = force;
        const self: *Sgp4Integrator = @ptrCast(@alignCast(ptr));
        const tsince_minutes = (t + dt) / 60.0;
        const result = try self.sgp4.propagate(tsince_minutes);
        return .{ result[0][0], result[0][1], result[0][2], result[1][0], result[1][1], result[1][2] };
    }
};

test "rk4 with twobody" {
    const TwoBody = @import("ForceModel.zig").TwoBody;
    var tb = TwoBody.init(398600.4418);
    var rk4 = Rk4{};
    const state = [6]f64{ 7000, 0, 0, 0, 7.5, 0 };
    const result = try rk4.integrator().step(state, 0, 60, tb.forceModel());
    try std.testing.expect(result[0] != state[0]);
}

test "sgp4 integrator" {
    const test_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;
    var tle = try Tle.parse(test_tle, std.testing.allocator);
    defer tle.deinit();

    var sgp4_int = try Sgp4Integrator.init(tle, constants.wgs84);
    const result = try sgp4_int.integrator().step([6]f64{ 0, 0, 0, 0, 0, 0 }, 0, 60, null);
    const r_mag = @sqrt(result[0] * result[0] + result[1] * result[1] + result[2] * result[2]);
    try std.testing.expect(r_mag > 6000.0 and r_mag < 8000.0);
}

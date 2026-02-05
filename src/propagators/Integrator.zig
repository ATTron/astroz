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

/// Dormand-Prince 8(7) adaptive step-size integrator
/// 13-stage embedded RK method with 8th order solution and 7th order error estimate
pub const DormandPrince87 = struct {
    rtol: f64 = 1e-9, // relative tolerance
    atol: f64 = 1e-12, // absolute tolerance
    hCurrent: f64 = 60.0, // current/suggested step size
    hMin: f64 = 0.001, // minimum step size (seconds)
    hMax: f64 = 3600.0, // maximum step size (seconds)
    safety: f64 = 0.9, // safety factor for step size control
    maxSubsteps: u32 = 10000, // prevent infinite loops

    // Dormand-Prince 8(7) Butcher tableau (Prince & Dormand 1981)
    // c_i coefficients (time fractions for each stage)
    const c = [13]f64{
        0.0,
        1.0 / 18.0,
        1.0 / 12.0,
        1.0 / 8.0,
        5.0 / 16.0,
        3.0 / 8.0,
        59.0 / 400.0,
        93.0 / 200.0,
        5490023248.0 / 9719169821.0,
        13.0 / 20.0,
        1201146811.0 / 1299019798.0,
        1.0,
        1.0,
    };

    // a_ij coefficients (stage coupling matrix)
    const a = [13][12]f64{
        .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        .{ 1.0 / 18.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        .{ 1.0 / 48.0, 1.0 / 16.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        .{ 1.0 / 32.0, 0, 3.0 / 32.0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        .{ 5.0 / 16.0, 0, -75.0 / 64.0, 75.0 / 64.0, 0, 0, 0, 0, 0, 0, 0, 0 },
        .{ 3.0 / 80.0, 0, 0, 3.0 / 16.0, 3.0 / 20.0, 0, 0, 0, 0, 0, 0, 0 },
        .{ 29443841.0 / 614563906.0, 0, 0, 77736538.0 / 692538347.0, -28693883.0 / 1125000000.0, 23124283.0 / 1800000000.0, 0, 0, 0, 0, 0, 0 },
        .{ 16016141.0 / 946692911.0, 0, 0, 61564180.0 / 158732637.0, 22789713.0 / 633445777.0, 545815736.0 / 2771057229.0, -180193667.0 / 1043307555.0, 0, 0, 0, 0, 0 },
        .{ 39632708.0 / 573591083.0, 0, 0, -433636366.0 / 683701615.0, -421739975.0 / 2616292301.0, 100302831.0 / 723423059.0, 790204164.0 / 839813087.0, 800635310.0 / 3783071287.0, 0, 0, 0, 0 },
        .{ 246121993.0 / 1340847787.0, 0, 0, -37695042795.0 / 15268766246.0, -309121744.0 / 1061227803.0, -12992083.0 / 490766935.0, 6005943493.0 / 2108947869.0, 393006217.0 / 1396673457.0, 123872331.0 / 1001029789.0, 0, 0, 0 },
        .{ -1028468189.0 / 846180014.0, 0, 0, 8478235783.0 / 508512852.0, 1311729495.0 / 1432422823.0, -10304129995.0 / 1701304382.0, -48777925059.0 / 3047939560.0, 15336726248.0 / 1032824649.0, -45442868181.0 / 3398467696.0, 3065993473.0 / 597172653.0, 0, 0 },
        .{ 185892177.0 / 718116043.0, 0, 0, -3185094517.0 / 667107341.0, -477755414.0 / 1098053517.0, -703635378.0 / 230739211.0, 5731566787.0 / 1027545527.0, 5232866602.0 / 850066563.0, -4093664535.0 / 808688257.0, 3962137247.0 / 1805957418.0, 65686358.0 / 487910083.0, 0 },
        .{ 403863854.0 / 491063109.0, 0, 0, -5068492393.0 / 434740067.0, -411421997.0 / 543043805.0, 652783627.0 / 914296604.0, 11173962825.0 / 925320556.0, -13158990841.0 / 6184727034.0, 3936647629.0 / 1978049680.0, -160528059.0 / 685178525.0, 248638103.0 / 1413531060.0, 0 },
    };

    // b_i coefficients for 8th order solution
    const b8 = [13]f64{
        14005451.0 / 335480064.0,   0,                           0,                         0,                            0,
        -59238493.0 / 1068277825.0, 181606767.0 / 758867731.0,   561292985.0 / 797845732.0, -1041891430.0 / 1371343529.0, 760417239.0 / 1151165299.0,
        118820643.0 / 751138087.0,  -528747749.0 / 2220607170.0, 1.0 / 4.0,
    };

    // b_i coefficients for 7th order error estimate
    const b7 = [13]f64{
        13451932.0 / 455176623.0,   0,                           0,                         0,                            0,
        -808719846.0 / 976000145.0, 1757004468.0 / 5645159321.0, 656045339.0 / 265891186.0, -3867574721.0 / 1518517206.0, 465885868.0 / 322736535.0,
        53011238.0 / 667516719.0,   2.0 / 45.0,                  0,
    };

    pub fn init() DormandPrince87 {
        return .{};
    }

    pub fn initWithTolerance(rtol: f64, atol: f64) DormandPrince87 {
        return .{ .rtol = rtol, .atol = atol };
    }

    pub fn integrator(self: *DormandPrince87) Integrator {
        return .{ .ptr = self, .vtable = &vtable };
    }

    const vtable = Integrator.VTable{ .step = step };

    fn step(ptr: *anyopaque, state: [6]f64, t: f64, dt: f64, force: ?ForceModel) anyerror![6]f64 {
        const self: *DormandPrince87 = @ptrCast(@alignCast(ptr));
        const fm = force orelse return error.ForceModelRequired;

        var currentState = state;
        var currentT = t;
        var remaining = dt;
        var substeps: u32 = 0;
        var h = @min(self.hCurrent, remaining);

        while (remaining > 1e-14 and substeps < self.maxSubsteps) {
            h = @min(h, remaining);
            h = @max(h, self.hMin);

            const result = adaptiveStep(self, currentState, currentT, h, fm);

            if (result.accepted) {
                currentState = result.yNew;
                currentT += h;
                remaining -= h;
                substeps += 1;
            }
            h = result.hNew;
        }

        self.hCurrent = h;
        return currentState;
    }

    const StepResult = struct { yNew: [6]f64, hNew: f64, accepted: bool };

    fn adaptiveStep(self: *const DormandPrince87, y: [6]f64, t: f64, h: f64, fm: ForceModel) StepResult {
        var k: [13][6]f64 = undefined;
        k[0] = derivative(y, t, fm);

        inline for (1..13) |i| {
            var yStage = y;
            inline for (0..i) |j| {
                if (a[i][j] != 0) {
                    for (0..6) |n| yStage[n] += a[i][j] * h * k[j][n];
                }
            }
            k[i] = derivative(yStage, t + c[i] * h, fm);
        }

        var y8: [6]f64 = y;
        var y7: [6]f64 = y;
        inline for (0..13) |i| {
            if (b8[i] != 0) {
                for (0..6) |n| {
                    y8[n] += b8[i] * h * k[i][n];
                }
            }
            if (b7[i] != 0) {
                for (0..6) |n| {
                    y7[n] += b7[i] * h * k[i][n];
                }
            }
        }

        var errNorm: f64 = 0;
        for (0..6) |i| {
            const scale = self.atol + self.rtol * @max(@abs(y[i]), @abs(y8[i]));
            const scaledErr = (y8[i] - y7[i]) / scale;
            errNorm += scaledErr * scaledErr;
        }
        errNorm = @sqrt(errNorm / 6.0);

        const accepted = errNorm <= 1.0;
        var hNew: f64 = if (errNorm < 1e-10) h * 5.0 else h * @min(5.0, @max(0.1, self.safety * std.math.pow(f64, 1.0 / errNorm, 1.0 / 8.0)));
        hNew = @min(hNew, self.hMax);
        hNew = @max(hNew, self.hMin);

        return .{ .yNew = y8, .hNew = hNew, .accepted = accepted };
    }

    fn derivative(state: [6]f64, t: f64, force: ForceModel) [6]f64 {
        const accel = force.acceleration(state, t);
        return .{ state[3], state[4], state[5], accel[0], accel[1], accel[2] };
    }

    pub fn getCurrentStepSize(self: *const DormandPrince87) f64 {
        return self.hCurrent;
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

    pub fn propagate(self: *const Sgp4Integrator, tsinceMinutes: f64) Error![2][3]f64 {
        return self.sgp4.propagate(tsinceMinutes);
    }

    const vtable = Integrator.VTable{ .step = step };

    fn step(ptr: *anyopaque, state: [6]f64, t: f64, dt: f64, force: ?ForceModel) anyerror![6]f64 {
        _ = state;
        _ = force;
        const self: *Sgp4Integrator = @ptrCast(@alignCast(ptr));
        const tsinceMinutes = (t + dt) / 60.0;
        const result = try self.sgp4.propagate(tsinceMinutes);
        return .{ result[0][0], result[0][1], result[0][2], result[1][0], result[1][1], result[1][2] };
    }
};

test "rk4 with twobody" {
    const fm = @import("ForceModel.zig");
    var tb = fm.TwoBody.init(398600.4418);
    var rk4 = Rk4{};
    const state = [6]f64{ 7000, 0, 0, 0, 7.5, 0 };
    const result = try rk4.integrator().step(state, 0, 60, fm.ForceModel.wrap(fm.TwoBody, &tb));
    try std.testing.expect(result[0] != state[0]);
}

test "sgp4 integrator" {
    const testTle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;
    var tle = try Tle.parse(testTle, std.testing.allocator);
    defer tle.deinit();

    var sgp4Int = try Sgp4Integrator.init(tle, constants.wgs84);
    const result = try sgp4Int.integrator().step([6]f64{ 0, 0, 0, 0, 0, 0 }, 0, 60, null);
    const rMag = @sqrt(result[0] * result[0] + result[1] * result[1] + result[2] * result[2]);
    try std.testing.expect(rMag > 6000.0 and rMag < 8000.0);
}

test "dormand-prince 8(7) basic" {
    const fm = @import("ForceModel.zig");
    var tb = fm.TwoBody.init(398600.4418);
    var dp87 = DormandPrince87.init();
    const state = [6]f64{ 7000, 0, 0, 0, 7.5, 0 };
    const result = try dp87.integrator().step(state, 0, 60, fm.ForceModel.wrap(fm.TwoBody, &tb));
    try std.testing.expect(result[0] != state[0]);
    try std.testing.expect(result[0] < state[0]);
}

test "dormand-prince 8(7) vs rk4 accuracy" {
    const fm = @import("ForceModel.zig");
    const mu = 398600.4418;
    var tb = fm.TwoBody.init(mu);
    const force = fm.ForceModel.wrap(fm.TwoBody, &tb);

    const r0 = 7000.0;
    const v0 = @sqrt(mu / r0);
    const state = [6]f64{ r0, 0, 0, 0, v0, 0 };
    const period = 2.0 * std.math.pi * @sqrt(r0 * r0 * r0 / mu);

    var rk4 = Rk4{};
    var rk4State = state;
    const rk4Dt = 10.0;
    var t: f64 = 0;
    while (t < period) {
        const stepDt = @min(rk4Dt, period - t);
        rk4State = try rk4.integrator().step(rk4State, t, stepDt, force);
        t += stepDt;
    }

    var dp87 = DormandPrince87.initWithTolerance(1e-12, 1e-14);
    const dp87State = try dp87.integrator().step(state, 0, period, force);

    const rk4R = @sqrt(rk4State[0] * rk4State[0] + rk4State[1] * rk4State[1] + rk4State[2] * rk4State[2]);
    const dp87R = @sqrt(dp87State[0] * dp87State[0] + dp87State[1] * dp87State[1] + dp87State[2] * dp87State[2]);

    const rk4Err = @abs(rk4R - r0);
    const dp87Err = @abs(dp87R - r0);

    try std.testing.expect(dp87Err < rk4Err);
}

test "dormand-prince 8(7) energy conservation" {
    const fm = @import("ForceModel.zig");
    const mu = 398600.4418;
    var tb = fm.TwoBody.init(mu);

    const state = [6]f64{ 7000, 0, 0, 0, 7.5, 0 };
    const r0 = @sqrt(state[0] * state[0] + state[1] * state[1] + state[2] * state[2]);
    const v0 = @sqrt(state[3] * state[3] + state[4] * state[4] + state[5] * state[5]);
    const energy0 = v0 * v0 / 2.0 - mu / r0;

    var dp87 = DormandPrince87.initWithTolerance(1e-10, 1e-12);
    const result = try dp87.integrator().step(state, 0, 3600, fm.ForceModel.wrap(fm.TwoBody, &tb));

    const r1 = @sqrt(result[0] * result[0] + result[1] * result[1] + result[2] * result[2]);
    const v1 = @sqrt(result[3] * result[3] + result[4] * result[4] + result[5] * result[5]);
    const energy1 = v1 * v1 / 2.0 - mu / r1;

    const relErr = @abs(energy1 - energy0) / @abs(energy0);
    try std.testing.expect(relErr < 1e-8);
}

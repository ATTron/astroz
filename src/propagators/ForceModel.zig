//! Force models for computing orbital accelerations
const std = @import("std");

pub const ForceModel = struct {
    ptr: *anyopaque,
    vtable: *const VTable,

    pub const VTable = struct {
        acceleration: *const fn (ptr: *anyopaque, state: [6]f64, t: f64) [3]f64,
    };

    pub fn acceleration(self: ForceModel, state: [6]f64, t: f64) [3]f64 {
        return self.vtable.acceleration(self.ptr, state, t);
    }
};

pub const TwoBody = struct {
    mu: f64,

    pub fn init(mu: f64) TwoBody {
        return .{ .mu = mu };
    }

    pub fn forceModel(self: *TwoBody) ForceModel {
        return .{ .ptr = self, .vtable = &vtable };
    }

    const vtable = ForceModel.VTable{ .acceleration = acceleration };

    fn acceleration(ptr: *anyopaque, state: [6]f64, t: f64) [3]f64 {
        _ = t;
        const self: *TwoBody = @ptrCast(@alignCast(ptr));
        const x, const y, const z = .{ state[0], state[1], state[2] };
        const r = @sqrt(x * x + y * y + z * z);
        const factor = -self.mu / (r * r * r);
        return .{ factor * x, factor * y, factor * z };
    }
};

pub const J2 = struct {
    mu: f64,
    j2: f64,
    r_eq: f64,

    pub fn init(mu: f64, j2: f64, r_eq: f64) J2 {
        return .{ .mu = mu, .j2 = j2, .r_eq = r_eq };
    }

    pub fn forceModel(self: *J2) ForceModel {
        return .{ .ptr = self, .vtable = &vtable };
    }

    const vtable = ForceModel.VTable{ .acceleration = acceleration };

    fn acceleration(ptr: *anyopaque, state: [6]f64, t: f64) [3]f64 {
        _ = t;
        const self: *J2 = @ptrCast(@alignCast(ptr));
        const x, const y, const z = .{ state[0], state[1], state[2] };
        const r2 = x * x + y * y + z * z;
        const r = @sqrt(r2);
        const factor = -1.5 * self.j2 * self.mu * self.r_eq * self.r_eq / (r2 * r2 * r);
        const z2_r2 = (z * z) / r2;
        return .{
            factor * x * (5.0 * z2_r2 - 1.0),
            factor * y * (5.0 * z2_r2 - 1.0),
            factor * z * (5.0 * z2_r2 - 3.0),
        };
    }
};

pub const Drag = struct {
    r_eq: f64,
    rho0: f64,
    H: f64,
    cd: f64,
    area: f64,
    mass: f64,
    max_altitude: f64,

    pub fn init(r_eq: f64, rho0: f64, H: f64, cd: f64, area: f64, mass: f64, max_altitude: f64) Drag {
        return .{ .r_eq = r_eq, .rho0 = rho0, .H = H, .cd = cd, .area = area, .mass = mass, .max_altitude = max_altitude };
    }

    pub fn forceModel(self: *Drag) ForceModel {
        return .{ .ptr = self, .vtable = &vtable };
    }

    const vtable = ForceModel.VTable{ .acceleration = acceleration };

    fn acceleration(ptr: *anyopaque, state: [6]f64, t: f64) [3]f64 {
        _ = t;
        const self: *Drag = @ptrCast(@alignCast(ptr));
        const x, const y, const z = .{ state[0], state[1], state[2] };
        const vx, const vy, const vz = .{ state[3], state[4], state[5] };

        const r = @sqrt(x * x + y * y + z * z);
        const altitude = r - self.r_eq;
        if (altitude > self.max_altitude) return .{ 0, 0, 0 };

        const v = @sqrt(vx * vx + vy * vy + vz * vz);
        if (v < 1e-10) return .{ 0, 0, 0 };

        const rho = self.rho0 * @exp(-altitude / self.H);
        const factor = -0.5 * self.cd * self.area * rho * v * 1e3 / self.mass;
        return .{ factor * vx / v, factor * vy / v, factor * vz / v };
    }
};

pub const Composite = struct {
    models: []ForceModel,
    allocator: std.mem.Allocator,

    pub fn init(allocator: std.mem.Allocator, models: []const ForceModel) !Composite {
        const owned = try allocator.alloc(ForceModel, models.len);
        @memcpy(owned, models);
        return .{ .models = owned, .allocator = allocator };
    }

    pub fn deinit(self: *Composite) void {
        self.allocator.free(self.models);
    }

    pub fn forceModel(self: *Composite) ForceModel {
        return .{ .ptr = self, .vtable = &vtable };
    }

    const vtable = ForceModel.VTable{ .acceleration = acceleration };

    fn acceleration(ptr: *anyopaque, state: [6]f64, t: f64) [3]f64 {
        const self: *Composite = @ptrCast(@alignCast(ptr));
        var total = [3]f64{ 0, 0, 0 };
        for (self.models) |model| {
            const a = model.acceleration(state, t);
            total[0] += a[0];
            total[1] += a[1];
            total[2] += a[2];
        }
        return total;
    }
};

test "twobody" {
    var tb = TwoBody.init(398600.4418);
    const fm = tb.forceModel();
    const state = [6]f64{ 7000, 0, 0, 0, 7.5, 0 };
    const accel = fm.acceleration(state, 0);
    try std.testing.expect(accel[0] < 0);
}

test "composite" {
    var tb = TwoBody.init(398600.4418);
    var j2 = J2.init(398600.4418, 0.00108263, 6378.137);
    const models = [_]ForceModel{ tb.forceModel(), j2.forceModel() };
    var comp = try Composite.init(std.testing.allocator, &models);
    defer comp.deinit();
    const state = [6]f64{ 7000, 0, 0, 0, 7.5, 0 };
    const accel = comp.forceModel().acceleration(state, 0);
    try std.testing.expect(accel[0] < 0);
}

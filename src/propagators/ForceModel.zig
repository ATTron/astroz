//! Force models for computing orbital accelerations
//!
//! All force models implement a common interface via `ForceModel.wrap()`.
//! To create a new force model, define a struct with:
//!   `pub fn acceleration(self: *const Self, state: [6]f64, t: f64) [3]f64`
//!
//! Usage:
//!   var tb = TwoBody.init(398600.4418);
//!   const fm = ForceModel.wrap(TwoBody, &tb);
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

    /// Wraps any type with an `acceleration(*const T, [6]f64, f64) [3]f64` method
    /// into a ForceModel interface
    pub fn wrap(comptime T: type, ptr: *T) ForceModel {
        const Gen = struct {
            const vtable = VTable{
                .acceleration = struct {
                    fn f(p: *anyopaque, state: [6]f64, t: f64) [3]f64 {
                        const self: *const T = @ptrCast(@alignCast(p));
                        return self.acceleration(state, t);
                    }
                }.f,
            };
        };
        return .{ .ptr = ptr, .vtable = &Gen.vtable };
    }
};

pub const TwoBody = struct {
    mu: f64,

    pub fn init(mu: f64) TwoBody {
        return .{ .mu = mu };
    }

    pub fn acceleration(self: *const TwoBody, state: [6]f64, t: f64) [3]f64 {
        _ = t;
        const x, const y, const z = .{ state[0], state[1], state[2] };
        const r = @sqrt(x * x + y * y + z * z);
        const factor = -self.mu / (r * r * r);
        return .{ factor * x, factor * y, factor * z };
    }
};

pub const J2 = struct {
    mu: f64,
    j2: f64,
    rEq: f64,

    pub fn init(mu: f64, j2: f64, rEq: f64) J2 {
        return .{ .mu = mu, .j2 = j2, .rEq = rEq };
    }

    pub fn acceleration(self: *const J2, state: [6]f64, t: f64) [3]f64 {
        _ = t;
        const x, const y, const z = .{ state[0], state[1], state[2] };
        const r2 = x * x + y * y + z * z;
        const r = @sqrt(r2);
        const factor = -1.5 * self.j2 * self.mu * self.rEq * self.rEq / (r2 * r2 * r);
        const z2R2 = (z * z) / r2;
        return .{
            factor * x * (5.0 * z2R2 - 1.0),
            factor * y * (5.0 * z2R2 - 1.0),
            factor * z * (5.0 * z2R2 - 3.0),
        };
    }
};

pub const Drag = struct {
    rEq: f64,
    rho0: f64,
    H: f64,
    cd: f64,
    area: f64,
    mass: f64,
    maxAltitude: f64,

    pub fn init(rEq: f64, rho0: f64, H: f64, cd: f64, area: f64, mass: f64, maxAltitude: f64) Drag {
        return .{ .rEq = rEq, .rho0 = rho0, .H = H, .cd = cd, .area = area, .mass = mass, .maxAltitude = maxAltitude };
    }

    pub fn acceleration(self: *const Drag, state: [6]f64, t: f64) [3]f64 {
        _ = t;
        const x, const y, const z = .{ state[0], state[1], state[2] };
        const vx, const vy, const vz = .{ state[3], state[4], state[5] };

        const r = @sqrt(x * x + y * y + z * z);
        const altitude = r - self.rEq;
        if (altitude > self.maxAltitude) return .{ 0, 0, 0 };

        const v = @sqrt(vx * vx + vy * vy + vz * vz);
        if (v < 1e-10) return .{ 0, 0, 0 };

        const rho = self.rho0 * @exp(-altitude / self.H);
        const factor = -0.5 * self.cd * self.area * rho * v * 1e3 / self.mass;
        return .{ factor * vx / v, factor * vy / v, factor * vz / v };
    }
};

pub const J3 = struct {
    mu: f64,
    j3: f64,
    rEq: f64,

    pub fn init(mu: f64, j3: f64, rEq: f64) J3 {
        return .{ .mu = mu, .j3 = j3, .rEq = rEq };
    }

    pub fn acceleration(self: *const J3, state: [6]f64, t: f64) [3]f64 {
        _ = t;
        const x, const y, const z = .{ state[0], state[1], state[2] };
        const r2 = x * x + y * y + z * z;
        const r = @sqrt(r2);
        const rEq3 = self.rEq * self.rEq * self.rEq;

        // J3 perturbation formula (Vallado)
        // factor = (5/2) * J3 * mu * R_e^3 / r^7
        const factor = 2.5 * self.j3 * self.mu * rEq3 / (r2 * r2 * r2 * r);
        const z2R2 = (z * z) / r2;

        const xyCoeff = 3.0 * z / r - 7.0 * z * z2R2 / r;
        const zCoeff = 6.0 * z * z - 7.0 * z * z * z2R2 - 0.6 * r2;

        return .{
            factor * x * xyCoeff,
            factor * y * xyCoeff,
            factor * zCoeff,
        };
    }
};

pub const J4 = struct {
    mu: f64,
    j4: f64,
    rEq: f64,

    pub fn init(mu: f64, j4: f64, rEq: f64) J4 {
        return .{ .mu = mu, .j4 = j4, .rEq = rEq };
    }

    pub fn acceleration(self: *const J4, state: [6]f64, t: f64) [3]f64 {
        _ = t;
        const x, const y, const z = .{ state[0], state[1], state[2] };
        const r2 = x * x + y * y + z * z;
        const r = @sqrt(r2);
        const r4 = r2 * r2;
        const z2 = z * z;
        const z4 = z2 * z2;
        const z2R2 = z2 / r2;
        const z4R4 = z4 / r4;

        const rEq4 = self.rEq * self.rEq * self.rEq * self.rEq;
        const factor = 1.875 * self.j4 * self.mu * rEq4 / (r4 * r4 * r);
        const xyTerm = 3.0 - 42.0 * z2R2 + 63.0 * z4R4;
        const zTerm = 15.0 - 70.0 * z2R2 + 63.0 * z4R4;

        return .{
            factor * x * xyTerm,
            factor * y * xyTerm,
            factor * z * zTerm,
        };
    }
};

pub const SolarRadiationPressure = struct {
    cr: f64, // reflectivity coefficient (1.0 = absorber, 2.0 = perfect reflector)
    area: f64, // cross-sectional area (m^2)
    mass: f64, // spacecraft mass (kg)
    rEq: f64, // Earth radius for shadow model (km)

    // Solar radiation pressure at 1 AU in N/m^2
    const pSr: f64 = 4.56e-6;

    pub fn init(cr: f64, area: f64, mass: f64, rEq: f64) SolarRadiationPressure {
        return .{ .cr = cr, .area = area, .mass = mass, .rEq = rEq };
    }

    pub fn acceleration(self: *const SolarRadiationPressure, state: [6]f64, t: f64) [3]f64 {
        _ = t; // Future: use for Sun ephemeris
        const x, const y, const z = .{ state[0], state[1], state[2] };

        // Simplified model: Sun at +X direction, 1 AU distance
        const sunDir = [3]f64{ 1.0, 0.0, 0.0 };

        // Check for cylindrical shadow (spacecraft behind Earth relative to Sun)
        // If x < 0 (behind Earth) and sqrt(y^2 + z^2) < rEq (within shadow cylinder)
        if (x < 0) {
            const rho = @sqrt(y * y + z * z);
            if (rho < self.rEq) {
                return .{ 0, 0, 0 }; // In shadow
            }
        }

        // SRP acceleration: a = -Cr * pSr * A / m * sunDirection
        // Convert km/s^2: multiply by 1e-3 (N to kN equivalent for km units)
        const factor = -self.cr * pSr * self.area / self.mass * 1e-3;
        return .{ factor * sunDir[0], factor * sunDir[1], factor * sunDir[2] };
    }
};

pub const ImprovedDrag = struct {
    rEq: f64,
    cd: f64,
    area: f64,
    mass: f64,
    maxAltitude: f64,
    f107: f64, // F10.7 solar flux proxy (typically 70-250)

    // Atmospheric density layers (US Standard Atmosphere 1976 + extensions)
    const DensityLayer = struct {
        baseAlt: f64, // km
        baseDensity: f64, // kg/m^3
        scaleHeight: f64, // km
    };

    const densityLayers = [_]DensityLayer{
        .{ .baseAlt = 100.0, .baseDensity = 5.297e-7, .scaleHeight = 5.877 },
        .{ .baseAlt = 200.0, .baseDensity = 2.789e-10, .scaleHeight = 37.105 },
        .{ .baseAlt = 400.0, .baseDensity = 3.725e-12, .scaleHeight = 62.822 },
        .{ .baseAlt = 600.0, .baseDensity = 2.418e-13, .scaleHeight = 79.864 },
        .{ .baseAlt = 1000.0, .baseDensity = 3.561e-15, .scaleHeight = 200.0 },
    };

    // Earth rotation rate (rad/s)
    const omega: f64 = 7.2921150e-5;

    pub fn init(rEq: f64, cd: f64, area: f64, mass: f64, maxAltitude: f64, f107: f64) ImprovedDrag {
        return .{ .rEq = rEq, .cd = cd, .area = area, .mass = mass, .maxAltitude = maxAltitude, .f107 = f107 };
    }

    fn getDensity(altitude: f64, f107Val: f64) f64 {
        // Find appropriate layer
        var layerIdx: usize = 0;
        for (densityLayers, 0..) |layer, i| {
            if (altitude >= layer.baseAlt) {
                layerIdx = i;
            }
        }

        const layer = densityLayers[layerIdx];
        const baseRho = layer.baseDensity;
        const H = layer.scaleHeight;
        const deltaH = altitude - layer.baseAlt;

        // Exponential decay from base
        var rho = baseRho * @exp(-deltaH / H);

        // F10.7 scaling: nominal is 150, scale linearly
        const f107Scale = f107Val / 150.0;
        rho *= f107Scale;

        return rho;
    }

    pub fn acceleration(self: *const ImprovedDrag, state: [6]f64, t: f64) [3]f64 {
        _ = t;
        const x, const y, const z = .{ state[0], state[1], state[2] };
        const vx, const vy, const vz = .{ state[3], state[4], state[5] };

        const r = @sqrt(x * x + y * y + z * z);
        const altitude = r - self.rEq;
        if (altitude > self.maxAltitude or altitude < 100.0) return .{ 0, 0, 0 };

        // Relative velocity (atmosphere co-rotates with Earth)
        const vrelX = vx + omega * y;
        const vrelY = vy - omega * x;
        const vrelZ = vz;

        const vrel = @sqrt(vrelX * vrelX + vrelY * vrelY + vrelZ * vrelZ);
        if (vrel < 1e-10) return .{ 0, 0, 0 };

        const rho = getDensity(altitude, self.f107);
        // factor includes 1e3 to convert km/s to m/s for consistent units
        const factor = -0.5 * self.cd * self.area * rho * vrel * 1e3 / self.mass;

        return .{
            factor * vrelX / vrel,
            factor * vrelY / vrel,
            factor * vrelZ / vrel,
        };
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

    pub fn acceleration(self: *const Composite, state: [6]f64, t: f64) [3]f64 {
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
    const fm = ForceModel.wrap(TwoBody, &tb);
    const state = [6]f64{ 7000, 0, 0, 0, 7.5, 0 };
    const accel = fm.acceleration(state, 0);
    try std.testing.expect(accel[0] < 0);
}

test "composite" {
    var tb = TwoBody.init(398600.4418);
    var j2 = J2.init(398600.4418, 0.00108263, 6378.137);
    const models = [_]ForceModel{ ForceModel.wrap(TwoBody, &tb), ForceModel.wrap(J2, &j2) };
    var comp = try Composite.init(std.testing.allocator, &models);
    defer comp.deinit();
    const state = [6]f64{ 7000, 0, 0, 0, 7.5, 0 };
    const accel = ForceModel.wrap(Composite, &comp).acceleration(state, 0);
    try std.testing.expect(accel[0] < 0);
}

test "j3 magnitude smaller than j2" {
    const mu = 398600.5;
    const rEq = 6378.137;
    const j2Val = 0.00108262998905;
    const j3Val = -0.00000253215306;

    var j2Model = J2.init(mu, j2Val, rEq);
    var j3Model = J3.init(mu, j3Val, rEq);

    // Test at inclined position (z != 0 to see J3 effect)
    const state = [6]f64{ 6000, 1000, 2000, 0, 7.5, 0 };
    const aJ2 = ForceModel.wrap(J2, &j2Model).acceleration(state, 0);
    const aJ3 = ForceModel.wrap(J3, &j3Model).acceleration(state, 0);

    const magJ2 = @sqrt(aJ2[0] * aJ2[0] + aJ2[1] * aJ2[1] + aJ2[2] * aJ2[2]);
    const magJ3 = @sqrt(aJ3[0] * aJ3[0] + aJ3[1] * aJ3[1] + aJ3[2] * aJ3[2]);

    // J3 should be smaller than J2 (roughly 1000x smaller based on coefficient ratio)
    try std.testing.expect(magJ3 < magJ2);
    try std.testing.expect(magJ3 > 0);
}

test "j4 magnitude smaller than j2" {
    const mu = 398600.5;
    const rEq = 6378.137;
    const j2Val = 0.00108262998905;
    const j4Val = -0.00000161098761;

    var j2Model = J2.init(mu, j2Val, rEq);
    var j4Model = J4.init(mu, j4Val, rEq);

    const state = [6]f64{ 6000, 1000, 2000, 0, 7.5, 0 };
    const aJ2 = ForceModel.wrap(J2, &j2Model).acceleration(state, 0);
    const aJ4 = ForceModel.wrap(J4, &j4Model).acceleration(state, 0);

    const magJ2 = @sqrt(aJ2[0] * aJ2[0] + aJ2[1] * aJ2[1] + aJ2[2] * aJ2[2]);
    const magJ4 = @sqrt(aJ4[0] * aJ4[0] + aJ4[1] * aJ4[1] + aJ4[2] * aJ4[2]);

    // J4 should be smaller than J2
    try std.testing.expect(magJ4 < magJ2);
    try std.testing.expect(magJ4 > 0);
}

test "srp shadow model" {
    var srp = SolarRadiationPressure.init(1.5, 10.0, 500.0, 6378.137);
    const fm = ForceModel.wrap(SolarRadiationPressure, &srp);

    // In sunlight (positive x)
    const litState = [6]f64{ 7000, 0, 0, 0, 7.5, 0 };
    const aLit = fm.acceleration(litState, 0);
    try std.testing.expect(aLit[0] != 0);

    // In shadow (negative x, within Earth radius)
    const shadowState = [6]f64{ -7000, 0, 0, 0, 7.5, 0 };
    const aShadow = fm.acceleration(shadowState, 0);
    try std.testing.expect(aShadow[0] == 0);
    try std.testing.expect(aShadow[1] == 0);
    try std.testing.expect(aShadow[2] == 0);
}

test "improved drag uses relative velocity" {
    var drag = ImprovedDrag.init(6378.137, 2.2, 10.0, 500.0, 1500.0, 150.0);
    const fm = ForceModel.wrap(ImprovedDrag, &drag);

    // At 400 km altitude, position at +X, velocity in +Y (prograde)
    // Earth rotation: vrel_x = vx + omega*y = 0 + 0 = 0
    //                 vrel_y = vy - omega*x = 7.67 - omega*6778.137
    // omega*6778.137 = 7.2921e-5 * 6778.137 ≈ 0.494 km/s
    // So vrel_y ≈ 7.67 - 0.494 = 7.176 km/s (prograde relative to atmosphere)
    const state = [6]f64{ 6778.137, 0, 0, 0, 7.67, 0 };
    const accel = fm.acceleration(state, 0);

    // Drag should oppose relative velocity direction (primarily -Y)
    try std.testing.expect(accel[1] < 0);
    // Should have magnitude > 0
    const mag = @sqrt(accel[0] * accel[0] + accel[1] * accel[1] + accel[2] * accel[2]);
    try std.testing.expect(mag > 0);
}

test "improved drag outside altitude range" {
    var drag = ImprovedDrag.init(6378.137, 2.2, 10.0, 500.0, 1000.0, 150.0);
    const fm = ForceModel.wrap(ImprovedDrag, &drag);

    // At 2000 km altitude (above maxAltitude)
    const highState = [6]f64{ 8378.137, 0, 0, 0, 6.9, 0 };
    const aHigh = fm.acceleration(highState, 0);
    try std.testing.expect(aHigh[0] == 0);
    try std.testing.expect(aHigh[1] == 0);
    try std.testing.expect(aHigh[2] == 0);

    // At 50 km altitude (below minimum)
    const lowState = [6]f64{ 6428.137, 0, 0, 0, 7.8, 0 };
    const aLow = fm.acceleration(lowState, 0);
    try std.testing.expect(aLow[0] == 0);
    try std.testing.expect(aLow[1] == 0);
    try std.testing.expect(aLow[2] == 0);
}

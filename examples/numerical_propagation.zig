//! Numerical propagation with composable force models: the same orbit over one
//! day with two-body gravity (RK4 and Dormand-Prince 8(7)), adding J2, and the
//! full model with J2, solar radiation pressure and Sun/Moon gravity. With CSPICE
//! enabled the Sun and Moon positions come from NAIF ephemeris.

const std = @import("std");
const astroz = @import("astroz");
const Spice = astroz.Spice;
const constants = astroz.constants;
const p = astroz.propagators;

fn radius(s: [6]f64) f64 {
    return @sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
}

fn distance(a: [6]f64, b: [6]f64) f64 {
    return @sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
}

pub fn main(init: std.process.Init) !void {
    const gpa = init.gpa;
    const mu = constants.earth.mu;
    const rEq = constants.earth.eqRadius.?;

    // ISS-like orbit: circular at 420 km, 51.6 degrees
    const r0 = rEq + 420.0;
    const v0 = @sqrt(mu / r0);
    const inc = 51.6 * constants.deg2rad;
    const initial = [6]f64{ r0, 0, 0, 0, v0 * @cos(inc), v0 * @sin(inc) };
    const oneDay = 86400.0;

    var twobody = p.TwoBody.init(mu);
    var j2 = p.J2.init(mu, constants.earth.j2Perturbation, rEq);
    var srp = p.SolarRadiationPressure.init(1.5, 20.0, 1000.0, rEq);
    var sun = p.ThirdBody.init(constants.sun.mu, .{ constants.auKm, 0, 0 });
    var moon = p.ThirdBody.init(constants.moon.mu, .{ 384400, 0, 0 });

    const gravity = p.ForceModel.wrap(p.TwoBody, &twobody);
    var withJ2 = try p.Composite.init(gpa, &.{ gravity, p.ForceModel.wrap(p.J2, &j2) });
    defer withJ2.deinit();
    var full = try p.Composite.init(gpa, &.{
        gravity,
        p.ForceModel.wrap(p.J2, &j2),
        p.ForceModel.wrap(p.SolarRadiationPressure, &srp),
        p.ForceModel.wrap(p.ThirdBody, &sun),
        p.ForceModel.wrap(p.ThirdBody, &moon),
    });
    defer full.deinit();

    var rk4 = p.Rk4{};
    var dp87 = p.DormandPrince87.initWithTolerance(1e-10, 1e-12);

    var rk4Prop = p.Propagator.init(gpa, rk4.integrator(), gravity);
    const twoBodyRk4 = try rk4Prop.propagateTo(initial, 0, oneDay, 10);
    const twoBodyDp87 = try dp87.integrator().step(initial, 0, oneDay, gravity);
    const j2State = try dp87.integrator().step(initial, 0, oneDay, p.ForceModel.wrap(p.Composite, &withJ2));

    // Full model: RK4 in 10 s steps, refreshing the Sun and Moon from SPICE every 10 minutes
    const spiceEnabled = if (Spice.loadDefaultKernels()) true else |err| blk: {
        if (err == Spice.SpiceError.NotEnabled)
            std.debug.print("CSPICE not enabled, using fixed Sun/Moon positions.\n", .{})
        else
            std.debug.print("Failed to load kernels ({}), using fixed Sun/Moon positions.\n", .{err});
        std.debug.print("For real ephemeris: zig build fetch-kernels && zig build example -Denable-cspice=true\n\n", .{});
        break :blk false;
    };
    const epochEt = if (spiceEnabled) Spice.utcToEt("2024-06-21T12:00:00") catch Spice.jdToEt(2460483.0) else Spice.jdToEt(2460483.0);
    const fullForce = p.ForceModel.wrap(p.Composite, &full);
    var fullState = initial;
    var t: f64 = 0;
    while (t < oneDay) : (t += 10) {
        if (spiceEnabled and @mod(t, 600) == 0) {
            if (Spice.getSunPosition(epochEt + t)) |pos| {
                srp.updateSunPos(pos.array());
                sun.updatePos(pos.array());
            } else |_| {}
            if (Spice.getMoonPosition(epochEt + t)) |pos| moon.updatePos(pos.array()) else |_| {}
        }
        fullState = try rk4.integrator().step(fullState, t, 10, fullForce);
    }

    std.debug.print("After one day ({s} Sun/Moon):\n", .{if (spiceEnabled) "SPICE" else "fixed"});
    std.debug.print("  {s:<26} {s:>10} {s:>18}\n", .{ "model", "r (km)", "vs two-body (km)" });
    const rows = [_]struct { []const u8, [6]f64 }{
        .{ "two-body, RK4 10 s", twoBodyRk4 },
        .{ "two-body, DP 8(7)", twoBodyDp87 },
        .{ "two-body + J2", j2State },
        .{ "J2 + SRP + Sun + Moon", fullState },
    };
    for (rows) |row| std.debug.print("  {s:<26} {d:10.3} {d:18.3}\n", .{ row[0], radius(row[1]), distance(row[1], twoBodyDp87) });
}

const std = @import("std");
const astroz = @import("astroz");
const Spice = astroz.Spice;
const constants = astroz.constants;
const propagators = astroz.propagators;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    // Load SPICE kernels (graceful fallback to fixed positions)
    const spice_enabled = if (Spice.loadDefaultKernels()) true else |err| blk: {
        if (err == Spice.SpiceError.NotEnabled)
            std.debug.print("CSPICE not enabled — using fixed Sun/Moon positions.\n", .{})
        else
            std.debug.print("Failed to load kernels ({}) — using fixed positions.\n", .{err});
        std.debug.print("For real ephemeris: zig build fetch-kernels && zig build example -Denable-cspice=true\n\n", .{});
        break :blk false;
    };

    const epoch_et: f64 = if (spice_enabled)
        Spice.utcToEt("2024-06-21T12:00:00") catch Spice.jdToEt(2460483.0)
    else
        Spice.jdToEt(2460483.0);

    // ISS-like orbit: 420 km altitude, 51.6° inclination
    const rEq = constants.earth.eqRadius.?;
    const r0 = rEq + 420.0;
    const v0 = @sqrt(constants.earth.mu / r0);
    const inc = 51.6 * constants.deg2rad;
    const initial = [6]f64{ r0, 0, 0, 0, v0 * @cos(inc), v0 * @sin(inc) };

    // Force models: gravity + J2 + SRP + Sun/Moon third-body
    var twobody = propagators.TwoBody.init(constants.earth.mu);
    var j2 = propagators.J2.init(constants.earth.mu, constants.earth.j2Perturbation, rEq);
    var srp = propagators.SolarRadiationPressure.init(1.5, 20.0, 1000.0, rEq);
    var sun_tb = propagators.ThirdBody.init(constants.sun.mu, .{ constants.auKm, 0, 0 });
    var moon_tb = propagators.ThirdBody.init(constants.moon.mu, .{ 384400, 0, 0 });

    const models = [_]propagators.ForceModel{
        propagators.ForceModel.wrap(propagators.TwoBody, &twobody),
        propagators.ForceModel.wrap(propagators.J2, &j2),
        propagators.ForceModel.wrap(propagators.SolarRadiationPressure, &srp),
        propagators.ForceModel.wrap(propagators.ThirdBody, &sun_tb),
        propagators.ForceModel.wrap(propagators.ThirdBody, &moon_tb),
    };
    var composite = try propagators.Composite.init(allocator, &models);
    defer composite.deinit();
    const force = propagators.ForceModel.wrap(propagators.Composite, &composite);

    // Propagate 1 day, updating Sun/Moon from SPICE every 10 minutes
    var rk4 = propagators.Rk4{};
    const integrator = rk4.integrator();
    const dt: f64 = 10.0;
    const duration: f64 = 86400.0;

    std.debug.print("=== SPICE-Driven LEO Propagation (1 day) ===\n", .{});
    std.debug.print("{s:>8} {s:>10} {s:>10} {s:>10}\n", .{ "hr", "r (km)", "alt (km)", "v (km/s)" });

    var state = initial;
    var t: f64 = 0;
    var next_print: f64 = 0;
    var next_update: f64 = 600.0;

    while (t < duration) {
        if (t >= next_print) {
            const r = @sqrt(state[0] * state[0] + state[1] * state[1] + state[2] * state[2]);
            const v = @sqrt(state[3] * state[3] + state[4] * state[4] + state[5] * state[5]);
            std.debug.print("{d:8.1} {d:10.2} {d:10.2} {d:10.4}\n", .{ t / 3600.0, r, r - rEq, v });
            next_print += 3600.0;
        }

        if (spice_enabled and t >= next_update) {
            const et = epoch_et + t;
            if (Spice.getSunPosition(et)) |p| {
                srp.updateSunPos(p.array());
                sun_tb.updatePos(p.array());
            } else |_| {}
            if (Spice.getMoonPosition(et)) |p| moon_tb.updatePos(p.array()) else |_| {}
            next_update += 600.0;
        }

        state = try integrator.step(state, t, @min(dt, duration - t), force);
        t += @min(dt, duration - t);
    }

    const r_f = @sqrt(state[0] * state[0] + state[1] * state[1] + state[2] * state[2]);
    const v_f = @sqrt(state[3] * state[3] + state[4] * state[4] + state[5] * state[5]);
    std.debug.print("{d:8.1} {d:10.2} {d:10.2} {d:10.4}\n", .{ t / 3600.0, r_f, r_f - 6378.137, v_f });
}

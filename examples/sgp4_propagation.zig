const std = @import("std");
const astroz = @import("astroz");
const Tle = astroz.Tle;
const Sgp4 = astroz.Sgp4;
const constants = astroz.constants;
const propagators = astroz.propagators;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const tle_str =
        \\1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9025
        \\2 25544  51.6400 208.9163 0006703  30.5502 329.5947 15.49560830484850
    ;

    var tle = try Tle.parse(tle_str, allocator);
    defer tle.deinit();

    // Example 1: Direct SGP4 usage
    std.debug.print("=== Direct SGP4 Propagation ===\n", .{});
    var sgp4 = try Sgp4.init(tle, constants.wgs84);

    // Propagate at different times (minutes from TLE epoch)
    const times = [_]f64{ 0, 30, 60, 90, 120 };
    for (times) |tsince| {
        const result = try sgp4.propagate(tsince);
        const pos = result[0];
        const vel = result[1];
        const r = @sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
        const v = @sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
        std.debug.print("t={d:6.1} min: r={d:8.1} km, v={d:.3} km/s\n", .{ tsince, r, v });
    }

    // Example 2: Using the modular propagator interface
    std.debug.print("\n=== SGP4 via Propagator Interface ===\n", .{});
    var sgp4_int = try propagators.Sgp4Integrator.init(tle, constants.wgs84);
    var prop = propagators.Propagator.init(allocator, sgp4_int.integrator(), null);

    // Propagate for 2 hours (7200 seconds) with 30-minute steps
    const initial = [6]f64{ 0, 0, 0, 0, 0, 0 }; // ignored by SGP4
    var trajectory = try prop.propagate(initial, 0, 7200, 1800);
    defer trajectory.deinit(allocator);

    std.debug.print("Generated {d} trajectory points\n", .{trajectory.items.len});
    for (trajectory.items[1..]) |point| { // skip first point (SGP4 ignores initial state)
        const r = @sqrt(point.state[0] * point.state[0] +
            point.state[1] * point.state[1] +
            point.state[2] * point.state[2]);
        std.debug.print("t={d:6.0}s: r={d:8.1} km\n", .{ point.time, r });
    }

    // Example 3: Compare with numerical RK4 propagation
    std.debug.print("\n=== Comparison: SGP4 vs RK4+TwoBody ===\n", .{});

    // Get initial state from SGP4 at epoch
    const sgp4_state = try sgp4.propagate(0);
    const init_state = [6]f64{
        sgp4_state[0][0], sgp4_state[0][1], sgp4_state[0][2],
        sgp4_state[1][0], sgp4_state[1][1], sgp4_state[1][2],
    };

    // Setup RK4 with two-body gravity
    var twobody = propagators.TwoBody.init(constants.earth.mu);
    var rk4 = propagators.Rk4{};
    var rk4_prop = propagators.Propagator.init(allocator, rk4.integrator(), twobody.forceModel());

    // Propagate both for 90 minutes (one orbit)
    const duration = 90.0 * 60.0;
    const sgp4_final = try sgp4.propagate(90);
    const rk4_final = try rk4_prop.propagateTo(init_state, 0, duration, 10);

    const sgp4_r = @sqrt(sgp4_final[0][0] * sgp4_final[0][0] +
        sgp4_final[0][1] * sgp4_final[0][1] +
        sgp4_final[0][2] * sgp4_final[0][2]);
    const rk4_r = @sqrt(rk4_final[0] * rk4_final[0] +
        rk4_final[1] * rk4_final[1] +
        rk4_final[2] * rk4_final[2]);

    std.debug.print("After 90 minutes:\n", .{});
    std.debug.print("  SGP4:    r = {d:.1} km\n", .{sgp4_r});
    std.debug.print("  RK4:     r = {d:.1} km\n", .{rk4_r});
    std.debug.print("  Diff:    {d:.1} km\n", .{@abs(sgp4_r - rk4_r)});

    std.debug.print("\nNote: SGP4 includes J2/drag perturbations analytically,\n", .{});
    std.debug.print("while RK4+TwoBody is pure Keplerian - expect divergence.\n", .{});
}

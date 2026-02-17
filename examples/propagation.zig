const std = @import("std");
const astroz = @import("astroz");
const Tle = astroz.Tle;
const Sgp4 = astroz.Sgp4;
const Satellite = astroz.Satellite;
const constants = astroz.constants;
const propagators = astroz.propagators;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const tleStr =
        \\1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9025
        \\2 25544  51.6400 208.9163 0006703  30.5502 329.5947 15.49560830484850
    ;

    var tle = try Tle.parse(tleStr, allocator);
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

    // Example 1.5: Using SIMD V4 (4 steps at once)
    std.debug.print("\n=== SIMD V4 Batch Propagation (4 times) ===\n", .{});
    const batchTimes4 = [4]f64{ 0, 30, 60, 90 };
    const batchResults4 = try sgp4.propagateN(4, batchTimes4);
    for (batchTimes4, batchResults4) |tsince, result| {
        const pos = result[0];
        const vel = result[1];
        const r = @sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
        const v = @sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
        std.debug.print("t={d:6.1} min: r={d:8.1} km, v={d:.3} km/s\n", .{ tsince, r, v });
    }

    // Example 1.6: Using SIMD V8 (8 steps at once, AVX512)
    std.debug.print("\n=== SIMD V8 Batch Propagation (8 times) ===\n", .{});
    const batchTimes8 = [8]f64{ 0, 15, 30, 45, 60, 75, 90, 105 };
    const batchResults8 = try sgp4.propagateN(8, batchTimes8);
    for (batchTimes8, batchResults8) |tsince, result| {
        const pos = result[0];
        const vel = result[1];
        const r = @sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
        const v = @sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
        std.debug.print("t={d:6.1} min: r={d:8.1} km, v={d:.3} km/s\n", .{ tsince, r, v });
    }

    // Example 2: Using the modular propagator interface
    std.debug.print("\n=== SGP4 via Propagator Interface ===\n", .{});
    var sgp4Int = try propagators.Sgp4Integrator.init(tle, constants.wgs84);
    var prop = propagators.Propagator.init(allocator, sgp4Int.integrator(), null);

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
    const sgp4State = try sgp4.propagate(0);
    const initState = [6]f64{
        sgp4State[0][0], sgp4State[0][1], sgp4State[0][2],
        sgp4State[1][0], sgp4State[1][1], sgp4State[1][2],
    };

    // Setup RK4 with two-body gravity
    var twobody = propagators.TwoBody.init(constants.earth.mu);
    var rk4 = propagators.Rk4{};
    var rk4Prop = propagators.Propagator.init(allocator, rk4.integrator(), propagators.ForceModel.wrap(propagators.TwoBody, &twobody));

    // Propagate both for 90 minutes (one orbit)
    const duration = 90.0 * 60.0;
    const sgp4Final = try sgp4.propagate(90);
    const rk4Final = try rk4Prop.propagateTo(initState, 0, duration, 10);

    const sgp4R = @sqrt(sgp4Final[0][0] * sgp4Final[0][0] +
        sgp4Final[0][1] * sgp4Final[0][1] +
        sgp4Final[0][2] * sgp4Final[0][2]);
    const rk4R = @sqrt(rk4Final[0] * rk4Final[0] +
        rk4Final[1] * rk4Final[1] +
        rk4Final[2] * rk4Final[2]);

    std.debug.print("After 90 minutes:\n", .{});
    std.debug.print("  SGP4:    r = {d:.1} km\n", .{sgp4R});
    std.debug.print("  RK4:     r = {d:.1} km\n", .{rk4R});
    std.debug.print("  Diff:    {d:.1} km\n", .{@abs(sgp4R - rk4R)});

    std.debug.print("\nNote: SGP4 includes J2/drag perturbations analytically,\n", .{});
    std.debug.print("while RK4+TwoBody is pure Keplerian - expect divergence.\n", .{});

    // Example 4: Unified Satellite type (auto-dispatches SGP4/SDP4)
    std.debug.print("\n=== Unified Satellite (auto SGP4/SDP4) ===\n", .{});

    // Near-earth: automatically uses SGP4
    const sat_leo = try Satellite.init(tle, constants.wgs84);
    std.debug.print("ISS (LEO): deep_space={}\n", .{sat_leo.isDeepSpace()});
    const leo_result = try sat_leo.propagate(60);
    const leo_r = @sqrt(leo_result[0][0] * leo_result[0][0] +
        leo_result[0][1] * leo_result[0][1] +
        leo_result[0][2] * leo_result[0][2]);
    std.debug.print("  t=60 min: r={d:.1} km\n", .{leo_r});

    // Deep-space: automatically uses SDP4
    const gpsTleStr =
        \\1 20413U 90005A   24186.00000000  .00000012  00000+0  10000-3 0  9992
        \\2 20413  55.4408  61.4858 0112981 129.5765 231.5553  2.00561730104446
    ;
    var gpsTle = try Tle.parse(gpsTleStr, allocator);
    defer gpsTle.deinit();

    const sat_gps = try Satellite.init(gpsTle, constants.wgs72);
    std.debug.print("GPS (MEO): deep_space={}\n", .{sat_gps.isDeepSpace()});
    const gps_result = try sat_gps.propagate(720);
    const gps_r = @sqrt(gps_result[0][0] * gps_result[0][0] +
        gps_result[0][1] * gps_result[0][1] +
        gps_result[0][2] * gps_result[0][2]);
    std.debug.print("  t=720 min: r={d:.1} km\n", .{gps_r});
}

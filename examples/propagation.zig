const std = @import("std");
const astroz = @import("astroz");
const Tle = astroz.Tle;
const Sgp4 = astroz.Sgp4;
const Satellite = astroz.Satellite;
const constants = astroz.constants;
const propagators = astroz.propagators;

fn mag(v: [3]f64) f64 {
    return @sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;

    const tleStr =
        \\1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9025
        \\2 25544  51.6400 208.9163 0006703  30.5502 329.5947 15.49560830484850
    ;

    var tle = try Tle.parse(tleStr, allocator);
    defer tle.deinit();

    std.debug.print("=== SGP4 ===\n", .{});
    var sgp4 = try Sgp4.init(tle, constants.wgs84);

    // times are minutes from the TLE epoch
    for ([_]f64{ 0, 30, 60, 90, 120 }) |tsince| {
        const pv = try sgp4.propagate(tsince);
        std.debug.print("t={d:6.1} min: r={d:8.1} km, v={d:.3} km/s\n", .{ tsince, mag(pv[0]), mag(pv[1]) });
    }

    // propagateN runs N times at once with SIMD
    std.debug.print("\n=== SGP4, 8 times per call ===\n", .{});
    const batchTimes = [8]f64{ 0, 15, 30, 45, 60, 75, 90, 105 };
    const batch = try sgp4.propagateN(8, batchTimes);
    for (batchTimes, batch) |tsince, pv| {
        std.debug.print("t={d:6.1} min: r={d:8.1} km, v={d:.3} km/s\n", .{ tsince, mag(pv[0]), mag(pv[1]) });
    }

    std.debug.print("\n=== SGP4 through the Propagator interface ===\n", .{});
    var sgp4Int = try propagators.Sgp4Integrator.init(tle, constants.wgs84);
    var prop = propagators.Propagator.init(allocator, sgp4Int.integrator(), null);

    // 2 hours in 30 minute steps; SGP4 ignores the initial state
    var trajectory = try prop.propagate(.{ 0, 0, 0, 0, 0, 0 }, 0, 7200, 1800);
    defer trajectory.deinit(allocator);

    for (trajectory.items[1..]) |point| {
        std.debug.print("t={d:6.0}s: r={d:8.1} km\n", .{ point.time, mag(point.state[0..3].*) });
    }

    std.debug.print("\n=== SGP4 vs RK4 with two-body gravity ===\n", .{});
    const epochState = try sgp4.propagate(0);
    const initState = epochState[0] ++ epochState[1];

    var twobody = propagators.TwoBody.init(constants.earth.mu);
    var rk4 = propagators.Rk4{};
    var rk4Prop = propagators.Propagator.init(allocator, rk4.integrator(), propagators.ForceModel.wrap(propagators.TwoBody, &twobody));

    // one orbit
    const sgp4R = mag((try sgp4.propagate(90))[0]);
    const rk4R = mag((try rk4Prop.propagateTo(initState, 0, 90.0 * 60.0, 10))[0..3].*);

    std.debug.print("After 90 minutes:\n", .{});
    std.debug.print("  SGP4:    r = {d:.1} km\n", .{sgp4R});
    std.debug.print("  RK4:     r = {d:.1} km\n", .{rk4R});
    std.debug.print("  Diff:    {d:.1} km\n", .{@abs(sgp4R - rk4R)});
    std.debug.print("SGP4 models J2 and drag, two-body doesn't, so they drift apart.\n", .{});

    std.debug.print("\n=== Satellite (picks SGP4 or SDP4) ===\n", .{});
    const leo = try Satellite.init(tle, constants.wgs84);
    std.debug.print("ISS (LEO): deep_space={}\n", .{leo.isDeepSpace()});
    std.debug.print("  t=60 min: r={d:.1} km\n", .{mag((try leo.propagate(60))[0])});

    const gpsTleStr =
        \\1 20413U 90005A   24186.00000000  .00000012  00000+0  10000-3 0  9992
        \\2 20413  55.4408  61.4858 0112981 129.5765 231.5553  2.00561730104446
    ;
    var gpsTle = try Tle.parse(gpsTleStr, allocator);
    defer gpsTle.deinit();

    const gps = try Satellite.init(gpsTle, constants.wgs72);
    std.debug.print("GPS (MEO): deep_space={}\n", .{gps.isDeepSpace()});
    std.debug.print("  t=720 min: r={d:.1} km\n", .{mag((try gps.propagate(720))[0])});
}

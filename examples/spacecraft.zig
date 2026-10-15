//! Spacecraft maneuvers and attitude: impulse burns, a plane change and a phase
//! change applied during numerical propagation, and quaternion attitude
//! propagation under a body rate.

const std = @import("std");
const astroz = @import("astroz");
const Tle = astroz.Tle;
const constants = astroz.constants;
const Spacecraft = astroz.Spacecraft;
const Impulse = Spacecraft.Impulse;
const pi = std.math.pi;

// Hide the library's per-maneuver debug logging
pub const std_options: std.Options = .{ .log_level = .info };

const tleText =
    \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
    \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
;

// Impulse times are seconds after the start of propagation
const runs = [_]struct { name: []const u8, impulses: []const Impulse }{
    .{ .name = "no maneuvers", .impulses = &.{} },
    .{ .name = "three burns (km/s)", .impulses = &.{
        .{ .time = 3600.0, .maneuver = .{ .absolute = .{ 0.05, 0.03, 0.01 } } },
        .{ .time = 7200.0, .maneuver = .{ .absolute = .{ 0.011, -0.05, 0.02 } } },
        .{ .time = 10800.0, .maneuver = .{ .absolute = .{ -0.03, 0.08, -0.01 } } },
    } },
    .{ .name = "plane change +10 deg", .impulses = &.{
        .{ .time = 11 * 3600.0, .maneuver = .{ .planeChange = .{ .deltaInclination = pi / 18.0, .deltaRaan = pi / 36.0 } } },
    } },
    .{ .name = "phase change 90 deg", .impulses = &.{
        .{ .time = 11 * 3600.0, .maneuver = .{ .phase = .{ .angle = pi / 2.0, .orbits = 1.0 } } },
    } },
};

pub fn main(init: std.process.Init) !void {
    const gpa = init.gpa;
    var tle = try Tle.parse(tleText, gpa);
    defer tle.deinit();

    std.debug.print("Orbit after one day:\n", .{});
    std.debug.print("  {s:<22} {s:>10} {s:>8} {s:>9} {s:>9}\n", .{ "run", "a (km)", "e", "inc (deg)", "RAAN (deg)" });
    for (runs) |run| {
        var sc = Spacecraft.init("sc", tle, 300.0, .Cube, constants.earth, gpa);
        defer sc.deinit();
        try sc.propagate(sc.tle.epoch, 1, 1, run.impulses); // 1 day, 1 s steps
        const last = sc.orbitPredictions.items[sc.orbitPredictions.items.len - 1].state;
        const el = astroz.calculations.stateVectorToOrbitalElements(last[0..3].*, last[3..6].*, constants.earth.mu);
        std.debug.print("  {s:<22} {d:10.1} {d:8.5} {d:9.3} {d:9.3}\n", .{ run.name, el.a, el.e, el.i * constants.rad2deg, el.raan * constants.rad2deg });
    }

    // Attitude: an illustrative slowly varying body rate, integrated for 6 hours
    var sc = Spacecraft.init("sc", tle, 300.0, .Cube, constants.earth, gpa);
    defer sc.deinit();
    const dt = 60.0;
    std.debug.print("\nAttitude (quaternion w, x, y, z):\n", .{});
    var t: f64 = 0;
    while (t <= 6 * 3600.0) : (t += dt) {
        sc.angularVelocity = .{ 0.001 * @sin(t / 3000.0), 0.0005 * @cos(t / 5000.0), 0.0002 };
        sc.updateAttitude();
        sc.propagateAttitude(dt);
        if (@mod(t, 3600.0) == 0) {
            const q = sc.quaternion;
            std.debug.print("  {d:2.0} h  [{d:7.4} {d:7.4} {d:7.4} {d:7.4}]\n", .{ t / 3600.0, q[0], q[1], q[2], q[3] });
        }
    }
}

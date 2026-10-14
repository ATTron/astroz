const std = @import("std");
const math = std.math;

const astroz = @import("astroz");
const Tle = astroz.Tle;
const constants = astroz.constants;
const Spacecraft = astroz.Spacecraft;
const Impulse = Spacecraft.Impulse;

const testTle =
    \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
    \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
;

// Impulse times are seconds after the start of propagation
const runs = [_]struct { name: []const u8, impulses: []const Impulse }{
    .{ .name = "No maneuvers", .impulses = &.{} },
    .{ .name = "Three absolute burns", .impulses = &.{
        .{ .time = 3600.0, .maneuver = .{ .absolute = .{ 0.05, 0.03, 0.01 } } },
        .{ .time = 7200.0, .maneuver = .{ .absolute = .{ 0.011, -0.05, 0.02 } } },
        .{ .time = 10800.0, .maneuver = .{ .absolute = .{ -0.03, 0.08, -0.01 } } },
    } },
    .{ .name = "Plane change: 10 deg inclination, 5 deg RAAN", .impulses = &.{
        .{ .time = 11 * 3600.0, .maneuver = .{ .planeChange = .{ .deltaInclination = math.pi / 18.0, .deltaRaan = math.pi / 36.0 } } },
    } },
    .{ .name = "Phase change: 90 deg over one transfer orbit", .impulses = &.{
        .{ .time = 11 * 3600.0, .maneuver = .{ .phase = .{ .angle = math.pi / 2.0, .orbits = 1.0 } } },
    } },
};

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;

    var tle = try Tle.parse(testTle, allocator);
    defer tle.deinit();

    for (runs) |run| {
        var sc = Spacecraft.init("sc", tle, 300.0, Spacecraft.SatelliteSize.Cube, constants.earth, allocator);
        defer sc.deinit();

        // 3 days in 1 second steps
        try sc.propagate(sc.tle.epoch, 3, 1, run.impulses);
        std.debug.print("{s}: {d} states\n", .{ run.name, sc.orbitPredictions.items.len });
    }
}

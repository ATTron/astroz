const std = @import("std");
const math = std.math;

const astroz = @import("astroz");
const Tle = astroz.Tle;
const constants = astroz.constants;
const Spacecraft = astroz.Spacecraft;
const Impulse = Spacecraft.Impulse;

pub fn main() !void {
    var dbga = std.heap.DebugAllocator(.{}).init;
    defer _ = dbga.deinit();
    const allocator = dbga.allocator();

    const testTle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;

    // Example 1: Basic Orbit Propagation (no maneuvers)
    std.debug.print("=== Example 1: Basic Orbit Propagation ===\n", .{});
    var tle1 = try Tle.parse(testTle, allocator);
    defer tle1.deinit();

    var sc1 = Spacecraft.init("basic_sc", tle1, 300.0, Spacecraft.SatelliteSize.Cube, constants.earth, allocator);
    defer sc1.deinit();

    try sc1.propagate(sc1.tle.firstLine.epoch, 3, 1, null);
    std.debug.print("Generated {d} orbit predictions\n", .{sc1.orbitPredictions.items.len});

    // Example 2: Orbit Propagation with Impulse Maneuvers
    std.debug.print("\n=== Example 2: Impulse Maneuvers ===\n", .{});
    var tle2 = try Tle.parse(testTle, allocator);
    defer tle2.deinit();

    var sc2 = Spacecraft.init("impulse_sc", tle2, 300.0, Spacecraft.SatelliteSize.Cube, constants.earth, allocator);
    defer sc2.deinit();

    const impulses = [_]Impulse{
        .{ .time = 3600.0, .maneuver = .{ .absolute = .{ 0.05, 0.03, 0.01 } } },
        .{ .time = 7200.0, .maneuver = .{ .absolute = .{ 1.1, -0.05, 0.02 } } },
        .{ .time = 10800.0, .maneuver = .{ .absolute = .{ -0.03, 0.08, -0.01 } } },
    };

    try sc2.propagate(sc2.tle.firstLine.epoch, 3, 1, &impulses);
    std.debug.print("Applied {d} impulse maneuvers\n", .{impulses.len});

    // Example 3: Plane Change Maneuver
    std.debug.print("\n=== Example 3: Plane Change Maneuver ===\n", .{});
    var tle3 = try Tle.parse(testTle, allocator);
    defer tle3.deinit();

    var sc3 = Spacecraft.init("plane_change_sc", tle3, 300.0, Spacecraft.SatelliteSize.Cube, constants.earth, allocator);
    defer sc3.deinit();

    const planeChangeManeuver = Impulse{
        .time = 2500000.0,
        .maneuver = .{
            .planeChange = .{
                .deltaInclination = math.pi / 18.0, // 10-degree inclination change
                .deltaRaan = math.pi / 36.0, // 5-degree RAAN change
            },
        },
    };

    const planeImpulses = [_]Impulse{planeChangeManeuver};
    try sc3.propagate(sc3.tle.firstLine.epoch, 3, 1, &planeImpulses);
    std.debug.print("Applied plane change: Δi={d:.1}°, ΔRAAN={d:.1}°\n", .{
        math.degreesToRadians(10.0),
        math.degreesToRadians(5.0),
    });

    // Example 4: Phase Change Maneuver
    std.debug.print("\n=== Example 4: Phase Change Maneuver ===\n", .{});
    var tle4 = try Tle.parse(testTle, allocator);
    defer tle4.deinit();

    var sc4 = Spacecraft.init("phase_change_sc", tle4, 300.0, Spacecraft.SatelliteSize.Cube, constants.earth, allocator);
    defer sc4.deinit();

    const phaseManeuver = Impulse{
        .time = 2500000.0,
        .maneuver = .{
            .phase = .{
                .angle = math.pi / 2.0, // 90-degree phase change
                .orbits = 1.0, // complete in 1 transfer orbit
            },
        },
    };

    const phaseImpulses = [_]Impulse{phaseManeuver};
    try sc4.propagate(sc4.tle.firstLine.epoch, 3, 1, &phaseImpulses);
    std.debug.print("Applied phase change: {d:.1}° (π/2 radians)\n", .{90.0});

    std.debug.print("\nAll maneuver examples completed successfully!\n", .{});
}

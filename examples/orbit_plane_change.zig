const std = @import("std");
const math = std.math;
const astroz = @import("astroz");
const Tle = astroz.Tle;
const constants = astroz.constants;
const Impulse = Spacecraft.Impulse;
const Spacecraft = astroz.Spacecraft;

pub fn main() !void {
    var dbga = std.heap.DebugAllocator(.{}).init;
    defer _ = dbga.deinit();
    const allocator = dbga.allocator();

    const testTle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;

    var tle = try Tle.parse(testTle, allocator);
    defer tle.deinit();

    var testSc = Spacecraft.init("dummy_sc", tle, 300.000, Spacecraft.SatelliteSize.Cube, constants.earth, allocator);
    defer testSc.deinit();

    const planeChangeManeuver = Impulse{
        .time = 2500000.0,
        .deltaV = .{ 0.0, 0.0, 0.0 }, // Not used for plane changes
        .mode = .PlaneChange,
        .planeChange = .{
            .deltaInclination = math.pi / 18.0, // 10-degree inclination change
            .deltaRaan = math.pi / 36.0, // 5-degree RAAN change
        },
    };

    const impulses = [_]Impulse{planeChangeManeuver};

    try testSc.propagate(
        testSc.tle.firstLine.epoch,
        3, // 3 days worth of orbit predictions
        1, // steps, i.e. repredict every simulated second
        &impulses,
    );

    for (testSc.orbitPredictions.items) |iter| {
        const r = math.sqrt(iter.state[0] * iter.state[0] + iter.state[1] * iter.state[1] + iter.state[2] * iter.state[2]);

        std.debug.print("Next Prediction is: {any}\n", .{r});
    }
}

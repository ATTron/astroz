const std = @import("std");
const math = std.math;
const astroz = @import("astroz");
const Tle = astroz.Tle;
const constants = astroz.constants;
const Spacecraft = astroz.Spacecraft;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const test_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;

    var tle = try Tle.parse(test_tle, allocator);
    defer tle.deinit();

    var test_sc = Spacecraft.init("dummy_sc", tle, 300.000, Spacecraft.SatelliteSize.Cube, constants.earth, allocator);
    defer test_sc.deinit();

    try test_sc.propagate(
        test_sc.tle.first_line.epoch,
        3, // 3 days worth of orbit predictions
        1, // steps, i.e. repredict every simulated second
        null,
    );

    for (test_sc.orbit_predictions.items) |iter| {
        const r = math.sqrt(iter.state[0] * iter.state[0] + iter.state[1] * iter.state[1] + iter.state[2] * iter.state[2]);

        std.debug.print("Next Prediction is: {any}\n", .{r});
    }
}

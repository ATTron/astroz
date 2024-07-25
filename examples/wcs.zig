const std = @import("std");
const astroz = @import("astroz");
const Tle = astroz.Tle;
const WorldCoordinateSystem = astroz.WorldCoordinateSystem;

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
    const wcs = WorldCoordinateSystem.fromTle(test_tle, 0.0);

    std.log.debug("WCS OUTPUT: {any}", wcs);

    tle.output();
}

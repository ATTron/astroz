const std = @import("std");
const astroz = @import("astroz");
const Tle = astroz.Tle;
const WorldCoordinateSystem = astroz.WorldCoordinateSystem;

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
    const wcs = WorldCoordinateSystem.fromTle(tle, 0.0, astroz.constants.earth);

    std.log.debug("WCS OUTPUT: {any}", .{wcs});

    tle.output();
}

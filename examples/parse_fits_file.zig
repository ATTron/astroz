const std = @import("std");
const astroz = @import("astroz");
const Fits = astroz.Fits;

pub fn main() !void {
    var dbga: std.heap.DebugAllocator(.{}) = .init;
    defer _ = dbga.deinit();
    const allocator = dbga.allocator();

    var fitsPng: Fits = try .open_and_parse("test/sample_fits.fits", allocator, .{ .createImages = true, .stretchOptions = .{ .stretch = 0.2 }});
    defer fitsPng.close();
}

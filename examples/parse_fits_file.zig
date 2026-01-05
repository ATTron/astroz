const std = @import("std");
const astroz = @import("astroz");
const Fits = astroz.Fits;

pub fn main() !void {
    var dbga: std.heap.DebugAllocator(.{}) = .init;
    defer _ = dbga.deinit();
    const allocator = dbga.allocator();
    const io = std.Io.Threaded.global_single_threaded.ioBasic();

    var fitsPng: Fits = try .open_and_parse("test/sample_fits.fits", io, allocator, .{ .createImages = true, .stretchOptions = .{ .stretch = 0.2 } });
    defer fitsPng.close();
}

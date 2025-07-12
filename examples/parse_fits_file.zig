const std = @import("std");
const astroz = @import("astroz");
const Fits = astroz.Fits;

pub fn main() !void {
    var dbga = std.heap.DebugAllocator(.{}).init;
    defer _ = dbga.deinit();
    const allocator = dbga.allocator();

    var fitsPng = try Fits.open_and_parse(allocator);
    defer fitsPng.close();
}

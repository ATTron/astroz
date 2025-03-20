const std = @import("std");
const astroz = @import("astroz");
const Fits = astroz.Fits;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    var fitsPng = try Fits.open_and_parse(allocator);
    defer fitsPng.close();
}

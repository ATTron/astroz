const std = @import("std");
const astroz = @import("astroz");
const Fits = astroz.Fits;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    var fits_png = Fits.init(allocator);
    try fits_png.toImage(
        "test/sample_fits.fits",
        "test/test.png",
        .{},
    );
}

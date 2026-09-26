const std = @import("std");
const astroz = @import("astroz");
const Fits = astroz.Fits;

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    const io = init.io;

    var fitsPng: Fits = try .open_and_parse("test/sample_fits.fits", io, allocator, .{ .createImages = true, .stretchOptions = .{ .stretch = 0.2 } });
    defer fitsPng.close();
}

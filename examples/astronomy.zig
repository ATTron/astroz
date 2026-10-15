//! Astronomy helpers: precess a star's equatorial coordinates to a new date, and
//! read a FITS file, writing its images out as PNGs (into sample_fits/).

const std = @import("std");
const astroz = @import("astroz");
const Ecs = astroz.EquatorialCoordinateSystem;

pub fn main(init: std.process.Init) !void {
    const star = Ecs.init(.init(40, 10, 10), .init(19, 52, 2));
    std.log.info("J2000:      {any}", .{star});
    std.log.info("2005-07-30: {any}", .{star.precess(astroz.Datetime.initDate(2005, 7, 30))});

    var fits: astroz.Fits = try .open_and_parse("test/sample_fits.fits", init.io, init.gpa, .{
        .createImages = true,
        .stretchOptions = .{ .stretch = 0.2 },
    });
    defer fits.close();
}

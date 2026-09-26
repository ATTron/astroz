const std = @import("std");
const astroz = @import("astroz");
const Ecs = astroz.EquatorialCoordinateSystem;

pub fn main() void {
    const star = Ecs.init(.init(40, 10, 10), .init(19, 52, 2));
    const moved = star.precess(astroz.Datetime.initDate(2005, 7, 30));

    std.log.info("J2000:      {any}", .{star});
    std.log.info("2005-07-30: {any}", .{moved});
}

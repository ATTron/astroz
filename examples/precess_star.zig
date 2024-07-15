const std = @import("std");
const astroz = @import("astroz");
const Ecs = astroz.Ecs;
const Datetime = astroz.Datetime;

pub fn main() !void {
    const declination = Ecs.Declination.init(40, 10, 10);
    const ra = Ecs.RightAscension.init(19, 52, 2);
    const j2000 = Ecs.init(declination, ra);

    std.debug.print("Precessed to July 30, 2005:\n{any}", .{j2000.precess(Datetime.initDate(2005, 7, 30))});
}

const std = @import("std");
const astroz = @import("astroz");
const Datetime = astroz.Datetime;
const EquatorialCoordinateSystem = astroz.EquatorialCoordinateSystem;

pub fn main() !void {
    const declination = EquatorialCoordinateSystem.Declination.init(40, 10, 10);
    const ra = EquatorialCoordinateSystem.RightAscension.init(19, 52, 2);
    const j2000 = EquatorialCoordinateSystem.init(declination, ra);

    std.debug.print("Precessed to July 30, 2005:\n{any}", .{j2000.precess(Datetime.initDate(2005, 7, 30))});
}

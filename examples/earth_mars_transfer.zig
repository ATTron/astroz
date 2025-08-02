const std = @import("std");
const math = std.math;
const astroz = @import("astroz");
const constants = astroz.constants;
const calculations = astroz.calculations;
const OrbitalMechanics = astroz.OrbitalMechanics;
const Mission = astroz.Mission;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    var orbitalMechanics = OrbitalMechanics.init(constants.sun.mu, constants.sun);

    _ = Mission{
        .allocator = allocator,
        .parameters = undefined,
        .orbitalMechanics = orbitalMechanics,
    };

    const earthRadius = constants.earth.semiMajorAxis;
    const marsRadius = constants.mars.semiMajorAxis;

    const transfer = try orbitalMechanics.hohmannTransfer(earthRadius, marsRadius);

    std.debug.print("Earth-Mars Hohmann Transfer:\n");
    std.debug.print("Total Delta-V: {d} km/s\n", .{transfer.totalDeltaV});
    std.debug.print("Transfer Time: {d} days\n", .{transfer.transferTimeDays});
    std.debug.print("Semi-Major Axis: {d} km\n", .{transfer.semiMajorAxis});

    const file = try std.fs.cwd().createFile("test/earth_mars_transfer.csv", .{});
    defer file.close();

    const writer = file.writer();
    try writer.print("body,time_days,x_km,y_km,z_km,label\n");

    const totalDays = transfer.transferTimeDays * 3.0;
    const timeStep = 1.0;

    var day: f64 = 0.0;
    while (day <= totalDays) : (day += timeStep) {
        const earthAngle = (day / constants.earth.period) * 2.0 * math.pi;
        const earthX = earthRadius * @cos(earthAngle);
        const earthY = earthRadius * @sin(earthAngle);
        try writer.print("Earth,{d:.1},{d:.1},{d:.1},0.0,planet\n", .{ day, earthX, earthY });

        const marsAngle = (day / constants.mars.period) * 2.0 * math.pi;
        const marsX = marsRadius * @cos(marsAngle);
        const marsY = marsRadius * @sin(marsAngle);
        try writer.print("Mars,{d:.1},{d:.1},{d:.1},0.0,planet\n", .{ day, marsX, marsY });

        if (day <= transfer.transferTimeDays) {
            const transferAngle = (day / transfer.transferTimeDays) * math.pi;
            const a = transfer.semiMajorAxis;
            const e = (marsRadius - earthRadius) / (marsRadius + earthRadius);
            const r = a * (1.0 - e * e) / (1.0 + e * @cos(transferAngle));

            const transferX = r * @cos(transferAngle);
            const transferY = r * @sin(transferAngle);
            try writer.print("Transfer,{d:.1},{d:.1},{d:.1},0.0,trajectory\n", .{ day, transferX, transferY });
        }

        day += timeStep;
    }

    try writer.print("Earth_Departure,0.0,{d:.1},0.0,0.0,waypoint\n", .{earthRadius});

    const arrivalAngle = (transfer.transferTimeDays / constants.mars.period) * 2.0 * math.pi;
    const marsArrivalX = marsRadius * @cos(arrivalAngle);
    const marsArrivalY = marsRadius * @sin(arrivalAngle);
    try writer.print("Mars_Arrival,{d:.1},{d:.1},{d:.1},0.0,waypoint\n", .{ transfer.transferTimeDays, marsArrivalX, marsArrivalY });
}

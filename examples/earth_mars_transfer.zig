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

    // Initialize orbital mechanics with Sun's gravitational parameter
    var orbitalMechanics = OrbitalMechanics.init(constants.sun.mu, constants.sun);

    _ = Mission{
        .allocator = allocator,
        .parameters = undefined,
        .orbitalMechanics = orbitalMechanics,
    };

    // Calculate Earth-Mars Hohmann transfer
    const earthRadius = constants.earth.semiMajorAxis; // km
    const marsRadius = constants.mars.semiMajorAxis; // km

    const transfer = try orbitalMechanics.hohmannTransfer(earthRadius, marsRadius);

    std.debug.print("Earth-Mars Hohmann Transfer:\n");
    std.debug.print("Total Delta-V: {d:.2} km/s\n", .{transfer.totalDeltaV});
    std.debug.print("Transfer Time: {d:.1} days\n", .{transfer.transferTimeDays});
    std.debug.print("Semi-Major Axis: {d:.0} km\n", .{transfer.semiMajorAxis});

    // Create CSV file for visualization
    const file = try std.fs.cwd().createFile("test/earth_mars_transfer.csv", .{});
    defer file.close();

    const writer = file.writer();
    try writer.print("body,time_days,x_km,y_km,z_km,label\n");

    // Time span: 3 times the transfer period to show full context
    const totalDays = transfer.transferTimeDays * 3.0;
    const timeStep = 1.0; // 1 day steps

    var day: f64 = 0.0;
    while (day <= totalDays) : (day += timeStep) {
        // Earth position (circular orbit)
        const earthAngle = (day / constants.earth.period) * 2.0 * math.pi;
        const earthX = earthRadius * @cos(earthAngle);
        const earthY = earthRadius * @sin(earthAngle);
        try writer.print("Earth,{d:.1},{d:.1},{d:.1},0.0,planet\n", .{ day, earthX, earthY });

        // Mars position (circular orbit)
        const marsAngle = (day / constants.mars.period) * 2.0 * math.pi;
        const marsX = marsRadius * @cos(marsAngle);
        const marsY = marsRadius * @sin(marsAngle);
        try writer.print("Mars,{d:.1},{d:.1},{d:.1},0.0,planet\n", .{ day, marsX, marsY });

        // Transfer trajectory (only during transfer period)
        if (day <= transfer.transferTimeDays) {
            // Parametric equations for elliptical transfer orbit
            const transferAngle = (day / transfer.transferTimeDays) * math.pi; // 0 to π
            const a = transfer.semiMajorAxis;
            const e = (marsRadius - earthRadius) / (marsRadius + earthRadius);
            const r = a * (1.0 - e * e) / (1.0 + e * @cos(transferAngle));

            const transferX = r * @cos(transferAngle);
            const transferY = r * @sin(transferAngle);
            try writer.print("Transfer,{d:.1},{d:.1},{d:.1},0.0,trajectory\n", .{ day, transferX, transferY });
        }

        day += timeStep;
    }

    // Add some key waypoints
    try writer.print("Earth_Departure,0.0,{d:.1},0.0,0.0,waypoint\n", .{earthRadius});

    const arrivalAngle = (transfer.transferTimeDays / constants.mars.period) * 2.0 * math.pi;
    const marsArrivalX = marsRadius * @cos(arrivalAngle);
    const marsArrivalY = marsRadius * @sin(arrivalAngle);
    try writer.print("Mars_Arrival,{d:.1},{d:.1},{d:.1},0.0,waypoint\n", .{ transfer.transferTimeDays, marsArrivalX, marsArrivalY });

    std.debug.print("\nVisualization data written to: test/earth_mars_transfer.csv\n");
    std.debug.print("Run Python visualization script to generate animation.\n");
}

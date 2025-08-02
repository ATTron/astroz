const std = @import("std");
const astroz = @import("astroz");
const Mission = astroz.Mission;
const OrbitalMechanics = astroz.OrbitalMechanics;
const constants = astroz.constants;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const orbitalMechanics = OrbitalMechanics.init(constants.sun.mu, constants.sun);
    std.debug.print("=== Earth to Mars Transfer ===\n", .{});
    const earthMarsParams = try Mission.MissionParameters.init(constants.earth, constants.mars, 0.0, null, "hohmann");

    var earthMarsMission = Mission.init(allocator, earthMarsParams, orbitalMechanics);
    defer earthMarsMission.deinit();

    try earthMarsMission.propagateTransfer(520.0, 5.0);
    std.debug.print("Generated {d} trajectory points\n", .{earthMarsMission.trajectoryPredictions.items.len});
    std.debug.print("Waypoints:\n", .{});
    for (earthMarsMission.trajectoryPredictions.items) |point| {
        if (std.mem.eql(u8, point.label, "waypoint")) {
            std.debug.print("Day {d:>6.1}: {s} at ({d:>8.0}, {d:>8.0}) km\n", .{
                point.time,
                point.body,
                point.position.x(),
                point.position.y(),
            });
        }
    }

    std.debug.print("\n=== Mars to Jupiter Transfer ===\n", .{});
    const marsJupiterParams = try Mission.MissionParameters.init(constants.mars, constants.jupiter, 0.0, null, "hohmann");
    var marsJupiterMission = Mission.init(allocator, marsJupiterParams, orbitalMechanics);
    defer marsJupiterMission.deinit();

    try marsJupiterMission.propagateTransfer(1200.0, 10.0);

    std.debug.print("Generated {d} trajectory points\n", .{marsJupiterMission.trajectoryPredictions.items.len});
    std.debug.print("Waypoints:\n", .{});
    for (marsJupiterMission.trajectoryPredictions.items) |point| {
        if (std.mem.eql(u8, point.label, "waypoint")) {
            std.debug.print("Day {d:>6.1}: {s} at ({d:>8.0}, {d:>8.0}) km\n", .{
                point.time,
                point.body,
                point.position.x(),
                point.position.y(),
            });
        }
    }

    std.debug.print("\n=== Venus to Earth Transfer ===\n", .{});
    const venusEarthParams = try Mission.MissionParameters.init(constants.venus, constants.earth, 0.0, null, "hohmann");

    var venusEarthMission = Mission.init(allocator, venusEarthParams, orbitalMechanics);
    defer venusEarthMission.deinit();

    try venusEarthMission.propagateTransfer(400.0, 2.0);
    std.debug.print("Generated {d} trajectory points\n", .{venusEarthMission.trajectoryPredictions.items.len});
    std.debug.print("Waypoints:\n", .{});
    for (venusEarthMission.trajectoryPredictions.items) |point| {
        if (std.mem.eql(u8, point.label, "waypoint")) {
            std.debug.print("Day {d:>6.1}: {s} at ({d:>8.0}, {d:>8.0}) km\n", .{
                point.time,
                point.body,
                point.position.x(),
                point.position.y(),
            });
        }
    }

    std.debug.print("\n=== Programmatic Access Example ===\n", .{});
    std.debug.print("First 5 Earth-Mars transfer points:\n", .{});
    std.debug.print("Time (days) | Body | X (km) | Y (km) | Label\n", .{});
    std.debug.print("-----------|------|--------|--------|----------\n", .{});

    var count: usize = 0;
    for (earthMarsMission.trajectoryPredictions.items) |point| {
        if (count < 5) {
            std.debug.print("{d:>10.1} | {s:<8} | {d:>8.0} | {d:>8.0} | {s}\n", .{
                point.time,
                point.body,
                point.position.x(),
                point.position.y(),
                point.label,
            });
            count += 1;
        }
    }

    std.debug.print("\nTransfer trajectory points (first 5):\n", .{});
    count = 0;
    for (earthMarsMission.trajectoryPredictions.items) |point| {
        if (std.mem.eql(u8, point.body, "Transfer") and count < 5) {
            std.debug.print("Day {d:>6.1}: Transfer at ({d:>8.0}, {d:>8.0}) km\n", .{
                point.time,
                point.position.x(),
                point.position.y(),
            });
            count += 1;
        }
    }
}

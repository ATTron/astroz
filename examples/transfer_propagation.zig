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

    // Example 1: Mission Planning - Compare Hohmann vs Bi-Elliptic
    std.debug.print("=== Mission Planning: Hohmann vs Bi-Elliptic ===\n", .{});
    const hohmannParams = try Mission.MissionParameters.init(constants.earth, constants.mars, 0.0, null, "hohmann");
    var hohmannMission = Mission.init(allocator, hohmannParams, orbitalMechanics);
    defer hohmannMission.deinit();

    const hohmannPlan = try hohmannMission.planMission(hohmannParams);
    std.debug.print("Hohmann Transfer:\n", .{});
    std.debug.print("  Delta-V: {d:.2} km/s\n", .{hohmannPlan.transfer.hohmann.totalDeltaV});
    std.debug.print("  Transfer Time: {d:.1} days\n", .{hohmannPlan.transfer.hohmann.transferTimeDays});
    std.debug.print("  Synodic Period: {d:.1} days\n", .{hohmannPlan.synodicPeriodDays});

    const biEllipticParams = try Mission.MissionParameters.init(constants.earth, constants.mars, 0.0, null, "bi_elliptic");
    var biEllipticMission = Mission.init(allocator, biEllipticParams, orbitalMechanics);
    defer biEllipticMission.deinit();

    const biEllipticPlan = try biEllipticMission.planMission(biEllipticParams);
    std.debug.print("\nBi-Elliptic Transfer:\n", .{});
    std.debug.print("  Delta-V: {d:.2} km/s\n", .{biEllipticPlan.transfer.bi_elliptic.totalDeltaV});
    std.debug.print("  Transfer Time: {d:.1} days\n", .{biEllipticPlan.transfer.bi_elliptic.totalTimeDays});
    std.debug.print("  Synodic Period: {d:.1} days\n", .{biEllipticPlan.synodicPeriodDays});

    // Example 2: Transfer Propagation - Earth to Mars
    std.debug.print("\n=== Earth to Mars Transfer Propagation ===\n", .{});
    const earthMarsParams = try Mission.MissionParameters.init(constants.earth, constants.mars, 0.0, null, "hohmann");
    var earthMarsMission = Mission.init(allocator, earthMarsParams, orbitalMechanics);
    defer earthMarsMission.deinit();

    try earthMarsMission.propagateTransfer(520.0, 5.0);
    std.debug.print("Generated {d} trajectory points\n", .{earthMarsMission.trajectoryPredictions.items.len});
    std.debug.print("Waypoints:\n", .{});
    earthMarsMission.printWaypoints(null);

    // Example 3: Programmatic Access to Trajectory Data
    std.debug.print("\n=== Programmatic Access Example ===\n", .{});
    std.debug.print("First 5 trajectory points:\n", .{});
    earthMarsMission.printTrajectories(5);
}

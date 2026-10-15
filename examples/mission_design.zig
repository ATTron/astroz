//! Interplanetary mission design: Hohmann vs bi-elliptic Earth-Mars transfers,
//! the transfer trajectory, and Monte Carlo dispersion of the Hohmann transfer.

const std = @import("std");
const astroz = @import("astroz");
const Mission = astroz.Mission;
const MonteCarlo = astroz.MonteCarlo;
const constants = astroz.constants;

// Hide the library's progress logging
pub const std_options: std.Options = .{ .log_level = .info };

pub fn main(init: std.process.Init) !void {
    const gpa = init.gpa;
    const sunGravity = astroz.OrbitalMechanics.init(constants.sun.mu);

    std.debug.print("Earth to Mars:\n", .{});
    for ([_][]const u8{ "hohmann", "biElliptic" }) |kind| {
        const params = try Mission.MissionParameters.init(constants.earth, constants.mars, 0.0, null, kind);
        var mission = Mission.init(gpa, params, sunGravity);
        defer mission.deinit();
        const plan = try mission.planMission(params);
        const dv, const days = switch (plan.transfer) {
            .hohmann => |h| .{ h.totalDeltaV, h.transferTimeDays },
            .biElliptic => |b| .{ b.totalDeltaV, b.totalTimeDays },
        };
        std.debug.print("  {s:<11} delta-v {d:.2} km/s, {d:.0} days (synodic period {d:.0} days)\n", .{ kind, dv, days, plan.synodicPeriodDays });
    }

    const params = try Mission.MissionParameters.init(constants.earth, constants.mars, 0.0, null, "hohmann");
    var transfer = Mission.init(gpa, params, sunGravity);
    defer transfer.deinit();
    try transfer.propagateTransfer(520.0, 5.0); // 520 days in 5-day steps
    std.debug.print("\nWaypoints:\n", .{});
    transfer.printWaypoints(null);
    std.debug.print("Transfer orbit, every 50 days:\n", .{});
    for (transfer.trajectoryPredictions.items) |point| {
        if (!std.mem.eql(u8, point.label, "trajectory") or @mod(point.time, 50.0) != 0) continue;
        const r = point.position.magnitude() / constants.auKm;
        std.debug.print("  day {d:>3.0}: {d:.3} AU from the Sun\n", .{ point.time, r });
    }

    var mc = MonteCarlo.init(gpa, .{
        .numSimulations = 1000,
        .uncertainty = .{
            .departureRadiusUncertainty = 0.01,
            .arrivalRadiusUncertainty = 0.01,
            .muUncertainty = 0.001,
            .launchWindowDays = 1.0,
            .measurementNoise = 0.001,
        },
        .seed = 42,
        .departureBody = constants.earth,
        .arrivalBody = constants.mars,
        .baseMu = constants.sun.mu,
    });
    defer mc.deinit();
    std.debug.print("\nHohmann transfer under 1% radius and 0.1% mu uncertainty:\n", .{});
    try mc.runHohmannSimulation();
    try mc.printStatistics();
}

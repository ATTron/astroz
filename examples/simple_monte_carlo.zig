const std = @import("std");
const astroz = @import("astroz");
const MonteCarlo = astroz.MonteCarlo;
const constants = astroz.constants;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    std.debug.print("=== Simple Monte Carlo Test ===\n", .{});

    // Basic uncertainty configuration
    const uncertainty = MonteCarlo.UncertaintyParams{
        .departureRadiusUncertainty = 0.01, // 1%
        .arrivalRadiusUncertainty = 0.01,   // 1%
        .muUncertainty = 0.001,              // 0.1%
        .launchWindowDays = 1.0,
        .measurementNoise = 0.001,
    };

    const config = MonteCarlo.MonteCarloConfig{
        .numSimulations = 1000,
        .uncertainty = uncertainty,
        .seed = 42,
        .departureBody = constants.earth,
        .arrivalBody = constants.mars,
        .baseMu = constants.sun.mu,
    };

    var mc = MonteCarlo.init(allocator, config);
    defer mc.deinit();

    std.debug.print("Running {d} Earth-Mars Hohmann transfer simulations...\n", .{config.numSimulations});

    try mc.runHohmannSimulation();

    std.debug.print("Simulation completed!\n", .{});

    try mc.printStatistics();

    std.debug.print("\nFirst 3 simulation results:\n", .{});
    for (mc.results.items[0..@min(3, mc.results.items.len)], 0..) |result, i| {
        std.debug.print("  {d}: deltaV={d:.3} km/s, time={d:.1} days\n", .{
            i + 1,
            result.deltaV,
            result.transferTime,
        });
    }

    std.debug.print("\n=== Monte Carlo Test Complete ===\n", .{});
}

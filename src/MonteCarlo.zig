//! Monte Carlo simulation framework

const std = @import("std");
const log = std.log;
const OrbitalMechanics = @import("OrbitalMechanics.zig");
const Mission = @import("Mission.zig");
const calculations = @import("calculations.zig");
const constants = @import("constants.zig");
const ValidationError = @import("OrbitalMechanics.zig").ValidationError;

const MonteCarlo = @This();

pub const SimulationResult = struct {
    deltaV: f64,
    transferTime: f64,
    semiMajorAxis: f64,
    success: bool,
    departureRadius: f64,
    arrivalRadius: f64,
    uncertainty: UncertaintyParams,
};

pub const StatisticalSummary = struct {
    meanDeltaV: f64,
    stdDeltaV: f64,
    minDeltaV: f64,
    maxDeltaV: f64,
    meanTransferTime: f64,
    stdTransferTime: f64,
    successRate: f64,
    confidence95DeltaV: [2]f64,
    confidence95TransferTime: [2]f64,
    totalSimulations: usize,
    successfulSimulations: usize,
};

pub const UncertaintyParams = struct {
    departureRadiusUncertainty: f64,
    arrivalRadiusUncertainty: f64,
    muUncertainty: f64,
    launchWindowDays: f64,
    measurementNoise: f64,
};

pub const MonteCarloConfig = struct {
    numSimulations: usize,
    uncertainty: UncertaintyParams,
    seed: u64,
    departureBody: constants.CelestialBody,
    arrivalBody: constants.CelestialBody,
    baseMu: f64,
};

allocator: std.mem.Allocator,
config: MonteCarloConfig,
results: std.ArrayList(SimulationResult),
random: std.Random.Xoshiro256,

pub fn init(allocator: std.mem.Allocator, config: MonteCarloConfig) MonteCarlo {
    const prng = std.Random.Xoshiro256.init(config.seed);

    return .{
        .allocator = allocator,
        .config = config,
        .results = std.ArrayList(SimulationResult){},
        .random = prng,
    };
}

pub fn deinit(self: *MonteCarlo) void {
    self.results.deinit(self.allocator);
}

fn generateNormal(self: *MonteCarlo, mean: f64, stdDev: f64) f64 {
    const rand = self.random.random();
    const uniform1 = rand.float(f64);
    const uniform2 = rand.float(f64);

    const z = @sqrt(-2.0 * @log(uniform1)) * @cos(2.0 * std.math.pi * uniform2);
    return mean + stdDev * z;
}

fn generateUniform(self: *MonteCarlo, min: f64, max: f64) f64 {
    const rand = self.random.random();
    return min + (max - min) * rand.float(f64);
}

fn applyUncertainty(self: *MonteCarlo, baseValue: f64, uncertaintyPercent: f64) f64 {
    const stdDev = baseValue * uncertaintyPercent;
    return self.generateNormal(baseValue, stdDev);
}

pub fn runHohmannSimulation(self: *MonteCarlo) !void {
    self.results.clearRetainingCapacity();

    const baseR1 = self.config.departureBody.semiMajorAxis;
    const baseR2 = self.config.arrivalBody.semiMajorAxis;
    const baseMu = self.config.baseMu;

    log.debug("Starting Monte Carlo simulation with {d} iterations", .{self.config.numSimulations});
    log.debug("Base departure radius: {d:.0} km", .{baseR1});
    log.debug("Base arrival radius: {d:.0} km", .{baseR2});
    log.debug("Uncertainties: departure {d:.1}%, arrival {d:.1}%, mu {d:.1}%", .{
        self.config.uncertainty.departureRadiusUncertainty * 100,
        self.config.uncertainty.arrivalRadiusUncertainty * 100,
        self.config.uncertainty.muUncertainty * 100,
    });

    for (0..self.config.numSimulations) |i| {
        const r1 = self.applyUncertainty(baseR1, self.config.uncertainty.departureRadiusUncertainty);
        const r2 = self.applyUncertainty(baseR2, self.config.uncertainty.arrivalRadiusUncertainty);
        const mu = self.applyUncertainty(baseMu, self.config.uncertainty.muUncertainty);

        const safeR1 = @max(r1, baseR1 * 0.5);
        const safeR2 = @max(r2, baseR2 * 0.5);

        var orbitalMechanics = OrbitalMechanics.init(mu, self.config.departureBody);

        const appliedUncertainty = UncertaintyParams{
            .departureRadiusUncertainty = (safeR1 - baseR1) / baseR1,
            .arrivalRadiusUncertainty = (safeR2 - baseR2) / baseR2,
            .muUncertainty = (mu - baseMu) / baseMu,
            .launchWindowDays = self.config.uncertainty.launchWindowDays,
            .measurementNoise = self.config.uncertainty.measurementNoise,
        };

        const transferResult = orbitalMechanics.hohmannTransfer(safeR1, safeR2) catch {
            try self.results.append(self.allocator, SimulationResult{
                .deltaV = 0.0,
                .transferTime = 0.0,
                .semiMajorAxis = 0.0,
                .success = false,
                .departureRadius = safeR1,
                .arrivalRadius = safeR2,
                .uncertainty = appliedUncertainty,
            });
            continue;
        };

        try self.results.append(self.allocator, SimulationResult{
            .deltaV = transferResult.totalDeltaV,
            .transferTime = transferResult.transferTimeDays,
            .semiMajorAxis = transferResult.semiMajorAxis,
            .success = true,
            .departureRadius = safeR1,
            .arrivalRadius = safeR2,
            .uncertainty = appliedUncertainty,
        });

        if ((i + 1) % (self.config.numSimulations / 10) == 0) {
            const progress = @as(f64, @floatFromInt(i + 1)) / @as(f64, @floatFromInt(self.config.numSimulations)) * 100.0;
            log.debug("Progress: {d:.0}% ({d}/{d} simulations)", .{ progress, i + 1, self.config.numSimulations });
        }
    }

    log.info("Monte Carlo simulation completed: {d} total simulations", .{self.results.items.len});
}

/// Calculate statistical summary
pub fn calculateStatistics(self: *MonteCarlo) !StatisticalSummary {
    if (self.results.items.len == 0) {
        return ValidationError.ValueError;
    }

    var successfulResults = std.ArrayList(SimulationResult){};
    defer successfulResults.deinit(self.allocator);

    for (self.results.items) |result| {
        if (result.success) {
            try successfulResults.append(self.allocator, result);
        }
    }

    const successfulCount = successfulResults.items.len;
    if (successfulCount == 0) {
        return ValidationError.ValueError;
    }

    var sumDeltaV: f64 = 0.0;
    var sumTransferTime: f64 = 0.0;
    var minDeltaV: f64 = std.math.inf(f64);
    var maxDeltaV: f64 = -std.math.inf(f64);

    for (successfulResults.items) |result| {
        sumDeltaV += result.deltaV;
        sumTransferTime += result.transferTime;
        minDeltaV = @min(minDeltaV, result.deltaV);
        maxDeltaV = @max(maxDeltaV, result.deltaV);
    }

    const meanDeltaV = sumDeltaV / @as(f64, @floatFromInt(successfulCount));
    const meanTransferTime = sumTransferTime / @as(f64, @floatFromInt(successfulCount));

    var sumSquaredDiffDeltaV: f64 = 0.0;
    var sumSquaredDiffTransferTime: f64 = 0.0;

    for (successfulResults.items) |result| {
        const diffDeltaV = result.deltaV - meanDeltaV;
        const diffTransferTime = result.transferTime - meanTransferTime;
        sumSquaredDiffDeltaV += diffDeltaV * diffDeltaV;
        sumSquaredDiffTransferTime += diffTransferTime * diffTransferTime;
    }

    const stdDeltaV = @sqrt(sumSquaredDiffDeltaV / @as(f64, @floatFromInt(successfulCount - 1)));
    const stdTransferTime = @sqrt(sumSquaredDiffTransferTime / @as(f64, @floatFromInt(successfulCount - 1)));

    const z95 = 1.96;
    const marginDeltaV = z95 * stdDeltaV / @sqrt(@as(f64, @floatFromInt(successfulCount)));
    const marginTransferTime = z95 * stdTransferTime / @sqrt(@as(f64, @floatFromInt(successfulCount)));

    const successRate = @as(f64, @floatFromInt(successfulCount)) / @as(f64, @floatFromInt(self.results.items.len));

    return StatisticalSummary{
        .meanDeltaV = meanDeltaV,
        .stdDeltaV = stdDeltaV,
        .minDeltaV = minDeltaV,
        .maxDeltaV = maxDeltaV,
        .meanTransferTime = meanTransferTime,
        .stdTransferTime = stdTransferTime,
        .successRate = successRate,
        .confidence95DeltaV = .{ meanDeltaV - marginDeltaV, meanDeltaV + marginDeltaV },
        .confidence95TransferTime = .{ meanTransferTime - marginTransferTime, meanTransferTime + marginTransferTime },
        .totalSimulations = self.results.items.len,
        .successfulSimulations = successfulCount,
    };
}

/// Print statistical summary to console
pub fn printStatistics(self: *MonteCarlo) !void {
    const stats = try self.calculateStatistics();

    log.info("\n=== Monte Carlo Statistical Summary ===", .{});
    log.info("Total simulations: {d}", .{stats.totalSimulations});
    log.info("Successful simulations: {d}", .{stats.successfulSimulations});
    log.info("Success rate: {d:.1}%", .{stats.successRate * 100.0});
    log.info("", .{});
    log.info("Delta-V Stats:", .{});
    log.info("  Mean: {d:.3} km/s", .{stats.meanDeltaV});
    log.info("  Std Dev: {d:.3} km/s", .{stats.stdDeltaV});
    log.info("  Min: {d:.3} km/s", .{stats.minDeltaV});
    log.info("  Max: {d:.3} km/s", .{stats.maxDeltaV});
    log.info("  95% Confidence Interval: [{d:.3}, {d:.3}] km/s", .{ stats.confidence95DeltaV[0], stats.confidence95DeltaV[1] });
    log.info("", .{});
    log.info("Transfer Time Statistics:", .{});
    log.info("  Mean: {d:.1} days", .{stats.meanTransferTime});
    log.info("  Std Dev: {d:.1} days", .{stats.stdTransferTime});
    log.info("  95% Confidence Interval: [{d:.1}, {d:.1}] days", .{ stats.confidence95TransferTime[0], stats.confidence95TransferTime[1] });
}

/// Export results to CSV format for Python analysis
pub fn exportToCSV(self: *MonteCarlo, filePath: []const u8) !void {
    const file = try std.fs.cwd().createFile(filePath, .{});
    defer file.close();

    try file.writeAll("simulation,success,deltaV_km_s,transferTime_days,semiMajorAxis_km,departureRadius_km,arrivalRadius_km,departure_uncertainty_pct,arrival_uncertainty_pct,mu_uncertainty_pct\n");
    for (self.results.items, 0..) |result, i| {
        var buffer: [512]u8 = undefined;
        const line = try std.fmt.bufPrint(buffer[0..], "{d},{},{d:.6},{d:.3},{d:.0},{d:.0},{d:.0},{d:.6},{d:.6},{d:.6}\n", .{
            i,
            if (result.success) @as(u8, 1) else @as(u8, 0),
            result.deltaV,
            result.transferTime,
            result.semiMajorAxis,
            result.departureRadius,
            result.arrivalRadius,
            result.uncertainty.departureRadiusUncertainty * 100.0,
            result.uncertainty.arrivalRadiusUncertainty * 100.0,
            result.uncertainty.muUncertainty * 100.0,
        });
        try file.writeAll(line);
    }

    log.info("Results exported to: {s}", .{filePath});
}

test "Monte Carlo basic functionality" {
    const testing = std.testing;
    const ta = std.testing.allocator;

    const uncertainty = UncertaintyParams{
        .departureRadiusUncertainty = 0.01, // 1%
        .arrivalRadiusUncertainty = 0.01,
        .muUncertainty = 0.005,
        .launchWindowDays = 1.0,
        .measurementNoise = 0.001,
    };

    const config = MonteCarloConfig{
        .numSimulations = 100,
        .uncertainty = uncertainty,
        .seed = 12345,
        .departureBody = constants.earth,
        .arrivalBody = constants.mars,
        .baseMu = constants.sun.mu,
    };

    var mc = MonteCarlo.init(ta, config);
    defer mc.deinit();

    try mc.runHohmannSimulation();
    try testing.expect(mc.results.items.len == 100);

    const stats = try mc.calculateStatistics();
    try testing.expect(stats.totalSimulations == 100);
    try testing.expect(stats.successRate > 0.8);
    try testing.expect(stats.meanDeltaV > 0.0);
    try testing.expect(stats.stdDeltaV > 0.0);
}

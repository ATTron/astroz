const std = @import("std");
const astroz = @import("astroz");

const Tle = astroz.Tle;
const Sgp4 = astroz.Sgp4;
const Sgp4Batch = astroz.Constellation.Sgp4Batch;
const Sgp4Constellation = astroz.Constellation;

// ISS TLE
const ISS_TLE =
    \\1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995
    \\2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123
;

// Batch size from the library (4 for AVX2, 8 for AVX512)
const BatchSize = Sgp4Batch.BatchSize;

// scenarios
const scenarios = [_]struct { name: []const u8, points: usize, step_min: f64 }{
    .{ .name = "1 day (minute)", .points = 1440, .step_min = 1.0 },
    .{ .name = "1 week (minute)", .points = 10080, .step_min = 1.0 },
    .{ .name = "2 weeks (minute)", .points = 20160, .step_min = 1.0 },
    .{ .name = "2 weeks (second)", .points = 1209600, .step_min = 1.0 / 60.0 },
    .{ .name = "1 month (minute)", .points = 43200, .step_min = 1.0 },
};

const iterations = 10;

fn nowNs() i64 {
    var ts: std.posix.timespec = undefined;
    _ = std.c.clock_gettime(std.c.CLOCK.MONOTONIC, &ts);
    return @as(i64, ts.sec) * 1_000_000_000 + ts.nsec;
}

// Prevent the optimizer from discarding computed results
var sink: f64 = 0;

fn doNotOptimize(val: f64) void {
    @as(*volatile f64, @ptrCast(@constCast(&sink))).* = val;
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    // Parse TLE
    var tle = try Tle.parse(ISS_TLE, allocator);
    defer tle.deinit();

    // Initialize SGP4
    var sgp4 = try Sgp4.init(tle, astroz.constants.wgs84);

    // Warmup
    for (0..100) |i| {
        const result = try sgp4.propagate(@as(f64, @floatFromInt(i)));
        doNotOptimize(result[0][0]);
    }

    std.debug.print("\nastroz SGP4 Benchmark\n", .{});
    std.debug.print("==================================================\n", .{});
    std.debug.print("SIMD Batch Size: {} ({s})\n", .{
        BatchSize,
        if (BatchSize == 8) "AVX512" else "AVX2/SSE",
    });
    std.debug.print("\n--- Scalar Propagation ---\n", .{});

    for (scenarios) |scenario| {
        const times = try allocator.alloc(f64, scenario.points);
        defer allocator.free(times);

        for (0..scenario.points) |i| {
            times[i] = @as(f64, @floatFromInt(i)) * scenario.step_min;
        }

        var total_ns: u64 = 0;

        for (0..iterations) |_| {
            const start = nowNs();

            for (times) |t| {
                if (sgp4.propagate(t)) |result| {
                    doNotOptimize(result[0][0]);
                } else |_| {}
            }

            const end = nowNs();
            total_ns += @intCast(end - start);
        }

        const avg_ms = @as(f64, @floatFromInt(total_ns)) / @as(f64, @floatFromInt(iterations)) / 1_000_000.0;
        const props_per_sec = @as(f64, @floatFromInt(scenario.points)) / (avg_ms / 1000.0);

        std.debug.print("{s:<25} {d:>10.3} ms  ({d:.2} prop/s)\n", .{
            scenario.name,
            avg_ms,
            props_per_sec,
        });
    }

    std.debug.print("\n--- SIMD Batch{} Propagation ---\n", .{BatchSize});

    for (scenarios) |scenario| {
        const times = try allocator.alloc(f64, scenario.points);
        defer allocator.free(times);

        for (0..scenario.points) |i| {
            times[i] = @as(f64, @floatFromInt(i)) * scenario.step_min;
        }

        var total_ns: u64 = 0;

        for (0..iterations) |_| {
            const start = nowNs();

            // Process 8 at a time (AVX512), then 4 (AVX2), then scalar
            var i: usize = 0;

            // AVX512: 8 times at once
            while (i + 8 <= scenario.points) : (i += 8) {
                if (sgp4.propagateN(8, .{
                    times[i],     times[i + 1], times[i + 2], times[i + 3],
                    times[i + 4], times[i + 5], times[i + 6], times[i + 7],
                })) |result| {
                    doNotOptimize(result[0][0][0]);
                } else |_| {}
            }

            // AVX2: 4 times at once
            while (i + 4 <= scenario.points) : (i += 4) {
                if (sgp4.propagateN(4, .{ times[i], times[i + 1], times[i + 2], times[i + 3] })) |result| {
                    doNotOptimize(result[0][0][0]);
                } else |_| {}
            }

            // Scalar remainder
            while (i < scenario.points) : (i += 1) {
                if (sgp4.propagate(times[i])) |result| {
                    doNotOptimize(result[0][0]);
                } else |_| {}
            }

            const end = nowNs();
            total_ns += @intCast(end - start);
        }

        const avg_ms = @as(f64, @floatFromInt(total_ns)) / @as(f64, @floatFromInt(iterations)) / 1_000_000.0;
        const props_per_sec = @as(f64, @floatFromInt(scenario.points)) / (avg_ms / 1000.0);

        std.debug.print("{s:<25} {d:>10.3} ms  ({d:.2} prop/s)\n", .{
            scenario.name,
            avg_ms,
            props_per_sec,
        });
    }

    // --- Multithreaded Constellation Benchmark ---
    // Simulate a constellation of N satellites (replicated ISS TLE)
    const constellationSizes = [_]usize{ 100, 1000, 5000 };
    const constellationTimes = 1440; // 1 day at 1-minute steps

    std.debug.print("\n--- Multithreaded Constellation (SIMD + Threading) ---\n", .{});
    std.debug.print("Threads: {}\n\n", .{Sgp4Constellation.getMaxThreads()});

    const el = try Sgp4.initElements(tle, astroz.constants.wgs84);

    for (constellationSizes) |numSats| {
        const numBatches = (numSats + BatchSize - 1) / BatchSize;

        // Build batch elements (all identical for benchmarking)
        const batches = try allocator.alloc(Sgp4Batch.BatchElements(BatchSize), numBatches);
        defer allocator.free(batches);

        var tleBatch: [BatchSize]Tle = undefined;
        for (0..BatchSize) |k| tleBatch[k] = tle;

        const batchEl = try Sgp4Batch.initBatchElements(BatchSize, tleBatch, astroz.constants.wgs84);
        for (batches) |*b| b.* = batchEl;

        // Epoch offsets (all zero since same TLE)
        const epochOffsets = try allocator.alloc(f64, numBatches * BatchSize);
        defer allocator.free(epochOffsets);
        @memset(epochOffsets, 0.0);

        // Time array
        const times = try allocator.alloc(f64, constellationTimes);
        defer allocator.free(times);
        for (0..constellationTimes) |k| {
            times[k] = @as(f64, @floatFromInt(k));
        }

        // Output buffer
        const resultsPos = try allocator.alloc(f64, numSats * constellationTimes * 3);
        defer allocator.free(resultsPos);

        // Warmup
        try Sgp4Constellation.propagateConstellation(
            batches, numSats, times, epochOffsets, resultsPos, null,
            .teme, el.epochJd, null, .timeMajor,
        );

        var total_ns: u64 = 0;
        const iters: usize = if (numSats >= 5000) 3 else iterations;

        for (0..iters) |_| {
            const start = nowNs();

            try Sgp4Constellation.propagateConstellation(
                batches, numSats, times, epochOffsets, resultsPos, null,
                .teme, el.epochJd, null, .timeMajor,
            );

            const end = nowNs();
            total_ns += @intCast(end - start);
        }

        doNotOptimize(resultsPos[0]);

        const avg_ms = @as(f64, @floatFromInt(total_ns)) / @as(f64, @floatFromInt(iters)) / 1_000_000.0;
        const total_props = @as(f64, @floatFromInt(numSats)) * @as(f64, @floatFromInt(constellationTimes));
        const props_per_sec = total_props / (avg_ms / 1000.0);

        std.debug.print("{d:>5} sats x {d} times   {d:>10.1} ms  ({d:.2} prop/s)\n", .{
            numSats,
            constellationTimes,
            avg_ms,
            props_per_sec,
        });
    }

    std.debug.print("\n", .{});
}

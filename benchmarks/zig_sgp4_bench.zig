const std = @import("std");
const astroz = @import("astroz");

const Tle = astroz.Tle;
const Sgp4 = astroz.Sgp4;
const Sgp4Batch = astroz.Constellation.Sgp4Batch;
const Sgp4Constellation = astroz.Constellation;
const dispatch = astroz.dispatch;

const ISS_TLE =
    \\1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995
    \\2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123
;

const ScenarioType = enum {
    Scalar,
    SIMD,
    Multithreaded
};

const Scenario = struct {
    name: []const u8,
    points: usize,
    stepMin: f64,
    sceneType: ScenarioType,
};

const BatchSize = Sgp4Batch.BatchSize;

const scenarios = [_]Scenario {
    .{ .name = "1 day (minute)", .points = 1440, .stepMin = 1.0, .sceneType = .Scalar },
    .{ .name = "1 week (minute)", .points = 10080, .stepMin = 1.0, .sceneType = .Scalar },
    .{ .name = "2 weeks (minute)", .points = 20160, .stepMin = 1.0, .sceneType = .Scalar },
    .{ .name = "2 weeks (second)", .points = 1209600, .stepMin = 1.0 / 60.0, .sceneType = .Scalar },
    .{ .name = "1 month (minute)", .points = 43200, .stepMin = 1.0, .sceneType = .Scalar },
};

const scenariosSIMD = [_]Scenario {
    .{ .name = "1 day (minute)", .points = 1440, .stepMin = 1.0, .sceneType = .SIMD },
    .{ .name = "1 week (minute)", .points = 10080, .stepMin = 1.0, .sceneType = .SIMD },
    .{ .name = "2 weeks (minute)", .points = 20160, .stepMin = 1.0, .sceneType = .SIMD },
    .{ .name = "2 weeks (second)", .points = 1209600, .stepMin = 1.0 / 60.0, .sceneType = .SIMD },
    .{ .name = "1 month (minute)", .points = 43200, .stepMin = 1.0, .sceneType = .SIMD },
};

const scenariosMT = [_]Scenario {
    .{ .name = "2 weeks (minute)", .points = 20160, .stepMin = 1.0, .sceneType = .Multithreaded },
    .{ .name = "1 month (second)", .points = 2592000, .stepMin = 1.0 / 60.0, .sceneType = .Multithreaded },
    .{ .name = "3 months (second)", .points = 7776000, .stepMin = 1.0 / 60.0, .sceneType = .Multithreaded },
    .{ .name = "1 year (minute)", .points = 525600, .stepMin = 1.0, .sceneType = .Multithreaded },
    .{ .name = "1 year (second)", .points = 31536000, .stepMin = 1.0 / 60.0, .sceneType = .Multithreaded },
};

const iterations = 10;

fn simdWorker(sgp4: *const Sgp4, times: []const f64) void {
    var i: usize = 0;
    while (i + BatchSize <= times.len) : (i += BatchSize) {
        const result = dispatch.sgp4Times8(sgp4, times[i..][0..BatchSize].*) catch continue;
        std.mem.doNotOptimizeAway(result);
    }
    while (i < times.len) : (i += 1) {
        const result = sgp4.propagate(times[i]) catch continue;
        std.mem.doNotOptimizeAway(result);
    }
}

pub fn main(init: std.process.Init) !void {
    const gpa = init.gpa;
    const io = init.io;

    var tle = try Tle.parse(ISS_TLE, gpa);
    defer tle.deinit();

    var sgp4 = try Sgp4.init(tle, astroz.constants.wgs84);

    for (0..100) |i| {
        _ = try sgp4.propagate(@floatFromInt(i));
    }

    std.debug.print("\nastroz SGP4 Benchmark\n", .{});
    std.debug.print("==================================================\n", .{});
    std.debug.print("SIMD Batch Size: {} (oma runtime dispatch)\n", .{BatchSize});
    std.debug.print("\n--- Scalar Propagation ---\n", .{});

    const num_threads = std.Thread.getCpuCount() catch 1;

    try runScenarios(sgp4, &scenarios, num_threads, gpa, io);

    std.debug.print("\n--- SIMD Batch{} Propagation ---\n", .{BatchSize});

    try runScenarios(sgp4, &scenariosSIMD, num_threads, gpa, io);

    std.debug.print("\n--- Multithreaded SIMD Batch{} Propagation ---\n", .{BatchSize});
    std.debug.print("Threads: {}\n", .{num_threads});


    try runScenarios(sgp4, &scenariosMT, num_threads, gpa, io);

    std.debug.print("\n", .{});
}

fn getResults(name: []const u8, total_ns: usize, iters: usize, points: usize, totalPropsPerSec: *f64) void {
        const avg_ms = @as(f64, @floatFromInt(total_ns / iters)) / 1_000_000.0;
        const props_per_sec = @as(f64, @floatFromInt(points)) / (avg_ms / 1000.0);
        totalPropsPerSec.* += props_per_sec;

        std.debug.print("{s:<25} {d:>10.3} ms  ({d:.2} prop/s)\n", .{
            name,
            avg_ms,
            props_per_sec,
        });
}

fn reportAvg(totalPropsPerSec: f64, scenarioSize: f64) void {
        std.debug.print("{s:<25} {d:>17.2} prop/s\n", .{
            "Average",
            totalPropsPerSec / scenarioSize,
        });
}


fn runScenarios(sgp4: Sgp4, scenarioArr: []const Scenario, num_threads: usize, gpa: std.mem.Allocator, io: std.Io) !void {
    var totalPropsPerSec: f64 = 0.0;

    for (scenarioArr) |scenario| {
        const times = try gpa.alloc(f64, scenario.points);
        defer gpa.free(times);

        for (0..scenario.points) |i| {
            times[i] = @as(f64, @floatFromInt(i)) * scenario.stepMin;
        }

        var total_ns: u64 = 0;

        for (0..iterations) |_| {
            const start = std.Io.Timestamp.now(io, .awake);

            switch (scenario.sceneType) {
                .Scalar => {
                    for (times) |t| {
                        const res = sgp4.propagate(t) catch continue;
                        std.mem.doNotOptimizeAway(res);
                    }
                },
                .SIMD => {
                    var i: usize = 0;
                    while (i + BatchSize <= scenario.points) : (i += BatchSize) {
                        const result = dispatch.sgp4Times8(&sgp4, times[i..][0..BatchSize].*) catch continue;
                        std.mem.doNotOptimizeAway(result);
                    }
                    // scalar remainder
                    while (i < scenario.points) : (i += 1) {
                        const res = sgp4.propagate(times[i]) catch continue;
                        std.mem.doNotOptimizeAway(res);
                    }
                },
                .Multithreaded => {
                    const threads = try gpa.alloc(std.Thread, num_threads);
                    defer gpa.free(threads);

                    const chunk_size = (scenario.points + num_threads - 1) / num_threads;
                    var spawned: usize = 0;

                    for (0..num_threads) |tid| {
                        const begin = tid * chunk_size;
                        const end_idx = @min(begin + chunk_size, scenario.points);
                        if (begin >= end_idx) break;
                        threads[spawned] = std.Thread.spawn(.{}, simdWorker, .{ &sgp4, times[begin..end_idx] }) catch break;
                        spawned += 1;
                    }

                    for (threads[0..spawned]) |thread| {
                        thread.join();
                    }
                },
            }

            const end = std.Io.Timestamp.now(io, .awake);
            total_ns += @intCast(start.durationTo(end).nanoseconds);
        }
        getResults(scenario.name, total_ns, iterations, scenario.points, &totalPropsPerSec);
    }
    reportAvg(totalPropsPerSec, @as(f64, @floatFromInt(scenarioArr.len)));
}

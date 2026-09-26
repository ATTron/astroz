const std = @import("std");
const astroz = @import("astroz");

const Sgp4 = astroz.Sgp4;
const lanes = astroz.Constellation.Sgp4Batch.batchSize;

comptime {
    // dispatch.sgp4Times8 is fixed at 8 lanes.
    std.debug.assert(lanes == 8);
}

const line1 = "1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995";
const line2 = "2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123";

const iterations = 10;
const warmup_count = 100;

const Scenario = struct {
    name: []const u8,
    points: usize,
    step: f64,

    fn time(s: Scenario, i: usize) f64 {
        return @as(f64, @floatFromInt(i)) * s.step;
    }
};

const standard = [_]Scenario{
    .{ .name = "1 day (minute)", .points = 1440, .step = 1.0 },
    .{ .name = "1 week (minute)", .points = 10080, .step = 1.0 },
    .{ .name = "2 weeks (minute)", .points = 20160, .step = 1.0 },
    .{ .name = "2 weeks (second)", .points = 1209600, .step = 1.0 / 60.0 },
    .{ .name = "1 month (minute)", .points = 43200, .step = 1.0 },
};

const large = [_]Scenario{
    .{ .name = "1 month (second)", .points = 2592000, .step = 1.0 / 60.0 },
    .{ .name = "3 months (second)", .points = 7776000, .step = 1.0 / 60.0 },
    .{ .name = "1 year (minute)", .points = 525600, .step = 1.0 },
    .{ .name = "1 year (second)", .points = 31536000, .step = 1.0 / 60.0 },
};

const Mode = enum { scalar, simd, threaded };

const Ctx = struct {
    io: std.Io,
    out: *std.Io.Writer,
    sgp4: *const Sgp4,
    threads: usize,
};

fn runScalar(sgp4: *const Sgp4, s: Scenario, start: usize, end: usize) !void {
    var i = start;
    while (i < end) : (i += 1) {
        const pv = try sgp4.propagate(s.time(i));
        std.mem.doNotOptimizeAway(pv);
    }
}

fn runSimd(sgp4: *const Sgp4, s: Scenario, start: usize, end: usize) !void {
    var i = start;
    while (i + lanes <= end) : (i += lanes) {
        var times: [lanes]f64 = undefined;
        inline for (0..lanes) |k| times[k] = s.time(i + k);
        const pv = try astroz.dispatch.sgp4Times8(sgp4, times);
        std.mem.doNotOptimizeAway(pv);
    }
    try runScalar(sgp4, s, i, end);
}

fn simdWorker(sgp4: *const Sgp4, s: Scenario, start: usize, end: usize) void {
    runSimd(sgp4, s, start, end) catch |err| std.debug.panic("propagation failed: {t}", .{err});
}

fn runThreaded(sgp4: *const Sgp4, s: Scenario, n_threads: usize, pool: []std.Thread) !void {
    const chunk = (s.points + n_threads - 1) / n_threads;
    var spawned: usize = 0;
    defer for (pool[0..spawned]) |t| t.join();
    for (0..n_threads) |t| {
        const start = @min(t * chunk, s.points);
        const end = @min(start + chunk, s.points);
        pool[t] = try std.Thread.spawn(.{}, simdWorker, .{ sgp4, s, start, end });
        spawned += 1;
    }
}

fn runOnce(ctx: Ctx, mode: Mode, s: Scenario, pool: []std.Thread) !void {
    switch (mode) {
        .scalar => try runScalar(ctx.sgp4, s, 0, s.points),
        .simd => try runSimd(ctx.sgp4, s, 0, s.points),
        .threaded => try runThreaded(ctx.sgp4, s, ctx.threads, pool),
    }
}

fn section(ctx: Ctx, title: []const u8, mode: Mode, scenarios: []const Scenario, pool: []std.Thread) !void {
    try ctx.out.print("\n--- {s} ---\n", .{title});
    try ctx.out.flush();

    var rate_sum: f64 = 0;
    for (scenarios) |s| {
        const begin = std.Io.Timestamp.now(ctx.io, .awake);
        for (0..iterations) |_| try runOnce(ctx, mode, s, pool);
        const elapsed_ns: f64 = @floatFromInt(begin.untilNow(ctx.io, .awake).toNanoseconds());

        const avg_s = elapsed_ns / iterations / 1e9;
        const rate = @as(f64, @floatFromInt(s.points)) / avg_s;
        rate_sum += rate;
        try ctx.out.print("{s:<25} {d:>10.3} ms  ({d:.2} prop/s)\n", .{ s.name, avg_s * 1e3, rate });
        try ctx.out.flush();
    }
    const mean = rate_sum / @as(f64, @floatFromInt(scenarios.len));
    try ctx.out.print("{s:<25} {d:>17.2} prop/s\n", .{ "Average", mean });
    try ctx.out.flush();
}

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const gpa = init.gpa;

    var buf: [4096]u8 = undefined;
    var stdout = std.Io.File.stdout().writer(io, &buf);
    const out = &stdout.interface;

    var tle = try astroz.Tle.parseLines(line1, line2, gpa);
    defer tle.deinit();
    const sgp4 = try Sgp4.init(tle, astroz.constants.wgs72);

    const n_threads = std.Thread.getCpuCount() catch 1;
    const pool = try gpa.alloc(std.Thread, n_threads);
    defer gpa.free(pool);

    const ctx: Ctx = .{ .io = io, .out = out, .sgp4 = &sgp4, .threads = n_threads };

    const warm: Scenario = .{ .name = "warmup", .points = warmup_count, .step = 1.0 };
    try runScalar(&sgp4, warm, 0, warm.points);
    try runSimd(&sgp4, warm, 0, warm.points);

    try out.print("\nastroz SGP4 Benchmark\n", .{});
    try out.print("{s}\n", .{"=" ** 50});
    try out.print("SIMD Batch Size: {d} (oma runtime dispatch)\n", .{lanes});

    try section(ctx, "Scalar Propagation", .scalar, &standard, pool);
    try section(ctx, std.fmt.comptimePrint("SIMD Batch{d} Propagation", .{lanes}), .simd, &standard, pool);

    try out.print("\nThreads: {d}\n", .{n_threads});
    try section(ctx, std.fmt.comptimePrint("Multithreaded SIMD Batch{d} Propagation", .{lanes}), .threaded, &large, pool);
}

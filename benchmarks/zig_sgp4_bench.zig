const std = @import("std");
const astroz = @import("astroz");

const Tle = astroz.Tle;
const Sgp4 = astroz.Sgp4;

// ISS TLE
const ISS_TLE =
    \\1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995
    \\2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123
;

// scenarios
const scenarios = [_]struct { name: []const u8, points: usize, step_min: f64 }{
    .{ .name = "1 day (minute)", .points = 1440, .step_min = 1.0 },
    .{ .name = "1 week (minute)", .points = 10080, .step_min = 1.0 },
    .{ .name = "2 weeks (minute)", .points = 20160, .step_min = 1.0 },
    .{ .name = "2 weeks (second)", .points = 1209600, .step_min = 1.0 / 60.0 },
    .{ .name = "1 month (minute)", .points = 43200, .step_min = 1.0 },
};

const iterations = 10;

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
        _ = try sgp4.propagate(@as(f64, @floatFromInt(i)));
    }

    std.debug.print("\nastroz SGP4 Benchmark\n", .{});
    std.debug.print("==================================================\n", .{});

    for (scenarios) |scenario| {
        const times = try allocator.alloc(f64, scenario.points);
        defer allocator.free(times);

        for (0..scenario.points) |i| {
            times[i] = @as(f64, @floatFromInt(i)) * scenario.step_min;
        }

        var total_ns: u64 = 0;

        for (0..iterations) |_| {
            const start = try std.time.Instant.now();

            for (times) |t| {
                _ = sgp4.propagate(t) catch continue;
            }

            const end = try std.time.Instant.now();
            total_ns += end.since(start);
        }

        const avg_ms = @as(f64, @floatFromInt(total_ns)) / @as(f64, @floatFromInt(iterations)) / 1_000_000.0;
        const props_per_sec = @as(f64, @floatFromInt(scenario.points)) / (avg_ms / 1000.0);

        std.debug.print("{s:<25} {d:>10.3} ms  ({d:.2} prop/s)\n", .{
            scenario.name,
            avg_ms,
            props_per_sec,
        });
    }

    std.debug.print("\n", .{});
}

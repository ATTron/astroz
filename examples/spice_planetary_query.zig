const std = @import("std");
const astroz = @import("astroz");
const Spice = astroz.Spice;
const constants = astroz.constants;

pub fn main() !void {
    // 1. Load default SPICE kernels (fetched via `zig build fetch-kernels`)
    Spice.loadDefaultKernels() catch |err| {
        if (err == Spice.SpiceError.NotEnabled) {
            std.debug.print(
                \\CSPICE is not enabled in this build.
                \\To run this example, rebuild with:
                \\  zig build fetch-kernels
                \\  zig build example -Denable-cspice=true
                \\
            , .{});
            return;
        }
        std.debug.print("Failed to load kernels: {}\n", .{err});
        std.debug.print("Make sure you've run: zig build fetch-kernels\n", .{});
        return;
    };

    // 2. Convert a UTC date to Ephemeris Time (ET)
    const utc = "2024-06-21T12:00:00";
    const et = Spice.utcToEt(utc) catch |err| {
        std.debug.print("Failed to convert UTC to ET: {}\n", .{err});
        return;
    };
    std.debug.print("=== SPICE Planetary Query ===\n", .{});
    std.debug.print("Epoch: {s} UTC\n", .{utc});
    std.debug.print("Ephemeris Time: {d:.3} s past J2000\n\n", .{et});

    // 3. Query positions of planets relative to Earth
    const au_km = constants.auKm;
    const planets = [_][]const u8{ "MERCURY", "VENUS", "MARS", "JUPITER BARYCENTER", "SATURN BARYCENTER" };
    const names = [_][]const u8{ "Mercury", "Venus", "Mars", "Jupiter", "Saturn" };

    std.debug.print("{s:<12} {s:>14} {s:>16}\n", .{ "Body", "Distance (AU)", "Distance (km)" });
    std.debug.print("{s:-<12} {s:->14} {s:->16}\n", .{ "", "", "" });

    for (planets, names) |planet, name| {
        const pos = Spice.getPlanetPosition(planet, et) catch |err| {
            std.debug.print("{s:<12} error: {}\n", .{ name, err });
            continue;
        };
        const dist_km = pos.magnitude();
        const dist_au = dist_km / au_km;
        std.debug.print("{s:<12} {d:>14.6} {d:>16.1}\n", .{ name, dist_au, dist_km });
    }

    // 4. Query Moon position relative to Earth
    const moon_pos = Spice.getMoonPosition(et) catch |err| {
        std.debug.print("\nMoon: error: {}\n", .{err});
        return;
    };
    const moon_dist = moon_pos.magnitude();
    std.debug.print("{s:<12} {d:>14.6} {d:>16.1}\n", .{ "Moon", moon_dist / au_km, moon_dist });

    // 5. Show Sun position
    const sun_pos = Spice.getSunPosition(et) catch |err| {
        std.debug.print("\nSun: error: {}\n", .{err});
        return;
    };
    const sun_dist = sun_pos.magnitude();
    std.debug.print("{s:<12} {d:>14.6} {d:>16.1}\n", .{ "Sun", sun_dist / au_km, sun_dist });

    std.debug.print("\nSun position (J2000, km): ({d:.1}, {d:.1}, {d:.1})\n", .{ sun_pos.x, sun_pos.y, sun_pos.z });
    std.debug.print("Moon position (J2000, km): ({d:.1}, {d:.1}, {d:.1})\n", .{ moon_pos.x, moon_pos.y, moon_pos.z });
}

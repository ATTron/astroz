//! Constellation Phasing Example
//!
//! Uses J2-induced RAAN drift for constellation plane separation
//! and sun-synchronous orbit design.

const std = @import("std");
const astroz = @import("astroz");
const constants = astroz.constants;

pub fn main() !void {
    const mu = constants.earth.mu;
    const rEarth = constants.wgs84.radiusEarthKm;
    const j2 = constants.wgs84.j2;

    // --- Sun-Synchronous Orbit Design ---
    std.debug.print("=== Sun-Synchronous Orbit Design ===\n", .{});
    std.debug.print("Target RAAN drift: +0.9856 deg/day\n\n", .{});

    const sunSyncTarget: f64 = 0.9856; // deg/day (360 / 365.25)
    const testAlts = [_]f64{ 400, 500, 600, 700, 800, 900 };

    std.debug.print("{s:>10} | {s:>10} | {s:>14}\n", .{ "Alt (km)", "Inc (deg)", "RAAN (deg/day)" });
    std.debug.print("{s:->10}-+-{s:->10}-+-{s:->14}\n", .{ "", "", "" });

    for (testAlts) |alt| {
        const a = rEarth + alt;
        const n = @sqrt(mu / (a * a * a));
        const targetRadS = sunSyncTarget * constants.deg2rad / constants.secondsPerDay;
        const denominator = -1.5 * j2 * (rEarth / a) * (rEarth / a) * n;
        const cosI = targetRadS / denominator;

        if (@abs(cosI) <= 1.0) {
            const incDeg = std.math.acos(cosI) * constants.rad2deg;
            const actualDrift = denominator * cosI * constants.rad2deg * constants.secondsPerDay;
            std.debug.print("{d:>10.0} | {d:>10.2} | {d:>14.4}\n", .{ alt, incDeg, actualDrift });
        } else {
            std.debug.print("{d:>10.0} | {s:>10} | {s:>14}\n", .{ alt, "N/A", "N/A" });
        }
    }

    // --- Constellation Deployment via Differential Drift ---
    std.debug.print("\n=== Constellation Deployment ===\n", .{});
    std.debug.print("6 planes, 60-deg separation, 53 deg inclination\n", .{});
    std.debug.print("Strategy: drift at parking altitude, then raise to operational\n\n", .{});

    const operationalAlt: f64 = 550.0;
    const parkingAlt: f64 = 520.0;
    const inc: f64 = 53.0;
    const iRad = inc * constants.deg2rad;

    const driftOp = raanDriftDegDay(mu, j2, rEarth, operationalAlt, iRad);
    const driftPark = raanDriftDegDay(mu, j2, rEarth, parkingAlt, iRad);
    const diffDrift = @abs(driftPark - driftOp);

    std.debug.print("Operational ({d:.0} km): {d:.4} deg/day\n", .{ operationalAlt, driftOp });
    std.debug.print("Parking ({d:.0} km):     {d:.4} deg/day\n", .{ parkingAlt, driftPark });
    std.debug.print("Differential:          {d:.4} deg/day\n\n", .{diffDrift});

    std.debug.print("{s:>6} | {s:>10} | {s:>12}\n", .{ "Plane", "RAAN (deg)", "Drift (days)" });
    std.debug.print("{s:->6}-+-{s:->10}-+-{s:->12}\n", .{ "", "", "" });

    for (0..6) |plane| {
        const targetRaan: f64 = @as(f64, @floatFromInt(plane)) * 60.0;
        const days = if (targetRaan > 0) targetRaan / diffDrift else 0;
        std.debug.print("{d:>6} | {d:>10.1} | {d:>12.1}\n", .{ plane + 1, targetRaan, days });
    }
}

fn raanDriftDegDay(mu: f64, j2Val: f64, rEarth: f64, alt: f64, iRad: f64) f64 {
    const a = rEarth + alt;
    const n = @sqrt(mu / (a * a * a));
    return -1.5 * j2Val * (rEarth / a) * (rEarth / a) * n * @cos(iRad) * constants.rad2deg * constants.secondsPerDay;
}

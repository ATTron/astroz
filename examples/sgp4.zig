//! SGP4/SDP4 from TLE and OMM input: single and SIMD propagation, automatic
//! near-earth/deep-space selection, and the Earth-fixed position.

const std = @import("std");
const astroz = @import("astroz");
const Tle = astroz.Tle;
const Sgp4 = astroz.Sgp4;
const Satellite = astroz.Satellite;
const constants = astroz.constants;

const issTle =
    \\1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995
    \\2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123
;
const gpsTle =
    \\1 20413U 90005A   24186.00000000  .00000012  00000+0  10000-3 0  9992
    \\2 20413  55.4408  61.4858 0112981 129.5765 231.5553  2.00561730104446
;
// CCSDS OMM as CelesTrak serves it; a TLE can't hold catalog numbers above 99999
const saramagoOmm =
    \\{"OBJECT_NAME":"SARAMAGO","OBJECT_ID":"2026-067CY","EPOCH":"2026-07-14T21:45:20.933856",
    \\ "MEAN_MOTION":15.20467281,"ECCENTRICITY":0.00055903,"INCLINATION":97.4593,"RA_OF_ASC_NODE":154.097,
    \\ "ARG_OF_PERICENTER":270.5113,"MEAN_ANOMALY":89.5482,"NORAD_CAT_ID":100000,"BSTAR":0.00022159168}
;

fn mag(v: [3]f64) f64 {
    return @sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

pub fn main(init: std.process.Init) !void {
    const gpa = init.gpa;

    var iss = try Tle.parse(issTle, gpa);
    defer iss.deinit();
    const sgp4 = try Sgp4.init(iss, constants.wgs72);

    std.debug.print("ISS, minutes after the TLE epoch:\n", .{});
    for ([_]f64{ 0, 30, 60, 90 }) |t| {
        const pv = try sgp4.propagate(t);
        std.debug.print("  t={d:5.0}  r={d:8.1} km  v={d:.3} km/s\n", .{ t, mag(pv[0]), mag(pv[1]) });
    }

    // propagateN evaluates N times in one SIMD pass
    const times = [8]f64{ 0, 15, 30, 45, 60, 75, 90, 105 };
    const batch = try sgp4.propagateN(8, times);
    std.debug.print("\nSame orbit, 8 times in one call:\n", .{});
    for (times, batch) |t, pv| std.debug.print("  t={d:5.0}  r={d:8.1} km\n", .{ t, mag(pv[0]) });

    var saramago = try Tle.parseOmm(saramagoOmm, gpa);
    defer saramago.deinit();
    const fromOmm = try Satellite.init(saramago, constants.wgs72);
    std.debug.print("\nFrom OMM: catalog number {d}, r={d:.1} km at epoch\n", .{ saramago.satelliteNumber, mag((try fromOmm.propagate(0))[0]) });

    // Satellite picks SGP4 or SDP4 from the orbital period
    var gps = try Tle.parse(gpsTle, gpa);
    defer gps.deinit();
    std.debug.print("\n", .{});
    for ([_]struct { []const u8, Tle }{ .{ "ISS", iss }, .{ "GPS", gps } }) |entry| {
        const sat = try Satellite.init(entry[1], constants.wgs72);
        const model = if (sat.isDeepSpace()) "SDP4" else "SGP4";
        std.debug.print("{s}: {s}, r={d:.1} km after 12 h\n", .{ entry[0], model, mag((try sat.propagate(720))[0]) });
    }

    const ecef = astroz.WorldCoordinateSystem.fromTle(iss, 0.0, constants.earth);
    std.debug.print("\nISS Earth-fixed position at epoch: ({d:.1}, {d:.1}, {d:.1}) km\n", .{ ecef.x, ecef.y, ecef.z });
}

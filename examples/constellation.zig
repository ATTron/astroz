//! A mixed near-earth/deep-space constellation: SIMD batch propagation to TEME
//! and Earth-fixed coordinates, then conjunction screening with exact times of
//! closest approach.

const std = @import("std");
const astroz = @import("astroz");
const Tle = astroz.Tle;
const Satellite = astroz.Satellite;
const Screening = astroz.Screening;
const constants = astroz.constants;

const tles = [_][]const u8{
    "1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998\n2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371",
    "1 55910U 23035C   24187.17717543  .00022897  00000+0  15695-2 0  9993\n2 55910  43.9975 313.1680 0011485 277.1866  82.7988 15.05748091 71356",
    "1 28626U 05004A   24186.00000000 -.00000098  00000+0  00000+0 0  9998\n2 28626   0.0163 279.8379 0003069  20.3251 343.1766  1.00270142 70992",
    "1 20413U 90005A   24186.00000000  .00000012  00000+0  10000-3 0  9992\n2 20413  55.4408  61.4858 0112981 129.5765 231.5553  2.00561730104446",
};
// Two made-up satellites at the same altitude in planes 10 degrees apart,
// a few km apart at the shared node, so they keep passing close to each other
const crossing = [_][]const u8{
    \\{"EPOCH":"2024-07-05T12:00:00","MEAN_MOTION":15.05,"ECCENTRICITY":0.0001,"INCLINATION":50,"RA_OF_ASC_NODE":0,"ARG_OF_PERICENTER":0,"MEAN_ANOMALY":0,"NORAD_CAT_ID":90001,"BSTAR":0}
    ,
    \\{"EPOCH":"2024-07-05T12:00:00","MEAN_MOTION":15.05,"ECCENTRICITY":0.0001,"INCLINATION":60,"RA_OF_ASC_NODE":0,"ARG_OF_PERICENTER":0,"MEAN_ANOMALY":0.05,"NORAD_CAT_ID":90002,"BSTAR":0}
    ,
};
const n = tles.len + crossing.len;

pub fn main(init: std.process.Init) !void {
    const gpa = init.gpa;

    var parsed: [n]Tle = undefined;
    for (tles, parsed[0..tles.len]) |s, *t| t.* = try Tle.parse(s, gpa);
    for (crossing, parsed[tles.len..]) |s, *t| t.* = try Tle.parseOmm(s, gpa);
    defer for (&parsed) |*t| t.deinit();

    var constellation = try astroz.Constellation.init(gpa, &parsed, constants.wgs72);
    defer constellation.deinit();
    std.debug.print("{d} satellites: {d} SGP4, {d} SDP4\n", .{ n, constellation.numSgp4, constellation.numSdp4 });

    // Grid: 3 hours in 1-minute steps from 10 minutes before the crossing pair's epoch
    const dt = 60.0;
    const nTimes = 180;
    const jd0 = parsed[tles.len].epochJd - 10.0 / constants.minutesPerDay;
    var jd: [nTimes]f64 = undefined;
    var fr: [nTimes]f64 = undefined;
    for (&jd, &fr, 0..) |*j, *f, k| {
        j.* = @floor(jd0);
        f.* = jd0 - @floor(jd0) + @as(f64, @floatFromInt(k)) * dt / 86400.0;
    }

    // Time-major output: pos[(k * n + s) * 3 ..]
    const pos = try gpa.alloc(f64, nTimes * n * 3);
    defer gpa.free(pos);
    const vel = try gpa.alloc(f64, nTimes * n * 3);
    defer gpa.free(vel);

    try constellation.propagate(jd[0..1], fr[0..1], pos[0 .. n * 3], null, .ecef, .timeMajor);
    std.debug.print("\nEarth-fixed positions at the grid start (km):\n", .{});
    for (0..n) |s| {
        const p = pos[s * 3 ..][0..3];
        std.debug.print("  {d:>6}  ({d:9.1}, {d:9.1}, {d:9.1})\n", .{ parsed[s].satelliteNumber, p[0], p[1], p[2] });
    }

    try constellation.propagate(&jd, &fr, pos, vel, .teme, .timeMajor);
    const found = try Screening.candidates(gpa, @ptrCast(pos), @ptrCast(vel), n, .{ .dt = dt, .threshold = 10.0 });
    defer gpa.free(found);

    // Refinement evaluates each satellite at arbitrary times with its own propagator
    var sats: [n]Satellite = undefined;
    var props: [n]Screening.Propagator = undefined;
    for (&sats, &props, parsed) |*sat, *p, t| {
        sat.* = try Satellite.init(t, constants.wgs72);
        p.* = .{
            .model = switch (sat.impl) {
                .sgp4 => |*m| .{ .sgp4 = m },
                .sdp4 => |*m| .{ .sdp4 = m },
            },
            .epochJd = sat.epochJd(),
        };
    }
    const events = try Screening.refine(gpa, &props, found, jd0, dt);
    defer gpa.free(events);
    std.mem.sort(Screening.Event, events, {}, struct {
        fn earlier(_: void, a: Screening.Event, b: Screening.Event) bool {
            return a.tSec < b.tSec;
        }
    }.earlier);

    std.debug.print("\nApproaches within 10 km over 3 hours:\n", .{});
    for (events) |ev| {
        if (!ev.converged or ev.distance > 10.0) continue;
        std.debug.print("  {d} and {d}: {d:6.3} km at t+{d:6.1} min, {d:.2} km/s\n", .{
            parsed[ev.i].satelliteNumber, parsed[ev.j].satelliteNumber, ev.distance, ev.tSec / 60.0, ev.speed,
        });
    }
}

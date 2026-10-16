//! Ground observer at a fixed geodetic site. Computes topocentric look angles
//! (azimuth, elevation, range, range rate) to a satellite.

const std = @import("std");
const constants = @import("constants.zig");
const Wcs = @import("WorldCoordinateSystem.zig");

const Observer = @This();

sinLat: f64,
cosLat: f64,
sinLon: f64,
cosLon: f64,
/// Site position in ECEF (km)
ecef: [3]f64,

pub const LookAngles = struct {
    /// Degrees clockwise from north, [0, 360)
    azimuth: f64,
    /// Degrees above the horizon
    elevation: f64,
    /// Slant range (km)
    range: f64,
    /// Range rate (km/s), positive when the satellite is moving away
    rangeRate: f64,
};

/// latitude and longitude in degrees, altitude above the WGS84 ellipsoid in km
pub fn init(latDeg: f64, lonDeg: f64, altKm: f64) Observer {
    const lat = latDeg * constants.deg2rad;
    const lon = lonDeg * constants.deg2rad;
    return .{
        .sinLat = @sin(lat),
        .cosLat = @cos(lat),
        .sinLon = @sin(lon),
        .cosLon = @cos(lon),
        .ecef = Wcs.geodeticToEcef(lat, lon, altKm),
    };
}

/// Look angles to a satellite given its ECEF position (km) and velocity (km/s)
pub fn lookAngles(self: Observer, pos: [3]f64, vel: [3]f64) LookAngles {
    const rho = [3]f64{ pos[0] - self.ecef[0], pos[1] - self.ecef[1], pos[2] - self.ecef[2] };

    // rotate the line of sight into the local south-east-zenith frame
    const south = self.sinLat * self.cosLon * rho[0] + self.sinLat * self.sinLon * rho[1] - self.cosLat * rho[2];
    const east = -self.sinLon * rho[0] + self.cosLon * rho[1];
    const zenith = self.cosLat * self.cosLon * rho[0] + self.cosLat * self.sinLon * rho[1] + self.sinLat * rho[2];

    const range = @sqrt(rho[0] * rho[0] + rho[1] * rho[1] + rho[2] * rho[2]);
    var az = std.math.atan2(east, -south);
    if (az < 0) az += 2.0 * std.math.pi;

    return .{
        .azimuth = az * constants.rad2deg,
        .elevation = std.math.asin(zenith / range) * constants.rad2deg,
        .range = range,
        .rangeRate = (rho[0] * vel[0] + rho[1] * vel[1] + rho[2] * vel[2]) / range,
    };
}

/// Look angles from a TEME state (SGP4/SDP4 output) at Julian date `jd`
pub fn lookAnglesTeme(self: Observer, pos: [3]f64, vel: [3]f64, jd: f64) LookAngles {
    const ecef = Wcs.temeToEcef(pos, vel, Wcs.julianToGmst(jd));
    return self.lookAngles(ecef[0], ecef[1]);
}

const testing = std.testing;

fn above(latDeg: f64, lonDeg: f64, altKm: f64) [3]f64 {
    return Wcs.geodeticToEcef(latDeg * constants.deg2rad, lonDeg * constants.deg2rad, altKm);
}

test "geodeticToEcef round trips through ecefToGeodeticDeg" {
    const ecef = Wcs.geodeticToEcef(40.7 * constants.deg2rad, -74.0 * constants.deg2rad, 0.05);
    const lla = Wcs.ecefToGeodeticDeg(ecef);
    try testing.expectApproxEqAbs(40.7, lla[0], 1e-9);
    try testing.expectApproxEqAbs(-74.0, lla[1], 1e-9);
    try testing.expectApproxEqAbs(0.05, lla[2], 1e-9);
}

test "satellite directly overhead" {
    const obs = Observer.init(40.7, -74.0, 0.0);
    const la = obs.lookAngles(above(40.7, -74.0, 500.0), .{ 0, 0, 0 });
    try testing.expectApproxEqAbs(90.0, la.elevation, 1e-6);
    try testing.expectApproxEqAbs(500.0, la.range, 1e-6);
}

test "azimuth points north and east" {
    const obs = Observer.init(0.0, 0.0, 0.0);
    // north: further along the meridian, east: further along the equator
    const north = obs.lookAngles(above(5.0, 0.0, 500.0), .{ 0, 0, 0 });
    const east = obs.lookAngles(above(0.0, 5.0, 500.0), .{ 0, 0, 0 });
    try testing.expectApproxEqAbs(0.0, north.azimuth, 1e-9);
    try testing.expectApproxEqAbs(90.0, east.azimuth, 1e-9);
    try testing.expect(north.elevation > 0 and north.elevation < 90);

    // antipode is below the horizon
    const below = obs.lookAngles(above(0.0, 180.0, 500.0), .{ 0, 0, 0 });
    try testing.expectApproxEqAbs(-90.0, below.elevation, 1e-9);
}

test "range rate is the line-of-sight velocity" {
    const obs = Observer.init(0.0, 0.0, 0.0);
    const pos = above(0.0, 0.0, 500.0);
    // moving straight up: receding at full speed
    try testing.expectApproxEqAbs(7.0, obs.lookAngles(pos, .{ 7.0, 0, 0 }).rangeRate, 1e-9);
    // moving horizontally: no range rate when overhead
    try testing.expectApproxEqAbs(0.0, obs.lookAngles(pos, .{ 0, 7.0, 0 }).rangeRate, 1e-9);
}

test "temeToEcef velocity is the derivative of the ECEF position" {
    const pos = [3]f64{ 6778.0, 100.0, 50.0 };
    const vel = [3]f64{ -0.1, 7.6, 0.5 };
    const jd = 2460500.5;
    const dt = 1.0; // s
    const a = Wcs.temeToEcef(pos, vel, Wcs.julianToGmst(jd));
    const pos2 = [3]f64{ pos[0] + vel[0] * dt, pos[1] + vel[1] * dt, pos[2] + vel[2] * dt };
    const b = Wcs.temeToEcef(pos2, vel, Wcs.julianToGmst(jd + dt / 86400.0));
    for (0..3) |i| try testing.expectApproxEqAbs(b[0][i] - a[0][i], a[1][i] * dt, 1e-3);
}

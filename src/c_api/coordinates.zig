//! Coordinate transformation C API exports

const astroz = @import("astroz");
const WCS = astroz.WorldCoordinateSystem;

pub fn eciToEcef(eci: *const [3]f64, gmst: f64, ecef: *[3]f64) void {
    ecef.* = WCS.eciToEcefGmst(eci.*, gmst);
}

pub fn ecefToGeodetic(ecef: *const [3]f64, lla: *[3]f64) void {
    lla.* = WCS.ecefToGeodeticDeg(ecef.*);
}

pub fn julianToGmst(jd: f64) f64 {
    return WCS.julianToGmst(jd);
}

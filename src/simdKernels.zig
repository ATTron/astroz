//! SIMD kernel functions with C calling convention for oma multi-versioning.
//! Each function wraps an existing batch propagation entry point.

const Sgp4 = @import("Sgp4.zig");
const Sdp4 = @import("Sdp4.zig");
const Sgp4Batch = @import("Sgp4Batch.zig");
const Sdp4Batch = @import("Sdp4Batch.zig");

pub fn sgp4Batch8(el: *const Sgp4Batch.BatchElements(8), t: *const [8]f64, out: *Sgp4.PosVelArray(8)) callconv(.c) c_int {
    const pv = Sgp4Batch.propagateBatchDirect(8, el, @as(@Vector(8, f64), t.*)) catch |e| return encode(e);
    out.* = Sgp4.pvToArrays(8, pv);
    return 0;
}

pub fn sdp4Batch8(el: *const Sdp4Batch.Sdp4BatchElements(8), t: *const [8]f64, carry: *Sdp4Batch.ResonanceCarryBatch(8), out: *Sgp4.PosVelArray(8)) callconv(.c) c_int {
    const pv = Sdp4Batch.propagateBatchDirect(8, el, @as(@Vector(8, f64), t.*), carry) catch |e| return encode(e);
    out.* = Sgp4.pvToArrays(8, pv);
    return 0;
}

pub fn sgp4Times8(sgp4: *const Sgp4, times: *const [8]f64, out: *Sgp4.PosVelArray(8)) callconv(.c) c_int {
    out.* = sgp4.propagateN(8, times.*) catch |e| return encode(e);
    return 0;
}

pub fn sdp4Times8(sdp4: *const Sdp4, times: *const [8]f64, out: *Sgp4.PosVelArray(8)) callconv(.c) c_int {
    out.* = sdp4.propagateN(8, times.*) catch |e| return encode(e);
    return 0;
}

fn encode(e: Sgp4.Error) c_int {
    return switch (e) {
        error.SatelliteDecayed => 1,
        error.InvalidEccentricity => 2,
        error.DeepSpaceNotSupported => 3,
        error.OutOfMemory => 4,
    };
}

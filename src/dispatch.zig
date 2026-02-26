//! Runtime SIMD dispatch via oma.
//! Resolves the best CPU-specific variant of each kernel function.
//! CPU detection is cached by oma, so resolve calls are cheap enough to call inline.

const oma = @import("oma");
const Sgp4 = @import("Sgp4.zig");
const Sdp4 = @import("Sdp4.zig");
const Sgp4Batch = @import("Sgp4Batch.zig");
const Sdp4Batch = @import("Sdp4Batch.zig");

const PosVelArray8 = Sgp4.PosVelArray(8);

const Sgp4Batch8Fn = *const fn (*const anyopaque, *const [8]f64, *anyopaque) callconv(.c) c_int;
const Sdp4Batch8Fn = *const fn (*const anyopaque, *const [8]f64, *anyopaque, *anyopaque) callconv(.c) c_int;
const Sgp4Times8Fn = *const fn (*const anyopaque, *const [8]f64, *anyopaque) callconv(.c) c_int;
const Sdp4Times8Fn = *const fn (*const anyopaque, *const [8]f64, *anyopaque) callconv(.c) c_int;

pub fn sgp4Batch8(el: *const Sgp4Batch.BatchElements(8), tsince: [8]f64) Sgp4.Error!PosVelArray8 {
    var out: PosVelArray8 = undefined;
    var t: [8]f64 = tsince;
    const rc = oma.resolveNoIo(Sgp4Batch8Fn, "sgp4Batch8")(@ptrCast(el), &t, @ptrCast(&out));
    return if (rc == 0) out else decode(rc);
}

pub fn sdp4Batch8(el: *const Sdp4Batch.Sdp4BatchElements(8), tsince: [8]f64, carry: *Sdp4Batch.ResonanceCarryBatch(8)) Sgp4.Error!PosVelArray8 {
    var out: PosVelArray8 = undefined;
    var t: [8]f64 = tsince;
    const rc = oma.resolveNoIo(Sdp4Batch8Fn, "sdp4Batch8")(@ptrCast(el), &t, @ptrCast(carry), @ptrCast(&out));
    return if (rc == 0) out else decode(rc);
}

pub fn sgp4Times8(sgp4: *const Sgp4, times: [8]f64) Sgp4.Error!PosVelArray8 {
    var out: PosVelArray8 = undefined;
    var t: [8]f64 = times;
    const rc = oma.resolveNoIo(Sgp4Times8Fn, "sgp4Times8")(@ptrCast(sgp4), &t, @ptrCast(&out));
    return if (rc == 0) out else decode(rc);
}

pub fn sdp4Times8(sdp4: *const Sdp4, times: [8]f64) Sgp4.Error!PosVelArray8 {
    var out: PosVelArray8 = undefined;
    var t: [8]f64 = times;
    const rc = oma.resolveNoIo(Sdp4Times8Fn, "sdp4Times8")(@ptrCast(sdp4), &t, @ptrCast(&out));
    return if (rc == 0) out else decode(rc);
}

fn decode(rc: c_int) Sgp4.Error {
    return switch (rc) {
        1 => error.SatelliteDecayed,
        2 => error.InvalidEccentricity,
        3 => error.DeepSpaceNotSupported,
        4 => error.OutOfMemory,
        else => unreachable,
    };
}

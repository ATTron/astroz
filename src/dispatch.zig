//! Runtime SIMD dispatch via oma.
//! Resolves the best CPU-specific variant of each kernel function.
//! CPU detection is cached by oma, so resolve calls are cheap enough to call inline.

const oma = @import("oma");
const Sgp4 = @import("Sgp4.zig");
const Sdp4 = @import("Sdp4.zig");
const Sgp4Batch = @import("Sgp4Batch.zig");
const Sdp4Batch = @import("Sdp4Batch.zig");

const PosVelArray8 = Sgp4.PosVelArray(8);

/// (propagator, times, out) -> status
const KernelFn = *const fn (*const anyopaque, *const [8]f64, *anyopaque) callconv(.c) c_int;
/// (batch, times, resonance carry, out) -> status
const CarryKernelFn = *const fn (*const anyopaque, *const [8]f64, *anyopaque, *anyopaque) callconv(.c) c_int;

pub fn sgp4Batch8(el: *const Sgp4Batch.BatchElements(8), tsince: [8]f64) Sgp4.Error!PosVelArray8 {
    var out: PosVelArray8 = undefined;
    const rc = oma.resolveNoIo(KernelFn, "sgp4Batch8")(@ptrCast(el), &tsince, @ptrCast(&out));
    return if (rc == 0) out else decode(rc);
}

pub fn sdp4Batch8(el: *const Sdp4Batch.Sdp4BatchElements(8), tsince: [8]f64, carry: *Sdp4Batch.ResonanceCarryBatch(8)) Sgp4.Error!PosVelArray8 {
    var out: PosVelArray8 = undefined;
    const rc = oma.resolveNoIo(CarryKernelFn, "sdp4Batch8")(@ptrCast(el), &tsince, @ptrCast(carry), @ptrCast(&out));
    return if (rc == 0) out else decode(rc);
}

pub fn sgp4Times8(sgp4: *const Sgp4, times: [8]f64) Sgp4.Error!PosVelArray8 {
    var out: PosVelArray8 = undefined;
    const rc = oma.resolveNoIo(KernelFn, "sgp4Times8")(@ptrCast(sgp4), &times, @ptrCast(&out));
    return if (rc == 0) out else decode(rc);
}

pub fn sdp4Times8(sdp4: *const Sdp4, times: [8]f64) Sgp4.Error!PosVelArray8 {
    var out: PosVelArray8 = undefined;
    const rc = oma.resolveNoIo(KernelFn, "sdp4Times8")(@ptrCast(sdp4), &times, @ptrCast(&out));
    return if (rc == 0) out else decode(rc);
}

fn decode(rc: c_int) Sgp4.Error {
    return switch (rc) {
        1 => error.SatelliteDecayed,
        2 => error.InvalidEccentricity,
        3 => error.DeepSpaceNotSupported,
        4 => error.OutOfMemory,
        5 => error.NegativeMeanMotion,
        6 => error.PerturbedEccentricity,
        7 => error.NegativeSemilatusRectum,
        else => unreachable,
    };
}

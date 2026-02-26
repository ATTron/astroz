//! SDP4 Cross-Satellite SIMD Batch Propagation
//! Processes N deep-space satellites simultaneously (8 AVX-512, 4 AVX2).
//! Unlike Sdp4.propagateN (1 sat × N times), this does N sats × 1 time.

const std = @import("std");
const constants = @import("constants.zig");
const simdMath = @import("simdMath.zig");
const Sgp4 = @import("Sgp4.zig");
const Sdp4 = @import("Sdp4.zig");
const Tle = @import("Tle.zig");

const Error = Sgp4.Error;

pub const BatchSize: usize = simdMath.BatchSize;

pub fn Sdp4BatchElements(comptime N: usize) type {
    const Vec = simdMath.VecN(N);

    return struct {
        const Self = @This();

        // Gravity constants (splatted)
        xke: Vec,
        j2: Vec,
        j3oj2: Vec,
        radiusEarthKm: Vec,
        vkmpersec: Vec,

        // SGP4 mean elements
        ecco: Vec,
        inclo: Vec,
        nodeo: Vec,
        argpo: Vec,
        mo: Vec,
        bstar: Vec,
        noUnkozai: Vec,

        // SGP4 precomputed trig
        sinio: Vec,
        cosio: Vec,
        con41: Vec,
        x1mth2: Vec,
        x7thm1: Vec,

        // SGP4 secular rates
        mdot: Vec,
        argpdot: Vec,
        nodedot: Vec,

        // SGP4 drag (deep-space always isimp=true, no higher-order)
        cc1: Vec,
        cc4: Vec,
        t2cof: Vec,
        xnodcf: Vec,
        xlcof: Vec,
        aycof: Vec,
        aBase: Vec,

        // Solar perturbation coefficients (12)
        solar_e2: Vec,
        solar_e3: Vec,
        solar_i2: Vec,
        solar_i3: Vec,
        solar_l2: Vec,
        solar_l3: Vec,
        solar_l4: Vec,
        solar_gh2: Vec,
        solar_gh3: Vec,
        solar_gh4: Vec,
        solar_h2: Vec,
        solar_h3: Vec,

        // Lunar perturbation coefficients (12)
        lunar_e2: Vec,
        lunar_e3: Vec,
        lunar_i2: Vec,
        lunar_i3: Vec,
        lunar_l2: Vec,
        lunar_l3: Vec,
        lunar_l4: Vec,
        lunar_gh2: Vec,
        lunar_gh3: Vec,
        lunar_gh4: Vec,
        lunar_h2: Vec,
        lunar_h3: Vec,

        // Deep-space secular rates
        zmol: Vec,
        zmos: Vec,
        dedt: Vec,
        didt: Vec,
        dmdt: Vec,
        domdt: Vec,
        dnodt: Vec,

        // Resonance type masks (f64: 1.0 = true, 0.0 = false)
        hasResonance: Vec, // irez != 0
        isHalfDay: Vec, // irez == 2

        // Half-day resonance coefficients (10)
        d2201: Vec,
        d2211: Vec,
        d3210: Vec,
        d3222: Vec,
        d4410: Vec,
        d4422: Vec,
        d5220: Vec,
        d5232: Vec,
        d5421: Vec,
        d5433: Vec,

        // GEO resonance coefficients (3)
        del1: Vec,
        del2: Vec,
        del3: Vec,

        // Resonance integration
        xlamo: Vec,
        xfact: Vec,
        gsto: Vec,

        // Epoch (per-satellite)
        epochJd: Vec,
    };
}

/// Resonance integration carry state for a batch of N satellites
pub fn ResonanceCarryBatch(comptime N: usize) type {
    const Vec = simdMath.VecN(N);
    return struct {
        atime: Vec,
        xli: Vec,
        xni: Vec,
    };
}

/// Initialize carry state from batch elements
pub fn initCarry(comptime N: usize, el: *const Sdp4BatchElements(N)) ResonanceCarryBatch(N) {
    var result: ResonanceCarryBatch(N) = undefined;
    simdMath.fillScalar(N, &result.atime, 0.0);
    simdMath.storeArray(N, &result.xli, simdMath.readArray(N, &el.xlamo));
    simdMath.storeArray(N, &result.xni, simdMath.readArray(N, &el.noUnkozai));
    return result;
}

/// Transpose N scalar Sdp4.Elements into Sdp4BatchElements(N).
pub fn initFromElements(comptime N: usize, els: [N]Sdp4.Elements, grav: constants.Sgp4GravityModel) Sdp4BatchElements(N) {
    var result: Sdp4BatchElements(N) = undefined;
    const Batch = Sdp4BatchElements(N);

    // Reflection based transpose: match batch field names to source paths
    @setEvalBranchQuota(10000);
    inline for (@typeInfo(Batch).@"struct".fields) |field| {
        const name = field.name;
        if (comptime std.mem.startsWith(u8, name, "solar_")) {
            const suffix = comptime name["solar_".len..];
            var arr: [N]f64 = undefined;
            inline for (0..N) |i| arr[i] = @field(els[i].solar, suffix);
            simdMath.storeArray(N, &@field(result, name), arr);
        } else if (comptime std.mem.startsWith(u8, name, "lunar_")) {
            const suffix = comptime name["lunar_".len..];
            var arr: [N]f64 = undefined;
            inline for (0..N) |i| arr[i] = @field(els[i].lunar, suffix);
            simdMath.storeArray(N, &@field(result, name), arr);
        } else if (comptime @hasField(Sgp4.Elements, name)) {
            var arr: [N]f64 = undefined;
            inline for (0..N) |i| arr[i] = @field(els[i].sgp4, name);
            simdMath.storeArray(N, &@field(result, name), arr);
        } else if (comptime @hasField(Sdp4.Elements, name)) {
            var arr: [N]f64 = undefined;
            inline for (0..N) |i| arr[i] = @field(els[i], name);
            simdMath.storeArray(N, &@field(result, name), arr);
        }
    }

    // Derived: resonance masks from irez
    var hasRes_arr: [N]f64 = undefined;
    var isHalf_arr: [N]f64 = undefined;
    inline for (0..N) |i| {
        hasRes_arr[i] = if (els[i].irez != 0) 1.0 else 0.0;
        isHalf_arr[i] = if (els[i].irez == 2) 1.0 else 0.0;
    }
    simdMath.storeArray(N, &result.hasResonance, hasRes_arr);
    simdMath.storeArray(N, &result.isHalfDay, isHalf_arr);

    // Splat gravity constants
    simdMath.fillScalar(N, &result.xke, grav.xke);
    simdMath.fillScalar(N, &result.j2, grav.j2);
    simdMath.fillScalar(N, &result.j3oj2, grav.j3oj2);
    simdMath.fillScalar(N, &result.radiusEarthKm, grav.radiusEarthKm);
    simdMath.fillScalar(N, &result.vkmpersec, grav.xke * grav.radiusEarthKm / 60.0);

    return result;
}

/// Core batch propagation: N satellites × 1 time.
/// tsince is a Vec with per-satellite tsince values (different epochs -> different tsince).
/// carry holds resonance integration state across sequential time calls.
pub fn propagateBatchDirect(comptime N: usize, el: *const Sdp4BatchElements(N), tsince: simdMath.VecN(N), carry: *ResonanceCarryBatch(N)) Error!Sgp4.PositionVelocity(N) {
    const Vec = simdMath.VecN(N);

    const one: Vec = @splat(1.0);
    const zero: Vec = @splat(0.0);
    const half: Vec = @splat(0.5);
    const quarter: Vec = @splat(0.25);
    const two: Vec = @splat(2.0);
    const three: Vec = @splat(3.0);
    const five: Vec = @splat(5.0);
    const seven: Vec = @splat(7.0);

    // Step 1: Simplified secular update (isimp=true always for deep space)
    const t2 = tsince * tsince;

    const tempa = one - el.cc1 * tsince;
    const tempe = el.bstar * el.cc4 * tsince;
    const templ = el.t2cof * t2;

    const xmdf = el.mo + el.mdot * tsince;
    const argpdf = el.argpo + el.argpdot * tsince;
    const nodedf = el.nodeo + el.nodedot * tsince;
    const nodem_init = nodedf + el.xnodcf * t2;

    // Step 2: Deep space secular + resonance (dspaceBatch)
    var em = el.ecco;
    var argpm = argpdf;
    var inclm = el.inclo;
    var mm = xmdf;
    var nodem = nodem_init;
    var nm = el.noUnkozai;

    // Apply lunar/solar secular perturbation rates
    em = em + el.dedt * tsince;
    inclm = inclm + el.didt * tsince;
    argpm = argpm + el.domdt * tsince;
    nodem = nodem + el.dnodt * tsince;
    mm = mm + el.dmdt * tsince;

    // Resonance effects (only for lanes with hasResonance)
    const hasRes = el.hasResonance != zero;
    if (@reduce(.Or, hasRes)) {
        // Check restart: tsince * atime <= 0 or |tsince| < |atime| or atime == 0
        const needsRestart = (carry.atime == zero) |
            ((tsince * carry.atime) <= zero) |
            (@abs(tsince) < @abs(carry.atime));
        // Restart lanes that need it (and have resonance)
        const doRestart = needsRestart & hasRes;
        carry.atime = @select(f64, doRestart, zero, carry.atime);
        carry.xni = @select(f64, doRestart, el.noUnkozai, carry.xni);
        carry.xli = @select(f64, doRestart, el.xlamo, carry.xli);

        const steppV: Vec = @splat(Sdp4.stepp);
        const step2V: Vec = @splat(Sdp4.step2);
        const positive = tsince > zero;
        const delt = @select(f64, positive, @as(Vec, @splat(Sdp4.stepp)), @as(Vec, @splat(Sdp4.stepn)));

        // Masked integration loop
        var active = hasRes & (@abs(tsince - carry.atime) >= steppV);
        while (@reduce(.Or, active)) {
            const accel = computeResonanceAccelBatch(N, el, carry.xli, carry.xni, carry.atime);
            const new_xli = carry.xli + accel.xldot * delt + accel.xndt * step2V;
            const new_xni = carry.xni + accel.xndt * delt + accel.xnddt * step2V;
            const new_atime = carry.atime + delt;
            carry.xli = @select(f64, active, new_xli, carry.xli);
            carry.xni = @select(f64, active, new_xni, carry.xni);
            carry.atime = @select(f64, active, new_atime, carry.atime);
            active = hasRes & (@abs(tsince - carry.atime) >= steppV);
        }

        // Final partial step
        const ft = tsince - carry.atime;
        const accel = computeResonanceAccelBatch(N, el, carry.xli, carry.xni, carry.atime);

        const nm_res = carry.xni + accel.xndt * ft + accel.xnddt * ft * ft * half;
        const xl = carry.xli + accel.xldot * ft + accel.xndt * ft * ft * half;
        const rptimV: Vec = @splat(Sdp4.rptim);
        const theta = simdMath.modTwoPiN(N, el.gsto + tsince * rptimV);

        // GEO formula: mm = xl - nodem - argpm + theta
        const mm_geo = xl - nodem - argpm + theta;
        // Half-day formula: mm = xl - 2*nodem + 2*theta
        const mm_half = xl - two * nodem + two * theta;

        const isHalf = el.isHalfDay != zero;
        const mm_res = @select(f64, isHalf, mm_half, mm_geo);

        // For resonant lanes, use resonance nm and mm
        nm = @select(f64, hasRes, el.noUnkozai + (nm_res - el.noUnkozai), nm);
        mm = @select(f64, hasRes, mm_res, mm);
    }

    // Step 3: Compute semi major axis and mean motion
    const nm_bad = nm <= zero;
    if (@reduce(.Or, nm_bad)) return Error.SatelliteDecayed;

    const am = simdMath.pow23N(N, el.xke / nm) * tempa * tempa;
    nm = el.xke / simdMath.pow15N(N, am);
    em = em - tempe;

    const eccFloor: Vec = @splat(1.0e-6);
    const em_bad = em >= one;
    if (@reduce(.Or, em_bad)) return Error.InvalidEccentricity;
    em = @max(em, eccFloor);
    const am_bad = am < @as(Vec, @splat(0.95));
    if (@reduce(.Or, am_bad)) return Error.SatelliteDecayed;

    mm = mm + el.noUnkozai * templ;
    const xlm = mm + argpm + nodem;
    nodem = simdMath.modTwoPiN(N, nodem);
    argpm = simdMath.modTwoPiN(N, argpm);
    mm = simdMath.modTwoPiN(N, xlm - argpm - nodem);

    // Step 4: Deep space periodic perturbations (dpperBatch)
    dpperBatch(N, el, tsince, &em, &inclm, &nodem, &argpm, &mm);

    // Negative inclination fix
    const neg_incl = inclm < zero;
    inclm = @select(f64, neg_incl, -inclm, inclm);
    const piV: Vec = @splat(std.math.pi);
    nodem = @select(f64, neg_incl, nodem + piV, nodem);
    argpm = @select(f64, neg_incl, argpm - piV, argpm);

    em = @max(em, eccFloor);
    const em_bad2 = em >= one;
    if (@reduce(.Or, em_bad2)) return Error.InvalidEccentricity;

    // Step 5: Recompute inclination dependent terms (per lane, dpper modifies inclination)
    const sinip = simdMath.sinN(N, inclm);
    const cosip = simdMath.cosN(N, inclm);
    const cosip2 = cosip * cosip;

    const aycof = -half * el.j3oj2 * sinip;
    const singTol: Vec = @splat(1.5e-12);
    const denom = cosip + one;
    const safe_denom = @select(f64, @abs(denom) > singTol, denom, singTol);
    const xlcof = -quarter * el.j3oj2 * sinip * (three + five * cosip) / safe_denom;

    const x1mth2 = one - cosip2;
    const con41 = three * cosip2 - one;
    const x7thm1 = seven * cosip2 - one;

    // Steps 6-8: Kepler solver + short-period corrections + position/velocity (shared with SGP4)
    return Sgp4.keplerAndPosVel(N, am, em, mm, argpm, nodem, inclm, aycof, xlcof, con41, x1mth2, x7thm1, sinip, cosip, el.xke, el.j2, el.radiusEarthKm, el.vkmpersec);
}

/// Compute resonance acceleration for a heterogeneous batch.
/// Both GEO and half-day paths computed, merged with masks.
fn computeResonanceAccelBatch(comptime N: usize, el: *const Sdp4BatchElements(N), xli: simdMath.VecN(N), xni: simdMath.VecN(N), atime: simdMath.VecN(N)) struct { xndt: simdMath.VecN(N), xnddt: simdMath.VecN(N), xldot: simdMath.VecN(N) } {
    const Vec = simdMath.VecN(N);
    const zero: Vec = @splat(0.0);
    const two: Vec = @splat(2.0);
    const three: Vec = @splat(3.0);

    const xldot = xni + el.xfact;

    // Half-day resonance path
    const xomi = el.argpo + el.argpdot * atime;
    const x2omi = xomi + xomi;
    const x2li = xli + xli;

    const g22V: Vec = @splat(Sdp4.g22);
    const g32V: Vec = @splat(Sdp4.g32);
    const g44V: Vec = @splat(Sdp4.g44);
    const g52V: Vec = @splat(Sdp4.g52);
    const g54V: Vec = @splat(Sdp4.g54);

    const h1 = x2omi + xli - g22V;
    const h2 = xli - g22V;
    const h3 = xomi + xli - g32V;
    const h4 = -xomi + xli - g32V;
    const h5 = x2omi + x2li - g44V;
    const h6 = x2li - g44V;
    const h7 = xomi + xli - g52V;
    const h8 = -xomi + xli - g52V;
    const h9 = xomi + x2li - g54V;
    const h10 = -xomi + x2li - g54V;

    const sch1 = simdMath.sincosN(N, h1);
    const sch2 = simdMath.sincosN(N, h2);
    const sch3 = simdMath.sincosN(N, h3);
    const sch4 = simdMath.sincosN(N, h4);
    const sch5 = simdMath.sincosN(N, h5);
    const sch6 = simdMath.sincosN(N, h6);
    const sch7 = simdMath.sincosN(N, h7);
    const sch8 = simdMath.sincosN(N, h8);
    const sch9 = simdMath.sincosN(N, h9);
    const sch10 = simdMath.sincosN(N, h10);

    const xndtHalf = el.d2201 * sch1.sin + el.d2211 * sch2.sin +
        el.d3210 * sch3.sin + el.d3222 * sch4.sin +
        el.d4410 * sch5.sin + el.d4422 * sch6.sin +
        el.d5220 * sch7.sin + el.d5232 * sch8.sin +
        el.d5421 * sch9.sin + el.d5433 * sch10.sin;

    const xnddtHalf = (el.d2201 * sch1.cos + el.d2211 * sch2.cos +
        el.d3210 * sch3.cos + el.d3222 * sch4.cos +
        el.d5220 * sch7.cos + el.d5232 * sch8.cos +
        two * (el.d4410 * sch5.cos + el.d4422 * sch6.cos +
            el.d5421 * sch9.cos + el.d5433 * sch10.cos)) * xldot;

    // GEO resonance path
    const fasx2V: Vec = @splat(Sdp4.fasx2);
    const fasx4V: Vec = @splat(Sdp4.fasx4);
    const fasx6V: Vec = @splat(Sdp4.fasx6);

    const p1 = xli - fasx2V;
    const p2 = two * (xli - fasx4V);
    const p3 = three * (xli - fasx6V);

    const scp1 = simdMath.sincosN(N, p1);
    const scp2 = simdMath.sincosN(N, p2);
    const scp3 = simdMath.sincosN(N, p3);

    const xndtGeo = el.del1 * scp1.sin + el.del2 * scp2.sin + el.del3 * scp3.sin;
    const xnddtGeo = (el.del1 * scp1.cos + two * el.del2 * scp2.cos + three * el.del3 * scp3.cos) * xldot;

    // Merge: half-day lanes get half-day, GEO lanes get GEO, non-resonant get 0
    const isHalf = el.isHalfDay != zero;
    const hasRes = el.hasResonance != zero;
    var xndt = @select(f64, isHalf, xndtHalf, xndtGeo);
    xndt = @select(f64, hasRes, xndt, zero);
    var xnddt = @select(f64, isHalf, xnddtHalf, xnddtGeo);
    xnddt = @select(f64, hasRes, xnddt, zero);

    return .{ .xndt = xndt, .xnddt = xnddt, .xldot = xldot };
}

/// Deep space periodic perturbations for a batch of N satellites
fn dpperBatch(
    comptime N: usize,
    el: *const Sdp4BatchElements(N),
    tsince: simdMath.VecN(N),
    ep: *simdMath.VecN(N),
    inclp: *simdMath.VecN(N),
    nodep: *simdMath.VecN(N),
    argpp: *simdMath.VecN(N),
    mp: *simdMath.VecN(N),
) void {
    const Vec = simdMath.VecN(N);
    const two: Vec = @splat(2.0);
    const half: Vec = @splat(0.5);
    const quarter: Vec = @splat(0.25);

    // Solar terms
    const znsV: Vec = @splat(Sdp4.zns);
    const zesV: Vec = @splat(Sdp4.zes);

    var zm = el.zmos + znsV * tsince;
    var zf = zm + two * zesV * simdMath.sinN(N, zm);
    var sinzf = simdMath.sinN(N, zf);
    var f2 = half * sinzf * sinzf - quarter;
    var f3 = -half * sinzf * simdMath.cosN(N, zf);

    const ses = el.solar_e2 * f2 + el.solar_e3 * f3;
    const sis = el.solar_i2 * f2 + el.solar_i3 * f3;
    const sls = el.solar_l2 * f2 + el.solar_l3 * f3 + el.solar_l4 * sinzf;
    const sghs = el.solar_gh2 * f2 + el.solar_gh3 * f3 + el.solar_gh4 * sinzf;
    const shs = el.solar_h2 * f2 + el.solar_h3 * f3;

    // Lunar terms
    const znlV: Vec = @splat(Sdp4.znl);
    const zelV: Vec = @splat(Sdp4.zel);

    zm = el.zmol + znlV * tsince;
    zf = zm + two * zelV * simdMath.sinN(N, zm);
    sinzf = simdMath.sinN(N, zf);
    f2 = half * sinzf * sinzf - quarter;
    f3 = -half * sinzf * simdMath.cosN(N, zf);

    const sel = el.lunar_e2 * f2 + el.lunar_e3 * f3;
    const sil = el.lunar_i2 * f2 + el.lunar_i3 * f3;
    const sll = el.lunar_l2 * f2 + el.lunar_l3 * f3 + el.lunar_l4 * sinzf;
    const sghl = el.lunar_gh2 * f2 + el.lunar_gh3 * f3 + el.lunar_gh4 * sinzf;
    const shl = el.lunar_h2 * f2 + el.lunar_h3 * f3;

    // Combined
    const pe = ses + sel;
    const pinc = sis + sil;
    const pl = sls + sll;
    const pgh = sghs + sghl;
    const ph = shs + shl;

    inclp.* = inclp.* + pinc;
    ep.* = ep.* + pe;

    const sinip = simdMath.sinN(N, inclp.*);
    const cosip = simdMath.cosN(N, inclp.*);

    // Normal path (inclp >= 0.2)
    const ph_norm = ph / sinip;
    const pgh_norm = pgh - cosip * ph_norm;
    const argpp_norm = argpp.* + pgh_norm;
    const nodep_norm = nodep.* + ph_norm;
    const mp_norm = mp.* + pl;

    // Lyddane path (inclp < 0.2)
    const sinop = simdMath.sinN(N, nodep.*);
    const cosop = simdMath.cosN(N, nodep.*);
    var alfdp = sinip * sinop;
    var betdp = sinip * cosop;
    const dalf = ph * cosop + pinc * cosip * sinop;
    const dbet = -ph * sinop + pinc * cosip * cosop;
    alfdp = alfdp + dalf;
    betdp = betdp + dbet;
    const nodep_mod = simdMath.modTwoPiN(N, nodep.*);
    const xls = mp.* + argpp.* + cosip * nodep_mod;
    const dls = pl + pgh - pinc * nodep_mod * sinip;
    const xnoh = nodep_mod;
    var nodep_lyd = simdMath.atan2N(N, alfdp, betdp);

    // Continuity fix
    const diff_abs = @abs(xnoh - nodep_lyd);
    const piV: Vec = @splat(std.math.pi);
    const twoPiV: Vec = @splat(constants.twoPi);
    const need_adjust = diff_abs > piV;
    const adjust_sign = @select(f64, nodep_lyd < xnoh, twoPiV, -twoPiV);
    nodep_lyd = @select(f64, need_adjust, nodep_lyd + adjust_sign, nodep_lyd);
    const mp_lyd = mp.* + pl;
    const argpp_lyd = xls + dls - mp_lyd - cosip * nodep_lyd;

    // Merge with @select
    const thresh: Vec = @splat(0.2);
    const useNormal = inclp.* >= thresh;
    argpp.* = @select(f64, useNormal, argpp_norm, argpp_lyd);
    nodep.* = @select(f64, useNormal, nodep_norm, nodep_lyd);
    mp.* = @select(f64, useNormal, mp_norm, mp_lyd);
}

const testing = std.testing;

fn expectBatchMatchesScalar(tle_str: []const u8) !void {
    var tle = try Tle.parse(tle_str, testing.allocator);
    defer tle.deinit();
    const sdp4 = try Sdp4.init(tle, constants.wgs72);

    const els = [4]Sdp4.Elements{ sdp4.elements, sdp4.elements, sdp4.elements, sdp4.elements };
    const batch = initFromElements(4, els, constants.wgs72);

    for ([_]f64{ 0.0, 360.0, 720.0, 1440.0 }) |t| {
        var carry = initCarry(4, &batch);
        const tsince: simdMath.VecN(4) = @splat(t);
        const pv = try propagateBatchDirect(4, &batch, tsince, &carry);
        const scalar = try sdp4.propagate(t);
        try testing.expectApproxEqAbs(scalar[0][0], pv.rx[0], 1e-3);
        try testing.expectApproxEqAbs(scalar[0][1], pv.ry[0], 1e-3);
        try testing.expectApproxEqAbs(scalar[0][2], pv.rz[0], 1e-3);
        try testing.expectApproxEqAbs(scalar[1][0], pv.vx[0], 1e-6);
        try testing.expectApproxEqAbs(scalar[1][1], pv.vy[0], 1e-6);
        try testing.expectApproxEqAbs(scalar[1][2], pv.vz[0], 1e-6);
    }
}

test "Sdp4Batch matches scalar - GPS (irez=0)" {
    try expectBatchMatchesScalar("1 20413U 90005A   24186.00000000  .00000012  00000+0  10000-3 0  9992\n2 20413  55.4408  61.4858 0112981 129.5765 231.5553  2.00561730104446");
}

test "Sdp4Batch matches scalar - GEO (irez=1)" {
    try expectBatchMatchesScalar("1 28626U 05004A   24186.00000000 -.00000098  00000+0  00000+0 0  9998\n2 28626   0.0163 279.8379 0003069  20.3251 343.1766  1.00270142 70992");
}

test "Sdp4Batch matches scalar - HEO (irez=2)" {
    try expectBatchMatchesScalar("1 09880U 77021B   24186.00000000  .00000023  00000+0  00000+0 0  9999\n2 09880  63.4300  75.8891 7318036 269.8735  16.7549  2.00611684 54321");
}

test "Sdp4Batch mixed irez batch" {
    // Mix GPS (irez=0), GEO (irez=1), HEO (irez=2), and GPS again
    const tleStrs = [_][]const u8{
        "1 20413U 90005A   24186.00000000  .00000012  00000+0  10000-3 0  9992\n2 20413  55.4408  61.4858 0112981 129.5765 231.5553  2.00561730104446",
        "1 28626U 05004A   24186.00000000 -.00000098  00000+0  00000+0 0  9998\n2 28626   0.0163 279.8379 0003069  20.3251 343.1766  1.00270142 70992",
        "1 09880U 77021B   24186.00000000  .00000023  00000+0  00000+0 0  9999\n2 09880  63.4300  75.8891 7318036 269.8735  16.7549  2.00611684 54321",
        "1 20413U 90005A   24186.00000000  .00000012  00000+0  10000-3 0  9992\n2 20413  55.4408  61.4858 0112981 129.5765 231.5553  2.00561730104446",
    };

    var tles: [4]Tle = undefined;
    defer for (&tles) |*t| t.deinit();
    for (0..4) |i| tles[i] = try Tle.parse(tleStrs[i], testing.allocator);

    var sdp4s: [4]Sdp4 = undefined;
    var batch_els: [4]Sdp4.Elements = undefined;
    for (0..4) |i| {
        sdp4s[i] = try Sdp4.init(tles[i], constants.wgs72);
        batch_els[i] = sdp4s[i].elements;
    }

    const batch = initFromElements(4, batch_els, constants.wgs72);

    for ([_]f64{ 0.0, 720.0, 1440.0 }) |t| {
        var carry = initCarry(4, &batch);
        const tsince: simdMath.VecN(4) = @splat(t);
        const pv = try propagateBatchDirect(4, &batch, tsince, &carry);

        inline for (0..4) |i| {
            const scalar = try sdp4s[i].propagate(t);
            try testing.expectApproxEqAbs(scalar[0][0], pv.rx[i], 1e-3);
            try testing.expectApproxEqAbs(scalar[0][1], pv.ry[i], 1e-3);
            try testing.expectApproxEqAbs(scalar[0][2], pv.rz[i], 1e-3);
            try testing.expectApproxEqAbs(scalar[1][0], pv.vx[i], 1e-6);
            try testing.expectApproxEqAbs(scalar[1][1], pv.vy[i], 1e-6);
            try testing.expectApproxEqAbs(scalar[1][2], pv.vz[i], 1e-6);
        }
    }
}

test "Sdp4Batch carry convergence" {
    // Test that sequential carry matches non-carry
    const tle_str = "1 28626U 05004A   24186.00000000 -.00000098  00000+0  00000+0 0  9998\n2 28626   0.0163 279.8379 0003069  20.3251 343.1766  1.00270142 70992";
    var tle = try Tle.parse(tle_str, testing.allocator);
    defer tle.deinit();
    const sdp4 = try Sdp4.init(tle, constants.wgs72);

    const els = [4]Sdp4.Elements{ sdp4.elements, sdp4.elements, sdp4.elements, sdp4.elements };
    const batch = initFromElements(4, els, constants.wgs72);

    // Propagate sequentially with carry
    var carry = initCarry(4, &batch);
    const times = [_]f64{ 0.0, 360.0, 720.0, 1080.0, 1440.0 };

    for (times) |t| {
        const tsince: simdMath.VecN(4) = @splat(t);
        const pv_carry = try propagateBatchDirect(4, &batch, tsince, &carry);

        // Fresh propagation (no carry)
        var fresh_carry = initCarry(4, &batch);
        const pv_fresh = try propagateBatchDirect(4, &batch, tsince, &fresh_carry);

        try testing.expectApproxEqAbs(pv_fresh.rx[0], pv_carry.rx[0], 1e-6);
        try testing.expectApproxEqAbs(pv_fresh.ry[0], pv_carry.ry[0], 1e-6);
        try testing.expectApproxEqAbs(pv_fresh.rz[0], pv_carry.rz[0], 1e-6);
    }
}

test "Sdp4Batch Vec8 matches Vec4" {
    const tle_str = "1 20413U 90005A   24186.00000000  .00000012  00000+0  10000-3 0  9992\n2 20413  55.4408  61.4858 0112981 129.5765 231.5553  2.00561730104446";
    var tle = try Tle.parse(tle_str, testing.allocator);
    defer tle.deinit();
    const sdp4 = try Sdp4.init(tle, constants.wgs72);

    const els4 = [4]Sdp4.Elements{ sdp4.elements, sdp4.elements, sdp4.elements, sdp4.elements };
    const els8 = [8]Sdp4.Elements{ sdp4.elements, sdp4.elements, sdp4.elements, sdp4.elements, sdp4.elements, sdp4.elements, sdp4.elements, sdp4.elements };

    const batch4 = initFromElements(4, els4, constants.wgs72);
    const batch8 = initFromElements(8, els8, constants.wgs72);

    for ([_]f64{ 0.0, 720.0, 1440.0 }) |t| {
        var carry4 = initCarry(4, &batch4);
        var carry8 = initCarry(8, &batch8);
        const pv4 = try propagateBatchDirect(4, &batch4, @as(simdMath.VecN(4), @splat(t)), &carry4);
        const pv8 = try propagateBatchDirect(8, &batch8, @as(simdMath.VecN(8), @splat(t)), &carry8);

        // First 4 lanes of Vec8 should match Vec4
        inline for (0..4) |i| {
            try testing.expectApproxEqAbs(pv4.rx[i], pv8.rx[i], 1e-10);
            try testing.expectApproxEqAbs(pv4.ry[i], pv8.ry[i], 1e-10);
            try testing.expectApproxEqAbs(pv4.rz[i], pv8.rz[i], 1e-10);
        }
    }
}

//! Multi-satellite SIMD batch propagation (4 satellites per batch)
//! Provides struct-of-arrays ElementsV4 and batch propagation functions

const std = @import("std");
const constants = @import("constants.zig");
const simdMath = @import("simdMath.zig");
const Sgp4 = @import("Sgp4.zig");
const Tle = @import("Tle.zig");

const Vec4 = simdMath.Vec4;
const Elements = Sgp4.Elements;
const SecularStateV4 = Sgp4.SecularStateV4;
const KeplerStateV4 = Sgp4.KeplerStateV4;
const CorrectedStateV4 = Sgp4.CorrectedStateV4;
const eccentricityFloor = Sgp4.eccentricityFloor;
const Error = Sgp4.Error;

pub const ElementsV4 = struct {
    // Gravity model (shared across all 4 satellites)
    grav: constants.Sgp4GravityModel,

    // Mean elements from TLE
    noKozai: Vec4,
    ecco: Vec4,
    inclo: Vec4,
    nodeo: Vec4,
    argpo: Vec4,
    mo: Vec4,
    bstar: Vec4,

    // Derived mean elements
    noUnkozai: Vec4,
    a: Vec4,

    // Trigonometric terms
    sinio: Vec4,
    cosio: Vec4,
    cosio2: Vec4,
    cosio4: Vec4,

    // Polynomial terms in cos^2(i)
    con41: Vec4, // 3cos^2(i) - 1
    con42: Vec4, // 1 - 5cos^2(i)
    x1mth2: Vec4, // sin^2(i)
    x7thm1: Vec4, // 7cos^2(i) - 1

    // Secular rate coefficients
    mdot: Vec4,
    argpdot: Vec4,
    nodedot: Vec4,

    // Drag coefficients
    cc1: Vec4,
    cc4: Vec4,
    cc5: Vec4,
    t2cof: Vec4,
    omgcof: Vec4,
    xnodcf: Vec4,
    xlcof: Vec4,
    xmcof: Vec4,
    aycof: Vec4,
    eta: Vec4,
    delmo: Vec4,
    sinmao: Vec4,

    // Higher-order drag terms
    d2: Vec4,
    d3: Vec4,
    d4: Vec4,
    t3cof: Vec4,
    t4cof: Vec4,
    t5cof: Vec4,

    // Precomputed for propagation
    aBase: Vec4, // cbrt((xke/noUnkozai)^2)
    vkmpersec: Vec4, // velocity conversion factor

    // Per-satellite simplified drag model flag (1.0 = simplified, 0.0 = full)
    isimpMask: Vec4,

    // Julian dates of TLE epochs for each satellite
    epochJd: Vec4,

    // Pre-splatted constants (computed once at init, eliminates repeated splats)
    xke: Vec4,
    j2: Vec4,
    radiusEarthKm: Vec4,
    one: Vec4,
    zero: Vec4,
    half: Vec4,
    quarter: Vec4,
    oneHalf: Vec4,
    two: Vec4,
    eccFloor: Vec4,
};

/// Initialize 4 satellites into a struct-of-arrays format for multi-satellite SIMD propagation
pub fn initElementsV4(tles: [4]Tle, grav: constants.Sgp4GravityModel) Error!ElementsV4 {
    var els: [4]Elements = undefined;
    for (0..4) |i| {
        els[i] = try Sgp4.initElements(tles[i], grav);
    }

    return ElementsV4{
        .grav = grav,
        .noKozai = Vec4{ els[0].noKozai, els[1].noKozai, els[2].noKozai, els[3].noKozai },
        .ecco = Vec4{ els[0].ecco, els[1].ecco, els[2].ecco, els[3].ecco },
        .inclo = Vec4{ els[0].inclo, els[1].inclo, els[2].inclo, els[3].inclo },
        .nodeo = Vec4{ els[0].nodeo, els[1].nodeo, els[2].nodeo, els[3].nodeo },
        .argpo = Vec4{ els[0].argpo, els[1].argpo, els[2].argpo, els[3].argpo },
        .mo = Vec4{ els[0].mo, els[1].mo, els[2].mo, els[3].mo },
        .bstar = Vec4{ els[0].bstar, els[1].bstar, els[2].bstar, els[3].bstar },
        .noUnkozai = Vec4{ els[0].noUnkozai, els[1].noUnkozai, els[2].noUnkozai, els[3].noUnkozai },
        .a = Vec4{ els[0].a, els[1].a, els[2].a, els[3].a },
        .sinio = Vec4{ els[0].sinio, els[1].sinio, els[2].sinio, els[3].sinio },
        .cosio = Vec4{ els[0].cosio, els[1].cosio, els[2].cosio, els[3].cosio },
        .cosio2 = Vec4{ els[0].cosio2, els[1].cosio2, els[2].cosio2, els[3].cosio2 },
        .cosio4 = Vec4{ els[0].cosio4, els[1].cosio4, els[2].cosio4, els[3].cosio4 },
        .con41 = Vec4{ els[0].con41, els[1].con41, els[2].con41, els[3].con41 },
        .con42 = Vec4{ els[0].con42, els[1].con42, els[2].con42, els[3].con42 },
        .x1mth2 = Vec4{ els[0].x1mth2, els[1].x1mth2, els[2].x1mth2, els[3].x1mth2 },
        .x7thm1 = Vec4{ els[0].x7thm1, els[1].x7thm1, els[2].x7thm1, els[3].x7thm1 },
        .mdot = Vec4{ els[0].mdot, els[1].mdot, els[2].mdot, els[3].mdot },
        .argpdot = Vec4{ els[0].argpdot, els[1].argpdot, els[2].argpdot, els[3].argpdot },
        .nodedot = Vec4{ els[0].nodedot, els[1].nodedot, els[2].nodedot, els[3].nodedot },
        .cc1 = Vec4{ els[0].cc1, els[1].cc1, els[2].cc1, els[3].cc1 },
        .cc4 = Vec4{ els[0].cc4, els[1].cc4, els[2].cc4, els[3].cc4 },
        .cc5 = Vec4{ els[0].cc5, els[1].cc5, els[2].cc5, els[3].cc5 },
        .t2cof = Vec4{ els[0].t2cof, els[1].t2cof, els[2].t2cof, els[3].t2cof },
        .omgcof = Vec4{ els[0].omgcof, els[1].omgcof, els[2].omgcof, els[3].omgcof },
        .xnodcf = Vec4{ els[0].xnodcf, els[1].xnodcf, els[2].xnodcf, els[3].xnodcf },
        .xlcof = Vec4{ els[0].xlcof, els[1].xlcof, els[2].xlcof, els[3].xlcof },
        .xmcof = Vec4{ els[0].xmcof, els[1].xmcof, els[2].xmcof, els[3].xmcof },
        .aycof = Vec4{ els[0].aycof, els[1].aycof, els[2].aycof, els[3].aycof },
        .eta = Vec4{ els[0].eta, els[1].eta, els[2].eta, els[3].eta },
        .delmo = Vec4{ els[0].delmo, els[1].delmo, els[2].delmo, els[3].delmo },
        .sinmao = Vec4{ els[0].sinmao, els[1].sinmao, els[2].sinmao, els[3].sinmao },
        .d2 = Vec4{ els[0].d2, els[1].d2, els[2].d2, els[3].d2 },
        .d3 = Vec4{ els[0].d3, els[1].d3, els[2].d3, els[3].d3 },
        .d4 = Vec4{ els[0].d4, els[1].d4, els[2].d4, els[3].d4 },
        .t3cof = Vec4{ els[0].t3cof, els[1].t3cof, els[2].t3cof, els[3].t3cof },
        .t4cof = Vec4{ els[0].t4cof, els[1].t4cof, els[2].t4cof, els[3].t4cof },
        .t5cof = Vec4{ els[0].t5cof, els[1].t5cof, els[2].t5cof, els[3].t5cof },
        .aBase = Vec4{ els[0].aBase, els[1].aBase, els[2].aBase, els[3].aBase },
        .vkmpersec = Vec4{ els[0].vkmpersec, els[1].vkmpersec, els[2].vkmpersec, els[3].vkmpersec },
        .isimpMask = Vec4{
            if (els[0].isimp) 1.0 else 0.0,
            if (els[1].isimp) 1.0 else 0.0,
            if (els[2].isimp) 1.0 else 0.0,
            if (els[3].isimp) 1.0 else 0.0,
        },
        .epochJd = Vec4{ els[0].epochJd, els[1].epochJd, els[2].epochJd, els[3].epochJd },
        .xke = @splat(grav.xke),
        .j2 = @splat(grav.j2),
        .radiusEarthKm = @splat(grav.radiusEarthKm),
        .one = @splat(1.0),
        .zero = @splat(0.0),
        .half = @splat(0.5),
        .quarter = @splat(0.25),
        .oneHalf = @splat(1.5),
        .two = @splat(2.0),
        .eccFloor = @splat(eccentricityFloor),
    };
}

/// Propagate 4 satellites at a single time point using SIMD
pub inline fn propagateSatellitesV4(el: *const ElementsV4, tsince: f64) Error![4][2][3]f64 {
    const tsinceVec: Vec4 = @splat(tsince);
    return propagateSatellitesV4Vec(el, tsinceVec);
}

/// Propagate 4 satellites with per-satellite time offsets (Vec4)
pub inline fn propagateSatellitesV4Vec(el: *const ElementsV4, tsinceVec: Vec4) Error![4][2][3]f64 {
    const secular = updateSecularSatV4(el, tsinceVec);

    const decayed = secular.em < el.eccFloor;
    if (@reduce(.Or, decayed)) {
        return Error.SatelliteDecayed;
    }

    const nm: Vec4 = el.xke / simdMath.pow15V4(secular.a);
    const kepler = solveKeplerSatV4(el, secular);
    const corrected = applyShortPeriodCorrectionsSatV4(el, kepler, nm);
    return computePositionVelocitySatV4(el, corrected);
}

inline fn updateSecularSatV4(el: *const ElementsV4, tsince: Vec4) SecularStateV4 {
    const t2 = tsince * tsince;

    var tempa = el.one - el.cc1 * tsince;
    var tempe = el.bstar * el.cc4 * tsince;
    var templ = el.t2cof * t2;

    const xmdf = el.mo + el.mdot * tsince;
    const argpdf = el.argpo + el.argpdot * tsince;
    const nodedf = el.nodeo + el.nodedot * tsince;
    var argpm = argpdf;
    var mm = xmdf;
    var nodem = nodedf + el.xnodcf * t2;

    const delomg = el.omgcof * tsince;
    const delmtemp = el.one + el.eta * @cos(xmdf);
    const delmHo = el.xmcof * (delmtemp * delmtemp * delmtemp - el.delmo);
    const tempHo = delomg + delmHo;

    const temp = @select(f64, el.isimpMask == el.zero, tempHo, el.zero);

    mm = xmdf + temp;
    argpm = argpdf - temp;

    const t3 = t2 * tsince;
    const t4 = t3 * tsince;

    const tempaHo = tempa - el.d2 * t2 - el.d3 * t3 - el.d4 * t4;
    const tempeHo = tempe + el.bstar * el.cc5 * (@sin(mm) - el.sinmao);
    const templHo = templ + el.t3cof * t3 + t4 * (el.t4cof + tsince * el.t5cof);

    tempa = @select(f64, el.isimpMask == el.zero, tempaHo, tempa);
    tempe = @select(f64, el.isimpMask == el.zero, tempeHo, tempe);
    templ = @select(f64, el.isimpMask == el.zero, templHo, templ);

    const am = el.aBase * tempa * tempa;
    var em = el.ecco - tempe;
    em = @max(em, el.eccFloor);

    mm = mm + el.noUnkozai * templ;
    const xlm = mm + argpm + nodem;

    nodem = @mod(nodem, simdMath.twoPiVec);
    argpm = @mod(argpm, simdMath.twoPiVec);
    mm = @mod(xlm - argpm - nodem, simdMath.twoPiVec);

    return .{
        .mm = mm,
        .argpm = argpm,
        .nodem = nodem,
        .em = em,
        .a = am,
    };
}

inline fn solveKeplerSatV4(el: *const ElementsV4, sec: SecularStateV4) KeplerStateV4 {
    const temp = el.one / (sec.a * (el.one - sec.em * sec.em));
    const axnl = sec.em * @cos(sec.argpm);
    const aynl = sec.em * @sin(sec.argpm) + temp * el.aycof;
    const xl = @mod(sec.mm + sec.argpm + sec.nodem + temp * el.xlcof * axnl, simdMath.twoPiVec);

    var u = @mod(xl - sec.nodem, simdMath.twoPiVec);
    var eo1 = u;
    var sineo1: Vec4 = el.zero;
    var coseo1: Vec4 = el.one;

    const tolerance: Vec4 = @splat(1.0e-12);
    const clampVal: Vec4 = @splat(0.95);

    var ktr: u32 = 0;
    while (ktr < 10) : (ktr += 1) {
        sineo1 = @sin(eo1);
        coseo1 = @cos(eo1);
        var tem5 = el.one - coseo1 * axnl - sineo1 * aynl;
        tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
        const positive = tem5 > clampVal;
        const negative = tem5 < -clampVal;
        tem5 = @select(f64, positive, clampVal, tem5);
        tem5 = @select(f64, negative, -clampVal, tem5);
        eo1 = eo1 + tem5;

        const converged = @abs(tem5) < tolerance;
        if (@reduce(.And, converged)) break;
    }

    const ecose = axnl * coseo1 + aynl * sineo1;
    const esine = axnl * sineo1 - aynl * coseo1;
    const el2 = axnl * axnl + aynl * aynl;
    const pl = sec.a * (el.one - el2);
    const betal = @sqrt(el.one - el2);
    const rl = sec.a * (el.one - ecose);
    const rdotl = @sqrt(sec.a) * esine / rl;
    const rvdotl = @sqrt(pl) / rl;

    const aOverR = sec.a / rl;
    const esineTerm = esine / (el.one + betal);
    const sinu = aOverR * (sineo1 - aynl - axnl * esineTerm);
    const cosu = aOverR * (coseo1 - axnl + aynl * esineTerm);
    u = simdMath.atan2SIMD(sinu, cosu);

    return .{
        .u = u,
        .r = rl,
        .rdot = rdotl,
        .rvdot = rvdotl,
        .betal = betal,
        .sin2u = el.two * sinu * cosu,
        .cos2u = el.one - el.two * sinu * sinu,
        .nodem = sec.nodem,
        .pl = pl,
    };
}

inline fn applyShortPeriodCorrectionsSatV4(el: *const ElementsV4, kep: KeplerStateV4, nm: Vec4) CorrectedStateV4 {
    const temp = el.one / kep.pl;
    const temp1 = el.half * el.j2 * temp;
    const temp2 = temp1 * temp;

    const mrt = kep.r * (el.one - el.oneHalf * temp2 * kep.betal * el.con41) + el.half * temp1 * el.x1mth2 * kep.cos2u;
    const su = kep.u - el.quarter * temp2 * el.x7thm1 * kep.sin2u;
    const xnode = kep.nodem + el.oneHalf * temp2 * el.cosio * kep.sin2u;
    const xinc = el.inclo + el.oneHalf * temp2 * el.cosio * el.sinio * kep.cos2u;
    const mvt = kep.rdot - nm * temp1 * el.x1mth2 * kep.sin2u / el.xke;
    const rvdot = kep.rvdot + nm * temp1 * (el.x1mth2 * kep.cos2u + el.half * el.con41) / el.xke;

    return .{ .r = mrt, .rdot = mvt, .rvdot = rvdot, .u = su, .xnode = xnode, .xinc = xinc };
}

inline fn computePositionVelocitySatV4(el: *const ElementsV4, state: CorrectedStateV4) [4][2][3]f64 {
    const sinsu = @sin(state.u);
    const cossu = @cos(state.u);
    const snod = @sin(state.xnode);
    const cnod = @cos(state.xnode);
    const sini = @sin(state.xinc);
    const cosi = @cos(state.xinc);

    const xmx = -snod * cosi;
    const xmy = cnod * cosi;
    const ux = xmx * sinsu + cnod * cossu;
    const uy = xmy * sinsu + snod * cossu;
    const uz = sini * sinsu;
    const vx = xmx * cossu - cnod * sinsu;
    const vy = xmy * cossu - snod * sinsu;
    const vz = sini * cossu;

    const rScaled = state.r * el.radiusEarthKm;
    const rx = rScaled * ux;
    const ry = rScaled * uy;
    const rz = rScaled * uz;

    const vxOut = (state.rdot * ux + state.rvdot * vx) * el.vkmpersec;
    const vyOut = (state.rdot * uy + state.rvdot * vy) * el.vkmpersec;
    const vzOut = (state.rdot * uz + state.rvdot * vz) * el.vkmpersec;

    var results: [4][2][3]f64 = undefined;
    inline for (0..4) |i| {
        results[i][0][0] = rx[i];
        results[i][0][1] = ry[i];
        results[i][0][2] = rz[i];
        results[i][1][0] = vxOut[i];
        results[i][1][1] = vyOut[i];
        results[i][1][2] = vzOut[i];
    }

    return results;
}

test "Multi-satellite SIMD matches scalar" {
    const allocator = std.testing.allocator;

    const tleStrs = [_][]const u8{
        \\1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995
        \\2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123
        ,
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
        ,
        \\1 55910U 23035C   24187.17717543  .00022897  00000+0  15695-2 0  9993
        \\2 55910  43.9975 313.1680 0011485 277.1866  82.7988 15.05748091 71356
        ,
        \\1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995
        \\2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123
        ,
    };

    var tles: [4]Tle = undefined;
    for (0..4) |i| {
        tles[i] = try Tle.parse(tleStrs[i], allocator);
    }
    defer for (0..4) |i| {
        tles[i].deinit();
    };

    const elsV4 = try initElementsV4(tles, constants.wgs84);

    const testTimes = [_]f64{ 0.0, 10.0, 60.0, 1440.0 };

    for (testTimes) |t| {
        const simdResults = try propagateSatellitesV4(&elsV4, t);

        for (0..4) |i| {
            var sgp4 = try Sgp4.init(tles[i], constants.wgs84);
            const scalarResult = try sgp4.propagate(t);

            const tol = 1e-4;
            for (0..3) |j| {
                try std.testing.expectApproxEqAbs(scalarResult[0][j], simdResults[i][0][j], tol);
                try std.testing.expectApproxEqAbs(scalarResult[1][j], simdResults[i][1][j], tol);
            }
        }
    }
}

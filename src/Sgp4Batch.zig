//! Multi-satellite SIMD batch propagation (4 satellites per batch)

const std = @import("std");
const constants = @import("constants.zig");
const simdMath = @import("simdMath.zig");
const Sgp4 = @import("Sgp4.zig");
const Tle = @import("Tle.zig");

const Vec4 = simdMath.Vec4;
const Elements = Sgp4.Elements;
const Error = Sgp4.Error;

const one: Vec4 = @splat(1.0);
const zero: Vec4 = @splat(0.0);
const half: Vec4 = @splat(0.5);
const quarter: Vec4 = @splat(0.25);
const oneHalf: Vec4 = @splat(1.5);
const two: Vec4 = @splat(2.0);
const eccFloor: Vec4 = @splat(Sgp4.eccentricityFloor);
const tolerance: Vec4 = @splat(1.0e-12);
const clampVal: Vec4 = @splat(0.95);

/// Struct-of-arrays for 4 satellites
pub const ElementsV4 = struct {
    // Gravity constants (splatted)
    xke: Vec4,
    j2: Vec4,
    radiusEarthKm: Vec4,

    // Mean elements
    ecco: Vec4,
    inclo: Vec4,
    nodeo: Vec4,
    argpo: Vec4,
    mo: Vec4,
    bstar: Vec4,
    noUnkozai: Vec4,

    // Trig terms
    sinio: Vec4,
    cosio: Vec4,
    con41: Vec4,
    x1mth2: Vec4,
    x7thm1: Vec4,

    // Secular rates
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

    // Higher-order drag
    d2: Vec4,
    d3: Vec4,
    d4: Vec4,
    t3cof: Vec4,
    t4cof: Vec4,
    t5cof: Vec4,

    // Precomputed
    aBase: Vec4,
    vkmpersec: Vec4,
    isimpMask: Vec4,
    epochJd: Vec4,
};

/// Transpose 4 scalar Elements into Vec4 ElementsV4
pub fn initElementsV4(tles: [4]Tle, grav: constants.Sgp4GravityModel) Error!ElementsV4 {
    var els: [4]Elements = undefined;
    for (0..4) |i| {
        els[i] = try Sgp4.initElements(tles[i], grav);
    }

    var result: ElementsV4 = undefined;

    // Transpose matching fields automatically
    inline for (@typeInfo(ElementsV4).@"struct".fields) |field| {
        if (@hasField(Elements, field.name)) {
            @field(result, field.name) = .{
                @field(els[0], field.name),
                @field(els[1], field.name),
                @field(els[2], field.name),
                @field(els[3], field.name),
            };
        }
    }

    // Special handling
    result.isimpMask = .{
        if (els[0].isimp) 1.0 else 0.0,
        if (els[1].isimp) 1.0 else 0.0,
        if (els[2].isimp) 1.0 else 0.0,
        if (els[3].isimp) 1.0 else 0.0,
    };

    // Splat gravity constants
    result.xke = @splat(grav.xke);
    result.j2 = @splat(grav.j2);
    result.radiusEarthKm = @splat(grav.radiusEarthKm);

    return result;
}

/// Position and velocity as Vec4 (SoA)
pub const PositionVelocityV4 = struct {
    rx: Vec4,
    ry: Vec4,
    rz: Vec4,
    vx: Vec4,
    vy: Vec4,
    vz: Vec4,
};

/// Propagate 4 satellites, returns [4][2][3]f64 (AoS)
pub inline fn propagateSatellitesV4(el: *const ElementsV4, tsince: f64) Error![4][2][3]f64 {
    const pv = try propagateBatchV4Direct(el, @as(Vec4, @splat(tsince)));
    var results: [4][2][3]f64 = undefined;
    inline for (0..4) |i| {
        results[i] = .{ .{ pv.rx[i], pv.ry[i], pv.rz[i] }, .{ pv.vx[i], pv.vy[i], pv.vz[i] } };
    }
    return results;
}

/// Core propagation - returns Vec4 directly (hot path)
pub inline fn propagateBatchV4Direct(el: *const ElementsV4, tsince: Vec4) Error!PositionVelocityV4 {
    // Secular update
    const t2 = tsince * tsince;
    const t3 = t2 * tsince;
    const t4 = t3 * tsince;

    var tempa = one - el.cc1 * tsince;
    var tempe = el.bstar * el.cc4 * tsince;
    var templ = el.t2cof * t2;

    const xmdf = el.mo + el.mdot * tsince;
    const argpdf = el.argpo + el.argpdot * tsince;
    var nodem = el.nodeo + el.nodedot * tsince + el.xnodcf * t2;

    // Higher-order drag (branchless)
    const delomg = el.omgcof * tsince;
    const delmtemp = one + el.eta * simdMath.cosV4(xmdf);
    const delmHo = el.xmcof * (delmtemp * delmtemp * delmtemp - el.delmo);
    const tempHo = delomg + delmHo;
    const hoMask = el.isimpMask == zero;

    var mm = xmdf + @select(f64, hoMask, tempHo, zero);
    var argpm = argpdf - @select(f64, hoMask, tempHo, zero);

    tempa = @select(f64, hoMask, tempa - el.d2 * t2 - el.d3 * t3 - el.d4 * t4, tempa);
    tempe = @select(f64, hoMask, tempe + el.bstar * el.cc5 * (simdMath.sinV4(mm) - el.sinmao), tempe);
    templ = @select(f64, hoMask, templ + el.t3cof * t3 + t4 * (el.t4cof + tsince * el.t5cof), templ);

    const am = el.aBase * tempa * tempa;
    const em = @max(el.ecco - tempe, eccFloor);

    if (@reduce(.Or, em < eccFloor)) return Error.SatelliteDecayed;

    mm = mm + el.noUnkozai * templ;
    nodem = simdMath.modTwoPiV4(nodem);
    argpm = simdMath.modTwoPiV4(argpm);
    mm = simdMath.modTwoPiV4(mm + argpm + nodem - argpm - nodem);

    // Kepler solve
    const temp = one / (am * (one - em * em));
    const argpmSC = simdMath.sincosV4(argpm);
    const axnl = em * argpmSC.cos;
    const aynl = em * argpmSC.sin + temp * el.aycof;
    const xl = simdMath.modTwoPiV4(mm + argpm + nodem + temp * el.xlcof * axnl);

    var u = simdMath.modTwoPiV4(xl - nodem);
    var eo1 = u;
    var sineo1: Vec4 = zero;
    var coseo1: Vec4 = one;

    for (0..10) |_| {
        const sc = simdMath.sincosV4(eo1);
        sineo1 = sc.sin;
        coseo1 = sc.cos;
        var delta = (u - aynl * coseo1 + axnl * sineo1 - eo1) / (one - coseo1 * axnl - sineo1 * aynl);
        delta = @select(f64, delta > clampVal, clampVal, delta);
        delta = @select(f64, delta < -clampVal, -clampVal, delta);
        eo1 = eo1 + delta;
        if (@reduce(.And, @abs(delta) < tolerance)) break;
    }

    const ecose = axnl * coseo1 + aynl * sineo1;
    const esine = axnl * sineo1 - aynl * coseo1;
    const el2 = axnl * axnl + aynl * aynl;
    const pl = am * (one - el2);
    const betal = @sqrt(one - el2);
    const rl = am * (one - ecose);
    const rdotl = @sqrt(am) * esine / rl;
    const rvdotl = @sqrt(pl) / rl;

    const aOverR = am / rl;
    const esineTerm = esine / (one + betal);
    const sinu = aOverR * (sineo1 - aynl - axnl * esineTerm);
    const cosu = aOverR * (coseo1 - axnl + aynl * esineTerm);
    u = simdMath.atan2SIMD(sinu, cosu);

    const sin2u = two * sinu * cosu;
    const cos2u = one - two * sinu * sinu;

    // Short period corrections
    const temp1 = half * el.j2 / pl;
    const temp2 = temp1 / pl;
    const nm = el.xke / simdMath.pow15V4(am);

    const mrt = rl * (one - oneHalf * temp2 * betal * el.con41) + half * temp1 * el.x1mth2 * cos2u;
    const su = u - quarter * temp2 * el.x7thm1 * sin2u;
    const xnode = nodem + oneHalf * temp2 * el.cosio * sin2u;
    const xinc = el.inclo + oneHalf * temp2 * el.cosio * el.sinio * cos2u;
    const mvt = rdotl - nm * temp1 * el.x1mth2 * sin2u / el.xke;
    const rvdot = rvdotl + nm * temp1 * (el.x1mth2 * cos2u + half * el.con41) / el.xke;

    // Position/velocity
    const suSC = simdMath.sincosV4(su);
    const nodeSC = simdMath.sincosV4(xnode);
    const incSC = simdMath.sincosV4(xinc);

    const xmx = -nodeSC.sin * incSC.cos;
    const xmy = nodeSC.cos * incSC.cos;
    const ux = xmx * suSC.sin + nodeSC.cos * suSC.cos;
    const uy = xmy * suSC.sin + nodeSC.sin * suSC.cos;
    const uz = incSC.sin * suSC.sin;
    const vx = xmx * suSC.cos - nodeSC.cos * suSC.sin;
    const vy = xmy * suSC.cos - nodeSC.sin * suSC.sin;
    const vz = incSC.sin * suSC.cos;

    const rScaled = mrt * el.radiusEarthKm;

    return .{
        .rx = rScaled * ux,
        .ry = rScaled * uy,
        .rz = rScaled * uz,
        .vx = (mvt * ux + rvdot * vx) * el.vkmpersec,
        .vy = (mvt * uy + rvdot * vy) * el.vkmpersec,
        .vz = (mvt * uz + rvdot * vz) * el.vkmpersec,
    };
}

test "SIMD matches scalar" {
    const allocator = std.testing.allocator;
    const tleStrs = [_][]const u8{
        "1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995\n2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123",
        "1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998\n2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371",
        "1 55910U 23035C   24187.17717543  .00022897  00000+0  15695-2 0  9993\n2 55910  43.9975 313.1680 0011485 277.1866  82.7988 15.05748091 71356",
        "1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995\n2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123",
    };

    var tles: [4]Tle = undefined;
    defer for (&tles) |*t| t.deinit();
    for (0..4) |i| tles[i] = try Tle.parse(tleStrs[i], allocator);

    const elsV4 = try initElementsV4(tles, constants.wgs84);

    for ([_]f64{ 0.0, 10.0, 60.0, 1440.0 }) |t| {
        const simd = try propagateSatellitesV4(&elsV4, t);
        for (0..4) |i| {
            const scalar = try (try Sgp4.init(tles[i], constants.wgs84)).propagate(t);
            for (0..3) |j| {
                try std.testing.expectApproxEqAbs(scalar[0][j], simd[i][0][j], 2e-3);
                try std.testing.expectApproxEqAbs(scalar[1][j], simd[i][1][j], 2e-3);
            }
        }
    }
}

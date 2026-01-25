//! Multi-satellite SIMD batch propagation (N satellites per batch, N=4 or N=8)

const std = @import("std");
const builtin = @import("builtin");
const constants = @import("constants.zig");
const simdMath = @import("simdMath.zig");
const Sgp4 = @import("Sgp4.zig");
const Tle = @import("Tle.zig");

const Elements = Sgp4.Elements;
const Error = Sgp4.Error;

/// Detect AVX512 at compile time and select optimal batch size
const has_avx512 = blk: {
    const target = builtin.cpu;
    // Check for AVX512F (foundation) which is required for 512-bit operations
    break :blk target.features.isEnabled(@intFromEnum(std.Target.x86.Feature.avx512f));
};

/// Batch size: 8 for AVX512, 4 for AVX2/SSE
pub const BatchSize: usize = if (has_avx512) 8 else 4;

/// Struct-of-arrays for N satellites (generic over batch size)
pub fn BatchElements(comptime N: usize) type {
    const Vec = simdMath.VecN(N);

    return struct {
        const Self = @This();

        // Gravity constants (splatted)
        xke: Vec,
        j2: Vec,
        radiusEarthKm: Vec,

        // Mean elements
        ecco: Vec,
        inclo: Vec,
        nodeo: Vec,
        argpo: Vec,
        mo: Vec,
        bstar: Vec,
        noUnkozai: Vec,

        // Trig terms
        sinio: Vec,
        cosio: Vec,
        con41: Vec,
        x1mth2: Vec,
        x7thm1: Vec,

        // Secular rates
        mdot: Vec,
        argpdot: Vec,
        nodedot: Vec,

        // Drag coefficients
        cc1: Vec,
        cc4: Vec,
        cc5: Vec,
        t2cof: Vec,
        omgcof: Vec,
        xnodcf: Vec,
        xlcof: Vec,
        xmcof: Vec,
        aycof: Vec,
        eta: Vec,
        delmo: Vec,
        sinmao: Vec,

        // Higher-order drag
        d2: Vec,
        d3: Vec,
        d4: Vec,
        t3cof: Vec,
        t4cof: Vec,
        t5cof: Vec,

        // Precomputed
        aBase: Vec,
        vkmpersec: Vec,
        isimpMask: Vec,
        epochJd: Vec,
    };
}

/// Position and velocity as VecN (SoA)
pub fn PositionVelocity(comptime N: usize) type {
    const Vec = simdMath.VecN(N);
    return struct {
        rx: Vec,
        ry: Vec,
        rz: Vec,
        vx: Vec,
        vy: Vec,
        vz: Vec,
    };
}

/// Transpose N scalar Elements into BatchElements(N)
pub fn initBatchElements(comptime N: usize, tles: [N]Tle, grav: constants.Sgp4GravityModel) Error!BatchElements(N) {
    var els: [N]Elements = undefined;
    for (0..N) |i| {
        els[i] = try Sgp4.initElements(tles[i], grav);
    }

    var result: BatchElements(N) = undefined;

    // Transpose matching fields automatically
    inline for (@typeInfo(BatchElements(N)).@"struct".fields) |field| {
        if (@hasField(Elements, field.name)) {
            var arr: [N]f64 = undefined;
            inline for (0..N) |i| {
                arr[i] = @field(els[i], field.name);
            }
            @field(result, field.name) = arr;
        }
    }

    // Special handling for isimp boolean -> mask
    var isimpArr: [N]f64 = undefined;
    inline for (0..N) |i| {
        isimpArr[i] = if (els[i].isimp) 1.0 else 0.0;
    }
    result.isimpMask = isimpArr;

    // Splat gravity constants
    result.xke = @splat(grav.xke);
    result.j2 = @splat(grav.j2);
    result.radiusEarthKm = @splat(grav.radiusEarthKm);

    return result;
}

/// Propagate N satellites, returns [N][2][3]f64 (AoS)
pub fn propagateSatellites(comptime N: usize, el: *const BatchElements(N), tsince: f64) Error![N][2][3]f64 {
    const Vec = simdMath.VecN(N);
    const pv = try propagateBatchDirect(N, el, @as(Vec, @splat(tsince)));
    var results: [N][2][3]f64 = undefined;
    inline for (0..N) |i| {
        results[i] = .{ .{ pv.rx[i], pv.ry[i], pv.rz[i] }, .{ pv.vx[i], pv.vy[i], pv.vz[i] } };
    }
    return results;
}

/// Core propagation, returns VecN directly (hot path)
pub fn propagateBatchDirect(comptime N: usize, el: *const BatchElements(N), tsince: simdMath.VecN(N)) Error!PositionVelocity(N) {
    const Vec = simdMath.VecN(N);

    const one: Vec = @splat(1.0);
    const zero: Vec = @splat(0.0);
    const half: Vec = @splat(0.5);
    const quarter: Vec = @splat(0.25);
    const oneHalf: Vec = @splat(1.5);
    const two: Vec = @splat(2.0);
    const eccFloor: Vec = @splat(Sgp4.eccentricityFloor);
    const tolerance: Vec = @splat(1.0e-12);
    const clampVal: Vec = @splat(0.95);

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
    const delmtemp = one + el.eta * simdMath.cosN(N, xmdf);
    const delmHo = el.xmcof * (delmtemp * delmtemp * delmtemp - el.delmo);
    const tempHo = delomg + delmHo;
    const hoMask = el.isimpMask == zero;

    var mm = xmdf + @select(f64, hoMask, tempHo, zero);
    var argpm = argpdf - @select(f64, hoMask, tempHo, zero);

    tempa = @select(f64, hoMask, tempa - el.d2 * t2 - el.d3 * t3 - el.d4 * t4, tempa);
    tempe = @select(f64, hoMask, tempe + el.bstar * el.cc5 * (simdMath.sinN(N, mm) - el.sinmao), tempe);
    templ = @select(f64, hoMask, templ + el.t3cof * t3 + t4 * (el.t4cof + tsince * el.t5cof), templ);

    const am = el.aBase * tempa * tempa;
    const em = @max(el.ecco - tempe, eccFloor);

    if (@reduce(.Or, em < eccFloor)) return Error.SatelliteDecayed;

    mm = simdMath.modTwoPiN(N, mm + el.noUnkozai * templ);
    nodem = simdMath.modTwoPiN(N, nodem);
    argpm = simdMath.modTwoPiN(N, argpm);

    // Kepler solve
    const temp = one / (am * (one - em * em));
    const argpmSC = simdMath.sincosN(N, argpm);
    const axnl = em * argpmSC.cos;
    const aynl = em * argpmSC.sin + temp * el.aycof;
    const xl = mm + argpm + nodem + temp * el.xlcof * axnl;

    var u = xl - nodem;
    var eo1 = u;
    var sineo1: Vec = zero;
    var coseo1: Vec = one;

    for (0..10) |_| {
        const sc = simdMath.sincosN(N, eo1);
        sineo1 = sc.sin;
        coseo1 = sc.cos;
        const delta = (u - aynl * coseo1 + axnl * sineo1 - eo1) / (one - coseo1 * axnl - sineo1 * aynl);
        eo1 = eo1 + @max(-clampVal, @min(clampVal, delta));
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
    u = simdMath.atan2N(N, sinu, cosu);

    const sin2u = two * sinu * cosu;
    const cos2u = one - two * sinu * sinu;

    // Short period corrections
    const temp1 = half * el.j2 / pl;
    const temp2 = temp1 / pl;
    const nm = el.xke / simdMath.pow15N(N, am);

    const mrt = rl * (one - oneHalf * temp2 * betal * el.con41) + half * temp1 * el.x1mth2 * cos2u;
    const su = u - quarter * temp2 * el.x7thm1 * sin2u;
    const xnode = nodem + oneHalf * temp2 * el.cosio * sin2u;
    const xinc = el.inclo + oneHalf * temp2 * el.cosio * el.sinio * cos2u;
    const mvt = rdotl - nm * temp1 * el.x1mth2 * sin2u / el.xke;
    const rvdot = rvdotl + nm * temp1 * (el.x1mth2 * cos2u + oneHalf * el.con41) / el.xke;

    // Position/velocity
    const suSC = simdMath.sincosN(N, su);
    const nodeSC = simdMath.sincosN(N, xnode);
    const incSC = simdMath.sincosN(N, xinc);

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

    const els = try initBatchElements(4, tles, constants.wgs84);

    // Tolerances: position < 1m, velocity < 1µm/s
    for ([_]f64{ 0.0, 10.0, 60.0, 1440.0 }) |t| {
        const simd = try propagateSatellites(4, &els, t);
        for (0..4) |i| {
            const scalar = try (try Sgp4.init(tles[i], constants.wgs84)).propagate(t);
            for (0..3) |j| {
                try std.testing.expectApproxEqAbs(scalar[0][j], simd[i][0][j], 1e-3);
                try std.testing.expectApproxEqAbs(scalar[1][j], simd[i][1][j], 1e-6);
            }
        }
    }
}

test "Vec8 matches Vec4 output" {
    const allocator = std.testing.allocator;
    const tleStrs = [_][]const u8{
        "1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995\n2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123",
        "1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998\n2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371",
        "1 55910U 23035C   24187.17717543  .00022897  00000+0  15695-2 0  9993\n2 55910  43.9975 313.1680 0011485 277.1866  82.7988 15.05748091 71356",
        "1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995\n2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123",
    };

    var tles4: [4]Tle = undefined;
    defer for (&tles4) |*t| t.deinit();
    for (0..4) |i| tles4[i] = try Tle.parse(tleStrs[i], allocator);

    // Create 8 TLEs by duplicating the 4
    var tles8: [8]Tle = undefined;
    defer for (&tles8) |*t| t.deinit();
    for (0..8) |i| tles8[i] = try Tle.parse(tleStrs[i % 4], allocator);

    const elsV4 = try initBatchElements(4, tles4, constants.wgs84);
    const elsV8 = try initBatchElements(8, tles8, constants.wgs84);

    // Test that first 4 elements of Vec8 batch match Vec4 batch exactly
    for ([_]f64{ 0.0, 10.0, 60.0, 1440.0 }) |t| {
        const r4 = try propagateSatellites(4, &elsV4, t);
        const r8 = try propagateSatellites(8, &elsV8, t);

        for (0..4) |i| {
            for (0..3) |j| {
                try std.testing.expectApproxEqAbs(r4[i][0][j], r8[i][0][j], 1e-10);
                try std.testing.expectApproxEqAbs(r4[i][1][j], r8[i][1][j], 1e-10);
            }
        }
        // Also verify the duplicated satellites match
        for (0..4) |i| {
            for (0..3) |j| {
                try std.testing.expectApproxEqAbs(r8[i][0][j], r8[i + 4][0][j], 1e-10);
                try std.testing.expectApproxEqAbs(r8[i][1][j], r8[i + 4][1][j], 1e-10);
            }
        }
    }
}

test "Vallado AIAA 2006-6753 near-earth vectors" {
    // Test against official Vallado reference vectors (AIAA 2006-6753, tcppver.out)
    const allocator = std.testing.allocator;

    // Satellite 00005 - Vanguard 1 (high eccentricity test case)
    const tle00005 =
        \\1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753
        \\2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667
    ;

    // Satellite 06251 - Delta 1 Deb (normal drag test case)
    const tle06251 =
        \\1 06251U 62025E   06176.82412014  .00008885  00000-0  12808-3 0  3985
        \\2 06251  58.0579  54.0425 0030035 139.1568 221.1854 15.56387291  6774
    ;

    var t1 = try Tle.parse(tle00005, allocator);
    defer t1.deinit();
    var t2 = try Tle.parse(tle06251, allocator);
    defer t2.deinit();

    const tles: [4]Tle = .{ t1, t2, t1, t2 };
    const els = try initBatchElements(4, tles, constants.wgs72);

    // Vallado reference: sat 00005 at t=0
    // r = (7022.46529266, -1400.08296755, 0.03995155) km
    // v = (1.893841015, 6.405893759, 4.534807250) km/s
    {
        const result = try propagateSatellites(4, &els, 0.0);
        try std.testing.expectApproxEqAbs(@as(f64, 7022.46529266), result[0][0][0], 0.01);
        try std.testing.expectApproxEqAbs(@as(f64, -1400.08296755), result[0][0][1], 0.01);
        try std.testing.expectApproxEqAbs(@as(f64, 0.03995155), result[0][0][2], 0.01);
        try std.testing.expectApproxEqAbs(@as(f64, 1.893841015), result[0][1][0], 1e-6);
        try std.testing.expectApproxEqAbs(@as(f64, 6.405893759), result[0][1][1], 1e-6);
        try std.testing.expectApproxEqAbs(@as(f64, 4.534807250), result[0][1][2], 1e-6);
    }

    // Vallado reference: sat 00005 at t=360 min
    // r = (-7154.03120202, -3783.17682504, -3536.19412294) km
    // v = (4.741887409, -4.151817765, -2.093935425) km/s
    {
        const result = try propagateSatellites(4, &els, 360.0);
        try std.testing.expectApproxEqAbs(@as(f64, -7154.03120202), result[0][0][0], 0.01);
        try std.testing.expectApproxEqAbs(@as(f64, -3783.17682504), result[0][0][1], 0.01);
        try std.testing.expectApproxEqAbs(@as(f64, -3536.19412294), result[0][0][2], 0.01);
        try std.testing.expectApproxEqAbs(@as(f64, 4.741887409), result[0][1][0], 1e-6);
        try std.testing.expectApproxEqAbs(@as(f64, -4.151817765), result[0][1][1], 1e-6);
        try std.testing.expectApproxEqAbs(@as(f64, -2.093935425), result[0][1][2], 1e-6);
    }

    // Vallado reference: sat 06251 at t=0
    // r = (3988.31022699, 5498.96657235, 0.90055879) km
    // v = (-3.290032738, 2.357652820, 6.496623475) km/s
    {
        const result = try propagateSatellites(4, &els, 0.0);
        try std.testing.expectApproxEqAbs(@as(f64, 3988.31022699), result[1][0][0], 0.01);
        try std.testing.expectApproxEqAbs(@as(f64, 5498.96657235), result[1][0][1], 0.01);
        try std.testing.expectApproxEqAbs(@as(f64, 0.90055879), result[1][0][2], 0.01);
        try std.testing.expectApproxEqAbs(@as(f64, -3.290032738), result[1][1][0], 1e-6);
        try std.testing.expectApproxEqAbs(@as(f64, 2.357652820), result[1][1][1], 1e-6);
        try std.testing.expectApproxEqAbs(@as(f64, 6.496623475), result[1][1][2], 1e-6);
    }
}

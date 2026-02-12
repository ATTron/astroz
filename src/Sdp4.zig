//! SDP4 Deep-Space Orbit Propagator
//! Extension of SGP4 for satellites with orbital period > 225 minutes.
//! Adds lunar/solar gravitational perturbations and orbital resonance effects.
//! Reference: Vallado, Crawford, Hujsak, Kelso - AIAA 2006-6753

const std = @import("std");
const constants = @import("constants.zig");
const Datetime = @import("Datetime.zig");
const Tle = @import("Tle.zig");
const Sgp4 = @import("Sgp4.zig");

const Sdp4 = @This();

// ── Deep-space constants ────────────────────────────────────────────────
const zes = 0.01675;
const zel = 0.05490;
const c1ss = 2.9864797e-6;
const c1l = 4.7968065e-7;
const zsinis = 0.39785416;
const zcosis = 0.91744867;
const zcosgs = 0.1945905;
const zsings = -0.98088458;
const zns = 1.19459e-5;
const znl = 1.5835218e-4;

const q22 = 1.7891679e-6;
const q31 = 2.1460748e-6;
const q33 = 2.2123015e-7;

const root22 = 1.7891679e-6;
const root32 = 3.7393792e-7;
const root44 = 7.3636953e-9;
const root52 = 1.1428639e-7;
const root54 = 2.1765803e-9;

const rptim = 4.37526908801129966e-3;

// GEO resonance phase constants
const fasx2 = 0.13130908;
const fasx4 = 2.8843198;
const fasx6 = 0.37448087;
const g22 = 5.7686396;
const g32 = 0.95240898;
const g44 = 1.8014998;
const g52 = 1.0508330;
const g54 = 4.4108898;

const stepp = 720.0;
const stepn = -720.0;
const step2 = 259200.0;

pub const Error = Sgp4.Error;

pub const Elements = struct {
    sgp4: Sgp4.Elements,

    // Solar perturbation coefficients
    se2: f64,
    se3: f64,
    si2: f64,
    si3: f64,
    sl2: f64,
    sl3: f64,
    sl4: f64,
    sgh2: f64,
    sgh3: f64,
    sgh4: f64,
    sh2: f64,
    sh3: f64,

    // Lunar perturbation coefficients
    ee2: f64,
    e3: f64,
    xi2: f64,
    xi3: f64,
    xl2: f64,
    xl3: f64,
    xl4: f64,
    xgh2: f64,
    xgh3: f64,
    xgh4: f64,
    xh2: f64,
    xh3: f64,

    zmol: f64,
    zmos: f64,
    dedt: f64,
    didt: f64,
    dmdt: f64,
    domdt: f64,
    dnodt: f64,

    // Resonance: 0=none, 1=GEO, 2=half-day
    irez: u8,

    // Half-day resonance coefficients
    d2201: f64,
    d2211: f64,
    d3210: f64,
    d3222: f64,
    d4410: f64,
    d4422: f64,
    d5220: f64,
    d5232: f64,
    d5421: f64,
    d5433: f64,

    // GEO resonance coefficients
    del1: f64,
    del2: f64,
    del3: f64,

    // Resonance integration state
    xlamo: f64,
    xfact: f64,

    gsto: f64,
};

elements: Elements,

pub fn init(tle: Tle, grav: constants.Sgp4GravityModel) Error!Sdp4 {
    return .{ .elements = try initElements(tle, grav) };
}

pub fn propagate(self: *const Sdp4, tsince: f64) Error![2][3]f64 {
    return propagateElements(&self.elements, tsince);
}

pub fn initElements(tle: Tle, grav: constants.Sgp4GravityModel) Error!Elements {
    const meanEl = Sgp4.extractMeanElements(tle);
    if (meanEl.ecco < 0.0 or meanEl.ecco >= 1.0) return Error.InvalidEccentricity;

    const rec = Sgp4.recoverMeanMotion(meanEl, grav);
    const rp = rec.a * (1.0 - meanEl.ecco);
    if (rp < 1.0) return Error.SatelliteDecayed;

    const trig = Sgp4.computeTrigTerms(meanEl.inclo);
    const poly = Sgp4.computePolyTerms(trig);
    const secRates = Sgp4.computeSecularRates(meanEl, rec, trig, poly, grav);
    const drag = Sgp4.computeDragCoefficients(meanEl, rec, trig, poly, grav);

    const fullYear: u16 = if (tle.firstLine.epochYear < 57)
        2000 + tle.firstLine.epochYear
    else
        1900 + tle.firstLine.epochYear;
    const epochJd = Datetime.yearDoyToJulianDate(fullYear, tle.firstLine.epochDay);
    const gsto = gstime(epochJd);

    // Deep space initialization
    const day = epochJd - 2415020.0;
    const dc = dscom(meanEl, rec, trig, day, gsto);
    const di = dsinit(meanEl, rec, trig, secRates, &dc, gsto, grav);

    const sgp4El = Sgp4.Elements{
        .grav = grav,
        .epochJd = epochJd,
        .noKozai = meanEl.noKozai,
        .ecco = meanEl.ecco,
        .inclo = meanEl.inclo,
        .nodeo = meanEl.nodeo,
        .argpo = meanEl.argpo,
        .mo = meanEl.mo,
        .bstar = meanEl.bstar,
        .noUnkozai = rec.noUnkozai,
        .a = rec.a,
        .sinio = trig.sinio,
        .cosio = trig.cosio,
        .cosio2 = trig.cosio2,
        .cosio4 = trig.cosio4,
        .con41 = poly.con41,
        .con42 = poly.con42,
        .x1mth2 = poly.x1mth2,
        .x7thm1 = poly.x7thm1,
        .mdot = secRates.mdot,
        .argpdot = secRates.argpdot,
        .nodedot = secRates.nodedot,
        .cc1 = drag.cc1,
        .cc4 = drag.cc4,
        .cc5 = drag.cc5,
        .t2cof = drag.t2cof,
        .omgcof = drag.omgcof,
        .xnodcf = drag.xnodcf,
        .xlcof = drag.xlcof,
        .xmcof = drag.xmcof,
        .aycof = drag.aycof,
        .eta = drag.eta,
        .delmo = drag.delmo,
        .sinmao = drag.sinmao,
        .d2 = 0,
        .d3 = 0,
        .d4 = 0,
        .t3cof = 0,
        .t4cof = 0,
        .t5cof = 0,
        .aBase = blk: {
            const ratio = grav.xke / rec.noUnkozai;
            break :blk std.math.cbrt(ratio * ratio);
        },
        .vkmpersec = grav.xke * grav.radiusEarthKm / 60.0,
        .isimp = true,
    };

    return Elements{
        .sgp4 = sgp4El,
        .se2 = dc.se2,
        .se3 = dc.se3,
        .si2 = dc.si2,
        .si3 = dc.si3,
        .sl2 = dc.sl2,
        .sl3 = dc.sl3,
        .sl4 = dc.sl4,
        .sgh2 = dc.sgh2,
        .sgh3 = dc.sgh3,
        .sgh4 = dc.sgh4,
        .sh2 = dc.sh2,
        .sh3 = dc.sh3,
        .ee2 = dc.ee2,
        .e3 = dc.e3,
        .xi2 = dc.xi2,
        .xi3 = dc.xi3,
        .xl2 = dc.xl2,
        .xl3 = dc.xl3,
        .xl4 = dc.xl4,
        .xgh2 = dc.xgh2,
        .xgh3 = dc.xgh3,
        .xgh4 = dc.xgh4,
        .xh2 = dc.xh2,
        .xh3 = dc.xh3,
        .zmol = dc.zmol,
        .zmos = dc.zmos,
        .dedt = di.dedt,
        .didt = di.didt,
        .dmdt = di.dmdt,
        .domdt = di.domdt,
        .dnodt = di.dnodt,
        .irez = di.irez,
        .d2201 = di.d2201,
        .d2211 = di.d2211,
        .d3210 = di.d3210,
        .d3222 = di.d3222,
        .d4410 = di.d4410,
        .d4422 = di.d4422,
        .d5220 = di.d5220,
        .d5232 = di.d5232,
        .d5421 = di.d5421,
        .d5433 = di.d5433,
        .del1 = di.del1,
        .del2 = di.del2,
        .del3 = di.del3,
        .xlamo = di.xlamo,
        .xfact = di.xfact,
        .gsto = gsto,
    };
}

// ── Greenwich Sidereal Time ─────────────────────────────────────────────
pub fn gstime(jdut1: f64) f64 {
    const tut1 = (jdut1 - 2451545.0) / 36525.0;
    var temp = -6.2e-6 * tut1 * tut1 * tut1 +
        0.093104 * tut1 * tut1 +
        (876600.0 * 3600.0 + 8640184.812866) * tut1 + 67310.54841;
    temp = @mod(temp * constants.deg2rad / 240.0, constants.twoPi);
    if (temp < 0.0) temp += constants.twoPi;
    return temp;
}

// ── dscom ───────────────────────────────────────────────────────────────
const DscomResult = struct {
    se2: f64,
    se3: f64,
    si2: f64,
    si3: f64,
    sl2: f64,
    sl3: f64,
    sl4: f64,
    sgh2: f64,
    sgh3: f64,
    sgh4: f64,
    sh2: f64,
    sh3: f64,
    ee2: f64,
    e3: f64,
    xi2: f64,
    xi3: f64,
    xl2: f64,
    xl3: f64,
    xl4: f64,
    xgh2: f64,
    xgh3: f64,
    xgh4: f64,
    xh2: f64,
    xh3: f64,
    zmol: f64,
    zmos: f64,
    // Intermediate values for dsinit
    snodm: f64,
    cnodm: f64,
    sinim: f64,
    cosim: f64,
    sinomm: f64,
    cosomm: f64,
    emsq: f64,
    rtemsq: f64,
    nm: f64,
    gam: f64,
    ss1: f64,
    ss2: f64,
    ss3: f64,
    ss4: f64,
    ss5: f64,
    ss6: f64,
    ss7: f64,
    s1: f64,
    s2: f64,
    s3: f64,
    s4: f64,
    s5: f64,
    s6: f64,
    s7: f64,
    sz1: f64,
    sz2: f64,
    sz3: f64,
    sz11: f64,
    sz12: f64,
    sz13: f64,
    sz21: f64,
    sz22: f64,
    sz23: f64,
    sz31: f64,
    sz32: f64,
    sz33: f64,
    z1: f64,
    z2: f64,
    z3: f64,
    z11: f64,
    z12: f64,
    z13: f64,
    z21: f64,
    z22: f64,
    z23: f64,
    z31: f64,
    z32: f64,
    z33: f64,
};

fn dscom(
    el: Sgp4.MeanElements,
    rec: Sgp4.RecoveredElements,
    trig: Sgp4.TrigTerms,
    day: f64,
    gsto: f64,
) DscomResult {
    var r: DscomResult = std.mem.zeroes(DscomResult);

    r.nm = rec.noUnkozai;
    r.snodm = @sin(el.nodeo);
    r.cnodm = @cos(el.nodeo);
    r.sinomm = @sin(el.argpo);
    r.cosomm = @cos(el.argpo);
    r.sinim = trig.sinio;
    r.cosim = trig.cosio;
    r.emsq = el.ecco * el.ecco;
    r.rtemsq = @sqrt(1.0 - r.emsq);

    // Lunar node and inclination
    const xnodce = @mod(4.5236020 - 9.2422029e-4 * day, constants.twoPi);
    const stem = @sin(xnodce);
    const ctem = @cos(xnodce);
    const zcosil = 0.91375164 - 0.03568096 * ctem;
    const zsinil = @sqrt(1.0 - zcosil * zcosil);
    const zsinhl = 0.089683511 * stem / zsinil;
    const zcoshl = @sqrt(1.0 - zsinhl * zsinhl);
    r.gam = 5.8351514 + 0.0019443680 * day;
    var zx = 0.39785416 * stem / zsinil;
    const zy = zcoshl * ctem + 0.91744867 * zsinhl * stem;
    zx = std.math.atan2(zx, zy);
    zx += r.gam - xnodce;
    const zcosgl = @cos(zx);
    const zsingl = @sin(zx);

    const xnoi = 1.0 / r.nm;
    const betasq = 1.0 - r.emsq;

    // Two-pass loop: pass 1 = solar, pass 2 = lunar
    var zcosg: f64 = zcosgs;
    var zsing: f64 = zsings;
    var zcosi: f64 = zcosis;
    var zsini: f64 = zsinis;
    var zcosh: f64 = r.cnodm;
    var zsinh: f64 = r.snodm;
    var cc: f64 = c1ss;

    var lsflg: u8 = 1;
    while (lsflg <= 2) : (lsflg += 1) {
        const a1 = zcosg * zcosh + zsing * zcosi * zsinh;
        const a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
        const a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
        const a8 = zsing * zsini;
        const a9 = zsing * zsinh + zcosg * zcosi * zcosh;
        const a10 = zcosg * zsini;
        const a2 = r.cosim * a7 + r.sinim * a8;
        const a4 = r.cosim * a9 + r.sinim * a10;
        const a5 = -r.sinim * a7 + r.cosim * a8;
        const a6 = -r.sinim * a9 + r.cosim * a10;

        const x1 = a1 * r.cosomm + a2 * r.sinomm;
        const x2 = a3 * r.cosomm + a4 * r.sinomm;
        const x3 = -a1 * r.sinomm + a2 * r.cosomm;
        const x4 = -a3 * r.sinomm + a4 * r.cosomm;
        const x5 = a5 * r.sinomm;
        const x6 = a6 * r.sinomm;
        const x7 = a5 * r.cosomm;
        const x8 = a6 * r.cosomm;

        const z31v = 12.0 * x1 * x1 - 3.0 * x3 * x3;
        const z32v = 24.0 * x1 * x2 - 6.0 * x3 * x4;
        const z33v = 12.0 * x2 * x2 - 3.0 * x4 * x4;
        const z1v = 3.0 * (a1 * a1 + a2 * a2) + z31v * r.emsq;
        const z2v = 6.0 * (a1 * a3 + a2 * a4) + z32v * r.emsq;
        const z3v = 3.0 * (a3 * a3 + a4 * a4) + z33v * r.emsq;
        const z11v = -6.0 * a1 * a5 + r.emsq * (-24.0 * x1 * x7 - 6.0 * x3 * x5);
        const z12v = -6.0 * (a1 * a6 + a3 * a5) + r.emsq * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
        const z13v = -6.0 * a3 * a6 + r.emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
        const z21v = 6.0 * a2 * a5 + r.emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
        const z22v = 6.0 * (a4 * a5 + a2 * a6) + r.emsq * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
        const z23v = 6.0 * a4 * a6 + r.emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);

        const z1t = z1v + z1v + betasq * z31v;
        const z2t = z2v + z2v + betasq * z32v;
        const z3t = z3v + z3v + betasq * z33v;

        const s3v = cc * xnoi;
        const s2v = -0.5 * s3v / r.rtemsq;
        const s4v = s3v * r.rtemsq;
        const s1v = -15.0 * el.ecco * s4v;
        const s5v = x1 * x3 + x2 * x4;
        const s6v = x2 * x3 + x1 * x4;
        const s7v = x2 * x4 - x1 * x3;

        if (lsflg == 1) {
            // Store solar terms
            r.ss1 = s1v;
            r.ss2 = s2v;
            r.ss3 = s3v;
            r.ss4 = s4v;
            r.ss5 = s5v;
            r.ss6 = s6v;
            r.ss7 = s7v;
            r.sz1 = z1t;
            r.sz2 = z2t;
            r.sz3 = z3t;
            r.sz11 = z11v;
            r.sz12 = z12v;
            r.sz13 = z13v;
            r.sz21 = z21v;
            r.sz22 = z22v;
            r.sz23 = z23v;
            r.sz31 = z31v;
            r.sz32 = z32v;
            r.sz33 = z33v;

            // Solar perturbation coefficients
            r.se2 = 2.0 * s1v * s6v;
            r.se3 = 2.0 * s1v * s7v;
            r.si2 = 2.0 * s2v * z12v;
            r.si3 = 2.0 * s2v * (z13v - z11v);
            r.sl2 = -2.0 * s3v * z2t;
            r.sl3 = -2.0 * s3v * (z3t - z1t);
            r.sl4 = -2.0 * s3v * (-21.0 - 9.0 * r.emsq) * zes;
            r.sgh2 = 2.0 * s4v * z32v;
            r.sgh3 = 2.0 * s4v * (z33v - z31v);
            r.sgh4 = -18.0 * s4v * zes;
            r.sh2 = -2.0 * s2v * z22v;
            r.sh3 = -2.0 * s2v * (z23v - z21v);

            // Switch to lunar terms for pass 2
            zcosg = zcosgl;
            zsing = zsingl;
            zcosi = zcosil;
            zsini = zsinil;
            zcosh = zcoshl * r.cnodm + zsinhl * r.snodm;
            zsinh = r.snodm * zcoshl - r.cnodm * zsinhl;
            cc = c1l;
        } else {
            // Store lunar terms
            r.s1 = s1v;
            r.s2 = s2v;
            r.s3 = s3v;
            r.s4 = s4v;
            r.s5 = s5v;
            r.s6 = s6v;
            r.s7 = s7v;
            r.z1 = z1t;
            r.z2 = z2t;
            r.z3 = z3t;
            r.z11 = z11v;
            r.z12 = z12v;
            r.z13 = z13v;
            r.z21 = z21v;
            r.z22 = z22v;
            r.z23 = z23v;
            r.z31 = z31v;
            r.z32 = z32v;
            r.z33 = z33v;

            // Lunar perturbation coefficients
            r.ee2 = 2.0 * s1v * s6v;
            r.e3 = 2.0 * s1v * s7v;
            r.xi2 = 2.0 * s2v * z12v;
            r.xi3 = 2.0 * s2v * (z13v - z11v);
            r.xl2 = -2.0 * s3v * z2t;
            r.xl3 = -2.0 * s3v * (z3t - z1t);
            r.xl4 = -2.0 * s3v * (-21.0 - 9.0 * r.emsq) * zel;
            r.xgh2 = 2.0 * s4v * z32v;
            r.xgh3 = 2.0 * s4v * (z33v - z31v);
            r.xgh4 = -18.0 * s4v * zel;
            r.xh2 = -2.0 * s2v * z22v;
            r.xh3 = -2.0 * s2v * (z23v - z21v);
        }
    }

    r.zmol = @mod(4.7199672 + 0.22997150 * day - r.gam, constants.twoPi);
    r.zmos = @mod(6.2565837 + 0.017201977 * day, constants.twoPi);
    _ = gsto;

    return r;
}

// ── dsinit ──────────────────────────────────────────────────────────────
const DsinitResult = struct {
    dedt: f64,
    didt: f64,
    dmdt: f64,
    domdt: f64,
    dnodt: f64,
    irez: u8,
    d2201: f64,
    d2211: f64,
    d3210: f64,
    d3222: f64,
    d4410: f64,
    d4422: f64,
    d5220: f64,
    d5232: f64,
    d5421: f64,
    d5433: f64,
    del1: f64,
    del2: f64,
    del3: f64,
    xlamo: f64,
    xfact: f64,
};

fn dsinit(
    el: Sgp4.MeanElements,
    rec: Sgp4.RecoveredElements,
    trig: Sgp4.TrigTerms,
    secRates: Sgp4.SecularRates,
    dc: *const DscomResult,
    gsto: f64,
    grav: constants.Sgp4GravityModel,
) DsinitResult {
    _ = grav;
    var di: DsinitResult = std.mem.zeroes(DsinitResult);

    const eosq = el.ecco * el.ecco;
    const cosisq = trig.cosio2;
    const sini2 = trig.sinio * trig.sinio;
    const xpidot = secRates.argpdot + secRates.nodedot;

    // Solar secular perturbation rates
    const ses = dc.ss1 * zns * dc.ss5;
    const sis = dc.ss2 * zns * (dc.sz11 + dc.sz13);
    const sls = -zns * dc.ss3 * (dc.sz1 + dc.sz3 - 14.0 - 6.0 * dc.emsq);
    const sghs = dc.ss4 * zns * (dc.sz31 + dc.sz33 - 6.0);
    var shs = -zns * dc.ss2 * (dc.sz21 + dc.sz23);

    // Near-equatorial fix (inclination < 3° or > 177°)
    const inclm = el.inclo;
    if (inclm < 5.2359877e-2 or inclm > std.math.pi - 5.2359877e-2)
        shs = 0.0;
    if (dc.sinim != 0.0)
        shs = shs / dc.sinim;
    const sgs = sghs - dc.cosim * shs;

    // Lunar secular perturbation rates (added to solar)
    di.dedt = ses + dc.s1 * znl * dc.s5;
    di.didt = sis + dc.s2 * znl * (dc.z11 + dc.z13);
    di.dmdt = sls - znl * dc.s3 * (dc.z1 + dc.z3 - 14.0 - 6.0 * dc.emsq);
    const sghl = dc.s4 * znl * (dc.z31 + dc.z33 - 6.0);
    var shll = -znl * dc.s2 * (dc.z21 + dc.z23);

    // Near-equatorial fix
    if (inclm < 5.2359877e-2 or inclm > std.math.pi - 5.2359877e-2)
        shll = 0.0;

    di.domdt = sgs + sghl;
    di.dnodt = shs;
    if (dc.sinim != 0.0) {
        di.domdt -= dc.cosim / dc.sinim * shll;
        di.dnodt += shll / dc.sinim;
    }

    // Resonance detection
    if (rec.noUnkozai >= 0.00826 and rec.noUnkozai <= 0.00924 and el.ecco >= 0.5) {
        di.irez = 2; // half-day
    } else if (rec.noUnkozai >= 0.0034906585 and rec.noUnkozai <= 0.0052359877) {
        di.irez = 1; // GEO synchronous
    } else {
        di.irez = 0;
    }

    if (di.irez == 1) {
        // Geosynchronous resonance
        const g200 = 1.0 + eosq * (-2.5 + 0.8125 * eosq);
        const g310 = 1.0 + 2.0 * eosq;
        const g300 = 1.0 + eosq * (-6.0 + 6.60937 * eosq);
        const f220 = 0.75 * (1.0 + trig.cosio) * (1.0 + trig.cosio);
        const f311 = 0.9375 * sini2 * (1.0 + 3.0 * trig.cosio) - 0.75 * (1.0 + trig.cosio);
        var f330 = 1.0 + trig.cosio;
        f330 = 1.875 * f330 * f330 * f330;
        const aonv = 1.0 / rec.a;
        const temp1_geo = 3.0 * dc.nm * dc.nm * aonv * aonv;
        // Vallado cascading: del1 uses temp*aonv, del2 uses temp, del3 uses temp*aonv
        di.del2 = 2.0 * temp1_geo * f220 * g200 * q22;
        di.del3 = 3.0 * temp1_geo * f330 * g300 * q33 * aonv;
        di.del1 = temp1_geo * f311 * g310 * q31 * aonv;
        di.xlamo = @mod(el.mo + el.nodeo + el.argpo - gsto, constants.twoPi);
        di.xfact = secRates.mdot + xpidot - rptim + di.dmdt + di.domdt + di.dnodt - rec.noUnkozai;
    } else if (di.irez == 2) {
        // Half-day resonance — follow Vallado C++ exactly
        const g201 = -0.306 - (el.ecco - 0.64) * 0.440;
        const g211 = eccoG(el.ecco, &[_]f64{ 3.616, -13.2470, 16.2900 }, &[_]f64{ -72.099, 331.819, -508.738, 266.724 });
        const g310 = eccoG(el.ecco, &[_]f64{ -19.302, 117.3900, -228.4190, 156.591 }, &[_]f64{ -346.844, 1582.851, -2415.925, 1246.113 });
        const g322 = eccoG(el.ecco, &[_]f64{ -18.9068, 109.7927, -214.6334, 146.5816 }, &[_]f64{ -342.585, 1554.908, -2366.899, 1215.972 });
        const g410 = eccoG(el.ecco, &[_]f64{ -41.122, 242.6940, -471.0940, 313.953 }, &[_]f64{ -1052.797, 4758.686, -7193.992, 3651.957 });
        const g422 = eccoG(el.ecco, &[_]f64{ -146.407, 841.8800, -1629.014, 1083.435 }, &[_]f64{ -3581.690, 16178.110, -24462.770, 12422.520 });
        const g520 = if (el.ecco <= 0.65)
            polyEval(el.ecco, &[_]f64{ -532.114, 3017.977, -5740.032, 3708.276 })
        else if (el.ecco > 0.715)
            polyEval(el.ecco, &[_]f64{ -5149.66, 29936.92, -54087.36, 31324.56 })
        else
            1464.74 - 4664.75 * el.ecco + 3763.64 * el.ecco * el.ecco;
        const g521 = eccoG7(el.ecco, &[_]f64{ -822.71072, 4568.6173, -8491.4146, 5337.524 }, &[_]f64{ -51752.104, 218913.95, -309468.16, 146349.42 });
        const g532 = eccoG7(el.ecco, &[_]f64{ -853.66600, 4690.2500, -8624.7700, 5341.400 }, &[_]f64{ -40023.880, 170470.89, -242699.48, 115605.82 });
        const g533 = eccoG7(el.ecco, &[_]f64{ -919.22770, 4988.6100, -9064.7700, 5542.21 }, &[_]f64{ -37995.780, 161616.52, -229838.20, 109377.94 });

        const f220 = 0.75 * (1.0 + 2.0 * trig.cosio + cosisq);
        const f221 = 1.5 * sini2;
        const f321 = 1.875 * trig.sinio * (1.0 - 2.0 * trig.cosio - 3.0 * cosisq);
        const f322 = -1.875 * trig.sinio * (1.0 + 2.0 * trig.cosio - 3.0 * cosisq);
        const f441 = 35.0 * sini2 * f220;
        const f442 = 39.3750 * sini2 * sini2;
        const f522 = 9.84375 * trig.sinio * (sini2 * (1.0 - 2.0 * trig.cosio - 5.0 * cosisq) +
            0.33333333 * (-2.0 + 4.0 * trig.cosio + 6.0 * cosisq));
        const f523 = trig.sinio * (4.92187512 * sini2 * (-2.0 - 4.0 * trig.cosio + 10.0 * cosisq) +
            6.56250012 * (1.0 + 2.0 * trig.cosio - 3.0 * cosisq));
        const f542 = 29.53125 * trig.sinio * (2.0 - 8.0 * trig.cosio + cosisq * (-12.0 + 8.0 * trig.cosio + 10.0 * cosisq));
        const f543 = 29.53125 * trig.sinio * (-2.0 - 8.0 * trig.cosio + cosisq * (12.0 + 8.0 * trig.cosio - 10.0 * cosisq));

        // Vallado cascading: temp1 *= aonv between each group
        const aonv = 1.0 / rec.a;
        var temp1 = 3.0 * dc.nm * dc.nm * aonv * aonv;
        var temp = temp1 * root22;
        di.d2201 = temp * f220 * g201;
        di.d2211 = temp * f221 * g211;
        temp1 = temp1 * aonv;
        temp = temp1 * root32;
        di.d3210 = temp * f321 * g310;
        di.d3222 = temp * f322 * g322;
        temp1 = temp1 * aonv;
        temp = 2.0 * temp1 * root44;
        di.d4410 = temp * f441 * g410;
        di.d4422 = temp * f442 * g422;
        temp1 = temp1 * aonv;
        temp = temp1 * root52;
        di.d5220 = temp * f522 * g520;
        di.d5232 = temp * f523 * g532;
        temp = 2.0 * temp1 * root54;
        di.d5421 = temp * f542 * g521;
        di.d5433 = temp * f543 * g533;

        di.xlamo = @mod(el.mo + el.nodeo + el.nodeo - gsto - gsto, constants.twoPi);
        di.xfact = secRates.mdot + di.dmdt + 2.0 * (secRates.nodedot + di.dnodt - rptim) - rec.noUnkozai;
    }

    return di;
}

/// Evaluate eccentricity-dependent G function with threshold at e=0.65
fn eccoG(ecco: f64, lo: []const f64, hi: []const f64) f64 {
    if (ecco <= 0.65) return polyEval(ecco, lo);
    return polyEval(ecco, hi);
}

/// Evaluate eccentricity-dependent G function with threshold at e=0.7
fn eccoG7(ecco: f64, lo: []const f64, hi: []const f64) f64 {
    if (ecco < 0.7) return polyEval(ecco, lo);
    return polyEval(ecco, hi);
}

fn polyEval(x: f64, c: []const f64) f64 {
    var result: f64 = 0;
    var xn: f64 = 1.0;
    for (c) |ci| {
        result += ci * xn;
        xn *= x;
    }
    return result;
}

// ── dpper ───────────────────────────────────────────────────────────────
fn dpper(
    el: *const Elements,
    tsince: f64,
    ep: *f64,
    inclp: *f64,
    nodep: *f64,
    argpp: *f64,
    mp: *f64,
) void {
    // Solar terms
    var zm = el.zmos + zns * tsince;
    var zf = zm + 2.0 * zes * @sin(zm);
    var sinzf = @sin(zf);
    var f2 = 0.5 * sinzf * sinzf - 0.25;
    var f3 = -0.5 * sinzf * @cos(zf);

    const ses = el.se2 * f2 + el.se3 * f3;
    const sis = el.si2 * f2 + el.si3 * f3;
    const sls = el.sl2 * f2 + el.sl3 * f3 + el.sl4 * sinzf;
    const sghs = el.sgh2 * f2 + el.sgh3 * f3 + el.sgh4 * sinzf;
    const shs = el.sh2 * f2 + el.sh3 * f3;

    // Lunar terms
    zm = el.zmol + znl * tsince;
    zf = zm + 2.0 * zel * @sin(zm);
    sinzf = @sin(zf);
    f2 = 0.5 * sinzf * sinzf - 0.25;
    f3 = -0.5 * sinzf * @cos(zf);

    const sel = el.ee2 * f2 + el.e3 * f3;
    const sil = el.xi2 * f2 + el.xi3 * f3;
    const sll = el.xl2 * f2 + el.xl3 * f3 + el.xl4 * sinzf;
    const sghl = el.xgh2 * f2 + el.xgh3 * f3 + el.xgh4 * sinzf;
    const shl = el.xh2 * f2 + el.xh3 * f3;

    // Combined
    const pe = ses + sel;
    const pinc = sis + sil;
    const pl = sls + sll;
    var pgh = sghs + sghl;
    var ph = shs + shl;

    inclp.* += pinc;
    ep.* += pe;

    const sinip = @sin(inclp.*);
    const cosip = @cos(inclp.*);

    if (inclp.* >= 0.2) {
        ph /= sinip;
        pgh -= cosip * ph;
        argpp.* += pgh;
        nodep.* += ph;
        mp.* += pl;
    } else {
        // Lyddane modification for near-equatorial orbits
        const sinop = @sin(nodep.*);
        const cosop = @cos(nodep.*);
        var alfdp = sinip * sinop;
        var betdp = sinip * cosop;
        const dalf = ph * cosop + pinc * cosip * sinop;
        const dbet = -ph * sinop + pinc * cosip * cosop;
        alfdp += dalf;
        betdp += dbet;
        nodep.* = @mod(nodep.*, constants.twoPi);
        const xls = mp.* + argpp.* + cosip * nodep.*;
        const dls = pl + pgh - pinc * nodep.* * sinip;
        const xnoh = nodep.*;
        nodep.* = std.math.atan2(alfdp, betdp);
        if (@abs(xnoh - nodep.*) > std.math.pi) {
            if (nodep.* < xnoh)
                nodep.* += constants.twoPi
            else
                nodep.* -= constants.twoPi;
        }
        mp.* += pl;
        argpp.* = xls + dls - mp.* - cosip * nodep.*;
    }
}

// ── dspace ──────────────────────────────────────────────────────────────
const DspaceState = struct {
    em: f64,
    argpm: f64,
    inclm: f64,
    mm: f64,
    nodem: f64,
    nm: f64,
    // Mutable integration state
    xli: f64,
    xni: f64,
    atime: f64,
};

fn dspace(el: *const Elements, tsince: f64, s: *DspaceState) void {

    // Apply lunar-solar secular perturbation rates
    s.em += el.dedt * tsince;
    s.inclm += el.didt * tsince;
    s.argpm += el.domdt * tsince;
    s.nodem += el.dnodt * tsince;
    s.mm += el.dmdt * tsince;

    // Resonance effects
    if (el.irez == 0) return;

    // Restart integration if needed
    if (s.atime == 0.0 or tsince * s.atime <= 0.0 or @abs(tsince) < @abs(s.atime)) {
        s.atime = 0.0;
        s.xni = el.sgp4.noUnkozai;
        s.xli = el.xlamo;
    }

    const delt: f64 = if (tsince > 0.0) stepp else stepn;

    // Integration loop
    while (@abs(tsince - s.atime) >= stepp) {
        var xndt: f64 = undefined;
        var xnddt: f64 = undefined;
        var xldot: f64 = undefined;
        computeResonanceAccel(el, s.xli, s.xni, s.atime, &xndt, &xnddt, &xldot);
        s.xli += xldot * delt + xndt * step2;
        s.xni += xndt * delt + xnddt * step2;
        s.atime += delt;
    }

    // Final partial step
    const ft = tsince - s.atime;
    var xndt: f64 = undefined;
    var xnddt: f64 = undefined;
    var xldot: f64 = undefined;
    computeResonanceAccel(el, s.xli, s.xni, s.atime, &xndt, &xnddt, &xldot);

    s.nm = s.xni + xndt * ft + xnddt * ft * ft * 0.5;
    const xl = s.xli + xldot * ft + xndt * ft * ft * 0.5;
    const theta = @mod(el.gsto + tsince * rptim, constants.twoPi);

    if (el.irez != 2) {
        // GEO
        s.mm = xl - s.nodem - s.argpm + theta;
    } else {
        // Half-day
        s.mm = xl - 2.0 * s.nodem + 2.0 * theta;
    }
    const dndt = s.nm - el.sgp4.noUnkozai;
    s.nm = el.sgp4.noUnkozai + dndt;
}

fn computeResonanceAccel(el: *const Elements, xli: f64, xni: f64, atime: f64, xndt: *f64, xnddt: *f64, xldot: *f64) void {
    if (el.irez == 2) {
        // Half-day resonance
        const xomi = el.sgp4.argpo + el.sgp4.argpdot * atime;
        const x2omi = xomi + xomi;
        const x2li = xli + xli;

        xndt.* = el.d2201 * @sin(x2omi + xli - g22) +
            el.d2211 * @sin(xli - g22) +
            el.d3210 * @sin(xomi + xli - g32) +
            el.d3222 * @sin(-xomi + xli - g32) +
            el.d4410 * @sin(x2omi + x2li - g44) +
            el.d4422 * @sin(x2li - g44) +
            el.d5220 * @sin(xomi + xli - g52) +
            el.d5232 * @sin(-xomi + xli - g52) +
            el.d5421 * @sin(xomi + x2li - g54) +
            el.d5433 * @sin(-xomi + x2li - g54);

        xldot.* = xni + el.xfact;

        xnddt.* = el.d2201 * @cos(x2omi + xli - g22) +
            el.d2211 * @cos(xli - g22) +
            el.d3210 * @cos(xomi + xli - g32) +
            el.d3222 * @cos(-xomi + xli - g32) +
            el.d5220 * @cos(xomi + xli - g52) +
            el.d5232 * @cos(-xomi + xli - g52) +
            2.0 * (el.d4410 * @cos(x2omi + x2li - g44) +
                el.d4422 * @cos(x2li - g44) +
                el.d5421 * @cos(xomi + x2li - g54) +
                el.d5433 * @cos(-xomi + x2li - g54));
        xnddt.* *= xldot.*;
    } else {
        // GEO synchronous resonance — three different phase constants
        xndt.* = el.del1 * @sin(xli - fasx2) +
            el.del2 * @sin(2.0 * (xli - fasx4)) +
            el.del3 * @sin(3.0 * (xli - fasx6));
        xldot.* = xni + el.xfact;
        xnddt.* = el.del1 * @cos(xli - fasx2) +
            2.0 * el.del2 * @cos(2.0 * (xli - fasx4)) +
            3.0 * el.del3 * @cos(3.0 * (xli - fasx6));
        xnddt.* *= xldot.*;
    }
}

// ── Propagation Pipeline ────────────────────────────────────────────────
fn propagateElements(el: *const Elements, tsince: f64) Error![2][3]f64 {
    const sgp4El = &el.sgp4;
    const t2 = tsince * tsince;

    // Step 1: Simplified secular update (isimp=true always for deep space)
    const tempa = 1.0 - sgp4El.cc1 * tsince;
    const tempe = sgp4El.bstar * sgp4El.cc4 * tsince;
    const templ = sgp4El.t2cof * t2;

    const xmdf = sgp4El.mo + sgp4El.mdot * tsince;
    const argpdf = sgp4El.argpo + sgp4El.argpdot * tsince;
    const nodedf = sgp4El.nodeo + sgp4El.nodedot * tsince;
    const nodem_init = nodedf + sgp4El.xnodcf * t2;

    // Step 2: Deep space secular + resonance
    var ds = DspaceState{
        .em = sgp4El.ecco,
        .argpm = argpdf,
        .inclm = sgp4El.inclo,
        .mm = xmdf,
        .nodem = nodem_init,
        .nm = sgp4El.noUnkozai,
        .xli = el.xlamo,
        .xni = sgp4El.noUnkozai,
        .atime = 0.0,
    };
    dspace(el, tsince, &ds);

    // Step 3: Compute semi-major axis and mean motion
    var nm = ds.nm;
    if (nm <= 0.0) return Error.SatelliteDecayed;
    const am = std.math.pow(f64, sgp4El.grav.xke / nm, 2.0 / 3.0) * tempa * tempa;
    nm = sgp4El.grav.xke / std.math.pow(f64, am, 1.5);
    var em = ds.em - tempe;

    if (em >= 1.0 or em < -0.001) return Error.InvalidEccentricity;
    if (em < 1.0e-6) em = 1.0e-6;
    if (am < 0.95) return Error.SatelliteDecayed;

    var mm = ds.mm + sgp4El.noUnkozai * templ;
    const xlm = mm + ds.argpm + ds.nodem;
    var nodem = @mod(ds.nodem, constants.twoPi);
    var argpm = @mod(ds.argpm, constants.twoPi);
    mm = @mod(xlm - argpm - nodem, constants.twoPi);
    var inclm = ds.inclm;

    // Step 4: Deep space periodic perturbations
    dpper(el, tsince, &em, &inclm, &nodem, &argpm, &mm);

    if (inclm < 0.0) {
        inclm = -inclm;
        nodem += std.math.pi;
        argpm -= std.math.pi;
    }
    if (em < 1.0e-6) em = 1.0e-6;
    if (em >= 1.0) return Error.InvalidEccentricity;

    // Step 5: Recompute inclination-dependent terms for deep space
    const sinip = @sin(inclm);
    const cosip = @cos(inclm);
    const cosip2 = cosip * cosip;

    var elCopy = sgp4El.*;
    elCopy.aycof = -0.5 * sgp4El.grav.j3oj2 * sinip;
    if (@abs(cosip + 1.0) > 1.5e-12)
        elCopy.xlcof = -0.25 * sgp4El.grav.j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip)
    else
        elCopy.xlcof = -0.25 * sgp4El.grav.j3oj2 * sinip * (3.0 + 5.0 * cosip) / 1.5e-12;
    elCopy.inclo = inclm;
    elCopy.cosio = cosip;
    elCopy.sinio = sinip;
    elCopy.cosio2 = cosip2;
    elCopy.x1mth2 = 1.0 - cosip2;
    elCopy.con41 = 3.0 * cosip2 - 1.0;
    elCopy.x7thm1 = 7.0 * cosip2 - 1.0;

    // Step 6: Kepler solver, short-period corrections, position/velocity (reuse SGP4)
    const secular = Sgp4.SecularState{
        .mm = mm,
        .argpm = argpm,
        .nodem = nodem,
        .em = em,
        .a = am,
    };

    const kepler = Sgp4.solveKepler(&elCopy, secular);
    const corrected = Sgp4.applyShortPeriodCorrections(&elCopy, kepler, nm);

    if (corrected.r < 1.0) return Error.SatelliteDecayed;

    return Sgp4.computePositionVelocity(&elCopy, corrected);
}

// ── Tests ───────────────────────────────────────────────────────────────
const testing = std.testing;

test "gstime J2000" {
    const gst = gstime(2451545.0);
    // GMST at J2000.0 ≈ 280.46° ≈ 4.8950 rad
    try testing.expectApproxEqAbs(@as(f64, 4.8949612), gst, 0.001);
}

test "sdp4 init coefficients GPS 20413" {
    const tle_str =
        \\1 20413U 90005A   24186.00000000  .00000012  00000+0  10000-3 0  9992
        \\2 20413  55.4408  61.4858 0112981 129.5765 231.5553  2.00561730104446
    ;
    var tle = try Tle.parse(tle_str, testing.allocator);
    defer tle.deinit();
    const sdp4 = try Sdp4.init(tle, constants.wgs72);
    const el = &sdp4.elements;

    // Reference from python-sgp4
    try testing.expectApproxEqAbs(@as(f64, 0.0087511706), el.sgp4.noUnkozai, 1e-8);
    try testing.expectApproxEqAbs(@as(f64, 4.1643039044), el.sgp4.a, 1e-6);
    try testing.expect(el.irez == 0); // non-resonant (e < 0.5)
    try testing.expectApproxEqAbs(@as(f64, 4.9305096469), el.gsto, 1e-6);
    try testing.expectApproxEqAbs(@as(f64, 3.1269253661), el.zmos, 1e-6);
    try testing.expectApproxEqAbs(@as(f64, 1.2769002197), el.zmol, 1e-6);

    // Solar coefficients
    try testing.expectApproxEqAbs(@as(f64, 7.4611141471e-05), el.se2, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, -2.6550152994e-05), el.se3, 1e-12);

    // Secular rates
    try testing.expectApproxEqAbs(@as(f64, -1.3083083111e-10), el.dedt, 1e-17);
    try testing.expectApproxEqAbs(@as(f64, -1.9461479899e-08), el.dnodt, 1e-15);
}

test "sdp4 init coefficients GEO 28626" {
    const tle_str =
        \\1 28626U 05004A   24186.00000000 -.00000098  00000+0  00000+0 0  9998
        \\2 28626   0.0163 279.8379 0003069  20.3251 343.1766  1.00270142 70992
    ;
    var tle = try Tle.parse(tle_str, testing.allocator);
    defer tle.deinit();
    const sdp4 = try Sdp4.init(tle, constants.wgs72);
    const el = &sdp4.elements;

    try testing.expect(el.irez == 1); // GEO resonance
    try testing.expectApproxEqAbs(@as(f64, 0.0043749477), el.sgp4.noUnkozai, 1e-8);
    try testing.expectApproxEqAbs(@as(f64, -6.3971905151e-13), el.del1, 1e-20);
    try testing.expectApproxEqAbs(@as(f64, 1.4103492854e-11), el.del2, 1e-18);
    try testing.expectApproxEqAbs(@as(f64, 1.9783801291e-12), el.del3, 1e-19);
}

test "sdp4 init coefficients HEO 09880" {
    const tle_str =
        \\1 09880U 77021B   24186.00000000  .00000023  00000+0  00000+0 0  9999
        \\2 09880  63.4300  75.8891 7318036 269.8735  16.7549  2.00611684 54321
    ;
    var tle = try Tle.parse(tle_str, testing.allocator);
    defer tle.deinit();
    const sdp4 = try Sdp4.init(tle, constants.wgs72);
    const el = &sdp4.elements;

    try testing.expect(el.irez == 2); // half-day resonance
    try testing.expectApproxEqAbs(@as(f64, 0.0087538538), el.sgp4.noUnkozai, 1e-8);
    try testing.expectApproxEqAbs(@as(f64, -1.2912381830e-11), el.d2201, 1e-18);
    try testing.expectApproxEqAbs(@as(f64, 8.0024701733e-11), el.d2211, 1e-18);
}

test "sdp4 propagate GPS 20413" {
    const tle_str =
        \\1 20413U 90005A   24186.00000000  .00000012  00000+0  10000-3 0  9992
        \\2 20413  55.4408  61.4858 0112981 129.5765 231.5553  2.00561730104446
    ;
    var tle = try Tle.parse(tle_str, testing.allocator);
    defer tle.deinit();
    const sdp4 = try Sdp4.init(tle, constants.wgs72);

    // t=0 reference: R=[12743.39105131, 23518.44806062, 28.23195962]
    //                V=[-1.93584727, 1.00575086, 3.16821707]
    const r0 = try sdp4.propagate(0.0);
    try testing.expectApproxEqAbs(@as(f64, 12743.39105131), r0[0][0], 0.01);
    try testing.expectApproxEqAbs(@as(f64, 23518.44806062), r0[0][1], 0.01);
    try testing.expectApproxEqAbs(@as(f64, 28.23195962), r0[0][2], 0.01);
    try testing.expectApproxEqAbs(@as(f64, -1.93584727), r0[1][0], 1e-5);
    try testing.expectApproxEqAbs(@as(f64, 1.00575086), r0[1][1], 1e-5);
    try testing.expectApproxEqAbs(@as(f64, 3.16821707), r0[1][2], 1e-5);

    // t=720 reference: R=[12513.21380970, 23633.30126287, 414.75154681]
    //                  V=[-1.96762329, 0.94648802, 3.16764596]
    const r720 = try sdp4.propagate(720.0);
    try testing.expectApproxEqAbs(@as(f64, 12513.21380970), r720[0][0], 0.01);
    try testing.expectApproxEqAbs(@as(f64, 23633.30126287), r720[0][1], 0.01);
    try testing.expectApproxEqAbs(@as(f64, 414.75154681), r720[0][2], 0.01);

    // t=1440 reference: R=[12279.27857794, 23740.96214221, 801.15912584]
    //                   V=[-1.99885060, 0.88692941, 3.16608577]
    const r1440 = try sdp4.propagate(1440.0);
    try testing.expectApproxEqAbs(@as(f64, 12279.27857794), r1440[0][0], 0.01);
    try testing.expectApproxEqAbs(@as(f64, 23740.96214221), r1440[0][1], 0.01);
    try testing.expectApproxEqAbs(@as(f64, 801.15912584), r1440[0][2], 0.01);
}

test "sdp4 propagate GEO 28626" {
    const tle_str =
        \\1 28626U 05004A   24186.00000000 -.00000098  00000+0  00000+0 0  9998
        \\2 28626   0.0163 279.8379 0003069  20.3251 343.1766  1.00270142 70992
    ;
    var tle = try Tle.parse(tle_str, testing.allocator);
    defer tle.deinit();
    const sdp4 = try Sdp4.init(tle, constants.wgs72);

    // t=0 reference: R=[9727.65902952, -41014.43193269, -9.39974955]
    //               V=[2.99252607, 0.71003488, 0.00039212]
    const r0 = try sdp4.propagate(0.0);
    try testing.expectApproxEqAbs(@as(f64, 9727.65902952), r0[0][0], 0.01);
    try testing.expectApproxEqAbs(@as(f64, -41014.43193269), r0[0][1], 0.01);
    try testing.expectApproxEqAbs(@as(f64, -9.39974955), r0[0][2], 0.01);

    // t=1440 reference: R=[10430.11028794, -40841.32649665, -6.90212577]
    //                   V=[2.97990947, 0.76127679, 0.00014975]
    const r1440 = try sdp4.propagate(1440.0);
    try testing.expectApproxEqAbs(@as(f64, 10430.11028794), r1440[0][0], 0.01);
    try testing.expectApproxEqAbs(@as(f64, -40841.32649665), r1440[0][1], 0.01);
}

test "sdp4 propagate HEO 09880" {
    const tle_str =
        \\1 09880U 77021B   24186.00000000  .00000023  00000+0  00000+0 0  9999
        \\2 09880  63.4300  75.8891 7318036 269.8735  16.7549  2.00611684 54321
    ;
    var tle = try Tle.parse(tle_str, testing.allocator);
    defer tle.deinit();
    const sdp4 = try Sdp4.init(tle, constants.wgs72);

    // t=0 reference: R=[2575.93297901, 13237.27540045, 1419.07364030]
    //               V=[-1.59368245, 3.99224937, 5.03078896]
    const r0 = try sdp4.propagate(0.0);
    try testing.expectApproxEqAbs(@as(f64, 2575.93297901), r0[0][0], 0.01);
    try testing.expectApproxEqAbs(@as(f64, 13237.27540045), r0[0][1], 0.01);
    try testing.expectApproxEqAbs(@as(f64, 1419.07364030), r0[0][2], 0.01);

    // t=1440 reference: R=[2175.00194988, 14214.96552375, 2741.44350110]
    const r1440 = try sdp4.propagate(1440.0);
    try testing.expectApproxEqAbs(@as(f64, 2175.00194988), r1440[0][0], 0.01);
    try testing.expectApproxEqAbs(@as(f64, 14214.96552375), r1440[0][1], 0.01);
    try testing.expectApproxEqAbs(@as(f64, 2741.44350110), r1440[0][2], 0.01);
}

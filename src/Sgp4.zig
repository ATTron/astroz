//! SGP4 Orbit Propagator
//! Reference: Revisiting Spacetrack Report #3

const std = @import("std");
const constants = @import("constants.zig");
const calculations = @import("calculations.zig");
const simdMath = @import("simdMath.zig");
const Datetime = @import("Datetime.zig");
const Tle = @import("Tle.zig");

const Sgp4 = @This();

// i am speed . . .
const Vec4 = @Vector(4, f64);

// Atmospheric drag model parameters (km)
const perigeeSAdjustment = 156.0; // adjust s parameter if below this value
const perigeeSMinimum = 98.0; // use minimum s, if below this value
const sParameterStandard = 78.0; // atmospheric parameter
const qParameter = 120.0;
const sParameterMinimum = 20.0; // s value minimum
const perigeeSimplified = 220.0; // use simplified drag model, if below this value

// Numerical tolerances
const eccentricityMin = 1.0e-4; // minimum for cc3/xmcof
const eccentricityFloor = 1.0e-6; // floor to prevent numerical issues
const singularityTolerance = 1.5e-12; // near-polar orbit singularity check

pub const Error = error{
    DeepSpaceNotSupported,
    InvalidEccentricity,
    SatelliteDecayed,
};

pub const Elements = struct {
    grav: constants.Sgp4GravityModel,
    epochJd: f64,

    // Mean elements from TLE
    noKozai: f64,
    ecco: f64,
    inclo: f64,
    nodeo: f64,
    argpo: f64,
    mo: f64,
    bstar: f64,

    // Derived mean elements
    noUnkozai: f64,
    a: f64,

    // Trigonometric terms
    sinio: f64,
    cosio: f64,
    cosio2: f64,
    cosio4: f64,

    // Polynomial terms in cos^2(i)
    con41: f64, // 3cos^2(i) - 1
    con42: f64, // 1 - 5cos^2(i)
    x1mth2: f64, // sin^2(i)
    x7thm1: f64, // 7cos^2(i) - 1

    // Secular rate coefficients
    mdot: f64,
    argpdot: f64,
    nodedot: f64,

    // Drag coefficients
    cc1: f64,
    cc4: f64,
    cc5: f64,
    t2cof: f64,
    omgcof: f64,
    xnodcf: f64,
    xlcof: f64,
    xmcof: f64,
    aycof: f64,
    eta: f64,
    delmo: f64,
    sinmao: f64,

    // Higher-order drag terms (only used when isimp=false)
    d2: f64,
    d3: f64,
    d4: f64,
    t3cof: f64,
    t4cof: f64,
    t5cof: f64,

    // Precomputed for propagation
    aBase: f64, // cbrt((xke/noUnkozai)^2)
    vkmpersec: f64, // velocity conversion factor

    isimp: bool, // simplified drag model flag
};

elements: Elements,

/// initialize from tle
pub fn init(tle: Tle, grav: constants.Sgp4GravityModel) Error!Sgp4 {
    return .{ .elements = try initElements(tle, grav) };
}

/// propagate to time tsince (minutes from TLE epoch)
pub fn propagate(self: *const Sgp4, tsince: f64) Error![2][3]f64 {
    return propagateElements(&self.elements, tsince);
}

fn initElements(tle: Tle, grav: constants.Sgp4GravityModel) Error!Elements {
    const meanElements = extractMeanElements(tle);

    if (meanElements.ecco < 0.0 or meanElements.ecco >= 1.0) {
        return Error.InvalidEccentricity;
    }

    const recovered = recoverMeanMotion(meanElements, grav);

    const rp = recovered.a * (1.0 - meanElements.ecco);
    if (rp < 1.0) return Error.SatelliteDecayed;

    const period = constants.twoPi / recovered.noUnkozai;
    if (period > constants.sgp4DeepSpaceThresholdMinutes) {
        return Error.DeepSpaceNotSupported;
    }

    const trig = computeTrigTerms(meanElements.inclo);
    const poly = computePolyTerms(trig);
    const secular = computeSecularRates(meanElements, recovered, trig, poly, grav);
    const drag = computeDragCoefficients(meanElements, recovered, trig, poly, grav);
    const higherOrder = computeHigherOrderDrag(recovered, drag, grav);

    const fullYear: u16 = if (tle.firstLine.epochYear < 57)
        2000 + tle.firstLine.epochYear
    else
        1900 + tle.firstLine.epochYear;
    const epochJd = Datetime.yearDoyToJulianDate(fullYear, tle.firstLine.epochDay);

    return Elements{
        .grav = grav,
        .epochJd = epochJd,
        .noKozai = meanElements.noKozai,
        .ecco = meanElements.ecco,
        .inclo = meanElements.inclo,
        .nodeo = meanElements.nodeo,
        .argpo = meanElements.argpo,
        .mo = meanElements.mo,
        .bstar = meanElements.bstar,
        .noUnkozai = recovered.noUnkozai,
        .a = recovered.a,
        .sinio = trig.sinio,
        .cosio = trig.cosio,
        .cosio2 = trig.cosio2,
        .cosio4 = trig.cosio4,
        .con41 = poly.con41,
        .con42 = poly.con42,
        .x1mth2 = poly.x1mth2,
        .x7thm1 = poly.x7thm1,
        .mdot = secular.mdot,
        .argpdot = secular.argpdot,
        .nodedot = secular.nodedot,
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
        .d2 = higherOrder.d2,
        .d3 = higherOrder.d3,
        .d4 = higherOrder.d4,
        .t3cof = higherOrder.t3cof,
        .t4cof = higherOrder.t4cof,
        .t5cof = higherOrder.t5cof,
        .aBase = blk: {
            const ratio = grav.xke / recovered.noUnkozai;
            break :blk std.math.cbrt(ratio * ratio);
        },
        .vkmpersec = grav.xke * grav.radiusEarthKm / 60.0,
        .isimp = higherOrder.isimp,
    };
}

const MeanElements = struct {
    noKozai: f64,
    ecco: f64,
    inclo: f64,
    nodeo: f64,
    argpo: f64,
    mo: f64,
    bstar: f64,
};

fn extractMeanElements(tle: Tle) MeanElements {
    return .{
        .noKozai = tle.secondLine.mMotion * constants.twoPi / constants.minutes_per_day,
        .ecco = tle.secondLine.eccentricity,
        .inclo = tle.secondLine.inclination * constants.deg2rad,
        .nodeo = tle.secondLine.rightAscension * constants.deg2rad,
        .argpo = tle.secondLine.perigee * constants.deg2rad,
        .mo = tle.secondLine.mAnomaly * constants.deg2rad,
        .bstar = tle.firstLine.bstarDrag,
    };
}

const RecoveredElements = struct { noUnkozai: f64, a: f64 };

fn recoverMeanMotion(el: MeanElements, grav: constants.Sgp4GravityModel) RecoveredElements {
    // unkozai the mean motion
    const cosio = @cos(el.inclo);
    const theta2 = cosio * cosio;
    const x3thm1 = 3.0 * theta2 - 1.0;
    const eosq = el.ecco * el.ecco;
    const betao2 = 1.0 - eosq;
    const betao = @sqrt(betao2);

    const a1 = std.math.pow(f64, grav.xke / el.noKozai, 2.0 / 3.0);

    // Kozai correction: del = 0.75 * j2 * x3thm1 / (a^2 * beta^3)
    const del1 = 0.75 * grav.j2 * x3thm1 / (a1 * a1 * betao * betao2);

    // Second-order refinement: ao = a1 * (1 - del*(1/3 + del*(1 + 134/81*del)))
    const ao = a1 * (1.0 - del1 * (1.0 / 3.0 + del1 * (1.0 + 134.0 / 81.0 * del1)));
    const delo = 0.75 * grav.j2 * x3thm1 / (ao * ao * betao * betao2);

    const noUnkozai = el.noKozai / (1.0 + delo);
    const a = std.math.pow(f64, grav.xke / noUnkozai, 2.0 / 3.0);

    return .{ .noUnkozai = noUnkozai, .a = a };
}

const TrigTerms = struct { sinio: f64, cosio: f64, cosio2: f64, cosio4: f64 };

fn computeTrigTerms(inclo: f64) TrigTerms {
    const sinio = @sin(inclo);
    const cosio = @cos(inclo);
    const cosio2 = cosio * cosio;
    return .{ .sinio = sinio, .cosio = cosio, .cosio2 = cosio2, .cosio4 = cosio2 * cosio2 };
}

const PolyTerms = struct { con41: f64, con42: f64, x1mth2: f64, x7thm1: f64 };

fn computePolyTerms(trig: TrigTerms) PolyTerms {
    const theta2 = trig.cosio2;
    return .{
        .con41 = 3.0 * theta2 - 1.0,
        .con42 = 1.0 - 5.0 * theta2,
        .x1mth2 = 1.0 - theta2,
        .x7thm1 = 7.0 * theta2 - 1.0,
    };
}

const SecularRates = struct { mdot: f64, argpdot: f64, nodedot: f64 };

fn computeSecularRates(
    el: MeanElements,
    rec: RecoveredElements,
    trig: TrigTerms,
    poly: PolyTerms,
    grav: constants.Sgp4GravityModel,
) SecularRates {
    const omeosq = 1.0 - el.ecco * el.ecco;
    const rteosq = @sqrt(omeosq);
    const pinvsq = 1.0 / std.math.pow(f64, rec.a * omeosq, 2.0);

    // J2/J4 perturbation terms
    const temp1 = 1.5 * grav.j2 * pinvsq * rec.noUnkozai;
    const temp2 = 0.5 * temp1 * grav.j2 * pinvsq;
    const temp3 = -0.46875 * grav.j4 * pinvsq * pinvsq * rec.noUnkozai;

    // mean motion rate (rad/min)
    const mdot = rec.noUnkozai + 0.5 * temp1 * rteosq * poly.con41 +
        0.0625 * temp2 * rteosq * (13.0 - 78.0 * trig.cosio2 + 137.0 * trig.cosio4);

    // argument of perigee rate
    const argpdot = -0.5 * temp1 * poly.con42 +
        0.0625 * temp2 * (7.0 - 114.0 * trig.cosio2 + 395.0 * trig.cosio4) +
        temp3 * (3.0 - 36.0 * trig.cosio2 + 49.0 * trig.cosio4);

    // RAAN rate
    const xhdot1 = -temp1 * trig.cosio;
    const nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * trig.cosio2) +
        2.0 * temp3 * (3.0 - 7.0 * trig.cosio2)) * trig.cosio;

    return .{ .mdot = mdot, .argpdot = argpdot, .nodedot = nodedot };
}

const DragCoefficients = struct {
    cc1: f64,
    cc4: f64,
    cc5: f64,
    t2cof: f64,
    omgcof: f64,
    xnodcf: f64,
    xlcof: f64,
    xmcof: f64,
    aycof: f64,
    eta: f64,
    delmo: f64,
    sinmao: f64,
    perige: f64,
};

fn computeDragCoefficients(
    el: MeanElements,
    rec: RecoveredElements,
    trig: TrigTerms,
    poly: PolyTerms,
    grav: constants.Sgp4GravityModel,
) DragCoefficients {
    const omeosq = 1.0 - el.ecco * el.ecco;
    const rp = rec.a * (1.0 - el.ecco);
    const perige = (rp - 1.0) * grav.radiusEarthKm;

    // atmospheric drag parameters adjusted for low perigee
    const sfour, const qzms24 = if (perige < perigeeSAdjustment) blk: {
        const s = if (perige < perigeeSMinimum) sParameterMinimum else perige - sParameterStandard;
        const qtemp = (qParameter - s) / grav.radiusEarthKm;
        break :blk .{ s / grav.radiusEarthKm + 1.0, qtemp * qtemp * qtemp * qtemp };
    } else blk: {
        const qtemp = (qParameter - sParameterStandard) / grav.radiusEarthKm;
        break :blk .{ sParameterStandard / grav.radiusEarthKm + 1.0, qtemp * qtemp * qtemp * qtemp };
    };

    const pinvsq = 1.0 / std.math.pow(f64, rec.a * omeosq, 2.0);
    const tsi = 1.0 / (rec.a - sfour);
    const eta = rec.a * el.ecco * tsi;
    const etasq = eta * eta;
    const eeta = el.ecco * eta;
    const psisq = @abs(1.0 - etasq);
    const coef = qzms24 * std.math.pow(f64, tsi, 4.0);
    const coef1 = coef / std.math.pow(f64, psisq, 3.5);

    // drag coefficients cc1-cc5
    const cc2 = coef1 * rec.noUnkozai * (rec.a * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) +
        0.375 * grav.j2 * tsi / psisq * poly.con41 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
    const cc1 = el.bstar * cc2;

    const cc3 = if (el.ecco > eccentricityMin)
        -2.0 * coef * tsi * grav.j3oj2 * rec.noUnkozai * trig.sinio / el.ecco
    else
        0.0;

    const cc4 = 2.0 * rec.noUnkozai * coef1 * rec.a * omeosq *
        (eta * (2.0 + 0.5 * etasq) + el.ecco * (0.5 + 2.0 * etasq) -
            grav.j2 * tsi / (rec.a * psisq) *
                (-3.0 * poly.con41 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) +
                    0.75 * poly.x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) * @cos(2.0 * el.argpo)));

    const cc5 = 2.0 * coef1 * rec.a * omeosq * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);

    // Other drag-related coefficients
    const temp1 = 1.5 * grav.j2 * pinvsq * rec.noUnkozai;
    const xhdot1 = -temp1 * trig.cosio;
    const xnodcf = 3.5 * omeosq * xhdot1 * cc1;
    const t2cof = 1.5 * cc1;

    const xlcof = if (@abs(trig.cosio + 1.0) > singularityTolerance)
        -0.25 * grav.j3oj2 * trig.sinio * (3.0 + 5.0 * trig.cosio) / (1.0 + trig.cosio)
    else
        -0.25 * grav.j3oj2 * trig.sinio * (3.0 + 5.0 * trig.cosio) / singularityTolerance;

    const aycof = -0.5 * grav.j3oj2 * trig.sinio;
    const delmotemp = 1.0 + eta * @cos(el.mo);
    const delmo = delmotemp * delmotemp * delmotemp;
    const sinmao = @sin(el.mo);

    const xmcof = if (el.ecco > eccentricityMin)
        -(2.0 / 3.0) * coef * el.bstar / eeta
    else
        0.0;

    const omgcof = el.bstar * cc3 * @cos(el.argpo);

    return .{
        .cc1 = cc1,
        .cc4 = cc4,
        .cc5 = cc5,
        .t2cof = t2cof,
        .omgcof = omgcof,
        .xnodcf = xnodcf,
        .xlcof = xlcof,
        .xmcof = xmcof,
        .aycof = aycof,
        .eta = eta,
        .delmo = delmo,
        .sinmao = sinmao,
        .perige = perige,
    };
}

const HigherOrderDrag = struct {
    d2: f64,
    d3: f64,
    d4: f64,
    t3cof: f64,
    t4cof: f64,
    t5cof: f64,
    isimp: bool,
};

fn computeHigherOrderDrag(
    rec: RecoveredElements,
    drag: DragCoefficients,
    grav: constants.Sgp4GravityModel,
) HigherOrderDrag {
    if (drag.perige < perigeeSimplified) {
        return .{ .d2 = 0, .d3 = 0, .d4 = 0, .t3cof = 0, .t4cof = 0, .t5cof = 0, .isimp = true };
    }

    const s = sParameterStandard / grav.radiusEarthKm + 1.0;
    const tsi = 1.0 / (rec.a - s);

    const cc1sq = drag.cc1 * drag.cc1;
    const d2 = 4.0 * rec.a * tsi * cc1sq;
    const temp = d2 * tsi * drag.cc1 / 3.0;
    const d3 = (17.0 * rec.a + s) * temp;
    const d4 = 0.5 * temp * rec.a * tsi * (221.0 * rec.a + 31.0 * s) * drag.cc1;
    const t3cof = d2 + 2.0 * cc1sq;
    const t4cof = 0.25 * (3.0 * d3 + drag.cc1 * (12.0 * d2 + 10.0 * cc1sq));
    const t5cof = 0.2 * (3.0 * d4 + 12.0 * drag.cc1 * d3 + 6.0 * d2 * d2 + 15.0 * cc1sq * (2.0 * d2 + cc1sq));

    return .{ .d2 = d2, .d3 = d3, .d4 = d4, .t3cof = t3cof, .t4cof = t4cof, .t5cof = t5cof, .isimp = false };
}

fn propagateElements(el: *const Elements, tsince: f64) Error![2][3]f64 {
    const secular = updateSecular(el, tsince);
    const nm = el.grav.xke / std.math.pow(f64, secular.a, 1.5);
    const kepler = solveKepler(el, secular);
    const corrected = applyShortPeriodCorrections(el, kepler, nm);
    return computePositionVelocity(el, corrected);
}

const SecularState = struct {
    mm: f64,
    argpm: f64,
    nodem: f64,
    em: f64,
    a: f64,
};

fn updateSecular(el: *const Elements, tsince: f64) SecularState {
    const t2 = tsince * tsince;

    var tempa = 1.0 - el.cc1 * tsince;
    var tempe = el.bstar * el.cc4 * tsince;
    var templ = el.t2cof * t2;

    const xmdf = el.mo + el.mdot * tsince;
    const argpdf = el.argpo + el.argpdot * tsince;
    const nodedf = el.nodeo + el.nodedot * tsince;
    var argpm = argpdf;
    var mm = xmdf;
    var nodem = nodedf + el.xnodcf * t2;

    if (!el.isimp) {
        // Higher-order drag corrections
        const delomg = el.omgcof * tsince;
        const delmtemp = 1.0 + el.eta * @cos(xmdf);
        const delm = el.xmcof * (delmtemp * delmtemp * delmtemp - el.delmo);
        const temp = delomg + delm;
        mm = xmdf + temp;
        argpm = argpdf - temp;

        const t3 = t2 * tsince;
        const t4 = t3 * tsince;
        tempa = tempa - el.d2 * t2 - el.d3 * t3 - el.d4 * t4;
        tempe = tempe + el.bstar * el.cc5 * (@sin(mm) - el.sinmao);
        templ = templ + el.t3cof * t3 + t4 * (el.t4cof + tsince * el.t5cof);
    }

    const am = el.aBase * tempa * tempa;
    var em = el.ecco - tempe;
    em = @max(em, eccentricityFloor);

    mm = mm + el.noUnkozai * templ;
    const xlm = mm + argpm + nodem;

    nodem = @mod(nodem, constants.twoPi);
    argpm = @mod(argpm, constants.twoPi);
    mm = @mod(xlm - argpm - nodem, constants.twoPi);

    return .{ .mm = mm, .argpm = argpm, .nodem = nodem, .em = em, .a = am };
}

const KeplerState = struct {
    u: f64,
    r: f64,
    rdot: f64,
    rvdot: f64,
    betal: f64,
    sin2u: f64,
    cos2u: f64,
    nodem: f64,
    pl: f64,
};

/// Solve Kepler's equation in equinoctial form and compute orbital state.
/// Note: This uses SGP4's equinoctial formulation (axnl, aynl) rather than
/// the standard E - e*sin(E) = M form. See calculations.solveKeplerEquation
/// for the standard formulation.
fn solveKepler(el: *const Elements, sec: SecularState) KeplerState {
    const temp = 1.0 / (sec.a * (1.0 - sec.em * sec.em));
    const axnl = sec.em * @cos(sec.argpm);
    const aynl = sec.em * @sin(sec.argpm) + temp * el.aycof;
    const xl = @mod(sec.mm + sec.argpm + sec.nodem + temp * el.xlcof * axnl, constants.twoPi);

    var u = @mod(xl - sec.nodem, constants.twoPi);
    var eo1 = u;
    var sineo1: f64 = 0.0;
    var coseo1: f64 = 1.0;

    var tem5: f64 = 9999.9;
    var ktr: u32 = 1;
    while (@abs(tem5) >= 1.0e-12 and ktr <= 10) {
        sineo1 = @sin(eo1);
        coseo1 = @cos(eo1);
        tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl;
        tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
        if (@abs(tem5) >= 0.95) {
            tem5 = if (tem5 > 0.0) 0.95 else -0.95;
        }
        eo1 = eo1 + tem5;
        ktr += 1;
    }

    const ecose = axnl * coseo1 + aynl * sineo1;
    const esine = axnl * sineo1 - aynl * coseo1;
    const el2 = axnl * axnl + aynl * aynl;
    const pl = sec.a * (1.0 - el2);
    const betal = @sqrt(1.0 - el2);
    const rl = sec.a * (1.0 - ecose);
    const rdotl = @sqrt(sec.a) * esine / rl;
    const rvdotl = @sqrt(pl) / rl;

    const aOverR = sec.a / rl;
    const esineTerm = esine / (1.0 + betal);
    const sinu = aOverR * (sineo1 - aynl - axnl * esineTerm);
    const cosu = aOverR * (coseo1 - axnl + aynl * esineTerm);
    u = std.math.atan2(sinu, cosu);

    return .{
        .u = u,
        .r = rl,
        .rdot = rdotl,
        .rvdot = rvdotl,
        .betal = betal,
        .sin2u = 2.0 * sinu * cosu,
        .cos2u = 1.0 - 2.0 * sinu * sinu,
        .nodem = sec.nodem,
        .pl = pl,
    };
}

const CorrectedState = struct {
    r: f64,
    rdot: f64,
    rvdot: f64,
    u: f64,
    xnode: f64,
    xinc: f64,
};

fn applyShortPeriodCorrections(el: *const Elements, kep: KeplerState, nm: f64) CorrectedState {
    const temp = 1.0 / kep.pl;
    const temp1 = 0.5 * el.grav.j2 * temp;
    const temp2 = temp1 * temp;

    // short period corrections to position, velocity, and orbital elements
    const mrt = kep.r * (1.0 - 1.5 * temp2 * kep.betal * el.con41) + 0.5 * temp1 * el.x1mth2 * kep.cos2u;
    const su = kep.u - 0.25 * temp2 * el.x7thm1 * kep.sin2u;
    const xnode = kep.nodem + 1.5 * temp2 * el.cosio * kep.sin2u;
    const xinc = el.inclo + 1.5 * temp2 * el.cosio * el.sinio * kep.cos2u;
    const mvt = kep.rdot - nm * temp1 * el.x1mth2 * kep.sin2u / el.grav.xke;
    const rvdot = kep.rvdot + nm * temp1 * (el.x1mth2 * kep.cos2u + 1.5 * el.con41) / el.grav.xke;

    return .{ .r = mrt, .rdot = mvt, .rvdot = rvdot, .u = su, .xnode = xnode, .xinc = xinc };
}

fn computePositionVelocity(el: *const Elements, state: CorrectedState) [2][3]f64 {
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

    const rScaled = state.r * el.grav.radiusEarthKm;
    const r: [3]f64 = .{
        rScaled * ux,
        rScaled * uy,
        rScaled * uz,
    };
    const v: [3]f64 = .{
        (state.rdot * ux + state.rvdot * vx) * el.vkmpersec,
        (state.rdot * uy + state.rvdot * vy) * el.vkmpersec,
        (state.rdot * uz + state.rvdot * vz) * el.vkmpersec,
    };

    return .{ r, v };
}

// SIMD FUNCTIONS
// keeping the scalar functions to compare as they are gonna be
// easier to debug if shit breaks

const SecularStateV4 = struct {
    mm: Vec4,
    argpm: Vec4,
    nodem: Vec4,
    em: Vec4,
    a: Vec4,
};

const KeplerStateV4 = struct {
    u: Vec4,
    r: Vec4,
    rdot: Vec4,
    rvdot: Vec4,
    betal: Vec4,
    sin2u: Vec4,
    cos2u: Vec4,
    nodem: Vec4,
    pl: Vec4,
};

const CorrectedStateV4 = struct {
    r: Vec4,
    rdot: Vec4,
    rvdot: Vec4,
    u: Vec4,
    xnode: Vec4,
    xinc: Vec4,
};

pub fn propagateV4(self: *const Sgp4, times: [4]f64) Error![4][2][3]f64 {
    const el = &self.elements;
    const timeV4 = Vec4{ times[0], times[1], times[2], times[3] };

    const secular = updateSecularV4(el, timeV4);

    const emFloor: Vec4 = @splat(eccentricityFloor);
    const decayed = secular.em < emFloor;
    if (@reduce(.Or, decayed)) {
        return Error.SatelliteDecayed;
    }

    const nm: Vec4 = @as(Vec4, @splat(el.grav.xke)) / simdMath.pow15V4(secular.a);
    const kepler = solveKeplerV4(el, secular);
    const corrected = applyShortPeriodCorrectionsV4(el, kepler, nm);
    return computePositionVelocityV4(el, corrected);
}

fn updateSecularV4(el: *const Elements, tsince: Vec4) SecularStateV4 {
    const t2 = tsince * tsince;

    const one: Vec4 = @splat(1.0);
    const cc1V4: Vec4 = @splat(el.cc1);
    const cc4V4: Vec4 = @splat(el.cc4);
    const bstarV4: Vec4 = @splat(el.bstar);
    const t2cofV4: Vec4 = @splat(el.t2cof);

    var tempa = one - cc1V4 * tsince;
    var tempe = bstarV4 * cc4V4 * tsince;
    var templ = t2cofV4 * t2;

    const moV4: Vec4 = @splat(el.mo);
    const mdotV4: Vec4 = @splat(el.mdot);
    const argpoV4: Vec4 = @splat(el.argpo);
    const argpdotV4: Vec4 = @splat(el.argpdot);
    const nodeoV4: Vec4 = @splat(el.nodeo);
    const nodedotV4: Vec4 = @splat(el.nodedot);
    const xnodcfV4: Vec4 = @splat(el.xnodcf);

    const xmdf = moV4 + mdotV4 * tsince;
    const argpdf = argpoV4 + argpdotV4 * tsince;
    const nodedf = nodeoV4 + nodedotV4 * tsince;
    var argpm = argpdf;
    var mm = xmdf;
    var nodem = nodedf + xnodcfV4 * t2;

    if (!el.isimp) {
        const omgcofV4: Vec4 = @splat(el.omgcof);
        const etaV4: Vec4 = @splat(el.eta);
        const xmcofV4: Vec4 = @splat(el.xmcof);
        const delmoV4: Vec4 = @splat(el.delmo);
        const d2V4: Vec4 = @splat(el.d2);
        const d3V4: Vec4 = @splat(el.d3);
        const d4V4: Vec4 = @splat(el.d4);
        const cc5V4: Vec4 = @splat(el.cc5);
        const t3cofV4: Vec4 = @splat(el.t3cof);
        const t4cofV4: Vec4 = @splat(el.t4cof);
        const t5cofV4: Vec4 = @splat(el.t5cof);
        const sinmaoV4: Vec4 = @splat(el.sinmao);

        const delomg = omgcofV4 * tsince;
        const delmtemp = one + etaV4 * simdMath.cosSIMD(xmdf);
        // seems weird but it seems like just std multiply instead of cubing is WAY faster
        const delm = xmcofV4 * (delmtemp * delmtemp * delmtemp - delmoV4);
        const temp = delomg + delm;
        mm = xmdf + temp;
        argpm = argpdf - temp;

        const t3 = t2 * tsince;
        const t4 = t3 * tsince;
        tempa = tempa - d2V4 * t2 - d3V4 * t3 - d4V4 * t4;
        tempe = tempe + bstarV4 * cc5V4 * (simdMath.sinSIMD(mm) - sinmaoV4);
        templ = templ + t3cofV4 * t3 + t4 * (t4cofV4 + tsince * t5cofV4);
    }

    const aBaseV4: Vec4 = @splat(el.aBase);
    const eccoV4: Vec4 = @splat(el.ecco);
    const noUnkozaiV4: Vec4 = @splat(el.noUnkozai);

    const am = aBaseV4 * tempa * tempa;
    const eccFloorV4: Vec4 = @splat(eccentricityFloor);
    var em = eccoV4 - tempe;
    em = @max(em, eccFloorV4);

    mm = mm + noUnkozaiV4 * templ;
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

fn solveKeplerV4(el: *const Elements, sec: SecularStateV4) KeplerStateV4 {
    const one: Vec4 = @splat(1.0);

    const aycofV4: Vec4 = @splat(el.aycof);
    const xlcofV4: Vec4 = @splat(el.xlcof);

    const temp = one / (sec.a * (one - sec.em * sec.em));
    const axnl = sec.em * simdMath.cosSIMD(sec.argpm);
    const aynl = sec.em * simdMath.sinSIMD(sec.argpm) + temp * aycofV4;
    const xl = @mod(sec.mm + sec.argpm + sec.nodem + temp * xlcofV4 * axnl, simdMath.twoPiVec);

    var u = @mod(xl - sec.nodem, simdMath.twoPiVec);
    var eo1 = u;
    var sineo1: Vec4 = @splat(0.0);
    var coseo1: Vec4 = @splat(1.0);

    const tolerance: Vec4 = @splat(1.0e-12);
    const clampVal: Vec4 = @splat(0.95);

    var ktr: u32 = 0;
    while (ktr < 10) : (ktr += 1) {
        sineo1 = simdMath.sinSIMD(eo1);
        coseo1 = simdMath.cosSIMD(eo1);
        var tem5 = one - coseo1 * axnl - sineo1 * aynl;
        tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
        // branch setup
        const positive = tem5 > clampVal;
        const negative = tem5 < -clampVal;
        tem5 = @select(f64, positive, clampVal, tem5);
        tem5 = @select(f64, negative, -clampVal, tem5);
        eo1 = eo1 + tem5;

        // simd convergence checking
        const converged = @abs(tem5) < tolerance;
        if (@reduce(.And, converged)) break;
    }

    const ecose = axnl * coseo1 + aynl * sineo1;
    const esine = axnl * sineo1 - aynl * coseo1;
    const el2 = axnl * axnl + aynl * aynl;
    const pl = sec.a * (one - el2);
    const betal = @sqrt(one - el2);
    const rl = sec.a * (one - ecose);
    const rdotl = @sqrt(sec.a) * esine / rl;
    const rvdotl = @sqrt(pl) / rl;

    const aOverR = sec.a / rl;
    const esineTerm = esine / (one + betal);
    const sinu = aOverR * (sineo1 - aynl - axnl * esineTerm);
    const cosu = aOverR * (coseo1 - axnl + aynl * esineTerm);
    u = simdMath.atan2SIMD(sinu, cosu);

    return .{
        .u = u,
        .r = rl,
        .rdot = rdotl,
        .rvdot = rvdotl,
        .betal = betal,
        .sin2u = @as(Vec4, @splat(2.0)) * sinu * cosu,
        .cos2u = one - @as(Vec4, @splat(2.0)) * sinu * sinu,
        .nodem = sec.nodem,
        .pl = pl,
    };
}

fn applyShortPeriodCorrectionsV4(el: *const Elements, kep: KeplerStateV4, nm: Vec4) CorrectedStateV4 {
    const quarter: Vec4 = @splat(0.25);
    const half: Vec4 = @splat(0.5);
    const one: Vec4 = @splat(1.0);
    const oneHalf: Vec4 = @splat(1.5);
    const sinio: Vec4 = @splat(el.sinio);
    const con41: Vec4 = @splat(el.con41);
    const x1mth2: Vec4 = @splat(el.x1mth2);
    const x7thm1: Vec4 = @splat(el.x7thm1);
    const cosio: Vec4 = @splat(el.cosio);
    const inclo: Vec4 = @splat(el.inclo);
    const j2: Vec4 = @splat(el.grav.j2);
    const xke: Vec4 = @splat(el.grav.xke);

    const temp = one / kep.pl;
    const temp1 = half * j2 * temp;
    const temp2 = temp1 * temp;

    const mrt = kep.r * (one - oneHalf * temp2 * kep.betal * con41) + half * temp1 * x1mth2 * kep.cos2u;
    const su = kep.u - quarter * temp2 * x7thm1 * kep.sin2u;
    const xnode = kep.nodem + oneHalf * temp2 * cosio * kep.sin2u;
    const xinc = inclo + oneHalf * temp2 * cosio * sinio * kep.cos2u;
    const mvt = kep.rdot - nm * temp1 * x1mth2 * kep.sin2u / xke;
    const rvdot = kep.rvdot + nm * temp1 * (x1mth2 * kep.cos2u + oneHalf * con41) / xke;

    return .{ .r = mrt, .rdot = mvt, .rvdot = rvdot, .u = su, .xnode = xnode, .xinc = xinc };
}

fn computePositionVelocityV4(el: *const Elements, state: CorrectedStateV4) [4][2][3]f64 {
    const radiusV4: Vec4 = @splat(el.grav.radiusEarthKm);
    const vkmpersecV4: Vec4 = @splat(el.vkmpersec);

    const sinsu = simdMath.sinSIMD(state.u);
    const cossu = simdMath.cosSIMD(state.u);
    const snod = simdMath.sinSIMD(state.xnode);
    const cnod = simdMath.cosSIMD(state.xnode);
    const sini = simdMath.sinSIMD(state.xinc);
    const cosi = simdMath.cosSIMD(state.xinc);

    const xmx = -snod * cosi;
    const xmy = cnod * cosi;
    const ux = xmx * sinsu + cnod * cossu;
    const uy = xmy * sinsu + snod * cossu;
    const uz = sini * sinsu;
    const vx = xmx * cossu - cnod * sinsu;
    const vy = xmy * cossu - snod * sinsu;
    const vz = sini * cossu;

    const rScaled = state.r * radiusV4;
    const rx = rScaled * ux;
    const ry = rScaled * uy;
    const rz = rScaled * uz;

    const vxOut = (state.rdot * ux + state.rvdot * vx) * vkmpersecV4;
    const vyOut = (state.rdot * uy + state.rvdot * vy) * vkmpersecV4;
    const vzOut = (state.rdot * uz + state.rvdot * vz) * vkmpersecV4;

    // need to extract the return from the SIMD vectors
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

test "sgp4 basic init" {
    const test_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;

    var tle = try Tle.parse(test_tle, std.testing.allocator);
    defer tle.deinit();

    const sgp4 = try Sgp4.init(tle, constants.wgs84);

    try std.testing.expect(sgp4.elements.ecco > 0.0);
    try std.testing.expect(sgp4.elements.ecco < 1.0);
    try std.testing.expect(sgp4.elements.noUnkozai > 0.0);
}

test "sgp4 propagate basic" {
    const test_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;

    var tle = try Tle.parse(test_tle, std.testing.allocator);
    defer tle.deinit();

    const sgp4 = try Sgp4.init(tle, constants.wgs84);
    const result = try sgp4.propagate(0.0);
    const r = result[0];
    const v = result[1];

    const r_mag = calculations.mag(r);
    const v_mag = calculations.mag(v);

    try std.testing.expect(r_mag > 6000.0 and r_mag < 8000.0);
    try std.testing.expect(v_mag > 6.0 and v_mag < 9.0);
}

test "sgp4 vallado reference ISS" {
    // ISS TLE
    const test_tle =
        \\1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995
        \\2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123
    ;

    var tle = try Tle.parse(test_tle, std.testing.allocator);
    defer tle.deinit();

    const sgp4 = try Sgp4.init(tle, constants.wgs84);
    const el = &sgp4.elements;

    // verify initialization matches python-sgp4
    try std.testing.expectApproxEqAbs(@as(f64, 0.06767329492593213), el.noKozai, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.064977141044385), el.a, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.067673302731475), el.mdot, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.000044767460455), el.argpdot, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, -6.016088837547746e-05), el.nodedot, 1e-12);

    const result = try sgp4.propagate(0.0);
    const r = result[0];
    const v = result[1];

    // reference from python-sgp4:
    // r = [-5887.061832, 3151.888264, -1263.887271]
    // v = [-3.250642, -3.745001, 5.837125]
    const ref_r = [3]f64{ -5887.061832, 3151.888264, -1263.887271 };
    const ref_v = [3]f64{ -3.250642, -3.745001, 5.837125 };

    // Check position error (should be < 1 meter = 0.001 km)
    const pos_err = @sqrt((r[0] - ref_r[0]) * (r[0] - ref_r[0]) +
        (r[1] - ref_r[1]) * (r[1] - ref_r[1]) +
        (r[2] - ref_r[2]) * (r[2] - ref_r[2]));
    const vel_err = @sqrt((v[0] - ref_v[0]) * (v[0] - ref_v[0]) +
        (v[1] - ref_v[1]) * (v[1] - ref_v[1]) +
        (v[2] - ref_v[2]) * (v[2] - ref_v[2]));

    // should be sub-meter accuracy (< 0.001 km = 1 meter)
    try std.testing.expect(pos_err < 0.001);
    // velocity error should be < 10 mm/s = 0.00001 km/s
    try std.testing.expect(vel_err < 0.00001);
}

test "SIMD matches scalar" {
    const allocator = std.testing.allocator;

    // ISS TLE
    const tle_str =
        \\1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995
        \\2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123
    ;

    var tle = try Tle.parse(tle_str, allocator);
    defer tle.deinit();

    var sgp4 = try Sgp4.init(tle, constants.wgs84);

    const times = [4]f64{ 0.0, 10.0, 60.0, 1440.0 };

    // scalar run
    var scalar_results: [4][2][3]f64 = undefined;
    for (0..4) |i| {
        const result = try sgp4.propagate(times[i]);
        scalar_results[i][0] = result[0];
        scalar_results[i][1] = result[1];
    }

    // SIMD run
    const simd_results = try sgp4.propagateV4(times);

    const tol = 1e-9;
    for (0..4) |i| {
        for (0..3) |j| {
            try std.testing.expectApproxEqAbs(scalar_results[i][0][j], simd_results[i][0][j], tol);
            try std.testing.expectApproxEqAbs(scalar_results[i][1][j], simd_results[i][1][j], tol);
        }
    }
}

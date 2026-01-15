//! SGP4 Orbit Propagator
//! The magic numbers in this implementation come from the analytical derivation
//! of the SGP4 equations and are documented in the original Spacetrack Report #3.

const std = @import("std");
const constants = @import("../constants.zig");
const calculations = @import("../calculations.zig");
const Datetime = @import("../Datetime.zig");
const Tle = @import("../Tle.zig");

const Sgp4 = @This();

// gave these names so my head doesn't explode
const kozai_j2_coeff = 1.5;
const kozai_second_order = 1.0 / 3.0;
const kozai_third_order = 134.0 / 81.0;

const secular_j2_coeff = 1.5;
const secular_j4_coeff = 0.46875;
const secular_poly_coeff = 0.0625;

const mdot_coeff_13 = 13.0;
const mdot_coeff_78 = 78.0;
const mdot_coeff_137 = 137.0;

const argpdot_coeff_7 = 7.0;
const argpdot_coeff_114 = 114.0;
const argpdot_coeff_395 = 395.0;
const argpdot_coeff_3 = 3.0;
const argpdot_coeff_36 = 36.0;
const argpdot_coeff_49 = 49.0;

const nodedot_coeff_4 = 4.0;
const nodedot_coeff_19 = 19.0;

const perigee_threshold_high = 156.0;
const perigee_threshold_low = 98.0;
const atm_param_s_standard = 78.0;
const atm_param_qo = 120.0;
const atm_param_s_min = 20.0;
const perigee_simple_threshold = 220.0;

const drag_cc2_etasq_coeff = 1.5;
const drag_cc2_j2_coeff = 0.375;
const drag_cc4_inner_coeff = 2.0;
const drag_xnodcf_coeff = 3.5;
const drag_t2cof_factor = 1.5;
const drag_xlcof_coeff = 0.25;
const drag_aycof_coeff = 0.5;
const drag_cc5_coeff = 2.75;

const ho_d3_coeff = 17.0;
const ho_d4_coeff_a = 221.0;
const ho_d4_coeff_s = 31.0;
const ho_t4_d2_coeff = 12.0;
const ho_t4_cc1_coeff = 10.0;
const ho_t5_coeff = 0.2;
const ho_t5_d3_coeff = 12.0;
const ho_t5_d2_coeff = 6.0;
const ho_t5_cc1_coeff = 15.0;

const sp_betal_coeff = 1.5;

const eccentricity_min = 1.0e-4;
const eccentricity_floor = 1.0e-6;
const kepler_tolerance = 1.0e-12;
const singularity_tolerance = 1.5e-12;

pub const Error = error{
    DeepSpaceNotSupported,
    InvalidEccentricity,
    SatelliteDecayed,
    InvalidMeanMotion,
};

pub const Elements = struct {
    grav: constants.Sgp4GravityModel,
    epoch_jd: f64,

    // Mean elements from TLE (radians where applicable)
    no_kozai: f64,
    ecco: f64,
    inclo: f64,
    nodeo: f64,
    argpo: f64,
    mo: f64,
    bstar: f64,

    // Derived mean elements
    no_unkozai: f64,
    a: f64,

    // Trigonometric terms
    sinio: f64,
    cosio: f64,
    cosio2: f64,
    cosio4: f64,

    // Polynomial terms
    con41: f64,
    con42: f64,
    x1mth2: f64,
    x3thm1: f64,
    x7thm1: f64,

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

    isimp: bool,
    method: enum { near_earth },
};

grav: constants.Sgp4GravityModel,
elements: Elements,

/// initialize from tle
pub fn init(tle: Tle, grav: constants.Sgp4GravityModel) Error!Sgp4 {
    const elements = try initElements(tle, grav);
    return .{ .grav = grav, .elements = elements };
}

/// propagate to time tsince (minutes from TLE epoch)
pub fn propagate(self: *const Sgp4, tsince: f64) Error![2][3]f64 {
    return propagateElements(&self.elements, tsince);
}

fn initElements(tle: Tle, grav: constants.Sgp4GravityModel) Error!Elements {
    const mean_elements = extractMeanElements(tle);

    if (mean_elements.ecco < 0.0 or mean_elements.ecco >= 1.0) {
        return Error.InvalidEccentricity;
    }

    const recovered = recoverMeanMotion(mean_elements, grav);

    const rp = recovered.a * (1.0 - mean_elements.ecco);
    if (rp < 1.0) return Error.SatelliteDecayed;

    const period = constants.twoPi / recovered.no_unkozai;
    if (period > constants.sgp4DeepSpaceThresholdMinutes) {
        return Error.DeepSpaceNotSupported;
    }

    const trig = computeTrigTerms(mean_elements.inclo);
    const poly = computePolyTerms(trig, grav);
    const secular = computeSecularRates(mean_elements, recovered, trig, poly, grav);
    const drag = computeDragCoefficients(mean_elements, recovered, trig, poly, grav);
    const higher_order = computeHigherOrderDrag(recovered, drag, grav);

    const fullYear: u16 = if (tle.firstLine.epochYear < 57)
        2000 + tle.firstLine.epochYear
    else
        1900 + tle.firstLine.epochYear;
    const epoch_jd = Datetime.yearDoyToJulianDate(fullYear, @floatCast(tle.firstLine.epochDay));

    return Elements{
        .grav = grav,
        .epoch_jd = epoch_jd,
        .no_kozai = mean_elements.no_kozai,
        .ecco = mean_elements.ecco,
        .inclo = mean_elements.inclo,
        .nodeo = mean_elements.nodeo,
        .argpo = mean_elements.argpo,
        .mo = mean_elements.mo,
        .bstar = mean_elements.bstar,
        .no_unkozai = recovered.no_unkozai,
        .a = recovered.a,
        .sinio = trig.sinio,
        .cosio = trig.cosio,
        .cosio2 = trig.cosio2,
        .cosio4 = trig.cosio4,
        .con41 = poly.con41,
        .con42 = poly.con42,
        .x1mth2 = poly.x1mth2,
        .x3thm1 = poly.x3thm1,
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
        .d2 = higher_order.d2,
        .d3 = higher_order.d3,
        .d4 = higher_order.d4,
        .t3cof = higher_order.t3cof,
        .t4cof = higher_order.t4cof,
        .t5cof = higher_order.t5cof,
        .isimp = higher_order.isimp,
        .method = .near_earth,
    };
}

const MeanElements = struct {
    no_kozai: f64,
    ecco: f64,
    inclo: f64,
    nodeo: f64,
    argpo: f64,
    mo: f64,
    bstar: f64,
};

fn extractMeanElements(tle: Tle) MeanElements {
    return .{
        .no_kozai = tle.secondLine.mMotion * constants.twoPi / constants.minutes_per_day,
        .ecco = @floatCast(tle.secondLine.eccentricity),
        .inclo = @as(f64, @floatCast(tle.secondLine.inclination)) * constants.deg2rad,
        .nodeo = @as(f64, @floatCast(tle.secondLine.rightAscension)) * constants.deg2rad,
        .argpo = @as(f64, @floatCast(tle.secondLine.perigee)) * constants.deg2rad,
        .mo = @as(f64, @floatCast(tle.secondLine.mAnomaly)) * constants.deg2rad,
        .bstar = @floatCast(tle.firstLine.bstarDrag),
    };
}

const RecoveredElements = struct { no_unkozai: f64, a: f64 };

fn recoverMeanMotion(el: MeanElements, grav: constants.Sgp4GravityModel) RecoveredElements {
    const cosio = @cos(el.inclo);
    const theta2 = cosio * cosio;
    const x3thm1 = 3.0 * theta2 - 1.0;
    const eosq = el.ecco * el.ecco;
    const betao2 = 1.0 - eosq;
    const betao = @sqrt(betao2);

    const a1 = std.math.pow(f64, grav.xke / el.no_kozai, 2.0 / 3.0);
    const del1 = kozai_j2_coeff * grav.j2 * x3thm1 / (a1 * a1 * betao * betao2);
    const ao = a1 * (1.0 - del1 * (kozai_second_order + del1 * (1.0 + kozai_third_order * del1)));
    const delo = kozai_j2_coeff * grav.j2 * x3thm1 / (ao * ao * betao * betao2);

    const no_unkozai = el.no_kozai / (1.0 + delo);
    const a = std.math.pow(f64, grav.xke / no_unkozai, 2.0 / 3.0);

    return .{ .no_unkozai = no_unkozai, .a = a };
}

const TrigTerms = struct { sinio: f64, cosio: f64, cosio2: f64, cosio4: f64 };

fn computeTrigTerms(inclo: f64) TrigTerms {
    const sinio = @sin(inclo);
    const cosio = @cos(inclo);
    const cosio2 = cosio * cosio;
    return .{ .sinio = sinio, .cosio = cosio, .cosio2 = cosio2, .cosio4 = cosio2 * cosio2 };
}

const PolyTerms = struct { con41: f64, con42: f64, x1mth2: f64, x3thm1: f64, x7thm1: f64 };

fn computePolyTerms(trig: TrigTerms, grav: constants.Sgp4GravityModel) PolyTerms {
    const theta2 = trig.cosio2;
    return .{
        .con41 = -grav.j2 * (1.0 - 3.0 * theta2),
        .con42 = -1.5 * grav.j2 * (1.0 - 5.0 * theta2),
        .x1mth2 = 1.0 - theta2,
        .x3thm1 = 3.0 * theta2 - 1.0,
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
    const eosq = el.ecco * el.ecco;
    const betao2 = 1.0 - eosq;
    const betao = @sqrt(betao2);

    const pinvsq = 1.0 / (rec.a * rec.a);
    const temp1 = secular_j2_coeff * grav.j2 * pinvsq * rec.no_unkozai / (betao * betao2);
    const temp2 = drag_aycof_coeff * temp1 * grav.j2 * pinvsq / betao2;
    const temp3 = -secular_j4_coeff * grav.j4 * pinvsq * pinvsq * rec.no_unkozai;

    const mdot = rec.no_unkozai + temp1 * (1.0 + drag_aycof_coeff * eosq) * poly.x3thm1 +
        secular_poly_coeff * temp2 * poly.x1mth2 * (mdot_coeff_13 - mdot_coeff_78 * trig.cosio2 + mdot_coeff_137 * trig.cosio4);

    const argpdot = -drag_aycof_coeff * temp1 * poly.con41 +
        secular_poly_coeff * temp2 * (argpdot_coeff_7 - argpdot_coeff_114 * trig.cosio2 + argpdot_coeff_395 * trig.cosio4) +
        temp3 * (argpdot_coeff_3 - argpdot_coeff_36 * trig.cosio2 + argpdot_coeff_49 * trig.cosio4);

    const nodedot = -temp1 * trig.cosio +
        drag_aycof_coeff * temp2 * trig.cosio * (nodedot_coeff_4 - nodedot_coeff_19 * trig.cosio2) +
        2.0 * temp3 * trig.cosio * (argpdot_coeff_3 - argpdot_coeff_7 * trig.cosio2);

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

// brain blast
fn computeDragCoefficients(
    el: MeanElements,
    rec: RecoveredElements,
    trig: TrigTerms,
    poly: PolyTerms,
    grav: constants.Sgp4GravityModel,
) DragCoefficients {
    const eosq = el.ecco * el.ecco;
    const betao2 = 1.0 - eosq;
    const rp = rec.a * (1.0 - el.ecco);
    const perige = (rp - 1.0) * grav.radiusEarthKm;

    var s: f64 = undefined;
    var qoms24: f64 = undefined;
    if (perige < perigee_threshold_high) {
        s = if (perige < perigee_threshold_low) atm_param_s_min else perige - atm_param_s_standard;
        qoms24 = std.math.pow(f64, (atm_param_qo - s) / grav.radiusEarthKm, 4.0);
        s = s / grav.radiusEarthKm + 1.0;
    } else {
        s = atm_param_s_standard / grav.radiusEarthKm + 1.0;
        qoms24 = std.math.pow(f64, (atm_param_qo - atm_param_s_standard) / grav.radiusEarthKm, 4.0);
    }

    const pinvsq = 1.0 / (rec.a * rec.a * betao2 * betao2);
    const tsi = 1.0 / (rec.a - s);
    const eta = rec.a * el.ecco * tsi;
    const etasq = eta * eta;
    const eeta = el.ecco * eta;
    const psisq = @abs(1.0 - etasq);
    const coef = qoms24 * std.math.pow(f64, tsi, 4.0);
    const coef1 = coef / std.math.pow(f64, psisq, drag_xnodcf_coeff);

    const cc2 = coef1 * rec.no_unkozai * (rec.a * (1.0 + drag_cc2_etasq_coeff * etasq + eeta * (4.0 + etasq)) +
        drag_cc2_j2_coeff * grav.j2 * tsi / psisq * poly.x3thm1 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
    const cc1 = el.bstar * cc2;
    const cc3 = if (el.ecco > eccentricity_min)
        -2.0 * coef * tsi * grav.j3oj2 * rec.no_unkozai * trig.sinio / el.ecco
    else
        0.0;

    // cc4: secular drag coefficient - broken into components for readability
    const cc4_scale = 2.0 * rec.no_unkozai * coef1 * rec.a * betao2;
    const cc4_eta_term = eta * (2.0 + drag_aycof_coeff * etasq) + el.ecco * (drag_aycof_coeff + 2.0 * etasq);
    const cc4_j2_scale = grav.j2 * tsi / (rec.a * psisq);
    const cc4_x3thm1_term = -3.0 * poly.x3thm1 * (1.0 - 2.0 * eeta + etasq * (drag_cc2_etasq_coeff - drag_aycof_coeff * eeta));
    const cc4_x1mth2_term = 0.75 * poly.x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) * @cos(2.0 * el.argpo);
    const cc4 = cc4_scale * (cc4_eta_term - cc4_j2_scale * (cc4_x3thm1_term + cc4_x1mth2_term));

    const cc5 = 2.0 * coef1 * rec.a * betao2 * (1.0 + drag_cc5_coeff * (etasq + eeta) + eeta * etasq);

    const temp1 = 3.0 * grav.j2 * pinvsq * rec.no_unkozai;
    const xnodedot = -drag_aycof_coeff * temp1 * (-drag_cc2_etasq_coeff * grav.j2 * (1.0 - 5.0 * trig.cosio2));
    const xnodcf = drag_xnodcf_coeff * betao2 * xnodedot * cc1;
    const t2cof = drag_t2cof_factor * cc1;

    const xlcof = if (@abs(trig.cosio + 1.0) > singularity_tolerance)
        -drag_xlcof_coeff * grav.j3oj2 * trig.sinio * (3.0 + 5.0 * trig.cosio) / (1.0 + trig.cosio)
    else
        -drag_xlcof_coeff * grav.j3oj2 * trig.sinio * (3.0 + 5.0 * trig.cosio) / singularity_tolerance;

    const aycof = -drag_aycof_coeff * grav.j3oj2 * trig.sinio;
    const delmo = std.math.pow(f64, 1.0 + eta * @cos(el.mo), 3.0);
    const sinmao = @sin(el.mo);
    const xmcof = if (@abs(el.ecco) > eccentricity_min)
        -constants.twoPi * coef * el.bstar / (3.0 * eeta)
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
    const isimp = drag.perige < perigee_simple_threshold;

    if (isimp) {
        return .{ .d2 = 0, .d3 = 0, .d4 = 0, .t3cof = 0, .t4cof = 0, .t5cof = 0, .isimp = true };
    }

    const s = atm_param_s_standard / grav.radiusEarthKm + 1.0;
    const tsi = 1.0 / (rec.a - s);

    const cc1sq = drag.cc1 * drag.cc1;
    const d2 = 4.0 * rec.a * tsi * cc1sq;
    const temp = d2 * tsi * drag.cc1 / 3.0;
    const d3 = (ho_d3_coeff * rec.a + s) * temp;
    const d4 = drag_aycof_coeff * temp * rec.a * tsi * (ho_d4_coeff_a * rec.a + ho_d4_coeff_s * s) * drag.cc1;
    const t3cof = d2 + 2.0 * cc1sq;
    const t4cof = drag_xlcof_coeff * (3.0 * d3 + drag.cc1 * (ho_t4_d2_coeff * d2 + ho_t4_cc1_coeff * cc1sq));
    const t5cof = ho_t5_coeff * (3.0 * d4 + ho_t5_d3_coeff * drag.cc1 * d3 + ho_t5_d2_coeff * d2 * d2 + ho_t5_cc1_coeff * cc1sq * (2.0 * d2 + cc1sq));

    return .{ .d2 = d2, .d3 = d3, .d4 = d4, .t3cof = t3cof, .t4cof = t4cof, .t5cof = t5cof, .isimp = false };
}

fn propagateElements(el: *const Elements, tsince: f64) Error![2][3]f64 {
    const secular = updateSecular(el, tsince);

    const kepler = solveKepler(el, secular);

    const corrected = applyShortPeriodCorrections(el, kepler);

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

    if (!el.isimp) {
        const t3 = t2 * tsince;
        const t4 = t3 * tsince;
        tempa = tempa - el.d2 * t2 - el.d3 * t3 - el.d4 * t4;
        tempe = tempe + el.bstar * el.cc5 * (@sin(xmdf) - el.sinmao);
        templ = templ + el.t3cof * t3 + t4 * (el.t4cof + tsince * el.t5cof);
    }

    const a = std.math.pow(f64, el.grav.xke / el.no_unkozai, 2.0 / 3.0) * tempa * tempa;
    var em = el.ecco - tempe;
    if (em < eccentricity_floor) em = eccentricity_floor;

    var mm = xmdf + el.no_unkozai * templ;
    const argpm = el.argpo + el.argpdot * tsince;
    var nodem = el.nodeo + el.nodedot * tsince + el.xnodcf * t2;

    mm = mm + el.no_unkozai * templ;
    const xlm = mm + argpm + nodem;
    nodem = @mod(nodem, constants.twoPi);
    const argpm_mod = @mod(argpm, constants.twoPi);
    mm = @mod(xlm - argpm_mod - nodem, constants.twoPi);

    return .{ .mm = mm, .argpm = argpm_mod, .nodem = nodem, .em = em, .a = a };
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
};

fn solveKepler(el: *const Elements, sec: SecularState) KeplerState {
    const axnl = sec.em * @cos(sec.argpm);
    const aynl = sec.em * @sin(sec.argpm) + el.aycof;
    const xl = @mod(sec.mm + sec.argpm + sec.nodem + el.xlcof * axnl, constants.twoPi);

    var u = @mod(xl - sec.nodem, constants.twoPi);
    var eo1 = u;
    var sineo1: f64 = 0.0;
    var coseo1: f64 = 1.0;

    for (0..10) |_| {
        sineo1 = @sin(eo1);
        coseo1 = @cos(eo1);
        const f_val = u + aynl * coseo1 - axnl * sineo1 - eo1;
        if (@abs(f_val) < kepler_tolerance) break;
        eo1 = eo1 + f_val / (1.0 - coseo1 * axnl - sineo1 * aynl);
    }

    const ecose = axnl * coseo1 + aynl * sineo1;
    const esine = axnl * sineo1 - aynl * coseo1;
    const el2 = axnl * axnl + aynl * aynl;
    const r = sec.a * (1.0 - ecose);
    const rdot = @sqrt(sec.a) * esine / r;
    const rvdot = @sqrt(sec.a * (1.0 - el2)) / r;
    const betal = @sqrt(1.0 - el2);

    const sinu = sec.a / r * (sineo1 - aynl - axnl * esine / (1.0 + betal));
    const cosu = sec.a / r * (coseo1 - axnl + aynl * esine / (1.0 + betal));
    u = std.math.atan2(sinu, cosu);

    return .{
        .u = u,
        .r = r,
        .rdot = rdot,
        .rvdot = rvdot,
        .betal = betal,
        .sin2u = 2.0 * sinu * cosu,
        .cos2u = 2.0 * cosu * cosu - 1.0,
        .nodem = sec.nodem,
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

fn applyShortPeriodCorrections(el: *const Elements, kep: KeplerState) CorrectedState {
    const temp1 = drag_aycof_coeff * el.grav.j2 / (kep.r * (1.0 - el.ecco * el.ecco));
    const temp2 = temp1 / (kep.r * (1.0 - el.ecco * el.ecco));

    const sinim = el.sinio;
    const cosim = el.cosio;

    const r = kep.r * (1.0 - sp_betal_coeff * temp2 * kep.betal * el.x3thm1) + drag_aycof_coeff * temp1 * el.x1mth2 * kep.cos2u;
    const u = kep.u - drag_xlcof_coeff * temp2 * el.x7thm1 * kep.sin2u;
    const xnode = kep.nodem + sp_betal_coeff * temp2 * cosim * kep.sin2u;
    const xinc = el.inclo + sp_betal_coeff * temp2 * cosim * sinim * kep.cos2u;

    const rdot = kep.rdot - el.no_unkozai * temp1 * el.x1mth2 * kep.sin2u / el.grav.xke;
    const rvdot = kep.rvdot + el.no_unkozai * temp1 * (el.x1mth2 * kep.cos2u + sp_betal_coeff * el.x3thm1) / el.grav.xke;

    return .{ .r = r, .rdot = rdot, .rvdot = rvdot, .u = u, .xnode = xnode, .xinc = xinc };
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

    const vkmpersec = el.grav.xke * el.grav.radiusEarthKm / 60.0;

    const r: [3]f64 = .{
        state.r * ux * el.grav.radiusEarthKm,
        state.r * uy * el.grav.radiusEarthKm,
        state.r * uz * el.grav.radiusEarthKm,
    };
    const v: [3]f64 = .{
        (state.rdot * ux + state.rvdot * vx) * vkmpersec,
        (state.rdot * uy + state.rvdot * vy) * vkmpersec,
        (state.rdot * uz + state.rvdot * vz) * vkmpersec,
    };

    return .{ r, v };
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
    try std.testing.expect(sgp4.elements.no_unkozai > 0.0);
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

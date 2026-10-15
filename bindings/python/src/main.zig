//! Python module entry point for astroz

const std = @import("std");
const py = @import("python.zig");
const c = py.c;
const astroz = @import("astroz");
const tle = @import("tle.zig");
const sgp4 = @import("sgp4.zig");
const satrec = @import("satrec.zig");
const screening = @import("screening.zig");
const orbital_mechanics = @import("orbital_mechanics.zig");
const propagator = @import("propagator.zig");

var module_methods = [_]c.PyMethodDef{
    .{ .ml_name = "screen_candidates", .ml_meth = @ptrCast(&screening.pyScreenCandidates), .ml_flags = c.METH_VARARGS | c.METH_KEYWORDS, .ml_doc = "screen_candidates(positions, velocities, num_sats, dt, threshold, primary=None) -> bytes\n\nCandidate (i, j, step) intervals, packed as uint32 triples, that may hold a closest approach within threshold km. positions/velocities are time-major (n_times, num_sats, 3) TEME on a grid dt seconds apart." },
    .{ .ml_name = "refine_candidates", .ml_meth = @ptrCast(&screening.pyRefineCandidates), .ml_flags = c.METH_VARARGS | c.METH_KEYWORDS, .ml_doc = "refine_candidates(satrecs, candidates, jd0, dt) -> bytes\n\nTime of closest approach in each (i, j, step) candidate from screen_candidates, with satrecs[s] propagating satellite s. Returns float64 quadruples (time s, distance km, speed km/s, converged)." },
    .{ .ml_name = "jday", .ml_meth = @ptrCast(&satrec.pyJday), .ml_flags = c.METH_VARARGS, .ml_doc = "jday(year, month, day, hour, minute, second) -> (jd, fr)\n\nConvert calendar date to Julian date (integer and fractional parts)." },
    .{ .ml_name = "days2mdhms", .ml_meth = @ptrCast(&satrec.pyDays2mdhms), .ml_flags = c.METH_VARARGS, .ml_doc = "days2mdhms(year, days) -> (month, day, hour, minute, second)\n\nConvert day of year to month, day, hour, minute, second." },
    .{ .ml_name = "sdp4_batch_propagate_into", .ml_meth = @ptrCast(&satrec.pySdp4BatchPropagateInto), .ml_flags = c.METH_VARARGS | c.METH_KEYWORDS, .ml_doc = "sdp4_batch_propagate_into(satrecs, jd, fr, positions, velocities, output_stride=-1, sat_offset=0) -> None\n\nBatch SDP4 propagation with threading." },
    .{ .ml_name = "hohmann_transfer", .ml_meth = @ptrCast(&orbital_mechanics.pyHohmannTransfer), .ml_flags = c.METH_VARARGS, .ml_doc = "hohmann_transfer(mu, r1, r2) -> dict\n\nCompute Hohmann transfer between two circular orbits." },
    .{ .ml_name = "bi_elliptic_transfer", .ml_meth = @ptrCast(&orbital_mechanics.pyBiEllipticTransfer), .ml_flags = c.METH_VARARGS, .ml_doc = "bi_elliptic_transfer(mu, r1, r2, r_intermediate) -> dict\n\nCompute bi-elliptic transfer between two circular orbits." },
    .{ .ml_name = "lambert", .ml_meth = @ptrCast(&orbital_mechanics.pyLambert), .ml_flags = c.METH_VARARGS, .ml_doc = "lambert(mu, r1, r2, tof) -> dict\n\nSolve Lambert's problem for two position vectors and time of flight." },
    .{ .ml_name = "orbital_velocity", .ml_meth = @ptrCast(&orbital_mechanics.pyOrbitalVelocity), .ml_flags = c.METH_VARARGS | c.METH_KEYWORDS, .ml_doc = "orbital_velocity(mu, radius, sma=None) -> float\n\nCompute orbital velocity (vis-viva). Circular if sma omitted." },
    .{ .ml_name = "orbital_period", .ml_meth = @ptrCast(&orbital_mechanics.pyOrbitalPeriod), .ml_flags = c.METH_VARARGS, .ml_doc = "orbital_period(mu, sma) -> float\n\nCompute orbital period in seconds from Kepler's third law." },
    .{ .ml_name = "escape_velocity", .ml_meth = @ptrCast(&orbital_mechanics.pyEscapeVelocity), .ml_flags = c.METH_VARARGS, .ml_doc = "escape_velocity(mu, radius) -> float\n\nCompute escape velocity at given radius." },
    .{ .ml_name = "propagate_numerical", .ml_meth = @ptrCast(&propagator.pyPropagateNumerical), .ml_flags = c.METH_VARARGS | c.METH_KEYWORDS, .ml_doc = "propagate_numerical(state, t0, duration, dt, mu, **kwargs) -> (times, states)\n\nNumerically propagate an orbit with optional J2 and drag perturbations." },
    .{ .ml_name = null, .ml_meth = null, .ml_flags = 0, .ml_doc = null },
};

var module_def = py.PyModuleDef{
    .m_base = .{ .ob_base = py.makeBaseObject(), .m_init = null, .m_index = 0, .m_copy = null },
    .m_name = "_astroz",
    .m_doc = "High-performance astrodynamics",
    .m_size = -1,
    .m_methods = &module_methods,
    .m_slots = null,
    .m_traverse = null,
    .m_clear = null,
    .m_free = null,
};

pub export fn PyInit__astroz() ?*c.PyObject {
    if (tle.ready() < 0 or sgp4.ready() < 0 or satrec.ready() < 0) return null;
    const m = py.moduleCreate(&module_def) orelse return null;

    py.incref(@as(*c.PyObject, @ptrCast(&tle.TleType)));
    if (c.PyModule_AddObject(m, "Tle", @ptrCast(&tle.TleType)) < 0) return null;

    py.incref(@as(*c.PyObject, @ptrCast(&sgp4.Sgp4ConstellationType)));
    if (c.PyModule_AddObject(m, "Sgp4Constellation", @ptrCast(&sgp4.Sgp4ConstellationType)) < 0) return null;

    py.incref(@as(*c.PyObject, @ptrCast(&satrec.SatrecType)));
    if (c.PyModule_AddObject(m, "Satrec", @ptrCast(&satrec.SatrecType)) < 0) return null;

    py.incref(@as(*c.PyObject, @ptrCast(&satrec.SatrecArrayType)));
    if (c.PyModule_AddObject(m, "SatrecArray", @ptrCast(&satrec.SatrecArrayType)) < 0) return null;

    if (c.PyModule_AddStringConstant(m, "__version__", astroz.version.ptr) < 0) return null;
    if (c.PyModule_AddIntConstant(m, "WGS72", satrec.WGS72) < 0) return null;
    if (c.PyModule_AddIntConstant(m, "WGS84", satrec.WGS84) < 0) return null;

    if (c.PyModule_AddObject(m, "EARTH_MU", c.PyFloat_FromDouble(astroz.constants.wgs84.mu)) < 0) return null;
    if (c.PyModule_AddObject(m, "EARTH_R_EQ", c.PyFloat_FromDouble(astroz.constants.wgs84.radiusEarthKm)) < 0) return null;
    if (c.PyModule_AddObject(m, "EARTH_J2", c.PyFloat_FromDouble(astroz.constants.wgs84.j2)) < 0) return null;
    if (c.PyModule_AddObject(m, "SUN_MU", c.PyFloat_FromDouble(astroz.constants.sun.mu)) < 0) return null;
    if (c.PyModule_AddObject(m, "MOON_MU", c.PyFloat_FromDouble(astroz.constants.moon.mu)) < 0) return null;

    return m;
}

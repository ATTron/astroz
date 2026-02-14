//! Python module entry point for astroz

const std = @import("std");
const py = @import("python.zig");
const c = py.c;
const tle = @import("tle.zig");
const sgp4 = @import("sgp4.zig");
const satrec = @import("satrec.zig");
const conjunction = @import("conjunction.zig");

fn astroz_version(_: [*c]c.PyObject, _: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    return c.PyUnicode_FromString("0.7.1");
}

var module_methods = [_]c.PyMethodDef{
    .{ .ml_name = "version", .ml_meth = @ptrCast(&astroz_version), .ml_flags = c.METH_NOARGS, .ml_doc = "Return version" },
    .{ .ml_name = "coarse_screen", .ml_meth = @ptrCast(&conjunction.py_coarse_screen), .ml_flags = c.METH_VARARGS, .ml_doc = "coarse_screen(positions, num_sats, threshold, [valid_mask]) -> (pairs, t_)\n\nFind all satellite pairs within threshold distance at any time step using cell-list spatial indexing." },
    .{ .ml_name = "jday", .ml_meth = @ptrCast(&satrec.pyJday), .ml_flags = c.METH_VARARGS, .ml_doc = "jday(year, month, day, hour, minute, second) -> (jd, fr)\n\nConvert calendar date to Julian date (integer and fractional parts)." },
    .{ .ml_name = "days2mdhms", .ml_meth = @ptrCast(&satrec.pyDays2mdhms), .ml_flags = c.METH_VARARGS, .ml_doc = "days2mdhms(year, days) -> (month, day, hour, minute, second)\n\nConvert day of year to month, day, hour, minute, second." },
    .{ .ml_name = "sdp4_batch_propagate_into", .ml_meth = @ptrCast(&satrec.pySdp4BatchPropagateInto), .ml_flags = c.METH_VARARGS | c.METH_KEYWORDS, .ml_doc = "sdp4_batch_propagate_into(satrecs, jd, fr, positions, velocities, output_stride=-1, sat_offset=0) -> None\n\nBatch SDP4 propagation with threading." },
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

    c.Py_INCREF(@as(*c.PyObject, @ptrCast(&tle.TleType)));
    if (c.PyModule_AddObject(m, "Tle", @ptrCast(&tle.TleType)) < 0) return null;

    c.Py_INCREF(@as(*c.PyObject, @ptrCast(&sgp4.Sgp4ConstellationType)));
    if (c.PyModule_AddObject(m, "Sgp4Constellation", @ptrCast(&sgp4.Sgp4ConstellationType)) < 0) return null;

    // python-sgp4 compatible API
    c.Py_INCREF(@as(*c.PyObject, @ptrCast(&satrec.SatrecType)));
    if (c.PyModule_AddObject(m, "Satrec", @ptrCast(&satrec.SatrecType)) < 0) return null;

    c.Py_INCREF(@as(*c.PyObject, @ptrCast(&satrec.SatrecArrayType)));
    if (c.PyModule_AddObject(m, "SatrecArray", @ptrCast(&satrec.SatrecArrayType)) < 0) return null;

    // Gravity model constants
    if (c.PyModule_AddIntConstant(m, "WGS72", satrec.WGS72) < 0) return null;
    if (c.PyModule_AddIntConstant(m, "WGS84", satrec.WGS84) < 0) return null;

    return m;
}

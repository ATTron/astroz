//! Satrec and SatrecArray: the python-sgp4 compatible API.

const std = @import("std");
const py = @import("python.zig");
const c = py.c;
const shared = @import("shared.zig");
const astroz = @import("astroz");
const Tle = astroz.Tle;
const Sgp4 = astroz.Sgp4;
const Sdp4 = astroz.Sdp4;
const Sdp4Batch = astroz.Constellation.Sdp4Batch;
const Datetime = astroz.Datetime;
const constants = astroz.constants;

const allocator = shared.allocator;

pub const WGS72 = shared.WGS72;
pub const WGS84 = shared.WGS84;

pub const SatrecObject = extern struct {
    ob_base: c.PyObject,
    tle: ?*Tle,
    sgp4: ?*Sgp4,
    sdp4: ?*anyopaque, // points to Sdp4, anyopaque for extern struct compat
    whichconst: c_int, // gravity model used
    errorCode: c_int, // last propagation error code
    t: f64, // last propagation time (minutes)
    jdsatepoch: f64, // epoch JD integer part
    jdsatepochF: f64, // epoch JD fractional part
};

pub var SatrecType: c.PyTypeObject = undefined;

fn initType() void {
    SatrecType = std.mem.zeroes(c.PyTypeObject);
    SatrecType.tp_name = "astroz.Satrec";
    SatrecType.tp_doc = "SGP4 satellite record (python-sgp4 compatible).\n\nUse Satrec.twoline2rv(line1, line2) to create from TLE lines.";
    SatrecType.tp_basicsize = @sizeOf(SatrecObject);
    SatrecType.tp_flags = c.Py_TPFLAGS_DEFAULT | c.Py_TPFLAGS_BASETYPE;
    SatrecType.tp_new = c.PyType_GenericNew;
    SatrecType.tp_dealloc = @ptrCast(&satrec_dealloc);
    SatrecType.tp_methods = @constCast(&satrec_methods);
    SatrecType.tp_getset = @constCast(&satrec_getset);
}

fn getSdp4(self: *const SatrecObject) ?*const Sdp4 {
    const ptr = self.sdp4 orelse return null;
    return @ptrCast(@alignCast(ptr));
}

fn satrec_dealloc(self_obj: [*c]c.PyObject) callconv(.c) void {
    const self: *SatrecObject = @ptrCast(@alignCast(self_obj));
    if (self.sgp4) |s| allocator.destroy(s);
    if (getSdp4(self)) |sdp4| allocator.destroy(sdp4);
    if (self.tle) |t| {
        var tle_copy = t.*;
        tle_copy.deinit();
        allocator.destroy(t);
    }
    if (py.pyType(self_obj)) |tp| if (tp.*.tp_free) |free| free(@as(?*anyopaque, @ptrCast(self)));
}

fn satrec_twoline2rv(cls: [*c]c.PyObject, args: [*c]c.PyObject, kwds: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    _ = cls;
    var line1_ptr: [*c]const u8 = null;
    var line2_ptr: [*c]const u8 = null;
    var whichconst: c_int = WGS72; // python-sgp4's default

    const kwlist = [_:null]?[*:0]const u8{ "line1", "line2", "whichconst", null };
    if (c.PyArg_ParseTupleAndKeywords(args, kwds, "ss|i", @ptrCast(@constCast(&kwlist)), &line1_ptr, &line2_ptr, &whichconst) == 0)
        return null;

    const tle = Tle.parseLines(std.mem.span(line1_ptr), std.mem.span(line2_ptr), allocator) catch |err| {
        py.raiseValue(if (err == error.CatalogNumberMismatch) "TLE lines 1 and 2 have different catalog numbers" else "Failed to parse TLE lines");
        return null;
    };
    return satrecFromTle(tle, whichconst);
}

/// Builds a Satrec that takes ownership of `tle`.
fn satrecFromTle(tle: Tle, whichconst: c_int) [*c]c.PyObject {
    const self: *SatrecObject = @ptrCast(@alignCast(c.PyType_GenericAlloc(&SatrecType, 0) orelse {
        var tle_copy = tle;
        tle_copy.deinit();
        return null;
    }));
    self.whichconst = whichconst;

    const tle_ptr = allocator.create(Tle) catch {
        py.raiseRuntime("Out of memory");
        var tle_copy = tle;
        tle_copy.deinit();
        py.decref(@as(*c.PyObject, @ptrCast(self)));
        return null;
    };
    tle_ptr.* = tle;
    self.tle = tle_ptr;

    // python-sgp4 splits the epoch into a half-integer JD and a fraction
    self.jdsatepoch = @floor(tle.epochJd - 0.5) + 0.5;
    self.jdsatepochF = tle.epochJd - self.jdsatepoch;

    const grav = shared.getGravity(whichconst);
    const sgp4 = Sgp4.init(tle, grav) catch |e| switch (e) {
        error.DeepSpaceNotSupported => {
            const sdp4 = Sdp4.init(tle, grav) catch |e2| {
                self.errorCode = shared.sgp4ErrorCode(e2);
                return @ptrCast(self);
            };
            const sdp4_ptr = allocator.create(Sdp4) catch {
                py.raiseRuntime("Out of memory");
                py.decref(@as(*c.PyObject, @ptrCast(self)));
                return null;
            };
            sdp4_ptr.* = sdp4;
            self.sdp4 = @ptrCast(sdp4_ptr);
            return @ptrCast(self);
        },
        else => {
            self.errorCode = shared.sgp4ErrorCode(e);
            return @ptrCast(self);
        },
    };

    const sgp4_ptr = allocator.create(Sgp4) catch {
        py.raiseRuntime("Out of memory");
        py.decref(@as(*c.PyObject, @ptrCast(self)));
        return null;
    };
    sgp4_ptr.* = sgp4;
    self.sgp4 = sgp4_ptr;

    return @ptrCast(self);
}

fn satrec_from_omm(cls: [*c]c.PyObject, args: [*c]c.PyObject, kwds: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    _ = cls;
    var text_ptr: [*c]const u8 = null;
    var text_len: c.Py_ssize_t = 0;
    var whichconst: c_int = WGS72;

    const kwlist = [_:null]?[*:0]const u8{ "json", "whichconst", null };
    if (c.PyArg_ParseTupleAndKeywords(args, kwds, "s#|i", @ptrCast(@constCast(&kwlist)), &text_ptr, &text_len, &whichconst) == 0)
        return null;

    const tles = shared.parseOmm(text_ptr[0..@intCast(text_len)]) orelse return null;
    defer allocator.free(tles);

    const list = c.PyList_New(@intCast(tles.len)) orelse {
        for (tles) |*t| t.deinit();
        return null;
    };
    for (tles, 0..) |tle, i| {
        const sat = satrecFromTle(tle, whichconst) orelse {
            for (tles[i + 1 ..]) |*t| t.deinit();
            py.decref(list);
            return null;
        };
        _ = c.PyList_SetItem(list, @intCast(i), sat);
    }
    return list;
}

fn satrec_sgp4(self_obj: [*c]c.PyObject, args: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    const self: *SatrecObject = @ptrCast(@alignCast(self_obj));
    var jd: f64 = undefined;
    var fr: f64 = undefined;

    if (c.PyArg_ParseTuple(args, "dd", &jd, &fr) == 0) return null;

    const epochJd = self.jdsatepoch + self.jdsatepochF;
    const tsince = (jd + fr - epochJd) * constants.minutesPerDay;
    self.t = tsince;

    const zero: [2][3]f64 = .{ @splat(0), @splat(0) };
    const propagated = if (self.sgp4) |sgp4|
        sgp4.propagate(tsince)
    else if (getSdp4(self)) |sdp4|
        sdp4.propagate(tsince)
    else {
        self.errorCode = 6;
        return sgp4Result(self.errorCode, zero);
    };
    const pv = propagated catch |e| {
        self.errorCode = shared.sgp4ErrorCode(e);
        return sgp4Result(self.errorCode, zero);
    };
    self.errorCode = 0;
    return sgp4Result(0, pv);
}

/// (error, (x, y, z), (vx, vy, vz)), as python-sgp4 returns it
fn sgp4Result(code: c_int, pv: [2][3]f64) [*c]c.PyObject {
    const r = pv[0];
    const v = pv[1];
    return c.Py_BuildValue("i(ddd)(ddd)", code, r[0], r[1], r[2], v[0], v[1], v[2]);
}

const SgpSimdN = astroz.Constellation.Sgp4Batch.batchSize;

/// One satellite (SGP4 or SDP4) over many times: SIMD chunks through oma
/// dispatch, then the remainder one at a time.
fn propagateArray(comptime T: type, propagator: *const T, n: usize, epochJd: f64, jd_ptr: [*]const f64, fr_ptr: [*]const f64, pos_ptr: [*]f64, vel_ptr: [*]f64) void {
    const PVA = Sgp4.PosVelArray(SgpSimdN);
    var i: usize = 0;
    while (i + SgpSimdN <= n) : (i += SgpSimdN) {
        var times: [SgpSimdN]f64 = undefined;
        for (0..SgpSimdN) |j| {
            times[j] = ((jd_ptr[i + j] + fr_ptr[i + j]) - epochJd) * constants.minutesPerDay;
        }
        const maybe_results: ?PVA = if (T == Sgp4)
            astroz.dispatch.sgp4Times8(propagator, times) catch null
        else
            astroz.dispatch.sdp4Times8(propagator, times) catch null;

        if (maybe_results) |results| {
            for (0..SgpSimdN) |j| {
                const base = (i + j) * 3;
                pos_ptr[base..][0..3].* = results[j][0];
                vel_ptr[base..][0..3].* = results[j][1];
            }
        } else {
            // A failure at any of the 8 times fails the chunk; redo each time exactly
            for (0..SgpSimdN) |j| propagateOne(T, propagator, epochJd, jd_ptr, fr_ptr, pos_ptr, vel_ptr, i + j);
        }
    }
    while (i < n) : (i += 1) propagateOne(T, propagator, epochJd, jd_ptr, fr_ptr, pos_ptr, vel_ptr, i);
}

/// Scalar propagation of one time; NaN on failure (reported as e=6 in Python).
fn propagateOne(comptime T: type, propagator: *const T, epochJd: f64, jd_ptr: [*]const f64, fr_ptr: [*]const f64, pos_ptr: [*]f64, vel_ptr: [*]f64, i: usize) void {
    const tsince = ((jd_ptr[i] + fr_ptr[i]) - epochJd) * constants.minutesPerDay;
    const base = i * 3;
    if (propagator.propagate(tsince)) |result| {
        pos_ptr[base..][0..3].* = result[0];
        vel_ptr[base..][0..3].* = result[1];
    } else |_| {
        pos_ptr[base..][0..3].* = @splat(std.math.nan(f64));
        vel_ptr[base..][0..3].* = @splat(std.math.nan(f64));
    }
}

/// The jd, fr, positions and velocities arrays shared by the *_into functions.
const TimeBuffers = struct {
    jd: []f64,
    fr: []f64,
    pos: []f64,
    vel: []f64,
    bufs: [4]c.Py_buffer,

    /// Raises unless all four are float64, fr matches jd, and positions/velocities
    /// hold at least `columns` satellites' worth of states for every time.
    fn get(jd_obj: [*c]c.PyObject, fr_obj: [*c]c.PyObject, pos_obj: [*c]c.PyObject, vel_obj: [*c]c.PyObject, columns: usize) ?TimeBuffers {
        var t: TimeBuffers = undefined;
        t.bufs = @splat(std.mem.zeroes(c.Py_buffer));
        const objs = [_][*c]c.PyObject{ jd_obj, fr_obj, pos_obj, vel_obj };
        var slices: [4][]f64 = undefined;
        const names = .{ "jd", "fr", "positions", "velocities" };
        inline for (objs, &slices, names, 0..) |obj, *slice, name, k| {
            // positions and velocities are written to
            slice.* = shared.getBuffer(f64, name, obj, &t.bufs[k], k >= 2) orelse {
                for (t.bufs[0..k]) |*b| c.PyBuffer_Release(b);
                return null;
            };
        }
        t.jd, t.fr, t.pos, t.vel = slices;
        const need = t.jd.len * columns * 3;
        if (t.fr.len != t.jd.len or t.pos.len < need or t.vel.len < need) {
            t.release();
            py.raiseValue(if (t.fr.len != t.jd.len) "jd and fr must have the same length" else "positions/velocities array too small");
            return null;
        }
        return t;
    }

    fn release(t: *TimeBuffers) void {
        for (&t.bufs) |*b| c.PyBuffer_Release(b);
    }
};

fn satrec_sgp4_array_into(self_obj: [*c]c.PyObject, args: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    const self: *SatrecObject = @ptrCast(@alignCast(self_obj));
    var jd_obj: [*c]c.PyObject = null;
    var fr_obj: [*c]c.PyObject = null;
    var pos_obj: [*c]c.PyObject = null;
    var vel_obj: [*c]c.PyObject = null;

    if (c.PyArg_ParseTuple(args, "OOOO", &jd_obj, &fr_obj, &pos_obj, &vel_obj) == 0) return null;

    var bufs = TimeBuffers.get(jd_obj, fr_obj, pos_obj, vel_obj, 1) orelse return null;
    defer bufs.release();

    const n = bufs.jd.len;
    const epochJd = self.jdsatepoch + self.jdsatepochF;
    if (self.sgp4) |sgp4| {
        propagateArray(Sgp4, sgp4, n, epochJd, bufs.jd.ptr, bufs.fr.ptr, bufs.pos.ptr, bufs.vel.ptr);
    } else if (getSdp4(self)) |sdp4| {
        propagateArraySdp4(sdp4, n, epochJd, bufs.jd.ptr, bufs.fr.ptr, bufs.pos.ptr, bufs.vel.ptr);
    } else {
        py.raiseRuntime("Satrec not initialized");
        return null;
    }

    return py.none();
}

const satrec_methods = [_]c.PyMethodDef{
    .{
        .ml_name = "twoline2rv",
        .ml_meth = @ptrCast(&satrec_twoline2rv),
        .ml_flags = c.METH_CLASS | c.METH_VARARGS | c.METH_KEYWORDS,
        .ml_doc = "twoline2rv(line1, line2, whichconst=WGS72) -> Satrec\n\nCreate Satrec from TLE lines.",
    },
    .{
        .ml_name = "from_omm",
        .ml_meth = @ptrCast(&satrec_from_omm),
        .ml_flags = c.METH_CLASS | c.METH_VARARGS | c.METH_KEYWORDS,
        .ml_doc = "from_omm(json, whichconst=WGS72) -> list[Satrec]\n\nOne Satrec per record in CCSDS OMM JSON (a single object or an array).",
    },
    .{
        .ml_name = "sgp4",
        .ml_meth = @ptrCast(&satrec_sgp4),
        .ml_flags = c.METH_VARARGS,
        .ml_doc = "sgp4(jd, fr) -> (error, (x,y,z), (vx,vy,vz))\n\nPropagate to given Julian date.",
    },
    .{
        .ml_name = "sgp4_array_into",
        .ml_meth = @ptrCast(&satrec_sgp4_array_into),
        .ml_flags = c.METH_VARARGS,
        .ml_doc = "sgp4_array_into(jd, fr, positions, velocities) -> None\n\nSIMD batch propagation into pre-allocated arrays.",
    },
    .{ .ml_name = null, .ml_meth = null, .ml_flags = 0, .ml_doc = null },
};

fn makeGetter(comptime getter: fn (*const SatrecObject) [*c]c.PyObject) fn ([*c]c.PyObject, ?*anyopaque) callconv(.c) [*c]c.PyObject {
    return struct {
        fn get(self_obj: [*c]c.PyObject, _: ?*anyopaque) callconv(.c) [*c]c.PyObject {
            const self: *SatrecObject = @ptrCast(@alignCast(self_obj));
            return getter(self);
        }
    }.get;
}

fn makeTleGetter(comptime getter: fn (*const Tle) [*c]c.PyObject) fn ([*c]c.PyObject, ?*anyopaque) callconv(.c) [*c]c.PyObject {
    return struct {
        fn get(self_obj: [*c]c.PyObject, _: ?*anyopaque) callconv(.c) [*c]c.PyObject {
            const self: *SatrecObject = @ptrCast(@alignCast(self_obj));
            const tle = self.tle orelse {
                py.raiseRuntime("Satrec not initialized");
                return null;
            };
            return getter(tle);
        }
    }.get;
}

fn getSatnum(tle: *const Tle) [*c]c.PyObject {
    return py.int(@intCast(tle.satelliteNumber));
}
fn getSatnumStr(tle: *const Tle) [*c]c.PyObject {
    var buf: [10]u8 = undefined;
    const text = Tle.satelliteNumberString(tle.satelliteNumber, &buf);
    return c.PyUnicode_FromStringAndSize(text.ptr, @intCast(text.len));
}
fn getIntldesg(tle: *const Tle) [*c]c.PyObject {
    // python-sgp4 writes it as in the TLE ("98067A"); Tle holds "98-067A" (TLE) or "1998-067A" (OMM),
    // so keep the two year digits before the dash and everything after it
    const d = tle.intlDesignator;
    var buf: [16]u8 = undefined;
    const dash = std.mem.indexOfScalar(u8, d, '-') orelse d.len;
    const text = if (dash >= 2 and dash < d.len)
        std.fmt.bufPrint(&buf, "{s}{s}", .{ d[dash - 2 .. dash], d[dash + 1 ..] }) catch d
    else
        d;
    return c.PyUnicode_FromStringAndSize(text.ptr, @intCast(text.len));
}
fn getClassification(tle: *const Tle) [*c]c.PyObject {
    return c.PyUnicode_FromStringAndSize(&tle.classification, 1);
}
fn getElnum(tle: *const Tle) [*c]c.PyObject {
    return py.int(@intCast(tle.elemNumber));
}
fn getRevnum(tle: *const Tle) [*c]c.PyObject {
    return py.int(@intCast(tle.revNum));
}
fn getEphtype(tle: *const Tle) [*c]c.PyObject {
    return py.int(@intCast(tle.ephemType));
}
fn getNddot(tle: *const Tle) [*c]c.PyObject {
    // rev/day^3 -> rad/min^3, keeping the TLE's divide-by-6
    return py.float(tle.secondDerMeanMotion * constants.twoPi / (constants.minutesPerDay * constants.minutesPerDay * constants.minutesPerDay));
}
fn getEpochyr(tle: *const Tle) [*c]c.PyObject {
    return py.int(@intCast(tle.epochYear));
}
fn getEpochdays(tle: *const Tle) [*c]c.PyObject {
    return py.float(tle.epochDay);
}
fn getEcco(tle: *const Tle) [*c]c.PyObject {
    return py.float(tle.eccentricity);
}
fn getInclo(tle: *const Tle) [*c]c.PyObject {
    return py.float(tle.inclination * constants.deg2rad);
}
fn getNodeo(tle: *const Tle) [*c]c.PyObject {
    return py.float(tle.rightAscension * constants.deg2rad);
}
fn getArgpo(tle: *const Tle) [*c]c.PyObject {
    return py.float(tle.perigee * constants.deg2rad);
}
fn getMo(tle: *const Tle) [*c]c.PyObject {
    return py.float(tle.mAnomaly * constants.deg2rad);
}
fn getNoKozai(tle: *const Tle) [*c]c.PyObject {
    // rev/day -> rad/min
    return py.float(tle.mMotion * constants.twoPi / constants.minutesPerDay);
}
fn getBstar(tle: *const Tle) [*c]c.PyObject {
    return py.float(tle.bstarDrag);
}
fn getNdot(tle: *const Tle) [*c]c.PyObject {
    // rev/day^2 -> rad/min^2, keeping the TLE's divide-by-2
    return py.float(tle.firstDerMeanMotion * constants.twoPi / (constants.minutesPerDay * constants.minutesPerDay));
}

fn getJdsatepoch(self: *const SatrecObject) [*c]c.PyObject {
    return py.float(self.jdsatepoch);
}
fn getJdsatepochF(self: *const SatrecObject) [*c]c.PyObject {
    return py.float(self.jdsatepochF);
}
fn getError(self: *const SatrecObject) [*c]c.PyObject {
    return py.int(@intCast(self.errorCode));
}
fn getT(self: *const SatrecObject) [*c]c.PyObject {
    return py.float(self.t);
}

fn getSgp4Elements(self: *const SatrecObject) ?*const Sgp4.Elements {
    if (self.sgp4) |sgp4| return &sgp4.elements;
    if (getSdp4(self)) |sdp4| return &sdp4.elements.sgp4;
    return null;
}

fn makeElementsGetter(comptime getter: fn (*const Sgp4.Elements) f64) fn ([*c]c.PyObject, ?*anyopaque) callconv(.c) [*c]c.PyObject {
    return struct {
        fn get(self_obj: [*c]c.PyObject, _: ?*anyopaque) callconv(.c) [*c]c.PyObject {
            const self: *SatrecObject = @ptrCast(@alignCast(self_obj));
            const el = getSgp4Elements(self) orelse {
                py.raiseRuntime("Satrec not initialized");
                return null;
            };
            return py.float(getter(el));
        }
    }.get;
}

fn getA(el: *const Sgp4.Elements) f64 {
    return el.a;
}
fn getAlta(el: *const Sgp4.Elements) f64 {
    return el.a * (1.0 + el.ecco) - 1.0;
}
fn getAltp(el: *const Sgp4.Elements) f64 {
    return el.a * (1.0 - el.ecco) - 1.0;
}
fn getIsDeepSpace(self: *const SatrecObject) [*c]c.PyObject {
    return c.PyBool_FromLong(if (self.sdp4 != null) @as(c_long, 1) else @as(c_long, 0));
}

const satrec_getset = [_]c.PyGetSetDef{
    .{ .name = "satnum", .get = @ptrCast(&makeTleGetter(getSatnum)), .set = null, .doc = "NORAD catalog number", .closure = null },
    .{ .name = "satnum_str", .get = @ptrCast(&makeTleGetter(getSatnumStr)), .set = null, .doc = "Catalog number as written in a TLE (Alpha-5 above 99999)", .closure = null },
    .{ .name = "intldesg", .get = @ptrCast(&makeTleGetter(getIntldesg)), .set = null, .doc = "International designator as written in a TLE (e.g. 98067A)", .closure = null },
    .{ .name = "classification", .get = @ptrCast(&makeTleGetter(getClassification)), .set = null, .doc = "Classification (U, C or S)", .closure = null },
    .{ .name = "elnum", .get = @ptrCast(&makeTleGetter(getElnum)), .set = null, .doc = "Element set number", .closure = null },
    .{ .name = "revnum", .get = @ptrCast(&makeTleGetter(getRevnum)), .set = null, .doc = "Revolution number at epoch", .closure = null },
    .{ .name = "ephtype", .get = @ptrCast(&makeTleGetter(getEphtype)), .set = null, .doc = "Ephemeris type", .closure = null },
    .{ .name = "nddot", .get = @ptrCast(&makeTleGetter(getNddot)), .set = null, .doc = "Second derivative of mean motion (rad/min^3)", .closure = null },
    .{ .name = "epochyr", .get = @ptrCast(&makeTleGetter(getEpochyr)), .set = null, .doc = "Epoch year (2-digit)", .closure = null },
    .{ .name = "epochdays", .get = @ptrCast(&makeTleGetter(getEpochdays)), .set = null, .doc = "Epoch day of year (fractional)", .closure = null },
    .{ .name = "ecco", .get = @ptrCast(&makeTleGetter(getEcco)), .set = null, .doc = "Eccentricity", .closure = null },
    .{ .name = "inclo", .get = @ptrCast(&makeTleGetter(getInclo)), .set = null, .doc = "Inclination (radians)", .closure = null },
    .{ .name = "nodeo", .get = @ptrCast(&makeTleGetter(getNodeo)), .set = null, .doc = "Right ascension of ascending node (radians)", .closure = null },
    .{ .name = "argpo", .get = @ptrCast(&makeTleGetter(getArgpo)), .set = null, .doc = "Argument of perigee (radians)", .closure = null },
    .{ .name = "mo", .get = @ptrCast(&makeTleGetter(getMo)), .set = null, .doc = "Mean anomaly (radians)", .closure = null },
    .{ .name = "no_kozai", .get = @ptrCast(&makeTleGetter(getNoKozai)), .set = null, .doc = "Mean motion (radians/minute, Kozai)", .closure = null },
    .{ .name = "bstar", .get = @ptrCast(&makeTleGetter(getBstar)), .set = null, .doc = "B* drag term", .closure = null },
    .{ .name = "ndot", .get = @ptrCast(&makeTleGetter(getNdot)), .set = null, .doc = "First derivative of mean motion (rad/min^2)", .closure = null },
    .{ .name = "jdsatepoch", .get = @ptrCast(&makeGetter(getJdsatepoch)), .set = null, .doc = "Epoch Julian date (integer part)", .closure = null },
    .{ .name = "jdsatepochF", .get = @ptrCast(&makeGetter(getJdsatepochF)), .set = null, .doc = "Epoch Julian date (fractional part)", .closure = null },
    .{ .name = "error", .get = @ptrCast(&makeGetter(getError)), .set = null, .doc = "Last error code (0 = success)", .closure = null },
    .{ .name = "t", .get = @ptrCast(&makeGetter(getT)), .set = null, .doc = "Last propagation time (minutes from epoch)", .closure = null },
    .{ .name = "a", .get = @ptrCast(&makeElementsGetter(getA)), .set = null, .doc = "Semi-major axis (Earth radii)", .closure = null },
    .{ .name = "alta", .get = @ptrCast(&makeElementsGetter(getAlta)), .set = null, .doc = "Apoapsis altitude (Earth radii)", .closure = null },
    .{ .name = "altp", .get = @ptrCast(&makeElementsGetter(getAltp)), .set = null, .doc = "Periapsis altitude (Earth radii)", .closure = null },
    .{ .name = "is_deep_space", .get = @ptrCast(&makeGetter(getIsDeepSpace)), .set = null, .doc = "True if using SDP4 deep-space propagator", .closure = null },
    .{ .name = null, .get = null, .set = null, .doc = null, .closure = null },
};

/// sdp4_batch_propagate_into(satrecs, jd, fr, positions, velocities, output_stride=-1, sat_offset=0)
/// Output is satellite-major (n_sdp4, n_times, 3) by default. With output_stride > 0
/// it's time-major, with satellite s at column sat_offset + s, so the caller can
/// write SDP4 results straight into a mixed constellation's array.
pub fn pySdp4BatchPropagateInto(_: [*c]c.PyObject, args: [*c]c.PyObject, kwds: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    var satrecs_obj: [*c]c.PyObject = null;
    var jd_obj: [*c]c.PyObject = null;
    var fr_obj: [*c]c.PyObject = null;
    var pos_obj: [*c]c.PyObject = null;
    var vel_obj: [*c]c.PyObject = null;
    var output_stride_arg: c_int = -1;
    var sat_offset_arg: c_int = 0;

    const kwlist = [_:null]?[*:0]const u8{ "satrecs", "jd", "fr", "positions", "velocities", "output_stride", "sat_offset", null };
    if (c.PyArg_ParseTupleAndKeywords(args, kwds, "OOOOO|ii", @ptrCast(@constCast(&kwlist)), &satrecs_obj, &jd_obj, &fr_obj, &pos_obj, &vel_obj, &output_stride_arg, &sat_offset_arg) == 0) return null;

    if (c.PySequence_Check(satrecs_obj) == 0) {
        py.raiseType("First argument must be a sequence of Satrec objects");
        return null;
    }

    const num_sats: usize = @intCast(c.PySequence_Size(satrecs_obj));
    if (num_sats == 0) return py.none();

    if (output_stride_arg > 0 and @as(usize, @intCast(sat_offset_arg)) + num_sats > @as(usize, @intCast(output_stride_arg))) {
        py.raiseValue("sat_offset + len(satrecs) must not exceed output_stride");
        return null;
    }
    if (sat_offset_arg < 0 or (output_stride_arg <= 0 and sat_offset_arg != 0)) {
        py.raiseValue("sat_offset needs output_stride and must be >= 0");
        return null;
    }
    const columns: usize = if (output_stride_arg > 0) @intCast(output_stride_arg) else num_sats;
    var bufs = TimeBuffers.get(jd_obj, fr_obj, pos_obj, vel_obj, columns) orelse return null;
    defer bufs.release();
    const n_times = bufs.jd.len;

    const sdp4_ptrs = allocator.alloc(*const Sdp4, num_sats) catch {
        py.raiseRuntime("Out of memory");
        return null;
    };
    defer allocator.free(sdp4_ptrs);

    const epochs = allocator.alloc(f64, num_sats) catch {
        py.raiseRuntime("Out of memory");
        return null;
    };
    defer allocator.free(epochs);

    for (0..num_sats) |i| {
        const item = c.PySequence_GetItem(satrecs_obj, @intCast(i)) orelse {
            py.raiseValue("Failed to get Satrec from sequence");
            return null;
        };
        defer py.decref(item);

        if (c.PyObject_TypeCheck(item, &SatrecType) == 0) {
            py.raiseType("All items must be Satrec objects");
            return null;
        }

        const satrec: *SatrecObject = @ptrCast(@alignCast(item));
        sdp4_ptrs[i] = getSdp4(satrec) orelse {
            py.raiseValue("Satrec at index is not a deep-space (SDP4) object");
            return null;
        };
        epochs[i] = satrec.jdsatepoch + satrec.jdsatepochF;
    }

    const time_major = output_stride_arg > 0;
    const stride: usize = if (time_major) @intCast(output_stride_arg) else num_sats;
    const offset: usize = @intCast(sat_offset_arg);
    const layout: astroz.Constellation.Layout = if (time_major) .timeMajor else .satelliteMajor;
    sdp4BatchPropagate(sdp4_ptrs, epochs, n_times, stride, offset, layout, bufs.jd.ptr, bufs.fr.ptr, bufs.pos.ptr, bufs.vel.ptr) catch |e| {
        py.raiseRuntime(switch (e) {
            error.OutputTooSmall => "Output arrays are too small",
            else => |err| shared.sgp4ErrorMsg(err),
        });
        return null;
    };

    return py.none();
}

const SdpSimdN = Sdp4Batch.batchSize;

fn sdp4BatchPropagate(
    sdp4_ptrs: []const *const Sdp4,
    epochs: []const f64,
    n_times: usize,
    output_stride: usize,
    sat_offset: usize,
    layout: astroz.Constellation.Layout,
    jd_ptr: [*]const f64,
    fr_ptr: [*]const f64,
    pos_ptr: [*]f64,
    vel_ptr: [*]f64,
) astroz.Constellation.Error!void {
    const num_sats = sdp4_ptrs.len;
    const num_simd_batches = @divCeil(num_sats, SdpSimdN);

    const batch_els = try allocator.alloc(Sdp4Batch.Sdp4BatchElements(SdpSimdN), num_simd_batches);
    defer allocator.free(batch_els);
    const batch_epochs = try allocator.alloc([SdpSimdN]f64, num_simd_batches);
    defer allocator.free(batch_epochs);

    const grav = sdp4_ptrs[0].elements.sgp4.grav;
    for (0..num_simd_batches) |b| {
        var els: [SdpSimdN]Sdp4.Elements = undefined;
        var ep: [SdpSimdN]f64 = undefined;
        for (0..SdpSimdN) |lane| {
            const idx = b * SdpSimdN + lane;
            const safe_idx = if (idx < num_sats) idx else num_sats - 1;
            els[lane] = sdp4_ptrs[safe_idx].elements;
            ep[lane] = epochs[safe_idx];
        }
        batch_els[b] = Sdp4Batch.initFromElements(SdpSimdN, els, grav);
        batch_epochs[b] = ep;
    }

    const padded = num_simd_batches * SdpSimdN;
    const origIndices = try allocator.alloc(u32, padded);
    defer allocator.free(origIndices);
    for (origIndices, 0..) |*v, i| v.* = @intCast(sat_offset + i);

    const out_size = output_stride * n_times * 3;
    try astroz.Constellation.propagateSdp4Constellation(
        batch_els,
        batch_epochs,
        num_sats,
        output_stride,
        origIndices,
        jd_ptr[0..n_times],
        fr_ptr[0..n_times],
        pos_ptr[0..out_size],
        vel_ptr[0..out_size],
        .teme,
        layout,
    );
}

/// Scalar SDP4 over many times. Carrying the resonance state between calls
/// avoids re-integrating from t=0 for every time.
fn propagateArraySdp4(
    propagator: *const Sdp4,
    n_times: usize,
    epochJd: f64,
    jd_ptr: [*]const f64,
    fr_ptr: [*]const f64,
    pos_ptr: [*]f64,
    vel_ptr: [*]f64,
) void {
    var carry = Sdp4.ResonanceCarry{
        .atime = 0.0,
        .xli = propagator.elements.xlamo,
        .xni = propagator.elements.sgp4.noUnkozai,
    };
    for (0..n_times) |i| {
        const base = i * 3;
        const tsince = ((jd_ptr[i] + fr_ptr[i]) - epochJd) * constants.minutesPerDay;
        if (propagator.propagateCarry(tsince, &carry)) |result| {
            pos_ptr[base..][0..3].* = result[0];
            vel_ptr[base..][0..3].* = result[1];
        } else |_| {
            pos_ptr[base..][0..3].* = @splat(std.math.nan(f64));
            vel_ptr[base..][0..3].* = @splat(std.math.nan(f64));
        }
    }
}

/// jday(year, month, day, hour, minute, second) -> (jd, fr)
pub fn pyJday(_: [*c]c.PyObject, args: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    var year: c_int = undefined;
    var month: c_int = undefined;
    var day: c_int = undefined;
    var hour: c_int = undefined;
    var minute: c_int = undefined;
    var second: f64 = undefined;

    if (c.PyArg_ParseTuple(args, "iiiiid", &year, &month, &day, &hour, &minute, &second) == 0)
        return null;

    const result = Datetime.jday(
        @intCast(year),
        @intCast(month),
        @intCast(day),
        @intCast(hour),
        @intCast(minute),
        second,
    );

    return c.Py_BuildValue("(dd)", result.jd, result.fr);
}

/// days2mdhms(year, days) -> (month, day, hour, minute, second)
pub fn pyDays2mdhms(_: [*c]c.PyObject, args: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    var year: c_int = undefined;
    var days: f64 = undefined;

    if (c.PyArg_ParseTuple(args, "id", &year, &days) == 0)
        return null;

    const result = Datetime.days2mdhms(@intCast(year), days);

    return c.Py_BuildValue(
        "(iiiid)",
        @as(c_int, result.month),
        @as(c_int, result.day),
        @as(c_int, result.hour),
        @as(c_int, result.minute),
        result.second,
    );
}

// SatrecArray holds SGP4 satellites only; the Python wrapper sends deep-space
// satellites through sdp4_batch_propagate_into instead.

pub const SatrecArrayObject = shared.BatchObjectFields;

pub var SatrecArrayType: c.PyTypeObject = undefined;

fn satrecarray_get_epochs(self_obj: [*c]c.PyObject, _: ?*anyopaque) callconv(.c) [*c]c.PyObject {
    const self: *SatrecArrayObject = @ptrCast(@alignCast(self_obj));
    const epoch_jds = self.epoch_jds orelse {
        py.raiseRuntime("SatrecArray not initialized");
        return null;
    };
    return py.listFromF64(epoch_jds[0..self.num_satellites]);
}

fn satrecarray_get_num_satellites(self_obj: [*c]c.PyObject, _: ?*anyopaque) callconv(.c) [*c]c.PyObject {
    const self: *SatrecArrayObject = @ptrCast(@alignCast(self_obj));
    return py.int(@intCast(self.num_satellites));
}

const satrecarray_getset = [_]c.PyGetSetDef{
    .{ .name = "epochs", .get = @ptrCast(&satrecarray_get_epochs), .set = null, .doc = "List of epoch Julian dates for each satellite", .closure = null },
    .{ .name = "num_satellites", .get = @ptrCast(&satrecarray_get_num_satellites), .set = null, .doc = "Number of satellites in the array", .closure = null },
    .{ .name = null, .get = null, .set = null, .doc = null, .closure = null },
};

fn initSatrecArrayType() void {
    SatrecArrayType = std.mem.zeroes(c.PyTypeObject);
    SatrecArrayType.tp_name = "astroz.SatrecArray";
    SatrecArrayType.tp_doc = "Batch SGP4 propagator for multiple satellites (python-sgp4 compatible).\n\nUse SatrecArray(satrecs) with a list of Satrec objects.";
    SatrecArrayType.tp_basicsize = @sizeOf(SatrecArrayObject);
    SatrecArrayType.tp_flags = c.Py_TPFLAGS_DEFAULT | c.Py_TPFLAGS_BASETYPE;
    SatrecArrayType.tp_new = c.PyType_GenericNew;
    SatrecArrayType.tp_init = @ptrCast(&satrecarray_init);
    SatrecArrayType.tp_dealloc = @ptrCast(&satrecarray_dealloc);
    SatrecArrayType.tp_methods = @constCast(&satrecarray_methods);
    SatrecArrayType.tp_getset = @constCast(&satrecarray_getset);
}

fn satrecarray_init(self_obj: [*c]c.PyObject, args: [*c]c.PyObject, _: [*c]c.PyObject) callconv(.c) c_int {
    const self: *SatrecArrayObject = @ptrCast(@alignCast(self_obj));
    var satrecs_obj: [*c]c.PyObject = null;

    if (c.PyArg_ParseTuple(args, "O", &satrecs_obj) == 0)
        return -1;

    if (c.PySequence_Check(satrecs_obj) == 0) {
        py.raiseType("Argument must be a sequence of Satrec objects");
        return -1;
    }

    const num_sats: usize = @intCast(c.PySequence_Size(satrecs_obj));
    if (num_sats == 0) {
        py.raiseValue("Must provide at least one Satrec object");
        return -1;
    }

    shared.freeBatchFields(self);

    const tles = allocator.alloc(Tle, num_sats) catch {
        py.raiseRuntime("Out of memory");
        return -1;
    };
    defer allocator.free(tles);

    var grav_model: c_int = WGS72;

    for (0..num_sats) |i| {
        const item = c.PySequence_GetItem(satrecs_obj, @intCast(i)) orelse {
            py.raiseValue("Failed to get Satrec from sequence");
            return -1;
        };
        defer py.decref(item);

        if (c.PyObject_TypeCheck(item, &SatrecType) == 0) {
            py.raiseType("All items must be Satrec objects");
            return -1;
        }

        const satrec = @as(*SatrecObject, @ptrCast(@alignCast(item)));
        tles[i] = (satrec.tle orelse {
            py.raiseValue("Satrec TLE not initialized");
            return -1;
        }).*;

        // The whole batch uses the first satellite's gravity model
        if (i == 0) grav_model = satrec.whichconst;
    }

    const result = shared.buildBatches(tles, shared.getGravity(grav_model)) orelse return -1;
    shared.storeBatchResult(self, result, num_sats);
    return 0;
}

fn satrecarray_dealloc(self_obj: [*c]c.PyObject) callconv(.c) void {
    const self: *SatrecArrayObject = @ptrCast(@alignCast(self_obj));
    shared.freeBatchFields(self);
    if (py.pyType(self_obj)) |tp| if (tp.*.tp_free) |free| free(@as(?*anyopaque, @ptrCast(self)));
}

fn satrecarray_propagate_into(self_obj: [*c]c.PyObject, args: [*c]c.PyObject, kwds: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    const self: *SatrecArrayObject = @ptrCast(@alignCast(self_obj));
    var times_obj: [*c]c.PyObject = null;
    var pos_obj: [*c]c.PyObject = null;
    var vel_obj: [*c]c.PyObject = null;
    var offsets_obj: [*c]c.PyObject = null;

    const kwlist = [_:null]?[*:0]const u8{ "times", "positions", "velocities", "epoch_offsets", null };
    if (c.PyArg_ParseTupleAndKeywords(args, kwds, "OO|OO", @ptrCast(@constCast(&kwlist)), &times_obj, &pos_obj, &vel_obj, &offsets_obj) == 0) return null;

    const batches = (self.batches orelse {
        py.raiseRuntime("SatrecArray not initialized");
        return null;
    })[0..self.num_batches];

    return shared.propagateInto(self, batches, times_obj, pos_obj, vel_obj, offsets_obj, null, .{});
}

const satrecarray_methods = [_]c.PyMethodDef{
    .{
        .ml_name = "propagate_into",
        .ml_meth = @ptrCast(&satrecarray_propagate_into),
        .ml_flags = c.METH_VARARGS | c.METH_KEYWORDS,
        .ml_doc = "propagate_into(times, positions, velocities=None, epoch_offsets=None) -> None\n\nHigh-performance SIMD batch propagation.",
    },
    .{ .ml_name = null, .ml_meth = null, .ml_flags = 0, .ml_doc = null },
};

pub fn ready() c_int {
    initType();
    initSatrecArrayType();
    if (c.PyType_Ready(&SatrecType) < 0) return -1;
    return c.PyType_Ready(&SatrecArrayType);
}

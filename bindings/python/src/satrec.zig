//! Python Satrec type - python-sgp4 compatible API
//! Provides drop-in replacement for sgp4.api.Satrec

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

// Re-export gravity constants for main.zig
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
    SatrecType.tp_new = @ptrCast(&satrec_new);
    SatrecType.tp_dealloc = @ptrCast(&satrec_dealloc);
    SatrecType.tp_methods = @constCast(&satrec_methods);
    SatrecType.tp_getset = @constCast(&satrec_getset);
}

fn satrec_new(typ: [*c]c.PyTypeObject, _: [*c]c.PyObject, _: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    const self: *SatrecObject = @ptrCast(@alignCast(c.PyType_GenericAlloc(typ, 0) orelse return null));
    self.tle = null;
    self.sgp4 = null;
    self.sdp4 = null;
    self.whichconst = WGS84;
    self.errorCode = 0;
    self.t = 0.0;
    self.jdsatepoch = 0.0;
    self.jdsatepochF = 0.0;
    return @ptrCast(self);
}

/// Cast the opaque SDP4 pointer back to a typed pointer.
fn getSdp4(self: *const SatrecObject) ?*const Sdp4 {
    const ptr = self.sdp4 orelse return null;
    return @ptrCast(@alignCast(ptr));
}

fn satrec_dealloc(self_obj: [*c]c.PyObject) callconv(.c) void {
    const self: *SatrecObject = @ptrCast(@alignCast(self_obj));
    if (self.sgp4) |s| allocator.destroy(s);
    if (getSdp4(self)) |sdp4| {
        allocator.destroy(sdp4);
    }
    if (self.tle) |t| {
        var tle_copy = t.*;
        tle_copy.deinit();
        allocator.destroy(t);
    }
    if (py.pyType(self_obj)) |tp| if (tp.*.tp_free) |free| free(@as(?*anyopaque, @ptrCast(self)));
}

/// Class method: twoline2rv(line1, line2, whichconst=WGS72) -> Satrec
fn satrec_twoline2rv(cls: [*c]c.PyObject, args: [*c]c.PyObject, kwds: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    _ = cls;
    var line1_ptr: [*c]const u8 = null;
    var line2_ptr: [*c]const u8 = null;
    var whichconst: c_int = WGS72; // default to WGS72 like python-sgp4

    const kwlist = [_:null]?[*:0]const u8{ "line1", "line2", "whichconst", null };
    if (c.PyArg_ParseTupleAndKeywords(args, kwds, "ss|i", @ptrCast(@constCast(&kwlist)), &line1_ptr, &line2_ptr, &whichconst) == 0)
        return null;

    if (line1_ptr == null or line2_ptr == null) {
        py.raiseValue("line1 and line2 are required");
        return null;
    }

    const line1 = std.mem.span(line1_ptr);
    const line2 = std.mem.span(line2_ptr);

    // Create new SatrecObject
    const self: *SatrecObject = @ptrCast(@alignCast(
        c.PyType_GenericAlloc(&SatrecType, 0) orelse return null,
    ));
    self.tle = null;
    self.sgp4 = null;
    self.sdp4 = null;
    self.whichconst = whichconst;
    self.errorCode = 0;
    self.t = 0.0;

    // Parse TLE
    const tle = Tle.parseLines(line1, line2, allocator) catch {
        py.raiseValue("Failed to parse TLE lines");
        c.Py_DECREF(@as([*c]c.PyObject, @ptrCast(self)));
        return null;
    };

    const tle_ptr = allocator.create(Tle) catch {
        py.raiseRuntime("Out of memory");
        var tle_copy = tle;
        tle_copy.deinit();
        c.Py_DECREF(@as([*c]c.PyObject, @ptrCast(self)));
        return null;
    };
    tle_ptr.* = tle;
    self.tle = tle_ptr;

    // Compute epoch JD
    const fullYear: u16 = if (tle.firstLine.epochYear < 57)
        2000 + tle.firstLine.epochYear
    else
        1900 + tle.firstLine.epochYear;
    const epochJd = Datetime.yearDoyToJulianDate(fullYear, tle.firstLine.epochDay);
    // Split into integer (at noon) and fractional parts like python-sgp4
    self.jdsatepoch = @floor(epochJd - 0.5) + 0.5;
    self.jdsatepochF = epochJd - self.jdsatepoch;

    // Initialize SGP4, fallback to SDP4 for deep-space orbits
    const grav = shared.getGravity(whichconst);
    const sgp4 = Sgp4.init(tle, grav) catch |e| switch (e) {
        error.DeepSpaceNotSupported => {
            // Deep-space orbit — use SDP4
            const sdp4 = Sdp4.init(tle, grav) catch |e2| {
                self.errorCode = shared.sgp4ErrorCode(e2);
                return @ptrCast(self);
            };
            const sdp4_ptr = allocator.create(Sdp4) catch {
                py.raiseRuntime("Out of memory");
                c.Py_DECREF(@as([*c]c.PyObject, @ptrCast(self)));
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
        c.Py_DECREF(@as([*c]c.PyObject, @ptrCast(self)));
        return null;
    };
    sgp4_ptr.* = sgp4;
    self.sgp4 = sgp4_ptr;

    return @ptrCast(self);
}

/// Instance method: sgp4(jd, fr) -> (error, (x,y,z), (vx,vy,vz))
fn satrec_sgp4(self_obj: [*c]c.PyObject, args: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    const self: *SatrecObject = @ptrCast(@alignCast(self_obj));
    var jd: f64 = undefined;
    var fr: f64 = undefined;

    if (c.PyArg_ParseTuple(args, "dd", &jd, &fr) == 0) return null;

    // Convert JD to minutes since epoch
    const epochJd = self.jdsatepoch + self.jdsatepochF;
    const totalJd = jd + fr;
    const tsince = (totalJd - epochJd) * constants.minutesPerDay;

    self.t = tsince;

    if (self.sgp4) |sgp4| {
        const result = sgp4.propagate(tsince) catch |e| {
            self.errorCode = shared.sgp4ErrorCode(e);
            return buildErrorResult(self.errorCode);
        };
        self.errorCode = 0;
        return buildSuccessResult(result);
    } else if (getSdp4(self)) |sdp4| {
        const result = sdp4.propagate(tsince) catch |e| {
            self.errorCode = shared.sgp4ErrorCode(e);
            return buildErrorResult(self.errorCode);
        };
        self.errorCode = 0;
        return buildSuccessResult(result);
    } else {
        self.errorCode = 6;
        return buildErrorResult(self.errorCode);
    }
}

fn buildErrorResult(errorCode: c_int) [*c]c.PyObject {
    // Return (error_code, (0,0,0), (0,0,0))
    const err = py.int(@intCast(errorCode)) orelse return null;
    const zero_pos = py.vec3Tuple(.{ 0.0, 0.0, 0.0 }) orelse {
        c.Py_DECREF(err);
        return null;
    };
    const zero_vel = py.vec3Tuple(.{ 0.0, 0.0, 0.0 }) orelse {
        c.Py_DECREF(err);
        c.Py_DECREF(zero_pos);
        return null;
    };
    const result = py.tuple(3) orelse {
        c.Py_DECREF(err);
        c.Py_DECREF(zero_pos);
        c.Py_DECREF(zero_vel);
        return null;
    };
    py.tupleSet(result, 0, err);
    py.tupleSet(result, 1, zero_pos);
    py.tupleSet(result, 2, zero_vel);
    return result;
}

fn buildSuccessResult(pv: [2][3]f64) [*c]c.PyObject {
    // Return (0, (x,y,z), (vx,vy,vz))
    const err = py.int(0) orelse return null;
    const pos = py.vec3Tuple(pv[0]) orelse {
        c.Py_DECREF(err);
        return null;
    };
    const vel = py.vec3Tuple(pv[1]) orelse {
        c.Py_DECREF(err);
        c.Py_DECREF(pos);
        return null;
    };
    const result = py.tuple(3) orelse {
        c.Py_DECREF(err);
        c.Py_DECREF(pos);
        c.Py_DECREF(vel);
        return null;
    };
    py.tupleSet(result, 0, err);
    py.tupleSet(result, 1, pos);
    py.tupleSet(result, 2, vel);
    return result;
}

/// Generic SIMD batch propagation for a single satellite (SGP4 or SDP4).
/// Processes times in chunks of 4 using propagateN, remainder with scalar propagate.
fn propagateArray(comptime T: type, propagator: *const T, n: usize, epochJd: f64, jd_ptr: [*]const f64, fr_ptr: [*]const f64, pos_ptr: [*]f64, vel_ptr: [*]f64) void {
    var i: usize = 0;
    // Process in SIMD chunks of 4
    while (i + 4 <= n) : (i += 4) {
        var times: [4]f64 = undefined;
        for (0..4) |j| {
            times[j] = ((jd_ptr[i + j] + fr_ptr[i + j]) - epochJd) * constants.minutesPerDay;
        }
        if (propagator.propagateN(4, times)) |results| {
            for (0..4) |j| {
                const base = (i + j) * 3;
                pos_ptr[base + 0] = results[j][0][0];
                pos_ptr[base + 1] = results[j][0][1];
                pos_ptr[base + 2] = results[j][0][2];
                vel_ptr[base + 0] = results[j][1][0];
                vel_ptr[base + 1] = results[j][1][1];
                vel_ptr[base + 2] = results[j][1][2];
            }
        } else |_| {
            for (0..4) |j| {
                const base = (i + j) * 3;
                pos_ptr[base + 0] = 0;
                pos_ptr[base + 1] = 0;
                pos_ptr[base + 2] = 0;
                vel_ptr[base + 0] = 0;
                vel_ptr[base + 1] = 0;
                vel_ptr[base + 2] = 0;
            }
        }
    }
    // Handle remainder with scalar propagation
    while (i < n) : (i += 1) {
        const tsince = ((jd_ptr[i] + fr_ptr[i]) - epochJd) * constants.minutesPerDay;
        if (propagator.propagate(tsince)) |result| {
            const base = i * 3;
            pos_ptr[base + 0] = result[0][0];
            pos_ptr[base + 1] = result[0][1];
            pos_ptr[base + 2] = result[0][2];
            vel_ptr[base + 0] = result[1][0];
            vel_ptr[base + 1] = result[1][1];
            vel_ptr[base + 2] = result[1][2];
        } else |_| {
            const base = i * 3;
            pos_ptr[base + 0] = 0;
            pos_ptr[base + 1] = 0;
            pos_ptr[base + 2] = 0;
            vel_ptr[base + 0] = 0;
            vel_ptr[base + 1] = 0;
            vel_ptr[base + 2] = 0;
        }
    }
}

/// Instance method: sgp4_array_into(jd, fr, positions, velocities) -> None
/// SIMD batch propagation writing directly into numpy buffers.
fn satrec_sgp4_array_into(self_obj: [*c]c.PyObject, args: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    const self: *SatrecObject = @ptrCast(@alignCast(self_obj));
    var jd_obj: [*c]c.PyObject = null;
    var fr_obj: [*c]c.PyObject = null;
    var pos_obj: [*c]c.PyObject = null;
    var vel_obj: [*c]c.PyObject = null;

    if (c.PyArg_ParseTuple(args, "OOOO", &jd_obj, &fr_obj, &pos_obj, &vel_obj) == 0) return null;

    // Get buffers
    var jd_buf = std.mem.zeroes(c.Py_buffer);
    var fr_buf = std.mem.zeroes(c.Py_buffer);
    var pos_buf = std.mem.zeroes(c.Py_buffer);
    var vel_buf = std.mem.zeroes(c.Py_buffer);

    if (c.PyObject_GetBuffer(jd_obj, &jd_buf, c.PyBUF_SIMPLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&jd_buf);
    if (c.PyObject_GetBuffer(fr_obj, &fr_buf, c.PyBUF_SIMPLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&fr_buf);
    if (c.PyObject_GetBuffer(pos_obj, &pos_buf, c.PyBUF_WRITABLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&pos_buf);
    if (c.PyObject_GetBuffer(vel_obj, &vel_buf, c.PyBUF_WRITABLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&vel_buf);

    const n = @as(usize, @intCast(jd_buf.len)) / @sizeOf(f64);
    const epochJd = self.jdsatepoch + self.jdsatepochF;

    const jd_ptr: [*]const f64 = @ptrCast(@alignCast(jd_buf.buf));
    const fr_ptr: [*]const f64 = @ptrCast(@alignCast(fr_buf.buf));
    const pos_ptr: [*]f64 = @ptrCast(@alignCast(pos_buf.buf));
    const vel_ptr: [*]f64 = @ptrCast(@alignCast(vel_buf.buf));

    if (self.sgp4) |sgp4| {
        propagateArray(Sgp4, sgp4, n, epochJd, jd_ptr, fr_ptr, pos_ptr, vel_ptr);
    } else if (getSdp4(self)) |sdp4| {
        // Use sorted carry for SDP4 — avoids re-integrating resonance from t=0
        // stride=1 sat, offset=0 for single-satellite contiguous output
        propagateArraySdp4SortedStrided(sdp4, n, epochJd, jd_ptr, fr_ptr, pos_ptr, vel_ptr, 3, 0);
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

// Property getters - using comptime-generated wrappers like tle.zig
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
    return py.int(@intCast(tle.firstLine.satelliteNumber));
}
fn getEpochyr(tle: *const Tle) [*c]c.PyObject {
    return py.int(@intCast(tle.firstLine.epochYear));
}
fn getEpochdays(tle: *const Tle) [*c]c.PyObject {
    return py.float(tle.firstLine.epochDay);
}
fn getEcco(tle: *const Tle) [*c]c.PyObject {
    return py.float(tle.secondLine.eccentricity);
}
fn getInclo(tle: *const Tle) [*c]c.PyObject {
    return py.float(tle.secondLine.inclination * constants.deg2rad);
}
fn getNodeo(tle: *const Tle) [*c]c.PyObject {
    return py.float(tle.secondLine.rightAscension * constants.deg2rad);
}
fn getArgpo(tle: *const Tle) [*c]c.PyObject {
    return py.float(tle.secondLine.perigee * constants.deg2rad);
}
fn getMo(tle: *const Tle) [*c]c.PyObject {
    return py.float(tle.secondLine.mAnomaly * constants.deg2rad);
}
fn getNoKozai(tle: *const Tle) [*c]c.PyObject {
    // Convert rev/day to rad/min: mMotion * 2π / 1440
    return py.float(tle.secondLine.mMotion * constants.twoPi / constants.minutesPerDay);
}
fn getBstar(tle: *const Tle) [*c]c.PyObject {
    return py.float(tle.firstLine.bstarDrag);
}
fn getNdot(tle: *const Tle) [*c]c.PyObject {
    // First derivative of mean motion (rev/day^2 -> rad/min^2)
    // The TLE stores it in rev/day^2 / 2
    return py.float(tle.firstLine.firstDerMeanMotion * constants.twoPi / (constants.minutesPerDay * constants.minutesPerDay));
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

fn getA(self: *const SatrecObject) [*c]c.PyObject {
    const el = getSgp4Elements(self) orelse {
        py.raiseRuntime("Satrec not initialized");
        return null;
    };
    return py.float(el.a);
}
fn getAlta(self: *const SatrecObject) [*c]c.PyObject {
    const el = getSgp4Elements(self) orelse {
        py.raiseRuntime("Satrec not initialized");
        return null;
    };
    return py.float(el.a * (1.0 + el.ecco) - 1.0);
}
fn getAltp(self: *const SatrecObject) [*c]c.PyObject {
    const el = getSgp4Elements(self) orelse {
        py.raiseRuntime("Satrec not initialized");
        return null;
    };
    return py.float(el.a * (1.0 - el.ecco) - 1.0);
}
fn getIsDeepSpace(self: *const SatrecObject) [*c]c.PyObject {
    return c.PyBool_FromLong(if (self.sdp4 != null) @as(c_long, 1) else @as(c_long, 0));
}

const satrec_getset = [_]c.PyGetSetDef{
    // TLE-derived attributes
    .{ .name = "satnum", .get = @ptrCast(&makeTleGetter(getSatnum)), .set = null, .doc = "NORAD catalog number", .closure = null },
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
    // Satrec state attributes
    .{ .name = "jdsatepoch", .get = @ptrCast(&makeGetter(getJdsatepoch)), .set = null, .doc = "Epoch Julian date (integer part)", .closure = null },
    .{ .name = "jdsatepochF", .get = @ptrCast(&makeGetter(getJdsatepochF)), .set = null, .doc = "Epoch Julian date (fractional part)", .closure = null },
    .{ .name = "error", .get = @ptrCast(&makeGetter(getError)), .set = null, .doc = "Last error code (0 = success)", .closure = null },
    .{ .name = "t", .get = @ptrCast(&makeGetter(getT)), .set = null, .doc = "Last propagation time (minutes from epoch)", .closure = null },
    // Derived orbital elements
    .{ .name = "a", .get = @ptrCast(&makeGetter(getA)), .set = null, .doc = "Semi-major axis (Earth radii)", .closure = null },
    .{ .name = "alta", .get = @ptrCast(&makeGetter(getAlta)), .set = null, .doc = "Apoapsis altitude (Earth radii)", .closure = null },
    .{ .name = "altp", .get = @ptrCast(&makeGetter(getAltp)), .set = null, .doc = "Periapsis altitude (Earth radii)", .closure = null },
    // Propagator type
    .{ .name = "is_deep_space", .get = @ptrCast(&makeGetter(getIsDeepSpace)), .set = null, .doc = "True if using SDP4 deep-space propagator", .closure = null },
    .{ .name = null, .get = null, .set = null, .doc = null, .closure = null },
};

// Module-level functions

/// Batch SDP4 propagation with threading.
/// sdp4_batch_propagate_into(satrecs, jd, fr, positions, velocities,
///     output_stride=-1, sat_offset=0) -> None
/// Propagates multiple SDP4 satellites in parallel using Zig threads.
/// Default (output_stride=-1): satellite-major output (n_sdp4, n_times, 3).
/// With output_stride/sat_offset: time-major output, writing satellite s at
/// pos[t * output_stride * 3 + (sat_offset + s) * 3] for direct mixed-constellation use.
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

    // Get buffers
    var jd_buf = std.mem.zeroes(c.Py_buffer);
    var fr_buf = std.mem.zeroes(c.Py_buffer);
    var pos_buf = std.mem.zeroes(c.Py_buffer);
    var vel_buf = std.mem.zeroes(c.Py_buffer);

    if (c.PyObject_GetBuffer(jd_obj, &jd_buf, c.PyBUF_SIMPLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&jd_buf);
    if (c.PyObject_GetBuffer(fr_obj, &fr_buf, c.PyBUF_SIMPLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&fr_buf);
    if (c.PyObject_GetBuffer(pos_obj, &pos_buf, c.PyBUF_WRITABLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&pos_buf);
    if (c.PyObject_GetBuffer(vel_obj, &vel_buf, c.PyBUF_WRITABLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&vel_buf);

    const n_times = @as(usize, @intCast(jd_buf.len)) / @sizeOf(f64);
    const jd_ptr: [*]const f64 = @ptrCast(@alignCast(jd_buf.buf));
    const fr_ptr: [*]const f64 = @ptrCast(@alignCast(fr_buf.buf));
    const pos_ptr: [*]f64 = @ptrCast(@alignCast(pos_buf.buf));
    const vel_ptr: [*]f64 = @ptrCast(@alignCast(vel_buf.buf));

    // Extract SDP4 propagators and epochs from Satrec objects
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
        defer c.Py_DECREF(item);

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

    // Propagate in parallel using threads
    const stride: usize = if (output_stride_arg > 0) @intCast(output_stride_arg) else num_sats;
    const offset: usize = @intCast(sat_offset_arg);
    sdp4BatchPropagate(sdp4_ptrs, epochs, n_times, stride, offset, jd_ptr, fr_ptr, pos_ptr, vel_ptr);

    return py.none();
}

const SdpSimdN = Sdp4Batch.BatchSize;

fn sdp4BatchPropagate(
    sdp4_ptrs: []const *const Sdp4,
    epochs: []const f64,
    n_times: usize,
    output_stride: usize,
    sat_offset: usize,
    jd_ptr: [*]const f64,
    fr_ptr: [*]const f64,
    pos_ptr: [*]f64,
    vel_ptr: [*]f64,
) void {
    const num_sats = sdp4_ptrs.len;
    const num_simd_batches = (num_sats + SdpSimdN - 1) / SdpSimdN;

    // Build SIMD batch elements by transposing groups of N Sdp4.Elements
    const batch_els = allocator.alloc(Sdp4Batch.Sdp4BatchElements(SdpSimdN), num_simd_batches) catch return;
    defer allocator.free(batch_els);
    const batch_epochs = allocator.alloc([SdpSimdN]f64, num_simd_batches) catch return;
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

    // Build origIndices: sequential from sat_offset
    const padded = num_simd_batches * SdpSimdN;
    const origIndices = allocator.alloc(u32, padded) catch return;
    defer allocator.free(origIndices);
    for (origIndices, 0..) |*v, i| v.* = @intCast(sat_offset + i);

    const out_size = output_stride * n_times * 3;
    astroz.Constellation.propagateSdp4Constellation(
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
        .timeMajor,
    ) catch return;
}

/// SDP4 sorted propagation writing to strided time-major output (scalar fallback).
/// For satellite at column `sat_col`, time `t` writes to:
///   pos_ptr[t * time_stride + sat_col * 3]
fn propagateArraySdp4SortedStrided(
    propagator: *const Sdp4,
    n_times: usize,
    epochJd: f64,
    jd_ptr: [*]const f64,
    fr_ptr: [*]const f64,
    pos_ptr: [*]f64,
    vel_ptr: [*]f64,
    time_stride: usize,
    sat_col: usize,
) void {
    var carry = Sdp4.ResonanceCarry{
        .atime = 0.0,
        .xli = propagator.elements.xlamo,
        .xni = propagator.elements.sgp4.noUnkozai,
    };
    const col_offset = sat_col * 3;
    for (0..n_times) |i| {
        const base = i * time_stride + col_offset;
        const tsince = ((jd_ptr[i] + fr_ptr[i]) - epochJd) * constants.minutesPerDay;
        if (propagator.propagateCarry(tsince, &carry)) |result| {
            pos_ptr[base + 0] = result[0][0];
            pos_ptr[base + 1] = result[0][1];
            pos_ptr[base + 2] = result[0][2];
            vel_ptr[base + 0] = result[1][0];
            vel_ptr[base + 1] = result[1][1];
            vel_ptr[base + 2] = result[1][2];
        } else |_| {
            pos_ptr[base + 0] = 0;
            pos_ptr[base + 1] = 0;
            pos_ptr[base + 2] = 0;
            vel_ptr[base + 0] = 0;
            vel_ptr[base + 1] = 0;
            vel_ptr[base + 2] = 0;
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

    const jd_obj = py.float(result.jd) orelse return null;
    const fr_obj = py.float(result.fr) orelse {
        c.Py_DECREF(jd_obj);
        return null;
    };
    const tup = py.tuple(2) orelse {
        c.Py_DECREF(jd_obj);
        c.Py_DECREF(fr_obj);
        return null;
    };
    py.tupleSet(tup, 0, jd_obj);
    py.tupleSet(tup, 1, fr_obj);
    return tup;
}

/// days2mdhms(year, days) -> (month, day, hour, minute, second)
pub fn pyDays2mdhms(_: [*c]c.PyObject, args: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    var year: c_int = undefined;
    var days: f64 = undefined;

    if (c.PyArg_ParseTuple(args, "id", &year, &days) == 0)
        return null;

    const result = Datetime.days2mdhms(@intCast(year), days);

    const month_obj = py.int(@intCast(result.month)) orelse return null;
    const day_obj = py.int(@intCast(result.day)) orelse {
        c.Py_DECREF(month_obj);
        return null;
    };
    const hour_obj = py.int(@intCast(result.hour)) orelse {
        c.Py_DECREF(month_obj);
        c.Py_DECREF(day_obj);
        return null;
    };
    const minute_obj = py.int(@intCast(result.minute)) orelse {
        c.Py_DECREF(month_obj);
        c.Py_DECREF(day_obj);
        c.Py_DECREF(hour_obj);
        return null;
    };
    const second_obj = py.float(result.second) orelse {
        c.Py_DECREF(month_obj);
        c.Py_DECREF(day_obj);
        c.Py_DECREF(hour_obj);
        c.Py_DECREF(minute_obj);
        return null;
    };

    const tup = py.tuple(5) orelse {
        c.Py_DECREF(month_obj);
        c.Py_DECREF(day_obj);
        c.Py_DECREF(hour_obj);
        c.Py_DECREF(minute_obj);
        c.Py_DECREF(second_obj);
        return null;
    };
    py.tupleSet(tup, 0, month_obj);
    py.tupleSet(tup, 1, day_obj);
    py.tupleSet(tup, 2, hour_obj);
    py.tupleSet(tup, 3, minute_obj);
    py.tupleSet(tup, 4, second_obj);
    return tup;
}

// ============================================================================
// SatrecArray - Batch propagator using Sgp4Constellation for SIMD
// ============================================================================

pub const SatrecArrayObject = shared.BatchObjectFields;

pub var SatrecArrayType: c.PyTypeObject = undefined;

fn satrecarray_get_epochs(self_obj: [*c]c.PyObject, _: ?*anyopaque) callconv(.c) [*c]c.PyObject {
    const self: *SatrecArrayObject = @ptrCast(@alignCast(self_obj));
    const epoch_jds = self.epoch_jds orelse {
        py.raiseRuntime("SatrecArray not initialized");
        return null;
    };
    const num_sats = self.num_satellites;

    const list = c.PyList_New(@intCast(num_sats)) orelse return null;
    for (0..num_sats) |i| {
        const jd_obj = py.float(epoch_jds[i]) orelse {
            c.Py_DECREF(list);
            return null;
        };
        _ = c.PyList_SetItem(list, @intCast(i), jd_obj);
    }
    return list;
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
    SatrecArrayType.tp_new = @ptrCast(&satrecarray_new);
    SatrecArrayType.tp_init = @ptrCast(&satrecarray_init);
    SatrecArrayType.tp_dealloc = @ptrCast(&satrecarray_dealloc);
    SatrecArrayType.tp_methods = @constCast(&satrecarray_methods);
    SatrecArrayType.tp_getset = @constCast(&satrecarray_getset);
}

fn satrecarray_new(typ: [*c]c.PyTypeObject, _: [*c]c.PyObject, _: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    const self: *SatrecArrayObject = @ptrCast(@alignCast(c.PyType_GenericAlloc(typ, 0) orelse return null));
    shared.initBatchFields(self);
    return @ptrCast(self);
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

    // Free old data if reinitializing
    shared.freeBatchFields(self);

    // Extract TLEs and gravity model from Satrec objects
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
        defer c.Py_DECREF(item);

        if (c.PyObject_TypeCheck(item, &SatrecType) == 0) {
            py.raiseType("All items must be Satrec objects");
            return -1;
        }

        const satrec = @as(*SatrecObject, @ptrCast(@alignCast(item)));
        tles[i] = (satrec.tle orelse {
            py.raiseValue("Satrec TLE not initialized");
            return -1;
        }).*;

        // Use first satellite's gravity model
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

/// propagate_into(times, positions, velocities, epoch_offsets) -> None
/// Fastest path - matches Sgp4Constellation.propagate_into exactly
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

    // Get buffers
    var times_buf = std.mem.zeroes(c.Py_buffer);
    var pos_buf = std.mem.zeroes(c.Py_buffer);
    var vel_buf = std.mem.zeroes(c.Py_buffer);
    var offsets_buf = std.mem.zeroes(c.Py_buffer);

    if (c.PyObject_GetBuffer(times_obj, &times_buf, c.PyBUF_SIMPLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&times_buf);

    if (c.PyObject_GetBuffer(pos_obj, &pos_buf, c.PyBUF_WRITABLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&pos_buf);

    const num_times = @as(usize, @intCast(times_buf.len)) / @sizeOf(f64);
    const num_sats = self.num_satellites;
    const out_size = num_sats * num_times * 3;

    if (@as(usize, @intCast(pos_buf.len)) < out_size * @sizeOf(f64)) {
        py.raiseValue("positions array too small");
        return null;
    }

    // Optional velocities
    var have_vel = false;
    if (vel_obj != null and !py.isNone(vel_obj)) {
        if (c.PyObject_GetBuffer(vel_obj, &vel_buf, c.PyBUF_WRITABLE | c.PyBUF_FORMAT) < 0) return null;
        have_vel = true;
        if (@as(usize, @intCast(vel_buf.len)) < out_size * @sizeOf(f64)) {
            c.PyBuffer_Release(&vel_buf);
            py.raiseValue("velocities array too small");
            return null;
        }
    }
    defer if (have_vel) c.PyBuffer_Release(&vel_buf);

    // Get epoch offsets (required for this fast path)
    var epoch_offsets: []f64 = undefined;
    var owned_offsets: ?[]f64 = null;
    defer if (owned_offsets) |o| allocator.free(o);

    if (offsets_obj != null and !py.isNone(offsets_obj)) {
        if (c.PyObject_GetBuffer(offsets_obj, &offsets_buf, c.PyBUF_SIMPLE | c.PyBUF_FORMAT) < 0) return null;
        defer c.PyBuffer_Release(&offsets_buf);
        const ptr: [*]f64 = @ptrCast(@alignCast(offsets_buf.buf));
        epoch_offsets = @constCast(ptr[0..self.padded_count]);
    } else {
        // Default: compute from stored epoch_jds assuming reference_jd = first epoch
        const alloc = allocator.alloc(f64, self.padded_count) catch {
            py.raiseRuntime("Out of memory");
            return null;
        };
        owned_offsets = alloc;
        for (0..self.padded_count) |i| {
            alloc[i] = 0.0; // Zero offset = propagate relative to each satellite's epoch
        }
        epoch_offsets = alloc;
    }

    const times: [*]const f64 = @ptrCast(@alignCast(times_buf.buf));
    const pos_out: [*]f64 = @ptrCast(@alignCast(pos_buf.buf));

    astroz.Constellation.propagateConstellation(
        batches,
        num_sats,
        times[0..num_times],
        epoch_offsets,
        pos_out[0..out_size],
        if (have_vel) @as([*]f64, @ptrCast(@alignCast(vel_buf.buf)))[0..out_size] else null,
        .teme,
        0.0, // reference_jd not used for TEME
        null, // no mask
        .timeMajor,
    ) catch {
        py.raiseValue("Propagation failed");
        return null;
    };

    return py.none();
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

//! Python Sgp4 type - wraps astroz.Sgp4

const std = @import("std");
const py = @import("python.zig");
const c = py.c;
const astroz = @import("astroz");
const Sgp4 = astroz.Sgp4;
const tle_mod = @import("tle.zig");

pub const Sgp4Object = extern struct {
    ob_base: c.PyObject,
    sgp4: ?*Sgp4,
};

const allocator = std.heap.c_allocator;
pub var Sgp4Type: c.PyTypeObject = undefined;

// Class constants
pub const WGS84: c_int = 0;
pub const WGS72: c_int = 1;

fn initType() void {
    Sgp4Type = std.mem.zeroes(c.PyTypeObject);
    Sgp4Type.tp_name = "astroz.Sgp4";
    Sgp4Type.tp_doc = "SGP4 propagator.\n\nArgs:\n    tle: Tle object\n    gravity_model: WGS84 (default) or WGS72";
    Sgp4Type.tp_basicsize = @sizeOf(Sgp4Object);
    Sgp4Type.tp_flags = c.Py_TPFLAGS_DEFAULT | c.Py_TPFLAGS_BASETYPE;
    Sgp4Type.tp_new = @ptrCast(&sgp4_new);
    Sgp4Type.tp_init = @ptrCast(&sgp4_init);
    Sgp4Type.tp_dealloc = @ptrCast(&sgp4_dealloc);
    Sgp4Type.tp_methods = @constCast(&sgp4_methods);
}

fn sgp4_new(typ: [*c]c.PyTypeObject, _: [*c]c.PyObject, _: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    const self: *Sgp4Object = @ptrCast(@alignCast(c.PyType_GenericAlloc(typ, 0) orelse return null));
    self.sgp4 = null;
    return @ptrCast(self);
}

fn sgp4_init(self_obj: [*c]c.PyObject, args: [*c]c.PyObject, kwds: [*c]c.PyObject) callconv(.c) c_int {
    const self: *Sgp4Object = @ptrCast(@alignCast(self_obj));
    var tle_obj: [*c]c.PyObject = null;
    var grav_model: c_int = 0;

    const kwlist = [_:null]?[*:0]const u8{ "tle", "gravity_model", null };
    if (c.PyArg_ParseTupleAndKeywords(args, kwds, "O|i", @ptrCast(@constCast(&kwlist)), &tle_obj, &grav_model) == 0)
        return -1;

    if (c.PyObject_TypeCheck(tle_obj, &tle_mod.TleType) == 0) {
        py.raiseType("First argument must be a Tle object");
        return -1;
    }

    const tle = @as(*tle_mod.TleObject, @ptrCast(@alignCast(tle_obj))).tle orelse {
        py.raiseValue("TLE not initialized");
        return -1;
    };

    if (self.sgp4) |old| allocator.destroy(old);

    const grav = if (grav_model == 1) astroz.constants.wgs72 else astroz.constants.wgs84;
    const sgp4 = Sgp4.init(tle.*, grav) catch |e| {
        py.raiseValue(switch (e) {
            Sgp4.Error.DeepSpaceNotSupported => "Deep space not supported",
            Sgp4.Error.InvalidEccentricity => "Invalid eccentricity",
            Sgp4.Error.SatelliteDecayed => "Satellite decayed",
        });
        return -1;
    };

    const ptr = allocator.create(Sgp4) catch {
        py.raiseRuntime("Out of memory");
        return -1;
    };
    ptr.* = sgp4;
    self.sgp4 = ptr;
    return 0;
}

fn sgp4_dealloc(self_obj: [*c]c.PyObject) callconv(.c) void {
    const self: *Sgp4Object = @ptrCast(@alignCast(self_obj));
    if (self.sgp4) |s| allocator.destroy(s);
    if (c.Py_TYPE(self_obj)) |tp| if (tp.*.tp_free) |free| free(@ptrCast(self));
}

fn sgp4_propagate(self_obj: [*c]c.PyObject, args: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    const self: *Sgp4Object = @ptrCast(@alignCast(self_obj));
    var tsince: f64 = undefined;
    if (c.PyArg_ParseTuple(args, "d", &tsince) == 0) return null;

    const sgp4 = self.sgp4 orelse {
        py.raiseRuntime("Sgp4 not initialized");
        return null;
    };

    const result = sgp4.propagate(tsince) catch {
        py.raiseValue("Propagation failed");
        return null;
    };

    const pos = py.vec3Tuple(result[0]) orelse return null;
    const vel = py.vec3Tuple(result[1]) orelse {
        c.Py_DECREF(pos);
        return null;
    };
    const outer = py.tuple(2) orelse {
        c.Py_DECREF(pos);
        c.Py_DECREF(vel);
        return null;
    };
    py.tupleSet(outer, 0, pos);
    py.tupleSet(outer, 1, vel);
    return outer;
}

fn sgp4_propagate_into(self_obj: [*c]c.PyObject, args: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    const self: *Sgp4Object = @ptrCast(@alignCast(self_obj));
    var times_obj: [*c]c.PyObject = null;
    var pos_obj: [*c]c.PyObject = null;
    var vel_obj: [*c]c.PyObject = null;
    if (c.PyArg_ParseTuple(args, "OOO", &times_obj, &pos_obj, &vel_obj) == 0) return null;

    const sgp4 = self.sgp4 orelse {
        py.raiseRuntime("Sgp4 not initialized");
        return null;
    };

    var times_buf = std.mem.zeroes(c.Py_buffer);
    var pos_buf = std.mem.zeroes(c.Py_buffer);
    var vel_buf = std.mem.zeroes(c.Py_buffer);

    if (c.PyObject_GetBuffer(times_obj, &times_buf, c.PyBUF_SIMPLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&times_buf);
    if (c.PyObject_GetBuffer(pos_obj, &pos_buf, c.PyBUF_WRITABLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&pos_buf);
    if (c.PyObject_GetBuffer(vel_obj, &vel_buf, c.PyBUF_WRITABLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&vel_buf);

    const n = @as(usize, @intCast(times_buf.len)) / @sizeOf(f64);
    if (@as(usize, @intCast(pos_buf.len)) < n * 3 * @sizeOf(f64) or
        @as(usize, @intCast(vel_buf.len)) < n * 3 * @sizeOf(f64))
    {
        py.raiseValue("Output arrays too small");
        return null;
    }

    const times: [*]f64 = @ptrCast(@alignCast(times_buf.buf));
    const pos_out: [*]f64 = @ptrCast(@alignCast(pos_buf.buf));
    const vel_out: [*]f64 = @ptrCast(@alignCast(vel_buf.buf));

    // SIMD batches of 4
    var i: usize = 0;
    while (i + 4 <= n) : (i += 4) {
        const results = sgp4.propagateV4(.{ times[i], times[i + 1], times[i + 2], times[i + 3] }) catch {
            py.raiseValue("Propagation failed");
            return null;
        };
        inline for (0..4) |j| {
            const idx = (i + j) * 3;
            pos_out[idx] = results[j][0][0];
            pos_out[idx + 1] = results[j][0][1];
            pos_out[idx + 2] = results[j][0][2];
            vel_out[idx] = results[j][1][0];
            vel_out[idx + 1] = results[j][1][1];
            vel_out[idx + 2] = results[j][1][2];
        }
    }
    // Remainder
    while (i < n) : (i += 1) {
        const result = sgp4.propagate(times[i]) catch {
            py.raiseValue("Propagation failed");
            return null;
        };
        const idx = i * 3;
        pos_out[idx] = result[0][0];
        pos_out[idx + 1] = result[0][1];
        pos_out[idx + 2] = result[0][2];
        vel_out[idx] = result[1][0];
        vel_out[idx + 1] = result[1][1];
        vel_out[idx + 2] = result[1][2];
    }
    return py.none();
}

const sgp4_methods = [_]c.PyMethodDef{
    .{ .ml_name = "propagate", .ml_meth = @ptrCast(&sgp4_propagate), .ml_flags = c.METH_VARARGS, .ml_doc = "propagate(tsince) -> ((x,y,z), (vx,vy,vz))" },
    .{ .ml_name = "propagate_into", .ml_meth = @ptrCast(&sgp4_propagate_into), .ml_flags = c.METH_VARARGS, .ml_doc = "propagate_into(times, positions, velocities) -> None" },
    .{ .ml_name = null, .ml_meth = null, .ml_flags = 0, .ml_doc = null },
};

pub fn ready() c_int {
    initType();
    return c.PyType_Ready(&Sgp4Type);
}

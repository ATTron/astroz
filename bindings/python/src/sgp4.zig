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

// =============================================================================
// Sgp4Batch: Multi-satellite batch propagator
// =============================================================================

pub const Sgp4BatchObject = extern struct {
    ob_base: c.PyObject,
    elements: ?*Sgp4.ElementsV4,
};

pub var Sgp4BatchType: c.PyTypeObject = undefined;

fn initBatchType() void {
    Sgp4BatchType = std.mem.zeroes(c.PyTypeObject);
    Sgp4BatchType.tp_name = "astroz.Sgp4Batch";
    Sgp4BatchType.tp_doc = "Batch SGP4 propagator for 4 satellites.\n\nArgs:\n    tles: tuple of 4 Tle objects\n    gravity_model: WGS84 (default) or WGS72";
    Sgp4BatchType.tp_basicsize = @sizeOf(Sgp4BatchObject);
    Sgp4BatchType.tp_flags = c.Py_TPFLAGS_DEFAULT | c.Py_TPFLAGS_BASETYPE;
    Sgp4BatchType.tp_new = @ptrCast(&sgp4batch_new);
    Sgp4BatchType.tp_init = @ptrCast(&sgp4batch_init);
    Sgp4BatchType.tp_dealloc = @ptrCast(&sgp4batch_dealloc);
    Sgp4BatchType.tp_methods = @constCast(&sgp4batch_methods);
}

fn sgp4batch_new(typ: [*c]c.PyTypeObject, _: [*c]c.PyObject, _: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    const self: *Sgp4BatchObject = @ptrCast(@alignCast(c.PyType_GenericAlloc(typ, 0) orelse return null));
    self.elements = null;
    return @ptrCast(self);
}

fn sgp4batch_init(self_obj: [*c]c.PyObject, args: [*c]c.PyObject, kwds: [*c]c.PyObject) callconv(.c) c_int {
    const self: *Sgp4BatchObject = @ptrCast(@alignCast(self_obj));
    var tles_obj: [*c]c.PyObject = null;
    var grav_model: c_int = 0;

    const kwlist = [_:null]?[*:0]const u8{ "tles", "gravity_model", null };
    if (c.PyArg_ParseTupleAndKeywords(args, kwds, "O|i", @ptrCast(@constCast(&kwlist)), &tles_obj, &grav_model) == 0)
        return -1;

    // Check that tles_obj is a tuple or list with exactly 4 elements
    if (c.PyTuple_Check(tles_obj) == 0 and c.PyList_Check(tles_obj) == 0) {
        py.raiseType("First argument must be a tuple or list of 4 Tle objects");
        return -1;
    }

    const size = c.PySequence_Size(tles_obj);
    if (size != 4) {
        py.raiseValue("Must provide exactly 4 Tle objects");
        return -1;
    }

    // Extract 4 TLEs
    var tles: [4]astroz.Tle = undefined;
    for (0..4) |i| {
        const item = c.PySequence_GetItem(tles_obj, @intCast(i)) orelse {
            py.raiseValue("Failed to get TLE from sequence");
            return -1;
        };
        defer c.Py_DECREF(item);

        if (c.PyObject_TypeCheck(item, &tle_mod.TleType) == 0) {
            py.raiseType("All items must be Tle objects");
            return -1;
        }

        const tle_obj_ptr = @as(*tle_mod.TleObject, @ptrCast(@alignCast(item)));
        tles[i] = (tle_obj_ptr.tle orelse {
            py.raiseValue("TLE not initialized");
            return -1;
        }).*;
    }

    if (self.elements) |old| allocator.destroy(old);

    const grav = if (grav_model == 1) astroz.constants.wgs72 else astroz.constants.wgs84;
    const elements = Sgp4.initElementsV4(tles, grav) catch |e| {
        py.raiseValue(switch (e) {
            Sgp4.Error.DeepSpaceNotSupported => "Deep space not supported",
            Sgp4.Error.InvalidEccentricity => "Invalid eccentricity",
            Sgp4.Error.SatelliteDecayed => "Satellite decayed",
        });
        return -1;
    };

    const ptr = allocator.create(Sgp4.ElementsV4) catch {
        py.raiseRuntime("Out of memory");
        return -1;
    };
    ptr.* = elements;
    self.elements = ptr;
    return 0;
}

fn sgp4batch_dealloc(self_obj: [*c]c.PyObject) callconv(.c) void {
    const self: *Sgp4BatchObject = @ptrCast(@alignCast(self_obj));
    if (self.elements) |e| allocator.destroy(e);
    if (c.Py_TYPE(self_obj)) |tp| if (tp.*.tp_free) |free| free(@ptrCast(self));
}

fn sgp4batch_propagate(self_obj: [*c]c.PyObject, args: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    const self: *Sgp4BatchObject = @ptrCast(@alignCast(self_obj));
    var tsince: f64 = undefined;
    if (c.PyArg_ParseTuple(args, "d", &tsince) == 0) return null;

    const elements = self.elements orelse {
        py.raiseRuntime("Sgp4Batch not initialized");
        return null;
    };

    const results = Sgp4.propagateSatellitesV4(elements, tsince) catch {
        py.raiseValue("Propagation failed");
        return null;
    };

    // Return a tuple of 4 ((pos), (vel)) tuples
    const outer = py.tuple(4) orelse return null;
    for (0..4) |i| {
        const pos = py.vec3Tuple(results[i][0]) orelse {
            c.Py_DECREF(outer);
            return null;
        };
        const vel = py.vec3Tuple(results[i][1]) orelse {
            c.Py_DECREF(pos);
            c.Py_DECREF(outer);
            return null;
        };
        const pv = py.tuple(2) orelse {
            c.Py_DECREF(pos);
            c.Py_DECREF(vel);
            c.Py_DECREF(outer);
            return null;
        };
        py.tupleSet(pv, 0, pos);
        py.tupleSet(pv, 1, vel);
        py.tupleSet(outer, @intCast(i), pv);
    }
    return outer;
}

fn sgp4batch_propagate_into(self_obj: [*c]c.PyObject, args: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    const self: *Sgp4BatchObject = @ptrCast(@alignCast(self_obj));
    var times_obj: [*c]c.PyObject = null;
    var pos_obj: [*c]c.PyObject = null;
    var vel_obj: [*c]c.PyObject = null;
    if (c.PyArg_ParseTuple(args, "OOO", &times_obj, &pos_obj, &vel_obj) == 0) return null;

    const elements = self.elements orelse {
        py.raiseRuntime("Sgp4Batch not initialized");
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
    // For batch: 4 satellites × 3 components × n times = 12*n values per array
    if (@as(usize, @intCast(pos_buf.len)) < n * 4 * 3 * @sizeOf(f64) or
        @as(usize, @intCast(vel_buf.len)) < n * 4 * 3 * @sizeOf(f64))
    {
        py.raiseValue("Output arrays too small (need n*4*3 elements each)");
        return null;
    }

    const times: [*]f64 = @ptrCast(@alignCast(times_buf.buf));
    const pos_out: [*]f64 = @ptrCast(@alignCast(pos_buf.buf));
    const vel_out: [*]f64 = @ptrCast(@alignCast(vel_buf.buf));

    for (0..n) |t| {
        const results = Sgp4.propagateSatellitesV4(elements, times[t]) catch {
            py.raiseValue("Propagation failed");
            return null;
        };
        // Layout: [t0_sat0_x, t0_sat0_y, t0_sat0_z, t0_sat1_x, ..., t0_sat3_z, t1_sat0_x, ...]
        const timeBase = t * 12; // 4 satellites × 3 components
        inline for (0..4) |s| {
            const idx = timeBase + s * 3;
            pos_out[idx] = results[s][0][0];
            pos_out[idx + 1] = results[s][0][1];
            pos_out[idx + 2] = results[s][0][2];
            vel_out[idx] = results[s][1][0];
            vel_out[idx + 1] = results[s][1][1];
            vel_out[idx + 2] = results[s][1][2];
        }
    }
    return py.none();
}

const sgp4batch_methods = [_]c.PyMethodDef{
    .{ .ml_name = "propagate", .ml_meth = @ptrCast(&sgp4batch_propagate), .ml_flags = c.METH_VARARGS, .ml_doc = "propagate(tsince) -> (((x,y,z), (vx,vy,vz)), ...) for 4 sats" },
    .{ .ml_name = "propagate_into", .ml_meth = @ptrCast(&sgp4batch_propagate_into), .ml_flags = c.METH_VARARGS, .ml_doc = "propagate_into(times, positions, velocities) -> None\n\nArrays have shape (n_times, 4, 3) flattened" },
    .{ .ml_name = null, .ml_meth = null, .ml_flags = 0, .ml_doc = null },
};

// =============================================================================
// Sgp4Constellation: Multi-batch propagator for large constellations
// Eliminates Python loop overhead by processing all batches in a single call
// =============================================================================

pub const Sgp4ConstellationObject = extern struct {
    ob_base: c.PyObject,
    batches: ?[*]Sgp4.ElementsV4,
    num_batches: usize,
    num_satellites: usize,
};

pub var Sgp4ConstellationType: c.PyTypeObject = undefined;

fn initConstellationType() void {
    Sgp4ConstellationType = std.mem.zeroes(c.PyTypeObject);
    Sgp4ConstellationType.tp_name = "astroz.Sgp4Constellation";
    Sgp4ConstellationType.tp_doc = "Constellation propagator for many satellites.\n\nArgs:\n    tles: list of Tle objects (will be grouped into batches of 4)\n    gravity_model: WGS84 (default) or WGS72";
    Sgp4ConstellationType.tp_basicsize = @sizeOf(Sgp4ConstellationObject);
    Sgp4ConstellationType.tp_flags = c.Py_TPFLAGS_DEFAULT | c.Py_TPFLAGS_BASETYPE;
    Sgp4ConstellationType.tp_new = @ptrCast(&constellation_new);
    Sgp4ConstellationType.tp_init = @ptrCast(&constellation_init);
    Sgp4ConstellationType.tp_dealloc = @ptrCast(&constellation_dealloc);
    Sgp4ConstellationType.tp_methods = @constCast(&constellation_methods);
}

fn constellation_new(typ: [*c]c.PyTypeObject, _: [*c]c.PyObject, _: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    const self: *Sgp4ConstellationObject = @ptrCast(@alignCast(c.PyType_GenericAlloc(typ, 0) orelse return null));
    self.batches = null;
    self.num_batches = 0;
    self.num_satellites = 0;
    return @ptrCast(self);
}

fn constellation_init(self_obj: [*c]c.PyObject, args: [*c]c.PyObject, kwds: [*c]c.PyObject) callconv(.c) c_int {
    const self: *Sgp4ConstellationObject = @ptrCast(@alignCast(self_obj));
    var tles_obj: [*c]c.PyObject = null;
    var grav_model: c_int = 0;

    const kwlist = [_:null]?[*:0]const u8{ "tles", "gravity_model", null };
    if (c.PyArg_ParseTupleAndKeywords(args, kwds, "O|i", @ptrCast(@constCast(&kwlist)), &tles_obj, &grav_model) == 0)
        return -1;

    if (c.PySequence_Check(tles_obj) == 0) {
        py.raiseType("First argument must be a sequence of Tle objects");
        return -1;
    }

    const num_tles: usize = @intCast(c.PySequence_Size(tles_obj));
    if (num_tles == 0) {
        py.raiseValue("Must provide at least one Tle object");
        return -1;
    }

    // Free old data if reinitializing
    if (self.batches) |old_batches| {
        allocator.free(old_batches[0..self.num_batches]);
    }

    // Pad to multiple of 4
    const padded_count = ((num_tles + 3) / 4) * 4;
    const num_batches = padded_count / 4;

    // Allocate batch array
    const batches = allocator.alloc(Sgp4.ElementsV4, num_batches) catch {
        py.raiseRuntime("Out of memory");
        return -1;
    };

    const grav = if (grav_model == 1) astroz.constants.wgs72 else astroz.constants.wgs84;

    // Process TLEs in groups of 4
    var batch_idx: usize = 0;
    while (batch_idx < num_batches) : (batch_idx += 1) {
        var batch_tles: [4]astroz.Tle = undefined;

        for (0..4) |i| {
            const tle_idx = batch_idx * 4 + i;
            // For padding, reuse the last TLE
            const actual_idx = if (tle_idx < num_tles) tle_idx else num_tles - 1;

            const item = c.PySequence_GetItem(tles_obj, @intCast(actual_idx)) orelse {
                allocator.free(batches);
                py.raiseValue("Failed to get TLE from sequence");
                return -1;
            };
            defer c.Py_DECREF(item);

            if (c.PyObject_TypeCheck(item, &tle_mod.TleType) == 0) {
                allocator.free(batches);
                py.raiseType("All items must be Tle objects");
                return -1;
            }

            const tle_obj_ptr = @as(*tle_mod.TleObject, @ptrCast(@alignCast(item)));
            batch_tles[i] = (tle_obj_ptr.tle orelse {
                allocator.free(batches);
                py.raiseValue("TLE not initialized");
                return -1;
            }).*;
        }

        batches[batch_idx] = Sgp4.initElementsV4(batch_tles, grav) catch |e| {
            allocator.free(batches);
            py.raiseValue(switch (e) {
                Sgp4.Error.DeepSpaceNotSupported => "Deep space not supported",
                Sgp4.Error.InvalidEccentricity => "Invalid eccentricity",
                Sgp4.Error.SatelliteDecayed => "Satellite decayed",
            });
            return -1;
        };
    }

    self.batches = batches.ptr;
    self.num_batches = num_batches;
    self.num_satellites = num_tles;
    return 0;
}

fn constellation_dealloc(self_obj: [*c]c.PyObject) callconv(.c) void {
    const self: *Sgp4ConstellationObject = @ptrCast(@alignCast(self_obj));
    if (self.batches) |batches| {
        allocator.free(batches[0..self.num_batches]);
    }
    if (c.Py_TYPE(self_obj)) |tp| if (tp.*.tp_free) |free| free(@ptrCast(self));
}

fn constellation_propagate_into(self_obj: [*c]c.PyObject, args: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    const self: *Sgp4ConstellationObject = @ptrCast(@alignCast(self_obj));
    var times_obj: [*c]c.PyObject = null;
    var results_obj: [*c]c.PyObject = null;

    if (c.PyArg_ParseTuple(args, "OO", &times_obj, &results_obj) == 0) return null;

    const batches_ptr = self.batches orelse {
        py.raiseRuntime("Constellation not initialized");
        return null;
    };
    const batches = batches_ptr[0..self.num_batches];

    var times_buf = std.mem.zeroes(c.Py_buffer);
    var results_buf = std.mem.zeroes(c.Py_buffer);

    if (c.PyObject_GetBuffer(times_obj, &times_buf, c.PyBUF_SIMPLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&times_buf);
    if (c.PyObject_GetBuffer(results_obj, &results_buf, c.PyBUF_WRITABLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&results_buf);

    const num_times = @as(usize, @intCast(times_buf.len)) / @sizeOf(f64);
    const required_size = self.num_batches * num_times * 4 * 6 * @sizeOf(f64);

    if (@as(usize, @intCast(results_buf.len)) < required_size) {
        py.raiseValue("Output array too small");
        return null;
    }

    const times: [*]const f64 = @ptrCast(@alignCast(times_buf.buf));
    const results: [*]f64 = @ptrCast(@alignCast(results_buf.buf));

    const results_slice = results[0 .. self.num_batches * num_times * 24];
    const times_slice = times[0..num_times];

    Sgp4.propagateConstellationV4(batches, times_slice, results_slice) catch {
        py.raiseValue("Propagation failed");
        return null;
    };

    return py.none();
}

fn constellation_get_num_satellites(self_obj: [*c]c.PyObject, _: ?*anyopaque) callconv(.c) [*c]c.PyObject {
    const self: *Sgp4ConstellationObject = @ptrCast(@alignCast(self_obj));
    return py.int(@intCast(self.num_satellites));
}

fn constellation_get_num_batches(self_obj: [*c]c.PyObject, _: ?*anyopaque) callconv(.c) [*c]c.PyObject {
    const self: *Sgp4ConstellationObject = @ptrCast(@alignCast(self_obj));
    return py.int(@intCast(self.num_batches));
}

const constellation_methods = [_]c.PyMethodDef{
    .{ .ml_name = "propagate_into", .ml_meth = @ptrCast(&constellation_propagate_into), .ml_flags = c.METH_VARARGS, .ml_doc = "propagate_into(times, results) -> None\n\nPropagates all satellites across all times in a single call.\nresults shape: (num_times, num_batches, 4, 6) flattened" },
    .{ .ml_name = null, .ml_meth = null, .ml_flags = 0, .ml_doc = null },
};

const constellation_getset = [_]c.PyGetSetDef{
    .{ .name = "num_satellites", .get = @ptrCast(&constellation_get_num_satellites), .set = null, .doc = "Number of satellites", .closure = null },
    .{ .name = "num_batches", .get = @ptrCast(&constellation_get_num_batches), .set = null, .doc = "Number of 4-satellite batches", .closure = null },
    .{ .name = null, .get = null, .set = null, .doc = null, .closure = null },
};

pub fn ready() c_int {
    initType();
    initBatchType();
    initConstellationType();
    Sgp4ConstellationType.tp_getset = @constCast(&constellation_getset);
    if (c.PyType_Ready(&Sgp4Type) < 0) return -1;
    if (c.PyType_Ready(&Sgp4BatchType) < 0) return -1;
    return c.PyType_Ready(&Sgp4ConstellationType);
}

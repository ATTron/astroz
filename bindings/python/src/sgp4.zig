//! Python Sgp4 type - wraps astroz.Sgp4

const std = @import("std");
const py = @import("python.zig");
const c = py.c;
const astroz = @import("astroz");
const Sgp4 = astroz.Sgp4;
const tle_mod = @import("tle.zig");

const allocator = std.heap.c_allocator;

/// Batch size from core library (4 for AVX2, 8 for AVX512)
const BatchSize = astroz.Sgp4Constellation.BatchSize;

/// Batch elements type for current batch size
const BatchElements = astroz.Sgp4Batch.BatchElements(BatchSize);

// Common error message for SGP4 init errors
fn sgp4ErrorMsg(e: Sgp4.Error) [:0]const u8 {
    return switch (e) {
        Sgp4.Error.DeepSpaceNotSupported => "Deep space not supported",
        Sgp4.Error.InvalidEccentricity => "Invalid eccentricity",
        Sgp4.Error.SatelliteDecayed => "Satellite decayed",
    };
}

fn getGravity(grav_model: c_int) @TypeOf(astroz.constants.wgs84) {
    return if (grav_model == 1) astroz.constants.wgs72 else astroz.constants.wgs84;
}

pub const Sgp4Object = extern struct {
    ob_base: c.PyObject,
    sgp4: ?*Sgp4,
};
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

    const sgp4 = Sgp4.init(tle.*, getGravity(grav_model)) catch |e| {
        py.raiseValue(sgp4ErrorMsg(e));
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

    // Process in batches of BatchSize, padding remainder with last valid time
    var i: usize = 0;
    while (i < n) : (i += BatchSize) {
        const remaining = n - i;
        var batch: [BatchSize]f64 = undefined;
        inline for (0..BatchSize) |j| batch[j] = times[i + @min(j, remaining - 1)];

        const results = sgp4.propagateN(BatchSize, batch) catch {
            py.raiseValue("Propagation failed");
            return null;
        };
        const count = @min(remaining, BatchSize);
        for (0..count) |j| {
            const idx = (i + j) * 3;
            pos_out[idx] = results[j][0][0];
            pos_out[idx + 1] = results[j][0][1];
            pos_out[idx + 2] = results[j][0][2];
            vel_out[idx] = results[j][1][0];
            vel_out[idx + 1] = results[j][1][1];
            vel_out[idx + 2] = results[j][1][2];
        }
    }
    return py.none();
}

const sgp4_methods = [_]c.PyMethodDef{
    .{ .ml_name = "propagate", .ml_meth = @ptrCast(&sgp4_propagate), .ml_flags = c.METH_VARARGS, .ml_doc = "propagate(tsince) -> ((x,y,z), (vx,vy,vz))" },
    .{ .ml_name = "propagate_into", .ml_meth = @ptrCast(&sgp4_propagate_into), .ml_flags = c.METH_VARARGS, .ml_doc = "propagate_into(times, positions, velocities) -> None" },
    .{ .ml_name = null, .ml_meth = null, .ml_flags = 0, .ml_doc = null },
};

// Sgp4Batch: Multi-satellite batch propagator (fixed at 4 satellites for Python API)

const Elements4 = astroz.Sgp4Batch.BatchElements(4);

pub const Sgp4BatchObject = extern struct {
    ob_base: c.PyObject,
    elements: ?*Elements4,
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

    const elements = astroz.Sgp4Batch.initBatchElements(4, tles, getGravity(grav_model)) catch |e| {
        py.raiseValue(sgp4ErrorMsg(e));
        return -1;
    };

    const ptr = allocator.create(Elements4) catch {
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

    const results = astroz.Sgp4Batch.propagateSatellites(4, elements, tsince) catch {
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
        const results = astroz.Sgp4Batch.propagateSatellites(4, elements, times[t]) catch {
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

// Sgp4Constellation: Multi-batch propagator for large constellations

pub const Sgp4ConstellationObject = extern struct {
    ob_base: c.PyObject,
    batches: ?[*]BatchElements,
    num_batches: usize,
    num_satellites: usize,
    epoch_jds: ?[*]f64, // per-satellite epoch Julian dates
    padded_count: usize, // padded to multiple of BatchSize
};

pub var Sgp4ConstellationType: c.PyTypeObject = undefined;

fn initConstellationType() void {
    Sgp4ConstellationType = std.mem.zeroes(c.PyTypeObject);
    Sgp4ConstellationType.tp_name = "astroz.Sgp4Constellation";
    Sgp4ConstellationType.tp_doc = "Constellation propagator for many satellites.\n\nArgs:\n    tles: list of Tle objects (will be grouped into batches)\n    gravity_model: WGS84 (default) or WGS72";
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
    self.epoch_jds = null;
    self.padded_count = 0;
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
    if (self.epoch_jds) |old_ejds| {
        allocator.free(old_ejds[0..self.padded_count]);
    }

    // Extract TLEs from Python sequence into temp slice
    const tles = allocator.alloc(astroz.Tle, num_tles) catch {
        py.raiseRuntime("Out of memory");
        return -1;
    };
    defer allocator.free(tles);

    for (0..num_tles) |i| {
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

    const grav = if (grav_model == 1) astroz.constants.wgs72 else astroz.constants.wgs84;
    const result = buildBatches(tles, grav) orelse return -1;

    self.batches = result.batches.ptr;
    self.num_batches = result.batches.len;
    self.num_satellites = num_tles;
    self.epoch_jds = result.epoch_jds.ptr;
    self.padded_count = result.padded_count;
    return 0;
}

const BatchResult = struct {
    batches: []BatchElements,
    epoch_jds: []f64,
    padded_count: usize,
};

/// Group TLEs into SIMD batches of BatchSize, pad last batch, and extract epoch JDs.
fn buildBatches(tles: []const astroz.Tle, grav: anytype) ?BatchResult {
    const num_tles = tles.len;
    const padded_count = ((num_tles + BatchSize - 1) / BatchSize) * BatchSize;
    const num_batches = padded_count / BatchSize;

    const batches = allocator.alloc(BatchElements, num_batches) catch {
        py.raiseRuntime("Out of memory");
        return null;
    };

    for (0..num_batches) |batch_idx| {
        var batch_tles: [BatchSize]astroz.Tle = undefined;
        inline for (0..BatchSize) |i| {
            const tle_idx = batch_idx * BatchSize + i;
            batch_tles[i] = tles[if (tle_idx < num_tles) tle_idx else num_tles - 1];
        }
        batches[batch_idx] = astroz.Sgp4Batch.initBatchElements(BatchSize, batch_tles, grav) catch |e| {
            allocator.free(batches);
            py.raiseValue(sgp4ErrorMsg(e));
            return null;
        };
    }

    const epoch_jds = allocator.alloc(f64, padded_count) catch {
        allocator.free(batches);
        py.raiseRuntime("Out of memory");
        return null;
    };
    for (0..num_batches) |bi| {
        inline for (0..BatchSize) |si| {
            epoch_jds[bi * BatchSize + si] = batches[bi].epochJd[si];
        }
    }

    return .{ .batches = batches, .epoch_jds = epoch_jds, .padded_count = padded_count };
}

fn constellation_dealloc(self_obj: [*c]c.PyObject) callconv(.c) void {
    const self: *Sgp4ConstellationObject = @ptrCast(@alignCast(self_obj));
    if (self.batches) |batches| {
        allocator.free(batches[0..self.num_batches]);
    }
    if (self.epoch_jds) |ejds| {
        allocator.free(ejds[0..self.padded_count]);
    }
    if (c.Py_TYPE(self_obj)) |tp| if (tp.*.tp_free) |free| free(@ptrCast(self));
}

fn parseOutputMode(output_str: [*c]const u8, output_str_len: c.Py_ssize_t) ?astroz.Sgp4Constellation.OutputMode {
    if (output_str == null or output_str_len <= 0) return .teme;
    const mode = output_str[0..@intCast(output_str_len)];
    if (std.mem.eql(u8, mode, "teme")) return .teme;
    if (std.mem.eql(u8, mode, "ecef")) return .ecef;
    if (std.mem.eql(u8, mode, "geodetic")) return .geodetic;
    py.raiseValue("output must be 'teme', 'ecef', or 'geodetic'");
    return null;
}

/// Result of preparing a padded array from Python input
const PaddedArray = struct {
    data: []f64,
    owned: bool,

    fn deinit(self: *@This()) void {
        if (self.owned) allocator.free(self.data);
    }
};

/// Prepare epoch offsets: pad to padded_count, or create zeros if not provided
fn prepareOffsets(obj: [*c]c.PyObject, buf: *c.Py_buffer, num_sats: usize, padded_count: usize) ?PaddedArray {
    if (py.isNone(obj)) {
        const zeros = allocator.alloc(f64, padded_count) catch {
            py.raiseRuntime("Out of memory");
            return null;
        };
        @memset(zeros, 0.0);
        return .{ .data = zeros, .owned = true };
    }

    if (c.PyObject_GetBuffer(obj, buf, c.PyBUF_SIMPLE | c.PyBUF_FORMAT) < 0) return null;
    const ptr: [*]const f64 = @ptrCast(@alignCast(buf.buf));
    const len = @as(usize, @intCast(buf.len)) / @sizeOf(f64);

    if (len >= padded_count) {
        return .{ .data = @constCast(ptr[0..padded_count]), .owned = false };
    } else if (len >= num_sats) {
        const padded = allocator.alloc(f64, padded_count) catch {
            py.raiseRuntime("Out of memory");
            return null;
        };
        @memcpy(padded[0..num_sats], ptr[0..num_sats]);
        @memset(padded[num_sats..padded_count], 0.0);
        return .{ .data = padded, .owned = true };
    }
    py.raiseValue("epoch_offsets must have at least num_satellites elements");
    return null;
}

/// Prepare satellite mask: pad to padded_count if provided
fn prepareMask(obj: [*c]c.PyObject, buf: *c.Py_buffer, num_sats: usize, padded_count: usize) struct { mask: ?[]const u8, owned: ?[]u8 } {
    if (py.isNone(obj)) return .{ .mask = null, .owned = null };

    if (c.PyObject_GetBuffer(obj, buf, c.PyBUF_SIMPLE | c.PyBUF_FORMAT) < 0) return .{ .mask = null, .owned = null };
    const ptr: [*]const u8 = @ptrCast(@alignCast(buf.buf));
    const len = @as(usize, @intCast(buf.len));

    if (len >= padded_count) {
        return .{ .mask = ptr[0..padded_count], .owned = null };
    } else if (len >= num_sats) {
        const padded = allocator.alloc(u8, padded_count) catch {
            py.raiseRuntime("Out of memory");
            return .{ .mask = null, .owned = null };
        };
        @memcpy(padded[0..num_sats], ptr[0..num_sats]);
        @memset(padded[num_sats..padded_count], 0);
        return .{ .mask = padded, .owned = padded };
    }
    py.raiseValue("satellite_mask must have at least num_satellites elements");
    return .{ .mask = null, .owned = null };
}

fn constellation_propagate_into(self_obj: [*c]c.PyObject, args: [*c]c.PyObject, kwds: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    const self: *Sgp4ConstellationObject = @ptrCast(@alignCast(self_obj));

    // Parse arguments
    var times_obj: [*c]c.PyObject = null;
    var pos_obj: [*c]c.PyObject = null;
    var vel_obj: [*c]c.PyObject = null;
    var offsets_obj: [*c]c.PyObject = null;
    var mask_obj: [*c]c.PyObject = null;
    var output_str: [*c]const u8 = null;
    var output_str_len: c.Py_ssize_t = 0;
    var reference_jd: f64 = 0.0;
    var time_major: c_int = 0;

    const kwlist = [_:null]?[*:0]const u8{
        "times", "positions", "velocities", "epoch_offsets", "satellite_mask", "output", "reference_jd", "time_major", null,
    };
    if (c.PyArg_ParseTupleAndKeywords(args, kwds, "OO|OOOs#di", @ptrCast(@constCast(&kwlist)), &times_obj, &pos_obj, &vel_obj, &offsets_obj, &mask_obj, &output_str, &output_str_len, &reference_jd, &time_major) == 0)
        return null;

    const batches = (self.batches orelse {
        py.raiseRuntime("Constellation not initialized");
        return null;
    })[0..self.num_batches];

    const output_mode = parseOutputMode(output_str, output_str_len) orelse return null;

    // Get buffers
    var times_buf = std.mem.zeroes(c.Py_buffer);
    var pos_buf = std.mem.zeroes(c.Py_buffer);
    var vel_buf = std.mem.zeroes(c.Py_buffer);
    var offsets_buf = std.mem.zeroes(c.Py_buffer);
    var mask_buf = std.mem.zeroes(c.Py_buffer);

    if (c.PyObject_GetBuffer(times_obj, &times_buf, c.PyBUF_SIMPLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&times_buf);

    if (c.PyObject_GetBuffer(pos_obj, &pos_buf, c.PyBUF_WRITABLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&pos_buf);

    const num_times = @as(usize, @intCast(times_buf.len)) / @sizeOf(f64);
    const required_size = self.num_satellites * num_times * 3 * @sizeOf(f64);

    if (@as(usize, @intCast(pos_buf.len)) < required_size) {
        py.raiseValue("positions array too small");
        return null;
    }

    // Optional velocities
    var have_vel = false;
    if (!py.isNone(vel_obj)) {
        if (c.PyObject_GetBuffer(vel_obj, &vel_buf, c.PyBUF_WRITABLE | c.PyBUF_FORMAT) < 0) return null;
        have_vel = true;
        if (@as(usize, @intCast(vel_buf.len)) < required_size) {
            c.PyBuffer_Release(&vel_buf);
            py.raiseValue("velocities array too small");
            return null;
        }
    }
    defer if (have_vel) c.PyBuffer_Release(&vel_buf);

    // Prepare padded arrays
    var offsets = prepareOffsets(offsets_obj, &offsets_buf, self.num_satellites, self.padded_count) orelse return null;
    defer offsets.deinit();
    defer if (!offsets.owned) c.PyBuffer_Release(&offsets_buf);

    const mask_result = prepareMask(mask_obj, &mask_buf, self.num_satellites, self.padded_count);
    defer if (mask_result.owned) |m| allocator.free(m);
    defer if (mask_result.mask != null and mask_result.owned == null) c.PyBuffer_Release(&mask_buf);

    // Build slices and propagate
    const times: [*]const f64 = @ptrCast(@alignCast(times_buf.buf));
    const pos_out: [*]f64 = @ptrCast(@alignCast(pos_buf.buf));
    const out_len = self.num_satellites * num_times * 3;

    astroz.Sgp4Constellation.propagateConstellation(
        batches,
        self.num_satellites,
        times[0..num_times],
        offsets.data,
        pos_out[0..out_len],
        if (have_vel) @as([*]f64, @ptrCast(@alignCast(vel_buf.buf)))[0..out_len] else null,
        output_mode,
        reference_jd,
        mask_result.mask,
        if (time_major != 0) .timeMajor else .satelliteMajor,
    ) catch {
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

fn constellation_get_batch_size(_: [*c]c.PyObject, _: ?*anyopaque) callconv(.c) [*c]c.PyObject {
    return py.int(@intCast(BatchSize));
}

fn constellation_get_epochs(self_obj: [*c]c.PyObject, _: ?*anyopaque) callconv(.c) [*c]c.PyObject {
    const self: *Sgp4ConstellationObject = @ptrCast(@alignCast(self_obj));
    const ejds = self.epoch_jds orelse {
        py.raiseRuntime("Constellation not initialized");
        return null;
    };
    return py.listFromF64(ejds[0..self.num_satellites]);
}

fn constellation_from_tle_text(cls: [*c]c.PyObject, args: [*c]c.PyObject, kwds: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    _ = cls;
    var text_ptr: [*c]const u8 = null;
    var text_len: c.Py_ssize_t = 0;
    var grav_model: c_int = 0;

    const kwlist = [_:null]?[*:0]const u8{ "text", "gravity_model", null };
    if (c.PyArg_ParseTupleAndKeywords(args, kwds, "s#|i", @ptrCast(@constCast(&kwlist)), &text_ptr, &text_len, &grav_model) == 0)
        return null;

    if (text_ptr == null or text_len <= 0) {
        py.raiseValue("TLE text must not be empty");
        return null;
    }

    const text = text_ptr[0..@intCast(text_len)];
    const grav = if (grav_model == 1) astroz.constants.wgs72 else astroz.constants.wgs84;

    // Parse TLEs from text using core MultiIterator
    var tles: std.ArrayList(astroz.Tle) = .empty;
    defer tles.deinit(allocator);

    var iter = astroz.Tle.MultiIterator.init(text);
    while (iter.next()) |pair| {
        const tle = astroz.Tle.parseLines(pair.line1, pair.line2, allocator) catch continue;

        // Filter: period <= 225 minutes (near-earth only, no deep space)
        const period = astroz.constants.minutesPerDay / tle.secondLine.mMotion;
        if (period <= astroz.constants.sgp4DeepSpaceThresholdMinutes) {
            tles.append(allocator, tle) catch {
                py.raiseRuntime("Out of memory");
                return null;
            };
        } else {
            var tle_mut = tle;
            tle_mut.deinit();
        }
    }

    if (tles.items.len == 0) {
        py.raiseValue("No valid near-earth TLEs found in text");
        return null;
    }

    const result = buildBatches(tles.items, grav) orelse return null;

    // Create new Sgp4ConstellationObject
    const self: *Sgp4ConstellationObject = @ptrCast(@alignCast(
        c.PyType_GenericAlloc(&Sgp4ConstellationType, 0) orelse {
            allocator.free(result.batches);
            allocator.free(result.epoch_jds);
            return null;
        },
    ));
    self.batches = result.batches.ptr;
    self.num_batches = result.batches.len;
    self.num_satellites = tles.items.len;
    self.epoch_jds = result.epoch_jds.ptr;
    self.padded_count = result.padded_count;

    return @ptrCast(self);
}

const constellation_methods = [_]c.PyMethodDef{
    .{
        .ml_name = "propagate_into",
        .ml_meth = @ptrCast(&constellation_propagate_into),
        .ml_flags = c.METH_VARARGS | c.METH_KEYWORDS,
        .ml_doc = "propagate_into(times, positions, velocities=None, *, epoch_offsets=None, output='teme', reference_jd=0.0) -> None\n\nPropagates all satellites. positions shape: (num_sats, num_times, 3)",
    },
    .{
        .ml_name = "from_tle_text",
        .ml_meth = @ptrCast(&constellation_from_tle_text),
        .ml_flags = c.METH_CLASS | c.METH_VARARGS | c.METH_KEYWORDS,
        .ml_doc = "from_tle_text(text, gravity_model=WGS84) -> Sgp4Constellation\n\nParse all TLEs from text in one call.",
    },
    .{
        .ml_name = null,
        .ml_meth = null,
        .ml_flags = 0,
        .ml_doc = null,
    },
};

const constellation_getset = [_]c.PyGetSetDef{
    .{
        .name = "num_satellites",
        .get = @ptrCast(&constellation_get_num_satellites),
        .set = null,
        .doc = "Number of satellites",
        .closure = null,
    },
    .{
        .name = "num_batches",
        .get = @ptrCast(&constellation_get_num_batches),
        .set = null,
        .doc = "Number of SIMD batches",
        .closure = null,
    },
    .{
        .name = "batch_size",
        .get = @ptrCast(&constellation_get_batch_size),
        .set = null,
        .doc = "SIMD batch size (4 for AVX2, 8 for AVX512)",
        .closure = null,
    },
    .{
        .name = "epochs",
        .get = @ptrCast(&constellation_get_epochs),
        .set = null,
        .doc = "Epoch Julian dates for each satellite",
        .closure = null,
    },
    .{
        .name = null,
        .get = null,
        .set = null,
        .doc = null,
        .closure = null,
    },
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

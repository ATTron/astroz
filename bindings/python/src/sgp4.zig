//! Sgp4Constellation: batched SGP4 propagation for many near-earth satellites.

const std = @import("std");
const py = @import("python.zig");
const c = py.c;
const shared = @import("shared.zig");
const astroz = @import("astroz");
const tle_mod = @import("tle.zig");
const satrec_mod = @import("satrec.zig");

const allocator = shared.allocator;
const BatchSize = shared.BatchSize;

pub const Sgp4ConstellationObject = shared.BatchObjectFields;

pub var Sgp4ConstellationType: c.PyTypeObject = undefined;

fn initConstellationType() void {
    Sgp4ConstellationType = std.mem.zeroes(c.PyTypeObject);
    Sgp4ConstellationType.tp_name = "astroz._astroz.Sgp4Constellation";
    Sgp4ConstellationType.tp_doc = "Constellation propagator for many satellites.\n\nArgs:\n    tles: list of Tle objects (will be grouped into batches)\n    gravity_model: WGS84 (default) or WGS72";
    Sgp4ConstellationType.tp_basicsize = @sizeOf(Sgp4ConstellationObject);
    Sgp4ConstellationType.tp_flags = c.Py_TPFLAGS_DEFAULT | c.Py_TPFLAGS_BASETYPE;
    Sgp4ConstellationType.tp_new = c.PyType_GenericNew;
    Sgp4ConstellationType.tp_init = @ptrCast(&constellation_init);
    Sgp4ConstellationType.tp_dealloc = @ptrCast(&constellation_dealloc);
    Sgp4ConstellationType.tp_methods = @constCast(&constellation_methods);
}

fn constellation_init(self_obj: [*c]c.PyObject, args: [*c]c.PyObject, kwds: [*c]c.PyObject) callconv(.c) c_int {
    const self: *Sgp4ConstellationObject = @ptrCast(@alignCast(self_obj));
    var tles_obj: [*c]c.PyObject = null;
    var grav_model: c_int = 0;

    const kwlist = [_:null]?[*:0]const u8{ "tles", "gravity_model", null };
    if (c.PyArg_ParseTupleAndKeywords(args, kwds, "O|i", @ptrCast(@constCast(&kwlist)), &tles_obj, &grav_model) == 0)
        return -1;

    if (c.PySequence_Check(tles_obj) == 0) {
        py.raiseType("First argument must be a sequence of Tle or Satrec objects");
        return -1;
    }

    const num_tles: usize = @intCast(c.PySequence_Size(tles_obj));
    if (num_tles == 0) {
        py.raiseValue("Must provide at least one Tle object");
        return -1;
    }

    shared.freeBatchFields(self);

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
        defer py.decref(item);

        const tle_ptr = if (c.PyObject_TypeCheck(item, &tle_mod.TleType) != 0)
            @as(*tle_mod.TleObject, @ptrCast(@alignCast(item))).tle
        else if (c.PyObject_TypeCheck(item, &satrec_mod.SatrecType) != 0)
            @as(*satrec_mod.SatrecObject, @ptrCast(@alignCast(item))).tle
        else {
            py.raiseType("All items must be Tle or Satrec objects");
            return -1;
        };
        tles[i] = (tle_ptr orelse {
            py.raiseValue("TLE not initialized");
            return -1;
        }).*;
    }

    const result = shared.buildBatches(tles, shared.getGravity(grav_model)) orelse return -1;
    shared.storeBatchResult(self, result, num_tles);
    return 0;
}

fn constellation_dealloc(self_obj: [*c]c.PyObject) callconv(.c) void {
    const self: *Sgp4ConstellationObject = @ptrCast(@alignCast(self_obj));
    shared.freeBatchFields(self);
    if (py.pyType(self_obj)) |tp| if (tp.*.tp_free) |free| free(@ptrCast(self));
}

fn constellation_propagate_into(self_obj: [*c]c.PyObject, args: [*c]c.PyObject, kwds: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    const self: *Sgp4ConstellationObject = @ptrCast(@alignCast(self_obj));

    var times_obj: [*c]c.PyObject = null;
    var pos_obj: [*c]c.PyObject = null;
    var vel_obj: [*c]c.PyObject = null;
    var offsets_obj: [*c]c.PyObject = null;
    var mask_obj: [*c]c.PyObject = null;
    var output_str: [*c]const u8 = null;
    var output_str_len: c.Py_ssize_t = 0;
    var reference_jd: f64 = 0.0;
    var time_major: c_int = 0;
    var output_stride: c_int = -1; // -1 = use num_satellites

    const kwlist = [_:null]?[*:0]const u8{
        "times", "positions", "velocities", "epoch_offsets", "satellite_mask", "output", "reference_jd", "time_major", "output_stride", null,
    };
    if (c.PyArg_ParseTupleAndKeywords(args, kwds, "OO|OOOs#dii", @ptrCast(@constCast(&kwlist)), &times_obj, &pos_obj, &vel_obj, &offsets_obj, &mask_obj, &output_str, &output_str_len, &reference_jd, &time_major, &output_stride) == 0)
        return null;

    const batches = (self.batches orelse {
        py.raiseRuntime("Constellation not initialized");
        return null;
    })[0..self.num_batches];

    const output_mode = shared.parseOutputMode(output_str, output_str_len) orelse return null;

    return shared.propagateInto(self, batches, times_obj, pos_obj, vel_obj, offsets_obj, mask_obj, .{
        .output = output_mode,
        .reference_jd = reference_jd,
        .layout = if (time_major != 0) .timeMajor else .satelliteMajor,
        .output_stride = if (output_stride > 0) @intCast(output_stride) else 0,
    });
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

    var tles: std.ArrayList(astroz.Tle) = .empty;
    defer deinitTles(&tles);

    var iter = astroz.Tle.MultiIterator.init(text_ptr[0..@intCast(text_len)]);
    while (iter.next()) |pair| {
        const tle = astroz.Tle.parseLines(pair.line1, pair.line2, allocator) catch continue;
        tles.append(allocator, tle) catch {
            var tle_mut = tle;
            tle_mut.deinit();
            py.raiseRuntime("Out of memory");
            return null;
        };
    }
    return constellationFromTles(&tles, grav_model);
}

fn deinitTles(tles: *std.ArrayList(astroz.Tle)) void {
    for (tles.items) |*t| t.deinit();
    tles.deinit(allocator);
}

/// Batches the near-earth TLEs; deep-space ones can't go into an SGP4 batch and are skipped.
fn constellationFromTles(tles: *std.ArrayList(astroz.Tle), grav_model: c_int) [*c]c.PyObject {
    var near: std.ArrayList(astroz.Tle) = .empty;
    defer near.deinit(allocator);
    for (tles.items) |tle| {
        const period = astroz.constants.minutesPerDay / tle.mMotion;
        if (period > astroz.constants.sgp4DeepSpaceThresholdMinutes) continue;
        near.append(allocator, tle) catch {
            py.raiseRuntime("Out of memory");
            return null;
        };
    }

    if (near.items.len == 0) {
        py.raiseValue("No valid near-earth TLEs found");
        return null;
    }

    const result = shared.buildBatches(near.items, shared.getGravity(grav_model)) orelse return null;

    const self: *Sgp4ConstellationObject = @ptrCast(@alignCast(
        c.PyType_GenericAlloc(&Sgp4ConstellationType, 0) orelse {
            var r = result;
            r.deinit();
            return null;
        },
    ));
    shared.storeBatchResult(self, result, near.items.len);

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
    .{ .ml_name = null, .ml_meth = null, .ml_flags = 0, .ml_doc = null },
};

const constellation_getset = [_]c.PyGetSetDef{
    .{ .name = "num_satellites", .get = @ptrCast(&constellation_get_num_satellites), .set = null, .doc = "Number of satellites", .closure = null },
    .{ .name = "num_batches", .get = @ptrCast(&constellation_get_num_batches), .set = null, .doc = "Number of SIMD batches", .closure = null },
    .{ .name = "batch_size", .get = @ptrCast(&constellation_get_batch_size), .set = null, .doc = "SIMD batch size", .closure = null },
    .{ .name = "epochs", .get = @ptrCast(&constellation_get_epochs), .set = null, .doc = "Epoch Julian dates for each satellite", .closure = null },
    .{ .name = null, .get = null, .set = null, .doc = null, .closure = null },
};

pub fn ready() c_int {
    initConstellationType();
    Sgp4ConstellationType.tp_getset = @constCast(&constellation_getset);
    return c.PyType_Ready(&Sgp4ConstellationType);
}

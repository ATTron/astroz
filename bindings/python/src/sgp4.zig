//! Python Sgp4Constellation type - multi-satellite batch propagator

const std = @import("std");
const py = @import("python.zig");
const c = py.c;
const shared = @import("shared.zig");
const astroz = @import("astroz");
const tle_mod = @import("tle.zig");

const allocator = shared.allocator;
const BatchSize = shared.BatchSize;
const BatchElements = shared.BatchElements;

// Sgp4Constellation: Multi-batch propagator for large constellations

pub const Sgp4ConstellationObject = shared.BatchObjectFields;

pub var Sgp4ConstellationType: c.PyTypeObject = undefined;

fn initConstellationType() void {
    Sgp4ConstellationType = std.mem.zeroes(c.PyTypeObject);
    Sgp4ConstellationType.tp_name = "astroz.Constellation";
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
    shared.initBatchFields(self);
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
    shared.freeBatchFields(self);

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

    const result = shared.buildBatches(tles, shared.getGravity(grav_model)) orelse return -1;
    shared.storeBatchResult(self, result, num_tles);
    return 0;
}

fn constellation_dealloc(self_obj: [*c]c.PyObject) callconv(.c) void {
    const self: *Sgp4ConstellationObject = @ptrCast(@alignCast(self_obj));
    shared.freeBatchFields(self);
    if (py.pyType(self_obj)) |tp| if (tp.*.tp_free) |free| free(@ptrCast(self));
}

fn parseOutputMode(output_str: [*c]const u8, output_str_len: c.Py_ssize_t) ?astroz.Constellation.OutputMode {
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

    // output_stride overrides num_satellites for the output layout stride
    const stride_sats: usize = if (output_stride > 0) @intCast(output_stride) else self.num_satellites;
    const required_size = stride_sats * num_times * 3 * @sizeOf(f64);

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
    const out_len = stride_sats * num_times * 3;

    astroz.Constellation.propagateConstellation(
        batches,
        stride_sats,
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
    const grav = shared.getGravity(grav_model);

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

    const result = shared.buildBatches(tles.items, grav) orelse return null;

    // Create new Sgp4ConstellationObject
    const self: *Sgp4ConstellationObject = @ptrCast(@alignCast(
        c.PyType_GenericAlloc(&Sgp4ConstellationType, 0) orelse {
            var r = result;
            r.deinit();
            return null;
        },
    ));
    shared.storeBatchResult(self, result, tles.items.len);

    return @ptrCast(self);
}

fn constellation_screen_conjunction(self_obj: [*c]c.PyObject, args: [*c]c.PyObject, kwds: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    const self: *Sgp4ConstellationObject = @ptrCast(@alignCast(self_obj));

    var times_obj: [*c]c.PyObject = null;
    var target_idx_long: c_long = 0;
    var threshold: f64 = 10.0;
    var offsets_obj: [*c]c.PyObject = null;
    var reference_jd: f64 = 0.0;

    const kwlist = [_:null]?[*:0]const u8{
        "times", "target_idx", "threshold", "epoch_offsets", "reference_jd", null,
    };
    if (c.PyArg_ParseTupleAndKeywords(args, kwds, "Ol|dOd", @ptrCast(@constCast(&kwlist)), &times_obj, &target_idx_long, &threshold, &offsets_obj, &reference_jd) == 0)
        return null;

    const batches = (self.batches orelse {
        py.raiseRuntime("Constellation not initialized");
        return null;
    })[0..self.num_batches];

    const target_idx: usize = @intCast(target_idx_long);
    if (target_idx >= self.num_satellites) {
        py.raiseValue("target_idx out of range");
        return null;
    }

    // Get times buffer
    var times_buf = std.mem.zeroes(c.Py_buffer);
    if (c.PyObject_GetBuffer(times_obj, &times_buf, c.PyBUF_SIMPLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&times_buf);

    const num_times = @as(usize, @intCast(times_buf.len)) / @sizeOf(f64);
    const times: [*]const f64 = @ptrCast(@alignCast(times_buf.buf));

    // Prepare epoch offsets
    var offsets_buf = std.mem.zeroes(c.Py_buffer);
    var offsets = prepareOffsets(offsets_obj, &offsets_buf, self.num_satellites, self.padded_count) orelse return null;
    defer offsets.deinit();
    defer if (!offsets.owned) c.PyBuffer_Release(&offsets_buf);

    // Get target's epoch offset
    const target_epoch_offset = offsets.data[target_idx];

    // Allocate output arrays
    const out_min_dists = allocator.alloc(f64, self.num_satellites) catch {
        py.raiseRuntime("Out of memory");
        return null;
    };
    defer allocator.free(out_min_dists);

    const out_min_t_indices = allocator.alloc(u32, self.num_satellites) catch {
        py.raiseRuntime("Out of memory");
        return null;
    };
    defer allocator.free(out_min_t_indices);

    // Run the kernel
    astroz.Constellation.screenConstellation(
        batches,
        self.num_satellites,
        times[0..num_times],
        offsets.data,
        target_idx,
        target_epoch_offset,
        threshold,
        reference_jd,
        out_min_dists,
        out_min_t_indices,
    ) catch {
        py.raiseValue("Screening failed");
        return null;
    };

    // Build Python result: (min_distances, min_t_indices)
    const dists_list = py.listFromF64(out_min_dists[0..self.num_satellites]) orelse return null;
    const indices_list = py.listFromU32(out_min_t_indices[0..self.num_satellites]) orelse {
        c.Py_DECREF(dists_list);
        return null;
    };

    const result = py.tuple(2) orelse {
        c.Py_DECREF(dists_list);
        c.Py_DECREF(indices_list);
        return null;
    };
    py.tupleSet(result, 0, dists_list);
    py.tupleSet(result, 1, indices_list);
    return result;
}

const constellation_methods = [_]c.PyMethodDef{
    .{
        .ml_name = "propagate_into",
        .ml_meth = @ptrCast(&constellation_propagate_into),
        .ml_flags = c.METH_VARARGS | c.METH_KEYWORDS,
        .ml_doc = "propagate_into(times, positions, velocities=None, *, epoch_offsets=None, output='teme', reference_jd=0.0) -> None\n\nPropagates all satellites. positions shape: (num_sats, num_times, 3)",
    },
    .{
        .ml_name = "screen_conjunction",
        .ml_meth = @ptrCast(&constellation_screen_conjunction),
        .ml_flags = c.METH_VARARGS | c.METH_KEYWORDS,
        .ml_doc = "screen_conjunction(times, target_idx, threshold=10.0, *, epoch_offsets=None, reference_jd=0.0) -> (min_distances, min_t_indices)\n\nFused propagate+screen: find minimum distance to target for each satellite.",
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
    .{ .name = "batch_size", .get = @ptrCast(&constellation_get_batch_size), .set = null, .doc = "SIMD batch size (4 for AVX2, 8 for AVX512)", .closure = null },
    .{ .name = "epochs", .get = @ptrCast(&constellation_get_epochs), .set = null, .doc = "Epoch Julian dates for each satellite", .closure = null },
    .{ .name = null, .get = null, .set = null, .doc = null, .closure = null },
};

pub fn ready() c_int {
    initConstellationType();
    Sgp4ConstellationType.tp_getset = @constCast(&constellation_getset);
    return c.PyType_Ready(&Sgp4ConstellationType);
}

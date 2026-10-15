//! Code shared by sgp4.zig (Sgp4Constellation) and satrec.zig (python-sgp4 API).

const std = @import("std");
const py = @import("python.zig");
const c = py.c;
const astroz = @import("astroz");
const Sgp4 = astroz.Sgp4;
const constants = astroz.constants;

pub const allocator = std.heap.c_allocator;

pub const BatchSize = astroz.Constellation.batchSize;

pub const BatchElements = astroz.Constellation.Sgp4Batch.BatchElements(BatchSize);

pub const WGS84: c_int = 0;
pub const WGS72: c_int = 1;

pub fn getGravity(grav_model: c_int) constants.Sgp4GravityModel {
    return if (grav_model == WGS72) constants.wgs72 else constants.wgs84;
}

pub fn sgp4ErrorMsg(e: Sgp4.Error) [:0]const u8 {
    return switch (e) {
        Sgp4.Error.DeepSpaceNotSupported => "Deep space not supported",
        Sgp4.Error.InvalidEccentricity => "Invalid eccentricity",
        Sgp4.Error.SatelliteDecayed => "Satellite decayed",
        Sgp4.Error.OutOfMemory => "Out of memory",
    };
}

/// python-sgp4 error codes
pub fn sgp4ErrorCode(e: Sgp4.Error) c_int {
    return switch (e) {
        Sgp4.Error.InvalidEccentricity => 1,
        Sgp4.Error.DeepSpaceNotSupported => 3,
        Sgp4.Error.OutOfMemory => 4,
        Sgp4.Error.SatelliteDecayed => 6,
    };
}

pub const BatchResult = struct {
    batches: []BatchElements,
    epoch_jds: []f64,
    padded_count: usize,

    pub fn deinit(self: *BatchResult) void {
        allocator.free(self.batches);
        allocator.free(self.epoch_jds);
    }
};

/// Group TLEs into SIMD batches of BatchSize, pad last batch, and extract epoch JDs.
pub fn buildBatches(tles: []const astroz.Tle, grav: constants.Sgp4GravityModel) ?BatchResult {
    const num_tles = tles.len;
    const padded_count = @divCeil(num_tles, BatchSize) * BatchSize;
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
        batches[batch_idx] = astroz.Constellation.Sgp4Batch.initBatchElements(BatchSize, batch_tles, grav) catch |e| {
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

/// Object layout shared by SatrecArray and Sgp4Constellation
pub const BatchObjectFields = extern struct {
    ob_base: c.PyObject,
    batches: ?[*]BatchElements,
    num_batches: usize,
    num_satellites: usize,
    epoch_jds: ?[*]f64,
    padded_count: usize,
};

pub fn freeBatchFields(self: *BatchObjectFields) void {
    if (self.batches) |batches| {
        allocator.free(batches[0..self.num_batches]);
    }
    if (self.epoch_jds) |ejds| {
        allocator.free(ejds[0..self.padded_count]);
    }
}

pub fn storeBatchResult(self: *BatchObjectFields, result: BatchResult, num_satellites: usize) void {
    self.batches = result.batches.ptr;
    self.num_batches = result.batches.len;
    self.num_satellites = num_satellites;
    self.epoch_jds = result.epoch_jds.ptr;
    self.padded_count = result.padded_count;
}

pub fn parseOutputMode(output_str: [*c]const u8, output_str_len: c.Py_ssize_t) ?astroz.Constellation.OutputMode {
    if (output_str == null or output_str_len <= 0) return .teme;
    const mode = output_str[0..@intCast(output_str_len)];
    if (std.mem.eql(u8, mode, "teme")) return .teme;
    if (std.mem.eql(u8, mode, "ecef")) return .ecef;
    if (std.mem.eql(u8, mode, "geodetic")) return .geodetic;
    py.raiseValue("output must be 'teme', 'ecef', or 'geodetic'");
    return null;
}

pub const PaddedArray = struct {
    data: []f64,
    owned: bool,

    pub fn deinit(self: *@This()) void {
        if (self.owned) allocator.free(self.data);
    }
};

/// Epoch offsets padded out to padded_count, or zeros if `obj` is None.
pub fn prepareOffsets(obj: [*c]c.PyObject, buf: *c.Py_buffer, num_sats: usize, padded_count: usize) ?PaddedArray {
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

/// Satellite mask padded out to padded_count, or null if `obj` is None.
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

pub const PropagateOptions = struct {
    output: astroz.Constellation.OutputMode = .teme,
    reference_jd: f64 = 0.0,
    layout: astroz.Constellation.Layout = .timeMajor,
    /// Satellite stride of the output arrays; 0 means num_satellites.
    output_stride: usize = 0,
};

/// Backs both SatrecArray.propagate_into and Sgp4Constellation.propagate_into.
pub fn propagateInto(
    self: *BatchObjectFields,
    batches: []const BatchElements,
    times_obj: [*c]c.PyObject,
    pos_obj: [*c]c.PyObject,
    vel_obj: [*c]c.PyObject,
    offsets_obj: [*c]c.PyObject,
    mask_obj: [*c]c.PyObject,
    opts: PropagateOptions,
) [*c]c.PyObject {
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
    const stride_sats = if (opts.output_stride > 0) opts.output_stride else self.num_satellites;
    const out_len = stride_sats * num_times * 3;

    if (@as(usize, @intCast(pos_buf.len)) < out_len * @sizeOf(f64)) {
        py.raiseValue("positions array too small");
        return null;
    }

    var have_vel = false;
    if (!py.isNone(vel_obj)) {
        if (c.PyObject_GetBuffer(vel_obj, &vel_buf, c.PyBUF_WRITABLE | c.PyBUF_FORMAT) < 0) return null;
        have_vel = true;
        if (@as(usize, @intCast(vel_buf.len)) < out_len * @sizeOf(f64)) {
            c.PyBuffer_Release(&vel_buf);
            py.raiseValue("velocities array too small");
            return null;
        }
    }
    defer if (have_vel) c.PyBuffer_Release(&vel_buf);

    var offsets = prepareOffsets(offsets_obj, &offsets_buf, self.num_satellites, self.padded_count) orelse return null;
    defer offsets.deinit();
    defer if (!offsets.owned) c.PyBuffer_Release(&offsets_buf);

    const mask_result = prepareMask(mask_obj, &mask_buf, self.num_satellites, self.padded_count);
    defer if (mask_result.owned) |m| allocator.free(m);
    defer if (mask_result.mask != null and mask_result.owned == null) c.PyBuffer_Release(&mask_buf);

    const times: [*]const f64 = @ptrCast(@alignCast(times_buf.buf));
    const pos_out: [*]f64 = @ptrCast(@alignCast(pos_buf.buf));

    astroz.Constellation.propagateConstellation(
        batches,
        stride_sats,
        times[0..num_times],
        offsets.data,
        pos_out[0..out_len],
        if (have_vel) @as([*]f64, @ptrCast(@alignCast(vel_buf.buf)))[0..out_len] else null,
        opts.output,
        opts.reference_jd,
        mask_result.mask,
        opts.layout,
    ) catch {
        py.raiseValue("Propagation failed");
        return null;
    };

    return py.none();
}

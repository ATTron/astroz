//! Code shared by sgp4.zig (Sgp4Constellation) and satrec.zig (python-sgp4 API).

const std = @import("std");
const py = @import("python.zig");
const c = py.c;
const astroz = @import("astroz");
const Sgp4 = astroz.Sgp4;
const constants = astroz.constants;

pub const allocator = std.heap.c_allocator;

/// Fills `buf` with a C-contiguous view of `obj` as a slice of `T` and returns it,
/// or raises (TypeError naming the `name` argument for the wrong dtype) and returns null. On success the
/// caller releases `buf`. `T` is f64 (float64), u32 (uint32), bool, or u8 (bool
/// or uint8, for masks).
pub fn getBuffer(comptime T: type, comptime name: []const u8, obj: [*c]c.PyObject, buf: *c.Py_buffer, writable: bool) ?[]T {
    const codes, const dtype = switch (T) {
        f64 => .{ "d", "float64" },
        u32 => .{ "IL", "uint32" },
        bool => .{ "?", "bool" },
        u8 => .{ "?B", "bool or uint8" },
        else => @compileError("unsupported buffer element type"),
    };
    const flags = c.PyBUF_FORMAT | c.PyBUF_C_CONTIGUOUS | @as(c_int, if (writable) c.PyBUF_WRITABLE else 0);
    if (c.PyObject_GetBuffer(obj, buf, flags) < 0) return null;
    if (buf.itemsize != @sizeOf(T) or !formatIs(buf.format, codes)) {
        c.PyBuffer_Release(buf);
        py.raiseType(name ++ " must be a C-contiguous " ++ dtype ++ " array");
        return null;
    }
    const ptr: [*]T = @ptrCast(@alignCast(buf.buf));
    return ptr[0..@intCast(@divExact(buf.len, @sizeOf(T)))];
}

/// Struct format string is one of `codes`, in native byte order.
fn formatIs(format: [*c]const u8, codes: []const u8) bool {
    var f: []const u8 = if (format == null) "B" else std.mem.span(format);
    const littleEndian = @import("builtin").target.cpu.arch.endian() == .little;
    if (f.len == 2 and (f[0] == '@' or f[0] == '=' or (f[0] == '<' and littleEndian))) f = f[1..];
    return f.len == 1 and std.mem.indexOfScalar(u8, codes, f[0]) != null;
}

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

/// Parses CCSDS OMM JSON, a single object or an array, into TLEs the caller owns.
/// Raises ValueError and returns null if the JSON is not valid OMM.
pub fn parseOmm(json: []const u8) ?[]astroz.Tle {
    const trimmed = std.mem.trimStart(u8, json, " \t\r\n");
    if (std.mem.startsWith(u8, trimmed, "{")) {
        const tles = allocator.alloc(astroz.Tle, 1) catch {
            py.raiseRuntime("Out of memory");
            return null;
        };
        tles[0] = astroz.Tle.parseOmm(trimmed, allocator) catch {
            allocator.free(tles);
            py.raiseValue("Failed to parse OMM JSON");
            return null;
        };
        return tles;
    }
    return astroz.Tle.parseOmmArray(trimmed, allocator) catch {
        py.raiseValue("Failed to parse OMM JSON");
        return null;
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

    const offsets = getBuffer(f64, "epoch_offsets", obj, buf, false) orelse return null;
    if (offsets.len >= padded_count) return .{ .data = offsets[0..padded_count], .owned = false };

    // Too short to use in place: copy into a padded array and let the buffer go
    defer c.PyBuffer_Release(buf);
    if (offsets.len < num_sats) {
        py.raiseValue("epoch_offsets must have at least num_satellites elements");
        return null;
    }
    const padded = allocator.alloc(f64, padded_count) catch {
        py.raiseRuntime("Out of memory");
        return null;
    };
    @memcpy(padded[0..num_sats], offsets[0..num_sats]);
    @memset(padded[num_sats..padded_count], 0.0);
    return .{ .data = padded, .owned = true };
}

/// Satellite mask (one byte per satellite) padded out to padded_count.
const Mask = struct {
    data: ?[]const u8 = null,
    /// Set when the mask was copied; otherwise `data` views the caller's buffer.
    owned: ?[]u8 = null,

    fn release(m: Mask, buf: *c.Py_buffer) void {
        if (m.owned) |o| allocator.free(o) else if (m.data != null) c.PyBuffer_Release(buf);
    }
};

/// Empty Mask if `obj` is None; null (with an exception set) on error.
fn prepareMask(obj: [*c]c.PyObject, buf: *c.Py_buffer, num_sats: usize, padded_count: usize) ?Mask {
    if (py.isNone(obj)) return .{};

    const mask = getBuffer(u8, "satellite_mask", obj, buf, false) orelse return null;
    if (mask.len >= padded_count) return .{ .data = mask[0..padded_count] };

    defer c.PyBuffer_Release(buf);
    if (mask.len < num_sats) {
        py.raiseValue("satellite_mask must have at least num_satellites elements");
        return null;
    }
    const padded = allocator.alloc(u8, padded_count) catch {
        py.raiseRuntime("Out of memory");
        return null;
    };
    @memcpy(padded[0..num_sats], mask[0..num_sats]);
    @memset(padded[num_sats..padded_count], 0);
    return .{ .data = padded, .owned = padded };
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

    const times = getBuffer(f64, "times", times_obj, &times_buf, false) orelse return null;
    defer c.PyBuffer_Release(&times_buf);

    const pos = getBuffer(f64, "positions", pos_obj, &pos_buf, true) orelse return null;
    defer c.PyBuffer_Release(&pos_buf);

    const stride_sats = if (opts.output_stride > 0) opts.output_stride else self.num_satellites;
    const out_len = stride_sats * times.len * 3;
    if (pos.len < out_len) {
        py.raiseValue("positions array too small");
        return null;
    }

    var vel: ?[]f64 = null;
    if (!py.isNone(vel_obj)) {
        vel = getBuffer(f64, "velocities", vel_obj, &vel_buf, true) orelse return null;
        if (vel.?.len < out_len) {
            c.PyBuffer_Release(&vel_buf);
            py.raiseValue("velocities array too small");
            return null;
        }
    }
    defer if (vel != null) c.PyBuffer_Release(&vel_buf);

    var offsets = prepareOffsets(offsets_obj, &offsets_buf, self.num_satellites, self.padded_count) orelse return null;
    defer offsets.deinit();
    defer if (!offsets.owned) c.PyBuffer_Release(&offsets_buf);

    const mask = prepareMask(mask_obj, &mask_buf, self.num_satellites, self.padded_count) orelse return null;
    defer mask.release(&mask_buf);

    astroz.Constellation.propagateConstellation(
        batches,
        stride_sats,
        times,
        offsets.data,
        pos[0..out_len],
        if (vel) |v| v[0..out_len] else null,
        opts.output,
        opts.reference_jd,
        mask.data,
        opts.layout,
    ) catch {
        py.raiseValue("Propagation failed");
        return null;
    };

    return py.none();
}

//! Shared utilities for Python bindings
//!
//! Common code used by both sgp4.zig (native API) and satrec.zig (python-sgp4 compatible API)

const std = @import("std");
const py = @import("python.zig");
const c = py.c;
const astroz = @import("astroz");
const Sgp4 = astroz.Sgp4;
const constants = astroz.constants;

pub const allocator = std.heap.c_allocator;

/// Batch size from core library (4 for AVX2, 8 for AVX512)
pub const BatchSize = astroz.Constellation.BatchSize;

/// Batch elements type for current batch size
pub const BatchElements = astroz.Constellation.Sgp4Batch.BatchElements(BatchSize);

// Gravity model constants for Python API
pub const WGS84: c_int = 0;
pub const WGS72: c_int = 1;

/// Convert Python gravity model constant to astroz gravity model
pub fn getGravity(grav_model: c_int) constants.Sgp4GravityModel {
    return if (grav_model == WGS72) constants.wgs72 else constants.wgs84;
}

/// SGP4 error to human-readable message
pub fn sgp4ErrorMsg(e: Sgp4.Error) [:0]const u8 {
    return switch (e) {
        Sgp4.Error.DeepSpaceNotSupported => "Deep space not supported",
        Sgp4.Error.InvalidEccentricity => "Invalid eccentricity",
        Sgp4.Error.SatelliteDecayed => "Satellite decayed",
        Sgp4.Error.OutOfMemory => "Out of memory",
    };
}

/// SGP4 error to python-sgp4 compatible error code
pub fn sgp4ErrorCode(e: Sgp4.Error) c_int {
    return switch (e) {
        Sgp4.Error.InvalidEccentricity => 1,
        Sgp4.Error.DeepSpaceNotSupported => 3,
        Sgp4.Error.OutOfMemory => 4,
        Sgp4.Error.SatelliteDecayed => 6,
    };
}

/// Result of building SIMD batches from TLEs
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

/// Common object fields for batch propagators (SatrecArray, Sgp4Constellation)
pub const BatchObjectFields = extern struct {
    ob_base: c.PyObject,
    batches: ?[*]BatchElements,
    num_batches: usize,
    num_satellites: usize,
    epoch_jds: ?[*]f64,
    padded_count: usize,
};

/// Initialize batch object fields to null/zero
pub fn initBatchFields(self: *BatchObjectFields) void {
    self.batches = null;
    self.num_batches = 0;
    self.num_satellites = 0;
    self.epoch_jds = null;
    self.padded_count = 0;
}

/// Free batch object resources
pub fn freeBatchFields(self: *BatchObjectFields) void {
    if (self.batches) |batches| {
        allocator.free(batches[0..self.num_batches]);
    }
    if (self.epoch_jds) |ejds| {
        allocator.free(ejds[0..self.padded_count]);
    }
}

/// Store batch result into object fields
pub fn storeBatchResult(self: *BatchObjectFields, result: BatchResult, num_satellites: usize) void {
    self.batches = result.batches.ptr;
    self.num_batches = result.batches.len;
    self.num_satellites = num_satellites;
    self.epoch_jds = result.epoch_jds.ptr;
    self.padded_count = result.padded_count;
}

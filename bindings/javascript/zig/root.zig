//! WebAssembly exports for the JavaScript bindings. A thin layer over the core
//! library: functions take plain numbers and pointers into linear memory.
//! Status codes must stay in sync with bindings/javascript/src/index.ts.
//!
//! Whole-catalog calls run the 8-wide SGP4/SDP4 batch kernels from
//! Constellation (single threaded, no runtime dispatch). Single-satellite
//! calls use the scalar propagator.

const std = @import("std");
const astroz = @import("astroz");
const Satellite = astroz.Satellite;
const Observer = astroz.Observer;
const Tle = astroz.Tle;
const Constellation = astroz.Constellation;
const Sgp4Batch = Constellation.Sgp4Batch;
const Sdp4Batch = Constellation.Sdp4Batch;
const Wcs = astroz.WorldCoordinateSystem;

const allocator = std.heap.wasm_allocator;
const lanes = Constellation.batchSize;

const Status = enum(u8) { ok = 0, badTle = 1, invalidEccentricity = 2, decayed = 3, outOfMemory = 4 };

fn toStatus(err: anyerror) Status {
    return switch (err) {
        error.InvalidEccentricity => .invalidEccentricity,
        error.SatelliteDecayed => .decayed,
        error.OutOfMemory => .outOfMemory,
        else => .badTle,
    };
}

const Frame = enum(u32) { teme = 0, ecef = 1, geodetic = 2 };

const Entry = struct { tle: Tle, sat: Satellite };

const Catalog = struct {
    grav: astroz.constants.Sgp4GravityModel,
    entries: std.ArrayList(Entry) = .empty,
    /// SIMD batches, rebuilt on the next whole-catalog call after an add
    batches: ?Constellation = null,

    fn batched(self: *Catalog) Status {
        if (self.batches == null and self.entries.items.len > 0) {
            const tles = allocator.alloc(Tle, self.entries.items.len) catch return .outOfMemory;
            defer allocator.free(tles);
            for (self.entries.items, tles) |e, *t| t.* = e.tle;
            self.batches = Constellation.init(allocator, tles, self.grav) catch |e| return toStatus(e);
        }
        return .ok;
    }

    fn invalidate(self: *Catalog) void {
        if (self.batches) |*b| b.deinit();
        self.batches = null;
    }
};

// 8-byte aligned so JS can view output buffers as Float64Array
const f64Align = std.mem.Alignment.of(f64);

export fn alloc(len: usize) ?[*]align(8) u8 {
    const buf = allocator.alignedAlloc(u8, f64Align, len) catch return null;
    return buf.ptr;
}

export fn free(ptr: [*]align(8) u8, len: usize) void {
    allocator.free(ptr[0..len]);
}

/// grav: 0 = WGS72, 1 = WGS84
export fn catalog_new(grav: u32) ?*Catalog {
    const cat = allocator.create(Catalog) catch return null;
    cat.* = .{ .grav = if (grav == 1) astroz.constants.wgs84 else astroz.constants.wgs72 };
    return cat;
}

export fn catalog_free(cat: *Catalog) void {
    cat.invalidate();
    for (cat.entries.items) |*e| e.tle.deinit();
    cat.entries.deinit(allocator);
    allocator.destroy(cat);
}

export fn catalog_len(cat: *const Catalog) usize {
    return cat.entries.items.len;
}

/// Parse a TLE and append it
export fn catalog_add(cat: *Catalog, line1: [*]const u8, len1: usize, line2: [*]const u8, len2: usize) Status {
    var tle = Tle.parseLines(line1[0..len1], line2[0..len2], allocator) catch |e| return toStatus(e);
    const sat = Satellite.init(tle, cat.grav) catch |e| {
        tle.deinit();
        return toStatus(e);
    };
    cat.entries.append(allocator, .{ .tle = tle, .sat = sat }) catch {
        tle.deinit();
        return .outOfMemory;
    };
    cat.invalidate();
    return .ok;
}

export fn catalog_norad_id(cat: *const Catalog, i: usize) u32 {
    return cat.entries.items[i].tle.satelliteNumber;
}

export fn catalog_epoch_jd(cat: *const Catalog, i: usize) f64 {
    return cat.entries.items[i].sat.epochJd();
}

export fn catalog_deep_space(cat: *const Catalog, i: usize) bool {
    return cat.entries.items[i].sat.isDeepSpace();
}

// ---------------------------------------------------------------------------
// Output writers. Each receives the TEME state of satellite `i`, or an error.

/// sin and cos of GMST at the requested time, shared by every satellite
const Earth = struct {
    sin: f64,
    cos: f64,

    fn at(jd: f64) Earth {
        const gmst = Wcs.julianToGmst(jd);
        return .{ .sin = @sin(gmst), .cos = @cos(gmst) };
    }

    fn toEcef(self: Earth, state: [2][3]f64) [2][3]f64 {
        return Wcs.temeToEcefSinCos(state[0], state[1], self.sin, self.cos);
    }
};

const StateWriter = struct {
    frame: Frame,
    earth: Earth,
    pos: [*]f64,
    vel: [*]f64,
    status: [*]u8,

    fn write(self: StateWriter, i: usize, result: anyerror![2][3]f64) void {
        const state = result catch |e| {
            @memset(self.pos[i * 3 ..][0..3], std.math.nan(f64));
            @memset(self.vel[i * 3 ..][0..3], std.math.nan(f64));
            self.status[i] = @intFromEnum(toStatus(e));
            return;
        };
        const out = switch (self.frame) {
            .teme => state,
            .ecef => self.earth.toEcef(state),
            .geodetic => blk: {
                const ecef = self.earth.toEcef(state);
                break :blk [2][3]f64{ Wcs.ecefToGeodeticDeg(ecef[0]), ecef[1] };
            },
        };
        self.pos[i * 3 ..][0..3].* = out[0];
        self.vel[i * 3 ..][0..3].* = out[1];
        self.status[i] = @intFromEnum(Status.ok);
    }
};

const LookWriter = struct {
    obs: Observer,
    earth: Earth,
    out: [*]f64,
    status: [*]u8,

    fn write(self: LookWriter, i: usize, result: anyerror![2][3]f64) void {
        const state = result catch |e| {
            @memset(self.out[i * 4 ..][0..4], std.math.nan(f64));
            self.status[i] = @intFromEnum(toStatus(e));
            return;
        };
        const ecef = self.earth.toEcef(state);
        const la = self.obs.lookAngles(ecef[0], ecef[1]);
        self.out[i * 4 ..][0..4].* = .{ la.azimuth, la.elevation, la.range, la.rangeRate };
        self.status[i] = @intFromEnum(Status.ok);
    }
};

fn propagateOne(cat: *const Catalog, i: usize, jd: f64) anyerror![2][3]f64 {
    const sat = &cat.entries.items[i].sat;
    return sat.propagate((jd - sat.epochJd()) * 1440.0);
}

/// Propagate every satellite to `jd` through the batch kernels, passing each
/// result to `writer` with its catalog index.
fn propagateAll(cat: *Catalog, jd: f64, writer: anytype) Status {
    const s = cat.batched();
    if (s != .ok) return s;
    const b = &(cat.batches orelse return .ok);

    const tBase = (jd - b.referenceEpochJd) * 1440.0;
    for (b.sgp4Batches, 0..) |*el, bi| {
        const base = bi * lanes;
        var t: [lanes]f64 = undefined;
        for (&t, 0..) |*v, lane| v.* = tBase + b.sgp4EpochOffsets[base + lane];
        writeBatch(cat, jd, writer, b.sgp4OrigIndices[base..][0..lanes], b.numSgp4 - base, Sgp4Batch.propagateBatchDirect(lanes, el, t));
    }
    for (b.sdp4Batches, b.sdp4BatchEpochs, b.sdp4Carries, 0..) |*el, epochs, *carry, bi| {
        const base = bi * lanes;
        var t: [lanes]f64 = undefined;
        for (&t, epochs) |*v, epoch| v.* = (jd - epoch) * 1440.0;
        writeBatch(cat, jd, writer, b.sdp4OrigIndices[base..][0..lanes], b.numSdp4 - base, Sdp4Batch.propagateBatchDirect(lanes, el, t, carry));
    }
    return .ok;
}

fn writeBatch(cat: *const Catalog, jd: f64, writer: anytype, indices: []const u32, remaining: usize, result: anyerror!astroz.Sgp4.PositionVelocity(lanes)) void {
    const real = @min(remaining, lanes);
    const pv = result catch {
        // one lane failed: redo the batch per satellite so the others still get results
        for (indices[0..real]) |i| writer.write(i, propagateOne(cat, i, jd));
        return;
    };
    const rx: [lanes]f64 = pv.rx;
    const ry: [lanes]f64 = pv.ry;
    const rz: [lanes]f64 = pv.rz;
    const vx: [lanes]f64 = pv.vx;
    const vy: [lanes]f64 = pv.vy;
    const vz: [lanes]f64 = pv.vz;
    for (indices[0..real], 0..) |i, lane| {
        // the batch kernels can leave a failed lane as NaN without an error;
        // the scalar propagator reports why
        if (!std.math.isFinite(rx[lane])) {
            writer.write(i, propagateOne(cat, i, jd));
            continue;
        }
        writer.write(i, .{ .{ rx[lane], ry[lane], rz[lane] }, .{ vx[lane], vy[lane], vz[lane] } });
    }
}

/// Propagate every satellite to Julian date `jd`. Writes 3 values per
/// satellite to `pos` and `vel` and a status byte to `status`. Geodetic
/// positions are (lat deg, lon deg, alt km) with ECEF velocity. Failed
/// satellites get NaN outputs. Returns a catalog-level status.
export fn catalog_propagate(cat: *Catalog, jd: f64, frame: Frame, pos: [*]f64, vel: [*]f64, status: [*]u8) Status {
    return propagateAll(cat, jd, StateWriter{ .frame = frame, .earth = Earth.at(jd), .pos = pos, .vel = vel, .status = status });
}

/// Single satellite `i`; output layout as catalog_propagate
export fn catalog_propagate_one(cat: *const Catalog, i: usize, jd: f64, frame: Frame, pos: [*]f64, vel: [*]f64, status: [*]u8) void {
    const w = StateWriter{ .frame = frame, .earth = Earth.at(jd), .pos = pos, .vel = vel, .status = status };
    w.write(0, propagateOne(cat, i, jd));
}

/// Look angles from an observer (deg, deg, km) to every satellite at Julian
/// date `jd`. Writes (azimuth deg, elevation deg, range km, range rate km/s)
/// per satellite to `out`.
export fn catalog_look_angles(cat: *Catalog, jd: f64, lat: f64, lon: f64, alt: f64, out: [*]f64, status: [*]u8) Status {
    return propagateAll(cat, jd, LookWriter{ .obs = Observer.init(lat, lon, alt), .earth = Earth.at(jd), .out = out, .status = status });
}

/// Single satellite `i`; output layout as catalog_look_angles
export fn catalog_look_angles_one(cat: *const Catalog, i: usize, jd: f64, lat: f64, lon: f64, alt: f64, out: [*]f64, status: [*]u8) void {
    const w = LookWriter{ .obs = Observer.init(lat, lon, alt), .earth = Earth.at(jd), .out = out, .status = status };
    w.write(0, propagateOne(cat, i, jd));
}

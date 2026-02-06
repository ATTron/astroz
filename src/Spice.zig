//! CSPICE bindings for NAIF SPICE toolkit
//!
//! Requires CSPICE to be installed. Build with:
//!   zig build -Denable-cspice=true -Dcspice-include=/path/to/cspice/include -Dcspice-lib=/path/to/cspice/lib
//!
//! Or if CSPICE is in standard system paths:
//!   zig build -Denable-cspice=true

const std = @import("std");
const buildOptions = @import("build_options");

/// Whether CSPICE support is enabled in this build
pub const enabled = buildOptions.enable_cspice;

/// Raw CSPICE C bindings (only available when enabled)
pub const c = if (enabled) @cImport({
    @cInclude("SpiceUsr.h");
}) else struct {
    pub const SpiceInt = c_int;
    pub const SpiceBoolean = c_int;
};

/// SPICE error information
pub const SpiceError = error{
    SpiceError,
    KernelNotFound,
    InsufficientData,
    InvalidTime,
    InvalidFrame,
    InvalidBody,
    NotEnabled,
};

/// 3D position vector (km)
pub const Position = struct {
    x: f64,
    y: f64,
    z: f64,

    pub fn array(self: Position) [3]f64 {
        return .{ self.x, self.y, self.z };
    }

    pub fn magnitude(self: Position) f64 {
        return @sqrt(self.x * self.x + self.y * self.y + self.z * self.z);
    }
};

/// 3D state vector (position in km, velocity in km/s)
pub const State = struct {
    position: Position,
    velocity: Position,
    lightTime: f64, // one-way light time (seconds)

    pub fn posArray(self: State) [3]f64 {
        return self.position.array();
    }

    pub fn velArray(self: State) [3]f64 {
        return self.velocity.array();
    }

    pub fn stateArray(self: State) [6]f64 {
        return .{
            self.position.x, self.position.y, self.position.z,
            self.velocity.x, self.velocity.y, self.velocity.z,
        };
    }
};

/// Common NAIF body IDs
pub const Body = struct {
    pub const sun: c_int = 10;
    pub const mercury: c_int = 199;
    pub const venus: c_int = 299;
    pub const earth: c_int = 399;
    pub const moon: c_int = 301;
    pub const mars: c_int = 499;
    pub const jupiter: c_int = 599;
    pub const saturn: c_int = 699;
    pub const uranus: c_int = 799;
    pub const neptune: c_int = 899;
    pub const pluto: c_int = 999;
    pub const earthBarycenter: c_int = 3;
    pub const solarSystemBarycenter: c_int = 0;
};

/// Common reference frames
pub const Frame = struct {
    pub const j2000 = "J2000";
    pub const eclipj2000 = "ECLIPJ2000";
    pub const iauEarth = "IAU_EARTH";
    pub const itrf93 = "ITRF93";
    pub const galactic = "GALACTIC";
};

/// Aberration correction types
pub const AberrationCorrection = enum {
    none,
    lightTime,
    lightTimeStellar,
    convergedNewtonian,
    convergedNewtonianStellar,

    pub fn toStr(self: AberrationCorrection) [*:0]const u8 {
        return switch (self) {
            .none => "NONE",
            .lightTime => "LT",
            .lightTimeStellar => "LT+S",
            .convergedNewtonian => "CN",
            .convergedNewtonianStellar => "CN+S",
        };
    }
};

/// Default directory for downloaded SPICE kernels (relative to project root)
pub const default_kernel_dir = "data/kernels/";

/// Standard kernel filenames downloaded by `zig build fetch-kernels`
pub const default_kernels = .{
    "naif0012.tls",
    "de440s.bsp",
    "pck00011.tpc",
    "gm_de440.tpc",
};

/// Load all standard kernels from the default directory (`data/kernels/`)
pub fn loadDefaultKernels() SpiceError!void {
    inline for (default_kernels) |kernel| {
        try loadKernel(default_kernel_dir ++ kernel);
    }
}

/// Load all standard kernels from a custom directory
pub fn loadDefaultKernelsFrom(dir: []const u8) SpiceError!void {
    var buf: [1024]u8 = undefined;
    inline for (default_kernels) |kernel| {
        if (dir.len + kernel.len >= buf.len) return SpiceError.KernelNotFound;
        @memcpy(buf[0..dir.len], dir);
        @memcpy(buf[dir.len .. dir.len + kernel.len], kernel);
        try loadKernel(buf[0 .. dir.len + kernel.len]);
    }
}

/// Load a SPICE kernel
pub fn loadKernel(path: []const u8) SpiceError!void {
    if (comptime !enabled) return SpiceError.NotEnabled;

    var buf: [1024]u8 = undefined;
    const c_path = toCString(path, &buf) orelse return SpiceError.KernelNotFound;

    c.furnsh_c(c_path);

    if (failed()) {
        reset();
        return SpiceError.KernelNotFound;
    }
}

/// Unload a specific kernel
pub fn unloadKernel(path: []const u8) void {
    if (comptime !enabled) return;

    var buf: [1024]u8 = undefined;
    const c_path = toCString(path, &buf) orelse return;
    c.unload_c(c_path);
    reset();
}

/// Unload all loaded kernels
pub fn clearKernels() void {
    if (comptime !enabled) return;
    c.kclear_c();
    reset();
}

/// Get the number of loaded kernels
pub fn kernelCount() usize {
    if (comptime !enabled) return 0;
    var count: c.SpiceInt = 0;
    c.ktotal_c("ALL", &count);
    return @intCast(count);
}

/// Convert UTC time string to Ephemeris Time (ET/TDB)
/// Accepts formats like: "2024-01-15T12:00:00", "2024-01-15 12:00:00 UTC", etc
pub fn utcToEt(utc: []const u8) SpiceError!f64 {
    if (comptime !enabled) return SpiceError.NotEnabled;

    var buf: [256]u8 = undefined;
    const cUtc = toCString(utc, &buf) orelse return SpiceError.InvalidTime;

    var et: f64 = 0;
    c.str2et_c(cUtc, &et);

    if (failed()) {
        reset();
        return SpiceError.InvalidTime;
    }
    return et;
}

/// Convert Ephemeris Time to UTC string
pub fn etToUtc(et: f64, precision: u8, buf: []u8) SpiceError![]const u8 {
    if (comptime !enabled) return SpiceError.NotEnabled;

    const len: c.SpiceInt = @intCast(buf.len);
    c.et2utc_c(et, "ISOC", @intCast(precision), len, buf.ptr);

    if (failed()) {
        reset();
        return SpiceError.InvalidTime;
    }

    const strLen = std.mem.indexOfScalar(u8, buf, 0) orelse buf.len;
    return buf[0..strLen];
}

/// Convert Julian Date (TDB) to Ephemeris Time (no CSPICE required)
pub fn jdToEt(jd: f64) f64 {
    return (jd - 2451545.0) * 86400.0;
}

/// Convert Ephemeris Time to Julian Date (TDB) (no CSPICE required)
pub fn etToJd(et: f64) f64 {
    return et / 86400.0 + 2451545.0;
}

/// Get the position of a target body relative to an observer
pub fn getPosition(
    target: []const u8,
    et: f64,
    frame: []const u8,
    abcorr: AberrationCorrection,
    observer: []const u8,
) SpiceError!Position {
    if (comptime !enabled) return SpiceError.NotEnabled;

    var targetBuf: [64]u8 = undefined;
    var frameBuf: [64]u8 = undefined;
    var observerBuf: [64]u8 = undefined;

    const cTarget = toCString(target, &targetBuf) orelse return SpiceError.InvalidBody;
    const cFrame = toCString(frame, &frameBuf) orelse return SpiceError.InvalidFrame;
    const cObserver = toCString(observer, &observerBuf) orelse return SpiceError.InvalidBody;

    var pos: [3]f64 = undefined;
    var lt: f64 = 0;

    c.spkpos_c(cTarget, et, cFrame, abcorr.toStr(), cObserver, &pos, &lt);

    if (failed()) {
        reset();
        return SpiceError.InsufficientData;
    }

    return Position{ .x = pos[0], .y = pos[1], .z = pos[2] };
}

/// Get the state (position + velocity) of a target body relative to an observer
pub fn getState(
    target: []const u8,
    et: f64,
    frame: []const u8,
    abcorr: AberrationCorrection,
    observer: []const u8,
) SpiceError!State {
    if (comptime !enabled) return SpiceError.NotEnabled;

    var targetBuf: [64]u8 = undefined;
    var frameBuf: [64]u8 = undefined;
    var observerBuf: [64]u8 = undefined;

    const cTarget = toCString(target, &targetBuf) orelse return SpiceError.InvalidBody;
    const cFrame = toCString(frame, &frameBuf) orelse return SpiceError.InvalidFrame;
    const cObserver = toCString(observer, &observerBuf) orelse return SpiceError.InvalidBody;

    var state: [6]f64 = undefined;
    var lt: f64 = 0;

    c.spkezr_c(cTarget, et, cFrame, abcorr.toStr(), cObserver, &state, &lt);

    if (failed()) {
        reset();
        return SpiceError.InsufficientData;
    }

    return State{
        .position = Position{ .x = state[0], .y = state[1], .z = state[2] },
        .velocity = Position{ .x = state[3], .y = state[4], .z = state[5] },
        .lightTime = lt,
    };
}

/// Get geometric position using NAIF integer IDs instead of names.
/// Note: spkgps_c always returns geometric (no aberration correction) positions.
pub fn getPositionById(
    target: c_int,
    et: f64,
    frame: []const u8,
    observer: c_int,
) SpiceError!Position {
    if (comptime !enabled) return SpiceError.NotEnabled;

    var frameBuf: [64]u8 = undefined;
    const cFrame = toCString(frame, &frameBuf) orelse return SpiceError.InvalidFrame;

    var pos: [3]f64 = undefined;
    var lt: f64 = 0;

    c.spkgps_c(target, et, cFrame, observer, &pos, &lt);

    if (failed()) {
        reset();
        return SpiceError.InsufficientData;
    }

    return Position{ .x = pos[0], .y = pos[1], .z = pos[2] };
}

/// Get the rotation matrix from one frame to another at a given time
pub fn getFrameRotation(from: []const u8, to: []const u8, et: f64) SpiceError![3][3]f64 {
    if (comptime !enabled) return SpiceError.NotEnabled;

    var fromBuf: [64]u8 = undefined;
    var toBuf: [64]u8 = undefined;

    const cFrom = toCString(from, &fromBuf) orelse return SpiceError.InvalidFrame;
    const cTo = toCString(to, &toBuf) orelse return SpiceError.InvalidFrame;

    var matrix: [3][3]f64 = undefined;
    c.pxform_c(cFrom, cTo, et, &matrix);

    if (failed()) {
        reset();
        return SpiceError.InvalidFrame;
    }

    return matrix;
}

/// Transform a position vector from one frame to another
pub fn transformPosition(pos: Position, from: []const u8, to: []const u8, et: f64) SpiceError!Position {
    const matrix = try getFrameRotation(from, to, et);

    const p = pos.array();
    return Position{
        .x = matrix[0][0] * p[0] + matrix[0][1] * p[1] + matrix[0][2] * p[2],
        .y = matrix[1][0] * p[0] + matrix[1][1] * p[1] + matrix[1][2] * p[2],
        .z = matrix[2][0] * p[0] + matrix[2][1] * p[1] + matrix[2][2] * p[2],
    };
}

/// Get the NAIF ID for a body name
pub fn getBodyId(name: []const u8) SpiceError!c_int {
    if (comptime !enabled) return SpiceError.NotEnabled;

    var buf: [64]u8 = undefined;
    const cName = toCString(name, &buf) orelse return SpiceError.InvalidBody;

    var id: c.SpiceInt = 0;
    var found: c.SpiceBoolean = 0;

    c.bodn2c_c(cName, &id, &found);

    if (found == 0) {
        return SpiceError.InvalidBody;
    }

    return id;
}

/// Get the name for a NAIF body ID
pub fn getBodyName(id: c_int, buf: []u8) SpiceError![]const u8 {
    if (comptime !enabled) return SpiceError.NotEnabled;

    var found: c.SpiceBoolean = 0;
    c.bodc2n_c(id, @intCast(buf.len), buf.ptr, &found);

    if (found == 0) {
        return SpiceError.InvalidBody;
    }

    const strLen = std.mem.indexOfScalar(u8, buf, 0) orelse buf.len;
    return buf[0..strLen];
}

fn failed() bool {
    if (comptime !enabled) return false;
    return c.failed_c() != 0;
}

fn reset() void {
    if (comptime !enabled) return;
    c.reset_c();
}

fn toCString(slice: []const u8, buf: []u8) ?[*:0]const u8 {
    if (slice.len >= buf.len) return null;
    @memcpy(buf[0..slice.len], slice);
    buf[slice.len] = 0;
    return @ptrCast(buf[0 .. slice.len + 1].ptr);
}

/// Get Sun position relative to Earth in J2000 frame
pub fn getSunPosition(et: f64) SpiceError!Position {
    return getPosition("SUN", et, Frame.j2000, .lightTimeStellar, "EARTH");
}

/// Get Moon position relative to Earth in J2000 frame
pub fn getMoonPosition(et: f64) SpiceError!Position {
    return getPosition("MOON", et, Frame.j2000, .lightTimeStellar, "EARTH");
}

/// Get position of a planet relative to Earth in J2000 frame
pub fn getPlanetPosition(planet: []const u8, et: f64) SpiceError!Position {
    return getPosition(planet, et, Frame.j2000, .lightTimeStellar, "EARTH");
}

test "jd to et conversion" {
    const et = jdToEt(2451545.0);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), et, 1e-10);

    const et1day = jdToEt(2451546.0);
    try std.testing.expectApproxEqAbs(@as(f64, 86400.0), et1day, 1e-10);
}

test "et to jd conversion" {
    const jd = etToJd(0.0);
    try std.testing.expectApproxEqAbs(@as(f64, 2451545.0), jd, 1e-10);

    const jd1day = etToJd(86400.0);
    try std.testing.expectApproxEqAbs(@as(f64, 2451546.0), jd1day, 1e-10);
}

test "position struct" {
    const pos = Position{ .x = 3.0, .y = 4.0, .z = 0.0 };
    try std.testing.expectApproxEqAbs(@as(f64, 5.0), pos.magnitude(), 1e-10);

    const arr = pos.array();
    try std.testing.expectEqual(@as(f64, 3.0), arr[0]);
    try std.testing.expectEqual(@as(f64, 4.0), arr[1]);
    try std.testing.expectEqual(@as(f64, 0.0), arr[2]);
}

test "state struct" {
    const state = State{
        .position = Position{ .x = 1.0, .y = 2.0, .z = 3.0 },
        .velocity = Position{ .x = 4.0, .y = 5.0, .z = 6.0 },
        .lightTime = 0.5,
    };

    const arr = state.stateArray();
    try std.testing.expectEqual(@as(f64, 1.0), arr[0]);
    try std.testing.expectEqual(@as(f64, 6.0), arr[5]);
}

test "enabled flag" {
    _ = enabled;
}

test "not enabled returns error" {
    if (!enabled) {
        try std.testing.expectError(SpiceError.NotEnabled, loadKernel("test.bsp"));
        try std.testing.expectError(SpiceError.NotEnabled, utcToEt("2024-01-01"));
        try std.testing.expectError(SpiceError.NotEnabled, getPosition("SUN", 0, "J2000", .none, "EARTH"));
    }
}

test "cspice body name/id lookups" {
    if (!enabled) return;

    // Built-in mappings work without loading any kernels
    const earthId = try getBodyId("EARTH");
    try std.testing.expectEqual(@as(c_int, 399), earthId);

    const moonId = try getBodyId("MOON");
    try std.testing.expectEqual(@as(c_int, 301), moonId);

    const sunId = try getBodyId("SUN");
    try std.testing.expectEqual(@as(c_int, 10), sunId);

    // Reverse lookup
    var nameBuf: [64]u8 = undefined;
    const name = try getBodyName(399, &nameBuf);
    try std.testing.expectEqualStrings("EARTH", name);
}

test "cspice kernel management" {
    if (!enabled) return;

    clearKernels();
    try std.testing.expectEqual(@as(usize, 0), kernelCount());
}

test "loadDefaultKernels returns NotEnabled when cspice is off" {
    if (!enabled) {
        try std.testing.expectError(SpiceError.NotEnabled, loadDefaultKernels());
        try std.testing.expectError(SpiceError.NotEnabled, loadDefaultKernelsFrom("/tmp/kernels/"));
    }
}

test "default kernel paths fit in buffer" {
    inline for (default_kernels) |kernel| {
        const path = default_kernel_dir ++ kernel;
        try std.testing.expect(path.len < 1024);
    }
}

//! Python wrapper for astroz.Screening.

const std = @import("std");
const py = @import("python.zig");
const c = py.c;
const astroz = @import("astroz");
const shared = @import("shared.zig");

const allocator = shared.allocator;

/// screen_candidates(positions, velocities, num_sats, dt, threshold, primary=None) -> bytes
pub fn pyScreenCandidates(_: [*c]c.PyObject, args: [*c]c.PyObject, kwds: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    var pos_obj: [*c]c.PyObject = null;
    var vel_obj: [*c]c.PyObject = null;
    var num_sats: c.Py_ssize_t = 0;
    var dt: f64 = 0;
    var threshold: f64 = 0;
    var primary_obj: [*c]c.PyObject = null;

    const kwlist = [_:null]?[*:0]const u8{ "positions", "velocities", "num_sats", "dt", "threshold", "primary", null };
    if (c.PyArg_ParseTupleAndKeywords(args, kwds, "OOndd|O", @ptrCast(@constCast(&kwlist)), &pos_obj, &vel_obj, &num_sats, &dt, &threshold, &primary_obj) == 0)
        return null;
    if (num_sats < 0 or !(dt > 0) or !(threshold >= 0)) {
        py.raiseValue("num_sats must be >= 0, dt > 0 and threshold >= 0");
        return null;
    }
    const n: usize = @intCast(num_sats);

    var pos_buf = std.mem.zeroes(c.Py_buffer);
    if (c.PyObject_GetBuffer(pos_obj, &pos_buf, c.PyBUF_C_CONTIGUOUS) < 0) return null;
    defer c.PyBuffer_Release(&pos_buf);
    var vel_buf = std.mem.zeroes(c.Py_buffer);
    if (c.PyObject_GetBuffer(vel_obj, &vel_buf, c.PyBUF_C_CONTIGUOUS) < 0) return null;
    defer c.PyBuffer_Release(&vel_buf);

    const row = @sizeOf([3]f64) * @max(n, 1);
    const len: usize = @intCast(pos_buf.len);
    if (len != @as(usize, @intCast(vel_buf.len)) or len % row != 0) {
        py.raiseValue("positions and velocities must both be float64 arrays of shape (n_times, num_sats, 3)");
        return null;
    }
    const states = len / @sizeOf([3]f64);
    const pos: []const [3]f64 = @as([*]const [3]f64, @ptrCast(@alignCast(pos_buf.buf)))[0..states];
    const vel: []const [3]f64 = @as([*]const [3]f64, @ptrCast(@alignCast(vel_buf.buf)))[0..states];

    var primary_buf = std.mem.zeroes(c.Py_buffer);
    var primary: ?[]const bool = null;
    if (primary_obj != null and !py.isNone(primary_obj)) {
        if (c.PyObject_GetBuffer(primary_obj, &primary_buf, c.PyBUF_C_CONTIGUOUS) < 0) return null;
        if (@as(usize, @intCast(primary_buf.len)) != n) {
            c.PyBuffer_Release(&primary_buf);
            py.raiseValue("primary must be a bool array of length num_sats");
            return null;
        }
        primary = @as([*]const bool, @ptrCast(primary_buf.buf))[0..n];
    }
    defer if (primary != null) c.PyBuffer_Release(&primary_buf);

    // The sweep can take a while on a full catalog; let other Python threads run.
    const save = c.PyEval_SaveThread();
    const result = astroz.Screening.candidates(allocator, pos, vel, n, .{ .dt = dt, .threshold = threshold, .primary = primary });
    c.PyEval_RestoreThread(save);
    const found = result catch {
        py.raiseRuntime("Out of memory");
        return null;
    };
    defer allocator.free(found);

    const out = c.PyBytes_FromStringAndSize(null, @intCast(found.len * 3 * @sizeOf(u32))) orelse return null;
    const words: [*]u32 = @ptrCast(@alignCast(c.PyBytes_AsString(out)));
    for (found, 0..) |cand, k| {
        words[3 * k] = cand.i;
        words[3 * k + 1] = cand.j;
        words[3 * k + 2] = cand.step;
    }
    return out;
}

//! Python wrapper for astroz.Screening.

const std = @import("std");
const py = @import("python.zig");
const c = py.c;
const astroz = @import("astroz");
const shared = @import("shared.zig");
const satrec_mod = @import("satrec.zig");

const allocator = shared.allocator;

/// Gets a C-contiguous buffer whose items are `itemsize` bytes with one of the
/// struct format codes in `codes` (byte-order prefixes are allowed). Raises
/// TypeError naming `what` otherwise.
fn getTypedBuffer(obj: [*c]c.PyObject, buf: *c.Py_buffer, comptime codes: []const u8, itemsize: usize, comptime what: [:0]const u8) bool {
    if (c.PyObject_GetBuffer(obj, buf, c.PyBUF_C_CONTIGUOUS | c.PyBUF_FORMAT) < 0) return false;
    const format: []const u8 = if (buf.format) |f| std.mem.span(f) else "B";
    const code = if (format.len > 0) format[format.len - 1] else 0;
    if (@as(usize, @intCast(buf.itemsize)) != itemsize or std.mem.indexOfScalar(u8, codes, code) == null) {
        c.PyBuffer_Release(buf);
        py.raiseType(what);
        return false;
    }
    return true;
}

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
    if (!getTypedBuffer(pos_obj, &pos_buf, "d", 8, "positions must be a contiguous float64 array")) return null;
    defer c.PyBuffer_Release(&pos_buf);
    var vel_buf = std.mem.zeroes(c.Py_buffer);
    if (!getTypedBuffer(vel_obj, &vel_buf, "d", 8, "velocities must be a contiguous float64 array")) return null;
    defer c.PyBuffer_Release(&vel_buf);

    const row = @sizeOf([3]f64) * @max(n, 1);
    const len: usize = @intCast(pos_buf.len);
    if (len != @as(usize, @intCast(vel_buf.len)) or len % row != 0) {
        py.raiseValue("positions and velocities must both have shape (n_times, num_sats, 3)");
        return null;
    }
    const states = len / @sizeOf([3]f64);
    const pos: []const [3]f64 = @as([*]const [3]f64, @ptrCast(@alignCast(pos_buf.buf)))[0..states];
    const vel: []const [3]f64 = @as([*]const [3]f64, @ptrCast(@alignCast(vel_buf.buf)))[0..states];

    var primary_buf = std.mem.zeroes(c.Py_buffer);
    var primary: ?[]const bool = null;
    if (primary_obj != null and !py.isNone(primary_obj)) {
        if (!getTypedBuffer(primary_obj, &primary_buf, "?", 1, "primary must be a contiguous bool array")) return null;
        if (@as(usize, @intCast(primary_buf.len)) != n) {
            c.PyBuffer_Release(&primary_buf);
            py.raiseValue("primary must have length num_sats");
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

/// refine_candidates(satrecs, candidates, jd0, dt) -> bytes
pub fn pyRefineCandidates(_: [*c]c.PyObject, args: [*c]c.PyObject, kwds: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    var satrecs_obj: [*c]c.PyObject = null;
    var cand_obj: [*c]c.PyObject = null;
    var jd0: f64 = 0;
    var dt: f64 = 0;

    const kwlist = [_:null]?[*:0]const u8{ "satrecs", "candidates", "jd0", "dt", null };
    if (c.PyArg_ParseTupleAndKeywords(args, kwds, "OOdd", @ptrCast(@constCast(&kwlist)), &satrecs_obj, &cand_obj, &jd0, &dt) == 0)
        return null;
    if (!(dt > 0)) {
        py.raiseValue("dt must be > 0");
        return null;
    }
    if (c.PySequence_Check(satrecs_obj) == 0) {
        py.raiseType("satrecs must be a sequence of Satrec objects");
        return null;
    }

    const num_sats: usize = @intCast(c.PySequence_Size(satrecs_obj));
    const propagators = allocator.alloc(astroz.Screening.Propagator, num_sats) catch {
        py.raiseRuntime("Out of memory");
        return null;
    };
    defer allocator.free(propagators);
    // Missing propagators (a Satrec whose init failed) never appear in
    // candidates, since their states are NaN; mark them so a bad index is caught.
    const usable = allocator.alloc(bool, num_sats) catch {
        py.raiseRuntime("Out of memory");
        return null;
    };
    defer allocator.free(usable);

    for (propagators, usable, 0..) |*p, *ok, s| {
        const item = c.PySequence_GetItem(satrecs_obj, @intCast(s)) orelse return null;
        defer py.decref(item);
        if (c.PyObject_TypeCheck(item, &satrec_mod.SatrecType) == 0) {
            py.raiseType("satrecs must be a sequence of Satrec objects");
            return null;
        }
        const sat: *satrec_mod.SatrecObject = @ptrCast(@alignCast(item));
        const epoch_jd = sat.jdsatepoch + sat.jdsatepochF;
        ok.* = true;
        if (sat.sgp4) |m| {
            p.* = .{ .model = .{ .sgp4 = m }, .epochJd = epoch_jd };
        } else if (sat.sdp4) |m| {
            p.* = .{ .model = .{ .sdp4 = @ptrCast(@alignCast(m)) }, .epochJd = epoch_jd };
        } else {
            ok.* = false;
        }
    }

    var cand_buf = std.mem.zeroes(c.Py_buffer);
    if (!getTypedBuffer(cand_obj, &cand_buf, "IL", 4, "candidates must be a contiguous uint32 array of shape (n, 3)")) return null;
    defer c.PyBuffer_Release(&cand_buf);
    const cand_len: usize = @intCast(cand_buf.len);
    if (cand_len % @sizeOf(astroz.Screening.Candidate) != 0) {
        py.raiseValue("candidates must have shape (n, 3)");
        return null;
    }
    const cands: []const astroz.Screening.Candidate = @as([*]const astroz.Screening.Candidate, @ptrCast(@alignCast(cand_buf.buf)))[0 .. cand_len / @sizeOf(astroz.Screening.Candidate)];
    for (cands) |cand| {
        if (cand.i >= num_sats or cand.j >= num_sats or !usable[cand.i] or !usable[cand.j]) {
            py.raiseValue("candidate refers to a missing or uninitialized satellite");
            return null;
        }
    }

    // Satrecs stay alive for the call: the caller holds the sequence.
    const save = c.PyEval_SaveThread();
    const result = astroz.Screening.refine(allocator, propagators, cands, jd0, dt);
    c.PyEval_RestoreThread(save);
    const events = result catch {
        py.raiseRuntime("Out of memory");
        return null;
    };
    defer allocator.free(events);

    // Per event: time (s), distance (km), speed (km/s), converged (1 or 0)
    const out = c.PyBytes_FromStringAndSize(null, @intCast(events.len * 4 * @sizeOf(f64))) orelse return null;
    const vals: [*]f64 = @ptrCast(@alignCast(c.PyBytes_AsString(out)));
    for (events, 0..) |ev, k| {
        vals[4 * k] = ev.tSec;
        vals[4 * k + 1] = ev.distance;
        vals[4 * k + 2] = ev.speed;
        vals[4 * k + 3] = if (ev.converged) 1.0 else 0.0;
    }
    return out;
}

//! All-on-all conjunction screening that cannot miss an encounter.
//!
//! Objects are sampled on a uniform grid t_k = k*dt. Every local minimum of the
//! range rho(t) between two objects lies in exactly one interval [t_k, t_k+1],
//! the one where f = rho . rho_dot changes sign from negative to non-negative.
//!
//! The nearer interval endpoint is at most dt/2 from the time of closest
//! approach (TCA). With |rho_ddot| <= A (twice the largest gravitational
//! acceleration), Taylor's theorem gives
//!
//!     |rho_k| <= d_tca + |rho_dot_k| * dt/2 + A/2 * (dt/2)^2
//!
//! so any encounter with d_tca <= D is caught by a grid search of radius
//! R = D + vrel_max * dt/2 + A/2 * (dt/2)^2 at the nearer endpoint. Each pair
//! inside R is then tested on both adjacent intervals: a sign change of f, and a
//! lower bound on the interval minimum (linear motion from each endpoint minus
//! A dt^2 / 2). Only intervals that provably cannot reach D are dropped; the
//! survivors go to `refine`, which finds the exact TCA in each.
//!
//! Assumption: at most one range minimum per interval. Two need f' < 0 inside,
//! i.e. |rho_dot|^2 < rho . rho_ddot; with real relative dynamics (tidal
//! acceleration ~3 n^2 |rho|) that only happens for co-moving pairs (under
//! ~10 m/s at 5 km): docked or formation objects whose range changes on orbit
//! timescales.
//!
//! Ported from the astrowatch pipeline.

const std = @import("std");
const Allocator = std.mem.Allocator;
const constants = @import("constants.zig");
const Constellation = @import("Constellation.zig");
const Sgp4 = @import("Sgp4.zig");
const Sdp4 = @import("Sdp4.zig");
const Satellite = @import("Satellite.zig");

/// WGS72, the gravity model SGP4 element sets are fitted with.
pub const mu = 398600.8;
const rEarth = constants.wgs84.radiusEarthKm;
/// Bound on relative acceleration, km/s^2: two objects at the surface. The 5%
/// covers SGP4's J2 and drag terms beyond point-mass gravity.
pub const accelMax = 1.05 * 2.0 * mu / (rEarth * rEarth);

/// A pair whose range has a local minimum that may be within the threshold,
/// somewhere in [step*dt, (step+1)*dt]. Always i < j.
pub const Candidate = extern struct {
    i: u32,
    j: u32,
    step: u32,

    fn lessThan(_: void, a: Candidate, b: Candidate) bool {
        if (a.i != b.i) return a.i < b.i;
        if (a.j != b.j) return a.j < b.j;
        return a.step < b.step;
    }

    fn eql(a: Candidate, b: Candidate) bool {
        return a.i == b.i and a.j == b.j and a.step == b.step;
    }
};

pub const Options = struct {
    /// Seconds between grid steps.
    dt: f64,
    /// Km.
    threshold: f64,
    /// If set, only pairs with at least one primary object are reported.
    primary: ?[]const bool = null,
};

/// Screens time-major states (`pos[k * nSats + s]`, TEME km and km/s) and returns
/// the sorted, unique candidate intervals. The caller owns the result.
pub fn candidates(gpa: Allocator, pos: []const [3]f64, vel: []const [3]f64, nSats: usize, opts: Options) ![]Candidate {
    std.debug.assert(pos.len == vel.len and pos.len % @max(nSats, 1) == 0);
    if (opts.primary) |p| std.debug.assert(p.len == nSats);
    const nTimes = if (nSats == 0) 0 else pos.len / nSats;
    if (nTimes == 0 or nSats < 2) return gpa.alloc(Candidate, 0);

    const valid = try gpa.alloc(bool, pos.len);
    defer gpa.free(valid);
    for (valid, pos, vel) |*ok, r, v| ok.* = validState(r, v);

    const grid: Grid = .{ .pos = pos, .vel = vel, .valid = valid, .nSats = nSats, .nTimes = nTimes, .opts = opts };

    const nt = @max(@min(Constellation.getMaxThreads(), nTimes), 1);
    const workers = try gpa.alloc(Worker, nt);
    defer gpa.free(workers);
    const per = @divCeil(nTimes, nt);
    for (workers, 0..) |*w, t| w.* = .{ .grid = &grid, .gpa = gpa, .start = @min(t * per, nTimes), .end = @min((t + 1) * per, nTimes) };
    defer for (workers) |*w| w.found.deinit(gpa);

    runParallel(Worker.run, workers);

    var total: usize = 0;
    for (workers) |w| {
        if (w.failed) return error.OutOfMemory;
        total += w.found.items.len;
    }
    var out = try gpa.alloc(Candidate, total);
    var n: usize = 0;
    for (workers) |w| {
        @memcpy(out[n..][0..w.found.items.len], w.found.items);
        n += w.found.items.len;
    }
    std.mem.sort(Candidate, out, {}, Candidate.lessThan);
    var unique: usize = 0;
    for (out) |c| {
        if (unique > 0 and c.eql(out[unique - 1])) continue;
        out[unique] = c;
        unique += 1;
    }
    return gpa.realloc(out, unique);
}

/// Calls `run(&jobs[k])` for every job, each on its own thread except the first,
/// which runs here. A job whose thread can't be spawned runs here too.
fn runParallel(comptime run: anytype, jobs: anytype) void {
    var threads: [Constellation.MaxThreads]?std.Thread = @splat(null);
    for (jobs[1..], 1..) |*job, t| threads[t] = std.Thread.spawn(.{}, run, .{job}) catch blk: {
        run(job);
        break :blk null;
    };
    run(&jobs[0]);
    for (threads[1..jobs.len]) |th| if (th) |t| t.join();
}

/// Usable state: finite, above the surface, and bound (v^2 < 2 mu / r). SGP4
/// element sets are bound orbits, so an unbound state is propagation garbage
/// (typically long after decay), and excluding it keeps one absurd velocity from
/// inflating the search radius for every other object.
fn validState(r: [3]f64, v: [3]f64) bool {
    const r2 = dot(r, r);
    const v2 = dot(v, v);
    if (!std.math.isFinite(r2) or !std.math.isFinite(v2)) return false;
    return r2 > rEarth * rEarth and v2 * @sqrt(r2) < 2.0 * mu;
}

const Grid = struct {
    pos: []const [3]f64,
    vel: []const [3]f64,
    valid: []const bool,
    nSats: usize,
    nTimes: usize,
    opts: Options,

    fn at(g: *const Grid, k: usize, s: usize) usize {
        return k * g.nSats + s;
    }

    /// Local minimum of range in [t_kk, t_kk+1] that could be within the threshold.
    /// With |rho_ddot| <= A, |rho(t)| >= |rho_e + rho_dot_e (t - t_e)| - A dt^2 / 2.
    fn intervalHit(g: *const Grid, i: usize, j: usize, kk: usize) bool {
        const a_i = g.at(kk, i);
        const a_j = g.at(kk, j);
        const b_i = g.at(kk + 1, i);
        const b_j = g.at(kk + 1, j);
        if (!(g.valid[a_i] and g.valid[a_j] and g.valid[b_i] and g.valid[b_j])) return false;
        const ra = sub(g.pos[a_i], g.pos[a_j]);
        const rb = sub(g.pos[b_i], g.pos[b_j]);
        const va = sub(g.vel[a_i], g.vel[a_j]);
        const vb = sub(g.vel[b_i], g.vel[b_j]);
        if (!(dot(ra, va) < 0.0 and 0.0 <= dot(rb, vb))) return false;
        const dt = g.opts.dt;
        const slack = 0.5 * accelMax * dt * dt;
        if (linearMin(ra, va, dt) - slack > g.opts.threshold) return false;
        if (linearMin(rb, neg(vb), dt) - slack > g.opts.threshold) return false;
        return true;
    }
};

/// min over tau in [0, span] of |r + v tau|
fn linearMin(r: [3]f64, v: [3]f64, span: f64) f64 {
    const vv = dot(v, v);
    const tau = if (vv > 0.0) std.math.clamp(-dot(r, v) / vv, 0.0, span) else 0.0;
    const p: [3]f64 = .{ r[0] + v[0] * tau, r[1] + v[1] * tau, r[2] + v[2] * tau };
    return @sqrt(dot(p, p));
}

const gridOff: i64 = 1 << 20;
const gridW: i64 = 1 << 21;

const Worker = struct {
    grid: *const Grid,
    gpa: Allocator,
    start: usize,
    end: usize,
    found: std.ArrayList(Candidate) = .empty,
    failed: bool = false,

    fn run(w: *Worker) void {
        w.sweep() catch {
            w.failed = true;
        };
    }

    const Entry = struct {
        key: i64,
        z: f64,
        id: u32,

        fn lessThan(_: void, a: Entry, b: Entry) bool {
            if (a.key != b.key) return a.key < b.key;
            return a.z < b.z;
        }
    };

    /// For each step: bin objects into x/y columns of width R, sort by z within a
    /// column, and compare each column against itself and 4 half-neighbor columns
    /// over a z-window of +-R. Every pair within R lands in one of those comparisons.
    fn sweep(w: *Worker) !void {
        const g = w.grid;
        const entries = try w.gpa.alloc(Entry, g.nSats);
        defer w.gpa.free(entries);
        const colStart = try w.gpa.alloc(usize, g.nSats + 1);
        defer w.gpa.free(colStart);
        const colKey = try w.gpa.alloc(i64, g.nSats);
        defer w.gpa.free(colKey);

        const half = 0.5 * g.opts.dt;
        const accelTerm = 0.5 * accelMax * half * half;

        for (w.start..w.end) |k| {
            var vmax: f64 = 0.0;
            for (0..g.nSats) |s| {
                const idx = g.at(k, s);
                if (g.valid[idx]) vmax = @max(vmax, @sqrt(dot(g.vel[idx], g.vel[idx])));
            }
            const r = g.opts.threshold + 2.0 * vmax * half + accelTerm;

            var m: usize = 0;
            for (0..g.nSats) |s| {
                const p = g.pos[g.at(k, s)];
                if (!g.valid[g.at(k, s)]) continue;
                const cx: i64 = @intFromFloat(@floor(p[0] / r));
                const cy: i64 = @intFromFloat(@floor(p[1] / r));
                entries[m] = .{ .key = (cx + gridOff) * gridW + (cy + gridOff), .z = p[2], .id = @intCast(s) };
                m += 1;
            }
            if (m == 0) continue;
            const es = entries[0..m];
            std.mem.sort(Entry, es, {}, Entry.lessThan);

            var nCols: usize = 0;
            for (es, 0..) |e, a| {
                if (a == 0 or e.key != es[a - 1].key) {
                    colStart[nCols] = a;
                    colKey[nCols] = e.key;
                    nCols += 1;
                }
            }
            colStart[nCols] = m;

            for (0..nCols) |c| {
                const s0 = colStart[c];
                const s1 = colStart[c + 1];
                for (s0..s1) |a| {
                    for (a + 1..s1) |b| {
                        if (es[b].z - es[a].z > r) break;
                        try w.checkPair(k, es[a].id, es[b].id, half, accelTerm);
                    }
                }
                // Half-neighborhood: (+1,-1), (+1,0), (+1,+1), (0,+1)
                const targets = [_]i64{ colKey[c] + gridW - 1, colKey[c] + gridW, colKey[c] + gridW + 1, colKey[c] + 1 };
                for (targets) |target| {
                    const t = lowerBound(colKey[0..nCols], target);
                    if (t >= nCols or colKey[t] != target) continue;
                    const t1 = colStart[t + 1];
                    var p = colStart[t];
                    for (s0..s1) |a| {
                        while (p < t1 and es[p].z < es[a].z - r) p += 1;
                        var q = p;
                        while (q < t1 and es[q].z <= es[a].z + r) : (q += 1) {
                            try w.checkPair(k, es[a].id, es[q].id, half, accelTerm);
                        }
                    }
                }
            }
        }
    }

    /// Tight per-pair bound at step k, then the interval test on both adjacent intervals.
    fn checkPair(w: *Worker, k: usize, a: u32, b: u32, half: f64, accelTerm: f64) !void {
        const g = w.grid;
        if (g.opts.primary) |p| if (!(p[a] or p[b])) return;
        const dr = sub(g.pos[g.at(k, b)], g.pos[g.at(k, a)]);
        const dv = sub(g.vel[g.at(k, b)], g.vel[g.at(k, a)]);
        const d = @sqrt(dot(dr, dr));
        const sv = @sqrt(dot(dv, dv));
        if (d > g.opts.threshold + sv * half + accelTerm) return;
        const lo = @min(a, b);
        const hi = @max(a, b);
        // The minimum may lie in the interval before or after step k
        for (k -| 1..@min(k + 1, g.nTimes - 1)) |kk| {
            if (g.intervalHit(lo, hi, kk)) try w.found.append(w.gpa, .{ .i = lo, .j = hi, .step = @intCast(kk) });
        }
    }
};

/// One satellite's scalar propagator and epoch.
pub const Propagator = struct {
    model: union(enum) { sgp4: *const Sgp4, sdp4: *const Sdp4 },
    epochJd: f64,

    /// The propagator `sat` already holds. `sat` must outlive the result.
    pub fn fromSatellite(sat: *const Satellite) Propagator {
        return .{
            .model = switch (sat.impl) {
                .sgp4 => |*m| .{ .sgp4 = m },
                .sdp4 => |*m| .{ .sdp4 = m },
            },
            .epochJd = sat.epochJd(),
        };
    }

    /// TEME state at `jd0` plus `tSec` seconds; NaN if propagation fails.
    fn state(p: Propagator, jd0: f64, tSec: f64) [2][3]f64 {
        const tsince = (jd0 - p.epochJd) * constants.minutesPerDay + tSec / 60.0;
        const pv = switch (p.model) {
            .sgp4 => |m| m.propagate(tsince),
            .sdp4 => |m| m.propagate(tsince),
        };
        const nan: [3]f64 = @splat(std.math.nan(f64));
        return pv catch .{ nan, nan };
    }
};

/// A refined closest approach. `tSec` is seconds after the `jd0` passed to
/// `refine`, `distance` is km and `speed` km/s.
pub const Event = struct {
    i: u32,
    j: u32,
    tSec: f64,
    distance: f64,
    speed: f64,
    converged: bool,
};

/// Seconds.
pub const tcaTolerance = 1e-4;
const tcaMaxIter = 60;

/// Finds the time of closest approach in each candidate interval, on a grid
/// that started at `jd0` with steps `dt` seconds apart. `propagators[s]`
/// propagates satellite `s`. The caller owns the result, which is in the same
/// order as `cands`.
pub fn refine(gpa: Allocator, propagators: []const Propagator, cands: []const Candidate, jd0: f64, dt: f64) ![]Event {
    const events = try gpa.alloc(Event, cands.len);
    if (cands.len == 0) return events;

    var jobs: [Constellation.MaxThreads]RefineJob = undefined;
    const nt = @max(@min(Constellation.getMaxThreads(), cands.len), 1);
    const per = @divCeil(cands.len, nt);
    for (jobs[0..nt], 0..) |*job, t| job.* = .{
        .propagators = propagators,
        .cands = cands[@min(t * per, cands.len)..@min((t + 1) * per, cands.len)],
        .events = events[@min(t * per, cands.len)..@min((t + 1) * per, cands.len)],
        .jd0 = jd0,
        .dt = dt,
    };
    runParallel(RefineJob.run, jobs[0..nt]);
    return events;
}

const RefineJob = struct {
    propagators: []const Propagator,
    cands: []const Candidate,
    events: []Event,
    jd0: f64,
    dt: f64,

    fn run(job: *RefineJob) void {
        for (job.cands, job.events) |c, *ev| {
            ev.* = refineOne(job.propagators[c.i], job.propagators[c.j], c, job.jd0, job.dt);
        }
    }
};

/// Root of f = rho . rho_dot in the candidate's interval, where f(a) < 0 <= f(b).
/// Safeguarded Newton with f' = |rho_dot|^2 + rho . rho_ddot (two-body
/// rho_ddot), bisecting whenever a step leaves the bracket or the bracket fails
/// to halve (flat f on slow pairs), so it converges within `tcaMaxIter`.
fn refineOne(pi: Propagator, pj: Propagator, c: Candidate, jd0: f64, dt: f64) Event {
    var a = @as(f64, @floatFromInt(c.step)) * dt;
    var b = a + dt;
    var t = a + 0.5 * dt;
    var width = b - a;
    var converged = false;
    for (0..tcaMaxIter) |_| {
        const si = pi.state(jd0, t);
        const sj = pj.state(jd0, t);
        const rho = sub(si[0], sj[0]);
        const rhoDot = sub(si[1], sj[1]);
        const ri3 = std.math.pow(f64, dot(si[0], si[0]), 1.5);
        const rj3 = std.math.pow(f64, dot(sj[0], sj[0]), 1.5);
        var acc: [3]f64 = undefined;
        for (&acc, si[0], sj[0]) |*x, ri, rj| x.* = -mu * (ri / ri3 - rj / rj3);
        const f = dot(rho, rhoDot);
        const fp = dot(rhoDot, rhoDot) + dot(rho, acc);

        if (f < 0) a = t else b = t;
        var tNew = t - f / fp;
        const stalled = (b - a) > 0.5 * width;
        width = b - a;
        if (!std.math.isFinite(tNew) or tNew <= a or tNew >= b or stalled) tNew = 0.5 * (a + b);
        const moved = @abs(tNew - t);
        t = tNew;
        if (!(moved > tcaTolerance and b - a > tcaTolerance)) {
            converged = true;
            break;
        }
    }
    const si = pi.state(jd0, t);
    const sj = pj.state(jd0, t);
    const dr = sub(si[0], sj[0]);
    const dv = sub(si[1], sj[1]);
    return .{ .i = c.i, .j = c.j, .tSec = t, .distance = @sqrt(dot(dr, dr)), .speed = @sqrt(dot(dv, dv)), .converged = converged };
}

fn lowerBound(keys: []const i64, target: i64) usize {
    var lo: usize = 0;
    var hi: usize = keys.len;
    while (lo < hi) {
        const mid = lo + (hi - lo) / 2;
        if (keys[mid] < target) lo = mid + 1 else hi = mid;
    }
    return lo;
}

fn sub(a: [3]f64, b: [3]f64) [3]f64 {
    return .{ a[0] - b[0], a[1] - b[1], a[2] - b[2] };
}

fn neg(a: [3]f64) [3]f64 {
    return .{ -a[0], -a[1], -a[2] };
}

fn dot(a: [3]f64, b: [3]f64) f64 {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/// Two-body states on a circular orbit, for tests.
fn circular(radius: f64, incl: f64, raan: f64, phase: f64, t: f64) [2][3]f64 {
    const n = @sqrt(mu / (radius * radius * radius));
    const u = phase + n * t;
    const cu, const su = .{ @cos(u), @sin(u) };
    const ci, const si = .{ @cos(incl), @sin(incl) };
    const co, const so = .{ @cos(raan), @sin(raan) };
    const v = n * radius;
    return .{
        .{ radius * (co * cu - so * ci * su), radius * (so * cu + co * ci * su), radius * si * su },
        .{ v * (-co * su - so * ci * cu), v * (-so * su + co * ci * cu), v * si * cu },
    };
}

fn sample(gpa: Allocator, orbits: []const [4]f64, nTimes: usize, dt: f64) !struct { pos: [][3]f64, vel: [][3]f64 } {
    const pos = try gpa.alloc([3]f64, nTimes * orbits.len);
    const vel = try gpa.alloc([3]f64, nTimes * orbits.len);
    for (0..nTimes) |k| for (orbits, 0..) |o, s| {
        const st = circular(o[0], o[1], o[2], o[3], @as(f64, @floatFromInt(k)) * dt);
        pos[k * orbits.len + s] = st[0];
        vel[k * orbits.len + s] = st[1];
    };
    return .{ .pos = pos, .vel = vel };
}

test "finds a crossing that falls between grid points" {
    const gpa = std.testing.allocator;
    // Two 7000 km circular orbits in different planes meet at the node at t=0
    // (phase 0); offsetting one phase by 0.5 km of arc puts the pass ~0.5 km
    // apart, and starting the grid 30 s early puts TCA between two samples.
    const dt = 60.0;
    const lead = 30.0;
    const n = @sqrt(mu / (7000.0 * 7000.0 * 7000.0));
    const orbits = [_][4]f64{
        .{ 7000.0, 0.0, 0.0, -n * lead },
        .{ 7000.0, 1.0, 0.0, -n * lead + 0.5 / 7000.0 },
        .{ 7000.0, 0.5, 2.0, 0.0 }, // far away
    };
    const st = try sample(gpa, &orbits, 4, dt);
    defer gpa.free(st.pos);
    defer gpa.free(st.vel);
    // At the samples the pair is ~200 km apart, nowhere near the 5 km threshold.
    try std.testing.expect(@sqrt(dot(sub(st.pos[0], st.pos[1]), sub(st.pos[0], st.pos[1]))) > 100.0);
    try std.testing.expect(@sqrt(dot(sub(st.pos[3], st.pos[4]), sub(st.pos[3], st.pos[4]))) > 100.0);

    const found = try candidates(gpa, st.pos, st.vel, orbits.len, .{ .dt = dt, .threshold = 5.0 });
    defer gpa.free(found);
    try std.testing.expectEqualSlices(Candidate, &.{.{ .i = 0, .j = 1, .step = 0 }}, found);
}

test "primary mask limits the pairs" {
    const gpa = std.testing.allocator;
    const n = @sqrt(mu / (7000.0 * 7000.0 * 7000.0));
    const orbits = [_][4]f64{
        .{ 7000.0, 0.0, 0.0, -n * 30.0 },
        .{ 7000.0, 1.0, 0.0, -n * 30.0 + 0.5 / 7000.0 },
        .{ 7000.0, 0.5, 2.0, 0.0 },
    };
    const st = try sample(gpa, &orbits, 4, 60.0);
    defer gpa.free(st.pos);
    defer gpa.free(st.vel);
    const none = try candidates(gpa, st.pos, st.vel, 3, .{ .dt = 60.0, .threshold = 5.0, .primary = &.{ false, false, true } });
    defer gpa.free(none);
    try std.testing.expectEqual(0, none.len);
    const some = try candidates(gpa, st.pos, st.vel, 3, .{ .dt = 60.0, .threshold = 5.0, .primary = &.{ false, true, false } });
    defer gpa.free(some);
    try std.testing.expectEqual(1, some.len);
}

test "invalid states are skipped" {
    const gpa = std.testing.allocator;
    const nan = std.math.nan(f64);
    const pos = [_][3]f64{ .{ 7000, 0, 0 }, .{ 7000, 0.1, 0 }, .{ nan, 0, 0 }, .{ 7000, 0.2, 0 } };
    const vel = [_][3]f64{ .{ 0, 7.5, 0 }, .{ 0, -7.5, 0 }, .{ 0, 7.5, 0 }, .{ 0, -7.5, 0 } };
    const found = try candidates(gpa, &pos, &vel, 2, .{ .dt = 1.0, .threshold = 5.0 });
    defer gpa.free(found);
    try std.testing.expectEqual(0, found.len);
}

test "refine matches a dense scan of each candidate interval" {
    const gpa = std.testing.allocator;
    const Tle = @import("Tle.zig");
    // Same altitude, planes 10 degrees apart, phases offset so they cross near the node.
    const omm = [_][]const u8{
        \\{"EPOCH":"2026-07-14T21:45:20","MEAN_MOTION":15.05,"ECCENTRICITY":0.0001,"INCLINATION":50,"RA_OF_ASC_NODE":0,"ARG_OF_PERICENTER":0,"MEAN_ANOMALY":0,"NORAD_CAT_ID":1,"BSTAR":0}
        ,
        \\{"EPOCH":"2026-07-14T21:45:20","MEAN_MOTION":15.05,"ECCENTRICITY":0.0001,"INCLINATION":60,"RA_OF_ASC_NODE":0,"ARG_OF_PERICENTER":0,"MEAN_ANOMALY":0.1,"NORAD_CAT_ID":2,"BSTAR":0}
        ,
        \\{"EPOCH":"2026-07-14T21:45:20","MEAN_MOTION":14.2,"ECCENTRICITY":0.001,"INCLINATION":98,"RA_OF_ASC_NODE":120,"ARG_OF_PERICENTER":0,"MEAN_ANOMALY":200,"NORAD_CAT_ID":3,"BSTAR":0}
        ,
    };
    var sgp4s: [omm.len]Sgp4 = undefined;
    var props: [omm.len]Propagator = undefined;
    for (omm, &sgp4s, &props) |json, *s, *p| {
        var tle = try Tle.parseOmm(json, gpa);
        defer tle.deinit();
        s.* = try Sgp4.init(tle, constants.wgs72);
        p.* = .{ .model = .{ .sgp4 = s }, .epochJd = tle.epochJd };
    }

    const jd0 = props[0].epochJd - 10.0 / constants.minutesPerDay; // grid starts 10 min before epoch
    const dt = 60.0;
    const nTimes = 200;
    var pos: [nTimes * omm.len][3]f64 = undefined;
    var vel: [nTimes * omm.len][3]f64 = undefined;
    for (0..nTimes) |k| for (props, 0..) |p, s| {
        const st = p.state(jd0, @as(f64, @floatFromInt(k)) * dt);
        pos[k * omm.len + s] = st[0];
        vel[k * omm.len + s] = st[1];
    };

    const found = try candidates(gpa, &pos, &vel, omm.len, .{ .dt = dt, .threshold = 50.0 });
    defer gpa.free(found);
    try std.testing.expect(found.len > 0);
    const events = try refine(gpa, &props, found, jd0, dt);
    defer gpa.free(events);

    for (found, events) |c, ev| {
        try std.testing.expect(ev.converged);
        try std.testing.expectEqual(c.i, ev.i);
        // Dense scan of the interval at 1 ms
        var best: f64 = std.math.inf(f64);
        var bestT: f64 = 0;
        const a = @as(f64, @floatFromInt(c.step)) * dt;
        var t = a;
        while (t <= a + dt) : (t += 1e-3) {
            const si = props[c.i].state(jd0, t);
            const sj = props[c.j].state(jd0, t);
            const dr = sub(si[0], sj[0]);
            const d = @sqrt(dot(dr, dr));
            if (d < best) {
                best = d;
                bestT = t;
            }
        }
        try std.testing.expectApproxEqAbs(bestT, ev.tSec, 2e-3);
        try std.testing.expect(ev.distance <= best + 1e-9);
        try std.testing.expectApproxEqAbs(best, ev.distance, 1e-6);
    }
}

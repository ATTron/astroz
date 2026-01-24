//! Conjunction screening: pairwise distance kernel for satellite close-approach detection.

const std = @import("std");
const py = @import("python.zig");
const c = py.c;

const allocator = std.heap.c_allocator;
const Vec4 = @Vector(4, f64);

/// Compute minimum distances between pairs of satellites across all time steps
pub fn minDistancesBatch(
    positions: [*]const f64,
    num_sats: usize,
    num_times: usize,
    pairs: [*]const u32,
    num_pairs: usize,
    out_dists: [*]f64,
    out_indices: [*]u32,
    velocities: ?[*]const f64,
    out_rel_vels: ?[*]f64,
) void {
    _ = num_sats;

    for (0..num_pairs) |p| {
        const sat_i = pairs[p * 2];
        const sat_j = pairs[p * 2 + 1];
        const sat_i_base = @as(usize, sat_i) * num_times * 3;
        const sat_j_base = @as(usize, sat_j) * num_times * 3;

        var min_dist_sq: f64 = std.math.floatMax(f64);
        var min_idx: usize = 0;

        // SIMD: process 4 time steps at once
        var t: usize = 0;
        while (t + 4 <= num_times) : (t += 4) {
            // Load positions for satellite i at 4 consecutive time steps
            var xi: Vec4 = undefined;
            var yi: Vec4 = undefined;
            var zi: Vec4 = undefined;
            var xj: Vec4 = undefined;
            var yj: Vec4 = undefined;
            var zj: Vec4 = undefined;

            inline for (0..4) |k| {
                const offset_i = sat_i_base + (t + k) * 3;
                const offset_j = sat_j_base + (t + k) * 3;
                xi[k] = positions[offset_i];
                yi[k] = positions[offset_i + 1];
                zi[k] = positions[offset_i + 2];
                xj[k] = positions[offset_j];
                yj[k] = positions[offset_j + 1];
                zj[k] = positions[offset_j + 2];
            }

            const dx = xi - xj;
            const dy = yi - yj;
            const dz = zi - zj;
            const dist_sq = dx * dx + dy * dy + dz * dz;

            // Find minimum in this SIMD group
            inline for (0..4) |k| {
                if (dist_sq[k] < min_dist_sq) {
                    min_dist_sq = dist_sq[k];
                    min_idx = t + k;
                }
            }
        }

        // Handle remainder
        while (t < num_times) : (t += 1) {
            const offset_i = sat_i_base + t * 3;
            const offset_j = sat_j_base + t * 3;
            const dx = positions[offset_i] - positions[offset_j];
            const dy = positions[offset_i + 1] - positions[offset_j + 1];
            const dz = positions[offset_i + 2] - positions[offset_j + 2];
            const dist_sq = dx * dx + dy * dy + dz * dz;
            if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
                min_idx = t;
            }
        }

        out_dists[p] = @sqrt(min_dist_sq);
        out_indices[p] = @intCast(min_idx);

        // Compute relative velocity at min_idx if requested
        if (velocities) |vels| {
            if (out_rel_vels) |rel_vels| {
                const vi_base = @as(usize, sat_i) * num_times * 3 + min_idx * 3;
                const vj_base = @as(usize, sat_j) * num_times * 3 + min_idx * 3;
                const dvx = vels[vi_base] - vels[vj_base];
                const dvy = vels[vi_base + 1] - vels[vj_base + 1];
                const dvz = vels[vi_base + 2] - vels[vj_base + 2];
                rel_vels[p] = @sqrt(dvx * dvx + dvy * dvy + dvz * dvz);
            }
        }
    }
}

/// Coarse screening: find all satellite pairs within threshold distance at any time step
pub fn coarseScreen(
    positions: [*]const f64,
    num_sats: usize,
    num_times: usize,
    threshold: f64,
    valid_mask: ?[*]const u8,
    out_pairs: [*]u32,
    out_t_indices: [*]u32,
    max_results: usize,
) usize {
    const inv_cell_size = 1.0 / threshold;
    const threshold_sq = threshold * threshold;
    var result_count: usize = 0;

    const TABLE_BITS = 16;
    const TABLE_SIZE: usize = 1 << TABLE_BITS;
    const TABLE_MASK: u32 = TABLE_SIZE - 1;

    // Per satellite cell coords (3 * num_sats i32s) + hashes + chain links
    const sat_cells = allocator.alloc(i32, num_sats * 3) catch return 0;
    defer allocator.free(sat_cells);
    const sat_cx = sat_cells[0..num_sats];
    const sat_cy = sat_cells[num_sats .. num_sats * 2];
    const sat_cz = sat_cells[num_sats * 2 .. num_sats * 3];

    const sat_hashes = allocator.alloc(u32, num_sats) catch return 0;
    defer allocator.free(sat_hashes);

    // Hash table: head[hash] → first satellite index, next[idx] -> next in chain
    const head = allocator.alloc(u32, TABLE_SIZE) catch return 0;
    defer allocator.free(head);
    const next = allocator.alloc(u32, num_sats) catch return 0;
    defer allocator.free(next);

    const EMPTY: u32 = std.math.maxInt(u32);

    for (0..num_times) |t| {
        // Reset hash table heads
        @memset(head, EMPTY);

        // Assign satellites to cells and build hash chains
        for (0..num_sats) |s| {
            if (valid_mask) |m| {
                if (m[s] == 0) {
                    sat_hashes[s] = EMPTY;
                    continue;
                }
            }
            const base = s * num_times * 3 + t * 3;
            const x = positions[base];
            if (!std.math.isFinite(x)) {
                sat_hashes[s] = EMPTY;
                continue;
            }
            const y = positions[base + 1];
            const z = positions[base + 2];

            const cx = @as(i32, @intFromFloat(@floor(x * inv_cell_size)));
            const cy = @as(i32, @intFromFloat(@floor(y * inv_cell_size)));
            const cz = @as(i32, @intFromFloat(@floor(z * inv_cell_size)));
            sat_cx[s] = cx;
            sat_cy[s] = cy;
            sat_cz[s] = cz;

            const h = spatialHash(cx, cy, cz) & TABLE_MASK;
            sat_hashes[s] = h;

            // Insert into chain
            next[s] = head[h];
            head[h] = @intCast(s);
        }

        // For each satellite, check 27 neighboring cells via their hashes
        for (0..num_sats) |s| {
            if (sat_hashes[s] == EMPTY) continue;

            const base_s = s * num_times * 3 + t * 3;
            const sx = positions[base_s];
            const sy = positions[base_s + 1];
            const sz = positions[base_s + 2];
            const scx = sat_cx[s];
            const scy = sat_cy[s];
            const scz = sat_cz[s];

            // Check 27 neighbor cells
            var dcx: i32 = -1;
            while (dcx <= 1) : (dcx += 1) {
                var dcy: i32 = -1;
                while (dcy <= 1) : (dcy += 1) {
                    var dcz: i32 = -1;
                    while (dcz <= 1) : (dcz += 1) {
                        const ncx = scx + dcx;
                        const ncy = scy + dcy;
                        const ncz = scz + dcz;
                        const nh = spatialHash(ncx, ncy, ncz) & TABLE_MASK;

                        // Walk the chain for this hash bucket
                        var idx = head[nh];
                        while (idx != EMPTY) {
                            const other = idx;
                            idx = next[idx];

                            if (other <= s) continue; // avoid duplicates

                            // Check actual cell match (hash collisions possible)
                            if (sat_cx[other] != ncx or sat_cy[other] != ncy or sat_cz[other] != ncz) continue;

                            const base_o = @as(usize, other) * num_times * 3 + t * 3;
                            const dx = sx - positions[base_o];
                            const dy = sy - positions[base_o + 1];
                            const dz = sz - positions[base_o + 2];
                            const dist_sq = dx * dx + dy * dy + dz * dz;

                            if (dist_sq < threshold_sq) {
                                if (result_count >= max_results) return result_count;
                                out_pairs[result_count * 2] = @intCast(s);
                                out_pairs[result_count * 2 + 1] = @intCast(other);
                                out_t_indices[result_count] = @intCast(t);
                                result_count += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    return result_count;
}

fn spatialHash(cx: i32, cy: i32, cz: i32) u32 {
    // Fast spatial hash combining cell coordinates
    var h: u32 = @bitCast(cx);
    h = h *% 2654435761; // Knuth's multiplicative hash
    h ^= @bitCast(cy);
    h = h *% 2654435761;
    h ^= @bitCast(cz);
    h = h *% 2654435761;
    return h;
}

/// Python wrapper for coarseScreen
pub fn py_coarse_screen(_: [*c]c.PyObject, args: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    var pos_obj: [*c]c.PyObject = null;
    var num_sats_long: c_long = 0;
    var threshold: f64 = 0;
    var mask_obj: [*c]c.PyObject = null;

    if (c.PyArg_ParseTuple(args, "Old|O", &pos_obj, &num_sats_long, &threshold, &mask_obj) == 0) return null;

    const num_sats: usize = @intCast(num_sats_long);

    // Get positions buffer
    var pos_buf = std.mem.zeroes(c.Py_buffer);
    if (c.PyObject_GetBuffer(pos_obj, &pos_buf, c.PyBUF_SIMPLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&pos_buf);

    const positions: [*]const f64 = @ptrCast(@alignCast(pos_buf.buf));
    const total_elements = @as(usize, @intCast(pos_buf.len)) / @sizeOf(f64);

    if (num_sats == 0 or total_elements % (num_sats * 3) != 0) {
        py.raiseValue("positions array size not consistent with num_sats");
        return null;
    }
    const num_times = total_elements / (num_sats * 3);

    // Get optional valid_mask
    var mask_buf = std.mem.zeroes(c.Py_buffer);
    var valid_mask: ?[*]const u8 = null;
    var have_mask = false;
    if (!py.isNone(mask_obj)) {
        if (c.PyObject_GetBuffer(mask_obj, &mask_buf, c.PyBUF_SIMPLE | c.PyBUF_FORMAT) < 0) return null;
        have_mask = true;
        valid_mask = @ptrCast(@alignCast(mask_buf.buf));
    }
    defer if (have_mask) c.PyBuffer_Release(&mask_buf);

    // Allocate output buffers (generous initial size)
    const max_results: usize = 10_000_000; // 10M pair-time entries max
    const out_pairs = allocator.alloc(u32, max_results * 2) catch {
        py.raiseRuntime("Out of memory for coarse screen results");
        return null;
    };
    defer allocator.free(out_pairs);

    const out_t_indices = allocator.alloc(u32, max_results) catch {
        py.raiseRuntime("Out of memory for coarse screen results");
        return null;
    };
    defer allocator.free(out_t_indices);

    // Run the kernel
    const count = coarseScreen(
        positions,
        num_sats,
        num_times,
        threshold,
        valid_mask,
        out_pairs.ptr,
        out_t_indices.ptr,
        max_results,
    );

    // Build Python result
    const pairs_list = c.PyList_New(@intCast(count)) orelse return null;
    for (0..count) |i| {
        const pair_tuple = py.tuple(2) orelse {
            c.Py_DECREF(pairs_list);
            return null;
        };
        py.tupleSet(pair_tuple, 0, c.PyLong_FromUnsignedLong(@as(c_ulong, out_pairs[i * 2])) orelse {
            c.Py_DECREF(pairs_list);
            return null;
        });
        py.tupleSet(pair_tuple, 1, c.PyLong_FromUnsignedLong(@as(c_ulong, out_pairs[i * 2 + 1])) orelse {
            c.Py_DECREF(pairs_list);
            return null;
        });
        _ = c.PyList_SetItem(pairs_list, @intCast(i), pair_tuple);
    }

    const times_list = py.listFromU32(out_t_indices[0..count]) orelse {
        c.Py_DECREF(pairs_list);
        return null;
    };

    const result = py.tuple(2) orelse {
        c.Py_DECREF(pairs_list);
        c.Py_DECREF(times_list);
        return null;
    };
    py.tupleSet(result, 0, pairs_list);
    py.tupleSet(result, 1, times_list);
    return result;
}

/// Python wrapper for minDistancesBatch
pub fn py_min_distances(_: [*c]c.PyObject, args: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    var pos_obj: [*c]c.PyObject = null;
    var pairs_obj: [*c]c.PyObject = null;
    var vel_obj: [*c]c.PyObject = null;
    var num_sats_obj: [*c]c.PyObject = null;

    if (c.PyArg_ParseTuple(args, "OO|OO", &pos_obj, &pairs_obj, &vel_obj, &num_sats_obj) == 0) return null;

    // Get positions buffer
    var pos_buf = std.mem.zeroes(c.Py_buffer);
    if (c.PyObject_GetBuffer(pos_obj, &pos_buf, c.PyBUF_SIMPLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&pos_buf);

    // Get pairs buffer (uint32)
    var pairs_buf = std.mem.zeroes(c.Py_buffer);
    if (c.PyObject_GetBuffer(pairs_obj, &pairs_buf, c.PyBUF_SIMPLE | c.PyBUF_FORMAT) < 0) return null;
    defer c.PyBuffer_Release(&pairs_buf);

    // Parse dimensions from positions buffer
    const pairs_count = @as(usize, @intCast(pairs_buf.len)) / (@sizeOf(u32) * 2);
    const pairs: [*]const u32 = @ptrCast(@alignCast(pairs_buf.buf));
    const positions: [*]const f64 = @ptrCast(@alignCast(pos_buf.buf));
    const total_pos_elements = @as(usize, @intCast(pos_buf.len)) / @sizeOf(f64);

    if (pairs_count == 0) {
        py.raiseValue("pairs array must not be empty");
        return null;
    }

    // Get num_sats: from explicit argument or infer from max pair index
    var num_sats: usize = 0;
    if (!py.isNone(num_sats_obj)) {
        const ns = c.PyLong_AsLong(num_sats_obj);
        if (ns <= 0) {
            py.raiseValue("num_sats must be a positive integer");
            return null;
        }
        num_sats = @intCast(ns);
    } else {
        // Infer from max index in pairs + 1
        var max_sat_idx: u32 = 0;
        for (0..pairs_count * 2) |i| {
            if (pairs[i] > max_sat_idx) max_sat_idx = pairs[i];
        }
        num_sats = @as(usize, max_sat_idx) + 1;
    }

    // num_times = total_elements / (num_sats * 3)
    if (num_sats == 0 or total_pos_elements % (num_sats * 3) != 0) {
        py.raiseValue("positions array size not consistent with satellite count");
        return null;
    }
    const num_times = total_pos_elements / (num_sats * 3);

    // Get optional velocities buffer
    var vel_buf = std.mem.zeroes(c.Py_buffer);
    var velocities: ?[*]const f64 = null;
    var have_vel = false;
    if (!py.isNone(vel_obj)) {
        if (c.PyObject_GetBuffer(vel_obj, &vel_buf, c.PyBUF_SIMPLE | c.PyBUF_FORMAT) < 0) return null;
        have_vel = true;
        velocities = @ptrCast(@alignCast(vel_buf.buf));
    }
    defer if (have_vel) c.PyBuffer_Release(&vel_buf);

    // Allocate output arrays
    const out_dists = allocator.alloc(f64, pairs_count) catch {
        py.raiseRuntime("Out of memory");
        return null;
    };
    defer allocator.free(out_dists);

    const out_indices = allocator.alloc(u32, pairs_count) catch {
        py.raiseRuntime("Out of memory");
        return null;
    };
    defer allocator.free(out_indices);

    var out_rel_vels: ?[]f64 = null;
    if (have_vel) {
        out_rel_vels = allocator.alloc(f64, pairs_count) catch {
            py.raiseRuntime("Out of memory");
            return null;
        };
    }
    defer if (out_rel_vels) |rv| allocator.free(rv);

    // Call the kernel
    minDistancesBatch(
        positions,
        num_sats,
        num_times,
        pairs,
        pairs_count,
        out_dists.ptr,
        out_indices.ptr,
        velocities,
        if (out_rel_vels) |rv| rv.ptr else null,
    );

    // Build Python result: (min_dists, min_indices[, rel_vels])
    const dists_list = py.listFromF64(out_dists[0..pairs_count]) orelse return null;
    const indices_list = py.listFromU32(out_indices[0..pairs_count]) orelse {
        c.Py_DECREF(dists_list);
        return null;
    };

    if (out_rel_vels) |rv| {
        const vels_list = py.listFromF64(rv[0..pairs_count]) orelse {
            c.Py_DECREF(dists_list);
            c.Py_DECREF(indices_list);
            return null;
        };
        const result = py.tuple(3) orelse {
            c.Py_DECREF(dists_list);
            c.Py_DECREF(indices_list);
            c.Py_DECREF(vels_list);
            return null;
        };
        py.tupleSet(result, 0, dists_list);
        py.tupleSet(result, 1, indices_list);
        py.tupleSet(result, 2, vels_list);
        return result;
    } else {
        const result = py.tuple(2) orelse {
            c.Py_DECREF(dists_list);
            c.Py_DECREF(indices_list);
            return null;
        };
        py.tupleSet(result, 0, dists_list);
        py.tupleSet(result, 1, indices_list);
        return result;
    }
}

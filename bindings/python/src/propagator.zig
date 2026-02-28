//! Python bindings for numerical orbit propagation

const std = @import("std");
const py = @import("python.zig");
const c = py.c;
const astroz = @import("astroz");
const propagators = astroz.propagators;
const constants = astroz.constants;

/// Maximum altitude (km) above which exponential atmosphere drag is negligible.
const drag_max_altitude_km = 1500.0;

pub fn pyPropagateNumerical(_: [*c]c.PyObject, args: [*c]c.PyObject, kwargs: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    var state_obj: [*c]c.PyObject = null;
    var t0: f64 = 0.0;
    var duration: f64 = undefined;
    var dt: f64 = 10.0;
    var mu: f64 = undefined;
    var j2_obj: [*c]c.PyObject = null;
    var r_eq_obj: [*c]c.PyObject = null;
    var drag_cd_obj: [*c]c.PyObject = null;
    var drag_area_obj: [*c]c.PyObject = null;
    var drag_mass_obj: [*c]c.PyObject = null;
    var integrator_str: [*c]const u8 = null;
    var rtol_obj: [*c]c.PyObject = null;
    var atol_obj: [*c]c.PyObject = null;

    var kwnames = [_:null]?[*:0]const u8{
        "state",      "t0",   "duration", "dt",        "mu",
        "j2",         "r_eq", "drag_cd",  "drag_area", "drag_mass",
        "integrator", "rtol", "atol",
    };

    if (c.PyArg_ParseTupleAndKeywords(
        args,
        kwargs,
        "Odddd|OOOOOsOO",
        @ptrCast(&kwnames),
        &state_obj,
        &t0,
        &duration,
        &dt,
        &mu,
        &j2_obj,
        &r_eq_obj,
        &drag_cd_obj,
        &drag_area_obj,
        &drag_mass_obj,
        &integrator_str,
        &rtol_obj,
        &atol_obj,
    ) == 0) return null;

    // Parse state vector (6 elements)
    if (c.PySequence_Size(state_obj) != 6) {
        py.raiseValue("state must have exactly 6 elements [x, y, z, vx, vy, vz]");
        return null;
    }
    var initial_state: [6]f64 = undefined;
    inline for (0..6) |i| {
        const item = c.PySequence_GetItem(state_obj, @intCast(i)) orelse return null;
        defer c.Py_DECREF(item);
        initial_state[i] = c.PyFloat_AsDouble(item);
        if (initial_state[i] == -1.0 and c.PyErr_Occurred() != null) return null;
    }

    // Parse optional parameters
    const has_j2 = !py.isNone(j2_obj);
    const has_drag = !py.isNone(drag_cd_obj);

    if ((has_j2 or has_drag) and py.isNone(r_eq_obj)) {
        py.raiseValue("r_eq is required when j2 or drag_cd is specified");
        return null;
    }

    const j2_val = py.optionalFloat(j2_obj, 0) orelse return null;
    const r_eq_val = py.optionalFloat(r_eq_obj, 0) orelse return null;
    const rtol = py.optionalFloat(rtol_obj, 1e-9) orelse return null;
    const atol = py.optionalFloat(atol_obj, 1e-12) orelse return null;

    var drag_cd: f64 = 0;
    var drag_area: f64 = 0;
    var drag_mass: f64 = 0;
    if (has_drag) {
        drag_cd = py.optionalFloat(drag_cd_obj, 0) orelse return null;
        if (py.isNone(drag_area_obj) or py.isNone(drag_mass_obj)) {
            py.raiseValue("drag_area and drag_mass are required when drag_cd is specified");
            return null;
        }
        drag_area = py.optionalFloat(drag_area_obj, 0) orelse return null;
        drag_mass = py.optionalFloat(drag_mass_obj, 0) orelse return null;
    }

    var use_dp87 = true;
    if (integrator_str != null) {
        const int_name = std.mem.span(integrator_str);
        if (std.mem.eql(u8, int_name, "rk4")) {
            use_dp87 = false;
        } else if (!std.mem.eql(u8, int_name, "dp87")) {
            py.raiseValue("integrator must be 'rk4' or 'dp87'");
            return null;
        }
    }

    // Build force models on the stack
    var twobody = propagators.TwoBody.init(mu);
    var j2_model = propagators.J2.init(mu, j2_val, r_eq_val);
    var drag_model = propagators.Drag.init(
        r_eq_val,
        constants.earth.seaLevelDensity,
        constants.earth.scaleHeight,
        drag_cd,
        drag_area,
        drag_mass,
        drag_max_altitude_km,
    );

    // Build force model array and composite
    var force_models_buf: [4]propagators.ForceModel = undefined;
    var n_forces: usize = 0;

    force_models_buf[n_forces] = propagators.ForceModel.wrap(propagators.TwoBody, &twobody);
    n_forces += 1;

    if (has_j2) {
        force_models_buf[n_forces] = propagators.ForceModel.wrap(propagators.J2, &j2_model);
        n_forces += 1;
    }
    if (has_drag) {
        force_models_buf[n_forces] = propagators.ForceModel.wrap(propagators.Drag, &drag_model);
        n_forces += 1;
    }

    const allocator = std.heap.c_allocator;
    var force: propagators.ForceModel = undefined;
    var composite: ?propagators.Composite = null;

    if (n_forces == 1) {
        force = force_models_buf[0];
    } else {
        composite = propagators.Composite.init(allocator, force_models_buf[0..n_forces]) catch {
            py.raiseRuntime("out of memory creating composite force model");
            return null;
        };
        force = propagators.ForceModel.wrap(propagators.Composite, &composite.?);
    }
    defer if (composite) |*comp| comp.deinit();

    // Build integrator and propagate
    var rk4 = propagators.Rk4{};
    var dp87 = propagators.DormandPrince87.initWithTolerance(rtol, atol);
    const integrator = if (use_dp87) dp87.integrator() else rk4.integrator();

    var prop = propagators.Propagator.init(allocator, integrator, force);
    var trajectory = prop.propagate(initial_state, t0, duration, dt) catch {
        py.raiseRuntime("propagation failed");
        return null;
    };
    defer trajectory.deinit(allocator);

    // Build Python result: (times_list, states_list)
    const n = trajectory.items.len;
    const times_list = c.PyList_New(@intCast(n)) orelse return null;
    const states_list = c.PyList_New(@intCast(n)) orelse {
        c.Py_DECREF(times_list);
        return null;
    };

    for (trajectory.items, 0..) |st, i| {
        const t_obj = c.PyFloat_FromDouble(st.time) orelse {
            c.Py_DECREF(times_list);
            c.Py_DECREF(states_list);
            return null;
        };
        _ = c.PyList_SetItem(times_list, @intCast(i), t_obj);

        const s_obj = py.vecTuple(st.state) orelse {
            c.Py_DECREF(times_list);
            c.Py_DECREF(states_list);
            return null;
        };
        _ = c.PyList_SetItem(states_list, @intCast(i), s_obj);
    }

    const result = py.tuple(2) orelse {
        c.Py_DECREF(times_list);
        c.Py_DECREF(states_list);
        return null;
    };
    py.tupleSet(result, 0, times_list);
    py.tupleSet(result, 1, states_list);
    return result;
}

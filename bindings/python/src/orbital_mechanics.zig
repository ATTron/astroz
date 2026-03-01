//! Python bindings for orbital mechanics functions

const py = @import("python.zig");
const c = py.c;
const astroz = @import("astroz");
const OrbitalMechanics = astroz.OrbitalMechanics;

pub fn pyHohmannTransfer(_: [*c]c.PyObject, args: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    var mu: f64 = undefined;
    var r1: f64 = undefined;
    var r2: f64 = undefined;
    if (c.PyArg_ParseTuple(args, "ddd", &mu, &r1, &r2) == 0) return null;

    var om = OrbitalMechanics.init(mu);
    const result = om.hohmannTransfer(r1, r2) catch {
        py.raiseValue("invalid transfer parameters (radii must be positive and differ by >1000 km)");
        return null;
    };
    return py.structToDict(result, .{ "sma", "dv1", "dv2", "total_dv", "transfer_time", "transfer_time_days" });
}

pub fn pyBiEllipticTransfer(_: [*c]c.PyObject, args: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    var mu: f64 = undefined;
    var r1: f64 = undefined;
    var r2: f64 = undefined;
    var r_aph: f64 = undefined;
    if (c.PyArg_ParseTuple(args, "dddd", &mu, &r1, &r2, &r_aph) == 0) return null;

    var om = OrbitalMechanics.init(mu);
    const result = om.biEllipicTransfer(r1, r2, r_aph) catch {
        py.raiseValue("invalid bi-elliptic transfer parameters");
        return null;
    };
    return py.structToDict(result, .{ "sma", "dv1", "dv2", "dv3", "total_dv", "total_time", "total_time_days" });
}

pub fn pyLambert(_: [*c]c.PyObject, args: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    var mu: f64 = undefined;
    var r1_obj: [*c]c.PyObject = null;
    var r2_obj: [*c]c.PyObject = null;
    var tof: f64 = undefined;
    if (c.PyArg_ParseTuple(args, "dOOd", &mu, &r1_obj, &r2_obj, &tof) == 0) return null;

    const r1 = py.parseVec3(r1_obj) orelse return null;
    const r2 = py.parseVec3(r2_obj) orelse return null;

    const Vector3D = astroz.calculations.Vector3D;
    var om = OrbitalMechanics.init(mu);
    const lr = om.lambertSolverSimple(
        Vector3D.new(r1[0], r1[1], r1[2]),
        Vector3D.new(r2[0], r2[1], r2[2]),
        tof,
    ) catch {
        py.raiseValue("Lambert solver failed (check inputs: tof>0, non-zero position vectors)");
        return null;
    };

    const dict = c.PyDict_New() orelse return null;
    const dv = py.vecTuple(lr.departureVelocity.array()) orelse {
        c.Py_DECREF(dict);
        return null;
    };
    const av = py.vecTuple(lr.arrivalVelocity.array()) orelse {
        c.Py_DECREF(dict);
        c.Py_DECREF(dv);
        return null;
    };

    if (c.PyDict_SetItemString(dict, "departure_velocity", dv) < 0 or
        c.PyDict_SetItemString(dict, "arrival_velocity", av) < 0)
    {
        c.Py_DECREF(dict);
        c.Py_DECREF(dv);
        c.Py_DECREF(av);
        return null;
    }
    c.Py_DECREF(dv);
    c.Py_DECREF(av);

    // Add scalar fields with proper refcount management
    inline for (.{
        .{ "transfer_angle", lr.transferAngle },
        .{ "sma", lr.semiMajorAxis },
        .{ "tof", lr.timeOfFlight },
    }) |entry| {
        const val = c.PyFloat_FromDouble(entry[1]) orelse {
            c.Py_DECREF(dict);
            return null;
        };
        defer c.Py_DECREF(val);
        if (c.PyDict_SetItemString(dict, entry[0], val) < 0) {
            c.Py_DECREF(dict);
            return null;
        }
    }
    return dict;
}

pub fn pyOrbitalVelocity(_: [*c]c.PyObject, args: [*c]c.PyObject, kwargs: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    var mu: f64 = undefined;
    var radius: f64 = undefined;
    var sma_obj: [*c]c.PyObject = null;
    var kwnames = [_:null]?[*:0]const u8{ "mu", "radius", "sma" };
    if (c.PyArg_ParseTupleAndKeywords(args, kwargs, "dd|O", @ptrCast(&kwnames), &mu, &radius, &sma_obj) == 0) return null;

    // optionalFloat returns the default (which we ignore) for None; we need ?f64 here
    var sma: ?f64 = null;
    if (!py.isNone(sma_obj)) {
        sma = c.PyFloat_AsDouble(sma_obj);
        if (sma.? == -1.0 and c.PyErr_Occurred() != null) return null;
    }

    var om = OrbitalMechanics.init(mu);
    const vel = om.orbitalVelocity(radius, sma) catch {
        py.raiseValue("invalid parameters (radius and sma must be positive)");
        return null;
    };
    return c.PyFloat_FromDouble(vel);
}

pub fn pyOrbitalPeriod(_: [*c]c.PyObject, args: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    var mu: f64 = undefined;
    var sma: f64 = undefined;
    if (c.PyArg_ParseTuple(args, "dd", &mu, &sma) == 0) return null;

    var om = OrbitalMechanics.init(mu);
    const period = om.orbitalPeriod(sma) catch {
        py.raiseValue("invalid semi-major axis (must be positive)");
        return null;
    };
    return c.PyFloat_FromDouble(period);
}

pub fn pyEscapeVelocity(_: [*c]c.PyObject, args: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    var mu: f64 = undefined;
    var radius: f64 = undefined;
    if (c.PyArg_ParseTuple(args, "dd", &mu, &radius) == 0) return null;

    var om = OrbitalMechanics.init(mu);
    const vel = om.escapeVelocity(radius) catch {
        py.raiseValue("invalid radius (must be positive)");
        return null;
    };
    return c.PyFloat_FromDouble(vel);
}

//! Python C API helpers for Zig

pub const c = @cImport({
    @cDefine("PY_SSIZE_T_CLEAN", {});
    @cInclude("Python.h");
});

const std = @import("std");

/// Get PyTypeObject from a PyObject - workaround for Py_TYPE being a static inline function in Python 3.10+
pub inline fn pyType(obj: [*c]c.PyObject) ?*c.PyTypeObject {
    if (obj == null) return null;
    // Access ob_type directly from the PyObject structure
    return obj.*.ob_type;
}

/// Manual PyModuleDef since C macro doesn't translate
pub const PyModuleDef_Base = extern struct {
    ob_base: c.PyObject,
    m_init: ?*const fn () callconv(.c) ?*c.PyObject,
    m_index: c.Py_ssize_t,
    m_copy: ?*c.PyObject,
};

pub const PyModuleDef = extern struct {
    m_base: PyModuleDef_Base,
    m_name: ?[*:0]const u8,
    m_doc: ?[*:0]const u8,
    m_size: c.Py_ssize_t,
    m_methods: ?[*]c.PyMethodDef,
    m_slots: ?*c.PyModuleDef_Slot,
    m_traverse: ?*const fn (?*c.PyObject, c.visitproc, ?*anyopaque) callconv(.c) c_int,
    m_clear: ?*const fn (?*c.PyObject) callconv(.c) c_int,
    m_free: ?*const fn (?*anyopaque) callconv(.c) void,
};

pub fn makeBaseObject() c.PyObject {
    // PyModule_Create will properly initialize the object header.
    return std.mem.zeroes(c.PyObject);
}

pub fn moduleCreate(def: *PyModuleDef) ?*c.PyObject {
    return c.PyModule_Create(@as(*c.PyModuleDef, @ptrCast(def)));
}

pub fn none() *c.PyObject {
    c.Py_INCREF(@constCast(&c._Py_NoneStruct));
    return @constCast(&c._Py_NoneStruct);
}

pub fn raiseValue(msg: [:0]const u8) void {
    c.PyErr_SetString(c.PyExc_ValueError, msg.ptr);
}

pub fn raiseType(msg: [:0]const u8) void {
    c.PyErr_SetString(c.PyExc_TypeError, msg.ptr);
}

pub fn raiseRuntime(msg: [:0]const u8) void {
    c.PyErr_SetString(c.PyExc_RuntimeError, msg.ptr);
}

pub fn float(val: f64) ?*c.PyObject {
    return c.PyFloat_FromDouble(val);
}

pub fn int(val: i64) ?*c.PyObject {
    return c.PyLong_FromLongLong(val);
}

pub fn tuple(size: c.Py_ssize_t) ?*c.PyObject {
    return c.PyTuple_New(size);
}

pub fn tupleSet(t: *c.PyObject, idx: c.Py_ssize_t, val: *c.PyObject) void {
    _ = c.PyTuple_SetItem(t, idx, val);
}

/// Check if a PyObject pointer is null or Python None.
pub fn isNone(obj: [*c]c.PyObject) bool {
    return obj == null or obj == @as([*c]c.PyObject, @ptrCast(@constCast(&c._Py_NoneStruct)));
}

pub fn listFromF64(data: []const f64) ?*c.PyObject {
    const list = c.PyList_New(@intCast(data.len)) orelse return null;
    for (data, 0..) |val, i| {
        const obj = c.PyFloat_FromDouble(val) orelse {
            c.Py_DECREF(list);
            return null;
        };
        _ = c.PyList_SetItem(list, @intCast(i), obj);
    }
    return list;
}

pub fn listFromU32(data: []const u32) ?*c.PyObject {
    const list = c.PyList_New(@intCast(data.len)) orelse return null;
    for (data, 0..) |val, i| {
        const obj = c.PyLong_FromUnsignedLong(@as(c_ulong, val)) orelse {
            c.Py_DECREF(list);
            return null;
        };
        _ = c.PyList_SetItem(list, @intCast(i), obj);
    }
    return list;
}

pub fn vecTuple(v: anytype) ?*c.PyObject {
    const N = @typeInfo(@TypeOf(v)).array.len;
    const t = tuple(N) orelse return null;
    inline for (0..N) |i| {
        tupleSet(t, @intCast(i), float(v[i]) orelse {
            c.Py_DECREF(t);
            return null;
        });
    }
    return t;
}

/// Parse an optional float from a PyObject. Returns `default` if None/null,
/// the parsed value on success, or `null` if a Python error occurred.
pub fn optionalFloat(obj: [*c]c.PyObject, default: f64) ?f64 {
    if (isNone(obj)) return default;
    const val = c.PyFloat_AsDouble(obj);
    if (val == -1.0 and c.PyErr_Occurred() != null) return null;
    return val;
}

/// Convert a Zig struct with all-f64 fields to a Python dict.
/// `keys` must be a comptime tuple of string literals, one per struct field (in order).
pub fn structToDict(result: anytype, comptime keys: anytype) ?*c.PyObject {
    const fields = @typeInfo(@TypeOf(result)).@"struct".fields;
    if (keys.len != fields.len) @compileError("key count must match struct field count");
    const dict = c.PyDict_New() orelse return null;
    inline for (fields, 0..) |field, i| {
        const val = c.PyFloat_FromDouble(@field(result, field.name)) orelse {
            c.Py_DECREF(dict);
            return null;
        };
        defer c.Py_DECREF(val);
        if (c.PyDict_SetItemString(dict, keys[i], val) < 0) {
            c.Py_DECREF(dict);
            return null;
        }
    }
    return dict;
}

pub fn parseVec3(obj: [*c]c.PyObject) ?[3]f64 {
    if (obj == null) {
        raiseType("expected a sequence of 3 floats");
        return null;
    }
    if (c.PySequence_Size(obj) != 3) {
        raiseValue("sequence must have exactly 3 elements");
        return null;
    }
    var result: [3]f64 = undefined;
    inline for (0..3) |i| {
        const item = c.PySequence_GetItem(obj, @intCast(i)) orelse return null;
        defer c.Py_DECREF(item);
        result[i] = c.PyFloat_AsDouble(item);
        if (result[i] == -1.0 and c.PyErr_Occurred() != null) return null;
    }
    return result;
}

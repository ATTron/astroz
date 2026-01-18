//! Python C API helpers for Zig

pub const c = @cImport({
    @cDefine("PY_SSIZE_T_CLEAN", {});
    @cInclude("Python.h");
});

const std = @import("std");

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
    var obj = std.mem.zeroes(c.PyObject);
    obj.unnamed_0.ob_refcnt = 1;
    return obj;
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

pub fn vec3Tuple(v: [3]f64) ?*c.PyObject {
    const t = tuple(3) orelse return null;
    tupleSet(t, 0, float(v[0]) orelse {
        c.Py_DECREF(t);
        return null;
    });
    tupleSet(t, 1, float(v[1]) orelse {
        c.Py_DECREF(t);
        return null;
    });
    tupleSet(t, 2, float(v[2]) orelse {
        c.Py_DECREF(t);
        return null;
    });
    return t;
}

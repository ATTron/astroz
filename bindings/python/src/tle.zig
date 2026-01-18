//! Python Tle type - wraps astroz.Tle

const std = @import("std");
const py = @import("python.zig");
const c = py.c;
const astroz = @import("astroz");
const Tle = astroz.Tle;

pub const TleObject = extern struct {
    ob_base: c.PyObject,
    tle: ?*Tle,
};

const allocator = std.heap.c_allocator;
pub var TleType: c.PyTypeObject = undefined;

fn initType() void {
    TleType = std.mem.zeroes(c.PyTypeObject);
    TleType.tp_name = "astroz.Tle";
    TleType.tp_doc = "Two-Line Element set.\n\nArgs:\n    tle_string: TLE string (both lines)";
    TleType.tp_basicsize = @sizeOf(TleObject);
    TleType.tp_flags = c.Py_TPFLAGS_DEFAULT | c.Py_TPFLAGS_BASETYPE;
    TleType.tp_new = @ptrCast(&tle_new);
    TleType.tp_init = @ptrCast(&tle_init);
    TleType.tp_dealloc = @ptrCast(&tle_dealloc);
    TleType.tp_getset = @constCast(&tle_getset);
}

fn tle_new(typ: [*c]c.PyTypeObject, _: [*c]c.PyObject, _: [*c]c.PyObject) callconv(.c) [*c]c.PyObject {
    const self: *TleObject = @ptrCast(@alignCast(c.PyType_GenericAlloc(typ, 0) orelse return null));
    self.tle = null;
    return @ptrCast(self);
}

fn tle_init(self_obj: [*c]c.PyObject, args: [*c]c.PyObject, _: [*c]c.PyObject) callconv(.c) c_int {
    const self: *TleObject = @ptrCast(@alignCast(self_obj));
    var tle_str: [*c]const u8 = undefined;
    if (c.PyArg_ParseTuple(args, "s", &tle_str) == 0) return -1;

    if (self.tle) |old| {
        old.deinit();
        allocator.destroy(old);
    }

    const tle = Tle.parse(std.mem.span(tle_str), allocator) catch {
        py.raiseValue("Failed to parse TLE");
        return -1;
    };
    const ptr = allocator.create(Tle) catch {
        py.raiseRuntime("Out of memory");
        return -1;
    };
    ptr.* = tle;
    self.tle = ptr;
    return 0;
}

fn tle_dealloc(self_obj: [*c]c.PyObject) callconv(.c) void {
    const self: *TleObject = @ptrCast(@alignCast(self_obj));
    if (self.tle) |tle| {
        tle.deinit();
        allocator.destroy(tle);
    }
    if (c.Py_TYPE(self_obj)) |tp| if (tp.*.tp_free) |free| free(@ptrCast(self));
}

// Comptime-generated property getters
fn makeGetter(comptime getter: fn (*const Tle) [*c]c.PyObject) fn ([*c]c.PyObject, ?*anyopaque) callconv(.c) [*c]c.PyObject {
    return struct {
        fn get(self_obj: [*c]c.PyObject, _: ?*anyopaque) callconv(.c) [*c]c.PyObject {
            const self: *TleObject = @ptrCast(@alignCast(self_obj));
            const tle = self.tle orelse {
                py.raiseRuntime("TLE not initialized");
                return null;
            };
            return getter(tle);
        }
    }.get;
}

fn getSatNum(tle: *const Tle) [*c]c.PyObject {
    return py.int(@intCast(tle.firstLine.satelliteNumber));
}
fn getEpoch(tle: *const Tle) [*c]c.PyObject {
    return py.float(tle.firstLine.epoch);
}
fn getInc(tle: *const Tle) [*c]c.PyObject {
    return py.float(tle.secondLine.inclination);
}
fn getEcc(tle: *const Tle) [*c]c.PyObject {
    return py.float(tle.secondLine.eccentricity);
}
fn getMM(tle: *const Tle) [*c]c.PyObject {
    return py.float(tle.secondLine.mMotion);
}

const tle_getset = [_]c.PyGetSetDef{
    .{ .name = "satellite_number", .get = @ptrCast(&makeGetter(getSatNum)), .set = null, .doc = "NORAD catalog number", .closure = null },
    .{ .name = "epoch", .get = @ptrCast(&makeGetter(getEpoch)), .set = null, .doc = "Epoch (J2000 seconds)", .closure = null },
    .{ .name = "inclination", .get = @ptrCast(&makeGetter(getInc)), .set = null, .doc = "Inclination (degrees)", .closure = null },
    .{ .name = "eccentricity", .get = @ptrCast(&makeGetter(getEcc)), .set = null, .doc = "Eccentricity", .closure = null },
    .{ .name = "mean_motion", .get = @ptrCast(&makeGetter(getMM)), .set = null, .doc = "Mean motion (rev/day)", .closure = null },
    .{ .name = null, .get = null, .set = null, .doc = null, .closure = null },
};

pub fn ready() c_int {
    initType();
    return c.PyType_Ready(&TleType);
}

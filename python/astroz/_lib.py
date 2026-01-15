"""Load the astroz shared library."""

import ctypes
import os
import platform


_LIB_NAMES = {
    "Linux": "libastroz_c.so",
    "Darwin": "libastroz_c.dylib",
    "Windows": "astroz_c.dll",
}


def _find_library():
    """Find libastroz_c shared library.

    Search order:
    1. ASTROZ_LIB env var (explicit override)
    2. Next to this package (installed via pip)
    3. zig-out/lib/ (development)
    """
    system = platform.system()
    name = _LIB_NAMES.get(system)
    if name is None:
        raise OSError(f"Unsupported platform: {system}")

    # Explicit override
    env_path = os.environ.get("ASTROZ_LIB")
    if env_path:
        if os.path.exists(env_path):
            return env_path
        raise OSError(f"ASTROZ_LIB set but file not found: {env_path}")

    # Standard search paths
    pkg_dir = os.path.dirname(__file__)
    search_paths = [
        os.path.join(pkg_dir, name),  # installed package
        os.path.join(pkg_dir, "..", "..", "zig-out", "lib", name),  # dev
    ]

    for path in search_paths:
        if os.path.exists(path):
            return path

    raise OSError(
        f"Cannot find {name}. Either:\n"
        f"  - Run 'zig build c-api' for development\n"
        f"  - Set ASTROZ_LIB=/path/to/{name}"
    )


_lib = ctypes.CDLL(_find_library())

# aliases
_ptr = ctypes.c_void_p
_i32 = ctypes.c_int32
_u32 = ctypes.c_uint32
_f32 = ctypes.c_float
_f64 = ctypes.c_double
_str = ctypes.c_char_p
_out_ptr = ctypes.POINTER(ctypes.c_void_p)
_out_vec3 = ctypes.POINTER(ctypes.c_double * 3)


def _sig(func, args, ret=_i32):
    """Set function signature: argtypes and restype."""
    func.argtypes = args
    func.restype = ret


_sig(_lib.astroz_init, [], None)
_sig(_lib.astroz_deinit, [], None)
_sig(_lib.astroz_version, [], _u32)

_sig(_lib.tle_parse, [_str, _out_ptr])
_sig(_lib.tle_free, [_ptr], None)
_sig(_lib.tle_get_satellite_number, [_ptr], _u32)
_sig(_lib.tle_get_epoch, [_ptr], _f64)
_sig(_lib.tle_get_inclination, [_ptr], _f32)
_sig(_lib.tle_get_eccentricity, [_ptr], _f32)
_sig(_lib.tle_get_mean_motion, [_ptr], _f64)

_sig(_lib.sgp4_init, [_ptr, _i32, _out_ptr])
_sig(_lib.sgp4_free, [_ptr], None)
_sig(_lib.sgp4_propagate, [_ptr, _f64, _out_vec3, _out_vec3])


class HohmannResult(ctypes.Structure):
    _fields_ = [
        ("semi_major_axis", _f64),
        ("delta_v1", _f64),
        ("delta_v2", _f64),
        ("total_delta_v", _f64),
        ("transfer_time", _f64),
        ("transfer_time_days", _f64),
    ]


_sig(_lib.orbital_hohmann, [_f64, _f64, _f64, ctypes.POINTER(HohmannResult)])
_sig(_lib.orbital_velocity, [_f64, _f64, _f64], _f64)
_sig(_lib.orbital_period, [_f64, _f64], _f64)
_sig(_lib.orbital_escape_velocity, [_f64, _f64], _f64)

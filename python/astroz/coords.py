"""Coordinate transformations."""

from ctypes import byref, c_double
from ._lib import _lib


def eci_to_ecef(eci: tuple, gmst: float) -> tuple:
    """ECI to ECEF. Returns (x, y, z) in km."""
    eci_arr = (c_double * 3)(*eci)
    ecef_arr = (c_double * 3)()
    _lib.coords_eci_to_ecef(byref(eci_arr), gmst, byref(ecef_arr))
    return (ecef_arr[0], ecef_arr[1], ecef_arr[2])


def ecef_to_geodetic(ecef: tuple) -> tuple:
    """ECEF to geodetic (WGS84). Returns (lat_deg, lon_deg, alt_km)."""
    ecef_arr = (c_double * 3)(*ecef)
    lla_arr = (c_double * 3)()
    _lib.coords_ecef_to_geodetic(byref(ecef_arr), byref(lla_arr))
    return (lla_arr[0], lla_arr[1], lla_arr[2])


def eci_to_geodetic(eci: tuple, gmst: float) -> tuple:
    """ECI to geodetic (WGS84). Returns (lat_deg, lon_deg, alt_km)."""
    return ecef_to_geodetic(eci_to_ecef(eci, gmst))


def julian_to_gmst(jd: float) -> float:
    """Julian date to GMST (radians)."""
    return _lib.coords_julian_to_gmst(jd)

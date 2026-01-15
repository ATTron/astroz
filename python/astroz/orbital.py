"""Orbital mechanics functions."""

from ctypes import byref
from ._lib import _lib, HohmannResult
from .exceptions import check


def hohmann(mu: float, r1: float, r2: float) -> dict:
    """Hohmann transfer between two circular orbits.

    Returns dict with: semi_major_axis, delta_v1, delta_v2,
    total_delta_v, transfer_time, transfer_time_days
    """
    result = HohmannResult()
    check(_lib.orbital_hohmann(mu, r1, r2, byref(result)))
    return {
        "semi_major_axis": result.semi_major_axis,
        "delta_v1": result.delta_v1,
        "delta_v2": result.delta_v2,
        "total_delta_v": result.total_delta_v,
        "transfer_time": result.transfer_time,
        "transfer_time_days": result.transfer_time_days,
    }


def velocity(mu: float, radius: float, sma: float = 0) -> float:
    """Orbital velocity (km/s). Use sma=0 for circular orbit."""
    return _lib.orbital_velocity(mu, radius, sma)


def period(mu: float, sma: float) -> float:
    """Orbital period (seconds)."""
    return _lib.orbital_period(mu, sma)


def escape_velocity(mu: float, radius: float) -> float:
    """Escape velocity (km/s)."""
    return _lib.orbital_escape_velocity(mu, radius)

"""TLE (Two-Line Element) wrapper."""

import ctypes
from ._lib import _lib
from .exceptions import check


class Tle:
    """Two-Line Element set for satellite orbit data."""

    def __init__(self, tle_string: str):
        """Parse a TLE string (both lines)."""
        self._handle = ctypes.c_void_p()
        check(_lib.tle_parse(tle_string.encode(), ctypes.byref(self._handle)))

    def __del__(self):
        if hasattr(self, "_handle") and self._handle:
            _lib.tle_free(self._handle)

    @property
    def satellite_number(self) -> int:
        """NORAD catalog number."""
        return _lib.tle_get_satellite_number(self._handle)

    @property
    def epoch(self) -> float:
        """Epoch in J2000 seconds."""
        return _lib.tle_get_epoch(self._handle)

    @property
    def inclination(self) -> float:
        """Orbital inclination in degrees."""
        return _lib.tle_get_inclination(self._handle)

    @property
    def eccentricity(self) -> float:
        """Orbital eccentricity."""
        return _lib.tle_get_eccentricity(self._handle)

    @property
    def mean_motion(self) -> float:
        """Mean motion in revolutions per day."""
        return _lib.tle_get_mean_motion(self._handle)

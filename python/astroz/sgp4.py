"""SGP4 orbit propagator wrapper."""

import ctypes
from ._lib import _lib
from .exceptions import check
from .tle import Tle


class Sgp4:
    """SGP4 orbit propagator."""

    WGS84 = 0
    WGS72 = 1

    def __init__(self, tle: Tle, gravity_model: int = WGS84):
        """Initialize SGP4 propagator from TLE.

        Args:
            tle: TLE object
            gravity_model: WGS84 (default) or WGS72
        """
        self._handle = ctypes.c_void_p()
        check(_lib.sgp4_init(tle._handle, gravity_model, ctypes.byref(self._handle)))

    def __del__(self):
        if hasattr(self, "_handle") and self._handle:
            _lib.sgp4_free(self._handle)

    def propagate(self, tsince: float) -> tuple:
        """Propagate to time since TLE epoch.

        Args:
            tsince: Minutes since TLE epoch

        Returns:
            (position, velocity) where each is (x, y, z) tuple
            Position in km, velocity in km/s (TEME frame)
        """
        pos = (ctypes.c_double * 3)()
        vel = (ctypes.c_double * 3)()
        check(
            _lib.sgp4_propagate(
                self._handle,
                tsince,
                ctypes.byref(pos),
                ctypes.byref(vel),
            )
        )
        return (pos[0], pos[1], pos[2]), (vel[0], vel[1], vel[2])

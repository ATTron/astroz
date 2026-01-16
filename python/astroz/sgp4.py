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
        """Initialize SGP4 propagator from TLE."""
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
            (position, velocity) tuples in km and km/s (TEME frame)
        """
        pos = (ctypes.c_double * 3)()
        vel = (ctypes.c_double * 3)()
        check(
            _lib.sgp4_propagate(
                self._handle, tsince, ctypes.byref(pos), ctypes.byref(vel)
            )
        )
        return (pos[0], pos[1], pos[2]), (vel[0], vel[1], vel[2])

    def propagate_batch(self, times: list) -> list:
        """Batch propagation - much faster for multiple time steps.

        Args:
            times: List of minutes since TLE epoch

        Returns:
            List of (position, velocity) tuples
        """
        n = len(times)
        times_arr = (ctypes.c_double * n)(*times)
        results_arr = (ctypes.c_double * (n * 6))()
        check(_lib.sgp4_propagate_batch(self._handle, times_arr, results_arr, n))
        out = []
        for i in range(n):
            base = i * 6
            pos = (results_arr[base], results_arr[base + 1], results_arr[base + 2])
            vel = (results_arr[base + 3], results_arr[base + 4], results_arr[base + 5])
            out.append((pos, vel))
        return out

    def propagate_batch_np(self, times):
        """Batch propagation returning numpy arrays (fastest).

        Args:
            times: Array-like of minutes since TLE epoch

        Returns:
            (positions, velocities) as numpy arrays of shape (n, 3)
        """
        import numpy as np

        times = np.asarray(times, dtype=np.float64)
        n = len(times)
        results = np.empty(n * 6, dtype=np.float64)
        check(
            _lib.sgp4_propagate_batch(
                self._handle,
                times.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                results.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                n,
            )
        )
        results = results.reshape(n, 6)
        return results[:, :3], results[:, 3:]

"""astroz - High-performance astrodynamics library.
Native Python bindings with zero FFI overhead.
"""

import numpy as np
from ._astroz import (
    Tle,
    Sgp4 as _Sgp4,
    Sgp4Batch,
    Sgp4Constellation,
    WGS84,
    WGS72,
    version,
)


class Sgp4(_Sgp4):
    """SGP4 orbit propagator with SIMD acceleration.
    Extends native Sgp4 with convenience methods.
    """

    def propagate_batch(self, times):
        """Batch propagation returning position and velocity arrays.
        Args:
            times: float64 array (n,) of minutes since TLE epoch
        Returns:
            (positions, velocities) tuple of float64 arrays (n, 3)
        """
        n = len(times)
        pos = np.empty((n, 3), dtype=np.float64)
        vel = np.empty((n, 3), dtype=np.float64)
        self.propagate_into(times, pos, vel)
        return pos, vel


__version__ = "0.3.0"

__all__ = [
    "__version__",
    "version",
    "Tle",
    "Sgp4",
    "Sgp4Batch",
    "Sgp4Constellation",
    "WGS84",
    "WGS72",
]

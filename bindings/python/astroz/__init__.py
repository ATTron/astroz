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
    min_distances as _min_distances,
    coarse_screen as _coarse_screen,
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
        times = np.asarray(times, dtype=np.float64)
        n = len(times)
        pos = np.empty((n, 3), dtype=np.float64)
        vel = np.empty((n, 3), dtype=np.float64)
        self.propagate_into(times, pos, vel)
        return pos, vel


def propagate_constellation(
    constellation,
    times,
    *,
    epoch_offsets=None,
    satellite_mask=None,
    output="teme",
    reference_jd=0.0,
    velocities=False,
):
    """Propagate constellation, allocating output arrays automatically.

    Args:
        constellation: Sgp4Constellation object
        times: float64 array of time offsets in minutes
        epoch_offsets: float64 array (num_satellites,) of per-satellite
            time offsets in minutes. Padding to SIMD width is handled
            internally.
        satellite_mask: uint8 array (num_satellites,) where non-zero means
            propagate. Padding is handled internally.
        output: coordinate frame - 'teme', 'ecef', or 'geodetic'
        reference_jd: reference Julian date for ECEF/geodetic conversion
        velocities: if True, also compute velocities

    Returns:
        positions array (num_sats, num_times, 3) if velocities=False
        (positions, velocities) tuple if velocities=True
    """
    times = np.ascontiguousarray(times, dtype=np.float64)
    n_sats = constellation.num_satellites
    shape = (n_sats, len(times), 3)

    pos = np.empty(shape, dtype=np.float64)
    vel = np.empty(shape, dtype=np.float64) if velocities else None

    constellation.propagate_into(
        times,
        pos,
        vel,
        epoch_offsets=epoch_offsets,
        satellite_mask=satellite_mask,
        output=output,
        reference_jd=reference_jd,
    )
    return (pos, vel) if velocities else pos


def min_distances(positions, pairs, velocities=None, num_sats=None):
    """Compute minimum pairwise distances between satellites across time steps.

    Args:
        positions: float64 array (num_sats, num_times, 3) - satellite positions
        pairs: uint32 array (num_pairs, 2) - satellite index pairs to check
        velocities: optional float64 array (num_sats, num_times, 3)
        num_sats: number of satellites (inferred from positions shape if 3D)

    Returns:
        Without velocities: (min_dists, min_indices) numpy arrays
        With velocities: (min_dists, min_indices, rel_vels) numpy arrays
    """
    positions = np.ascontiguousarray(positions, dtype=np.float64)
    pairs = np.ascontiguousarray(pairs, dtype=np.uint32)

    if num_sats is None and positions.ndim == 3:
        num_sats = positions.shape[0]

    if velocities is not None:
        velocities = np.ascontiguousarray(velocities, dtype=np.float64)
        result = _min_distances(positions, pairs, velocities, num_sats)
        return (
            np.array(result[0]),
            np.array(result[1], dtype=np.uint32),
            np.array(result[2]),
        )
    else:
        result = _min_distances(positions, pairs, None, num_sats)
        return np.array(result[0]), np.array(result[1], dtype=np.uint32)


def coarse_screen(positions, num_sats=None, threshold=10.0, valid_mask=None):
    """Find all satellite pairs within threshold distance at any time step.

    Uses a spatial hash grid for O(n) neighbor search per time step.

    Args:
        positions: float64 array (num_sats, num_times, 3) - satellite positions
        num_sats: number of satellites (inferred from positions shape if 3D)
        threshold: distance threshold in km (default: 10.0)
        valid_mask: optional uint8 array (num_sats,) - 0 means skip satellite

    Returns:
        (pairs, t_indices) where:
            pairs: list of (i, j) tuples - satellite index pairs
            t_indices: list of int - time step index for each pair detection
    """
    positions = np.ascontiguousarray(positions, dtype=np.float64)

    if num_sats is None and positions.ndim == 3:
        num_sats = positions.shape[0]
    elif num_sats is None:
        raise ValueError("num_sats must be provided for non-3D position arrays")

    if valid_mask is not None:
        valid_mask = np.ascontiguousarray(valid_mask, dtype=np.uint8)

    return _coarse_screen(positions, num_sats, threshold, valid_mask)


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
    "propagate_constellation",
    "min_distances",
    "coarse_screen",
]

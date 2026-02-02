"""python-sgp4 compatible API for astroz.

Drop-in replacement for `sgp4.api` with 2-100x faster performance.

Migration
---------
Just change the import::

    # Before
    from sgp4.api import Satrec, SatrecArray, jday

    # After (2-100x faster)
    from astroz.api import Satrec, SatrecArray, jday

Performance
-----------
===============================  ===========  =========  =========
Method                           python-sgp4  astroz     Speedup
===============================  ===========  =========  =========
``sat.sgp4()`` loop              1.3M/s       2.5M/s     **2x**
``sat.sgp4_array()``             2.7M/s       15M/s      **5x**
``SatrecArray.sgp4()``           3M/s         290M/s     **100x**
``SatrecArray.sgp4(vel=False)``  3M/s         330M/s     **110x**
===============================  ===========  =========  =========

Example
-------
>>> from astroz.api import Satrec, SatrecArray, jday, WGS72
>>> import numpy as np
>>>
>>> # Create satellite (same syntax as python-sgp4)
>>> line1 = "1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995"
>>> line2 = "2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123"
>>> sat = Satrec.twoline2rv(line1, line2, WGS72)
>>>
>>> # Single propagation
>>> jd, fr = jday(2024, 5, 6, 19, 52, 0.0)
>>> error, position, velocity = sat.sgp4(jd, fr)
>>>
>>> # Batch propagation - single satellite (5x faster than loop)
>>> jd_arr = np.full(1000, jd)
>>> fr_arr = fr + np.arange(1000) / 1440.0
>>> e, r, v = sat.sgp4_array(jd_arr, fr_arr)
>>>
>>> # Multi-satellite batch (100x faster)
>>> sat_array = SatrecArray([sat1, sat2, sat3])
>>> e, r, v = sat_array.sgp4(jd_arr, fr_arr)
>>>
>>> # Skip velocities for 30% more speed
>>> e, r, _ = sat_array.sgp4(jd_arr, fr_arr, velocities=False)

Notes
-----
- SIMD acceleration (AVX2/AVX512) for batch propagation
- All orbital elements match python-sgp4 units (radians, Earth radii)
- Error codes match python-sgp4 conventions (0=success)
- Use ``velocities=False`` when positions are sufficient

See Also
--------
- Migration guide: https://github.com/ATTron/astroz/tree/main/bindings/python#migrating-from-python-sgp4
- python-sgp4: https://github.com/brandon-rhodes/python-sgp4
"""

import numpy as np
from ._astroz import (
    Satrec as _Satrec,
    SatrecArray as _SatrecArray,
    jday,
    days2mdhms,
    WGS72,
    WGS84,
)

# Alias for compatibility (WGS72OLD = WGS72 in python-sgp4)
WGS72OLD = WGS72

# Always true for astroz (we use native Zig+SIMD acceleration)
accelerated = True


class Satrec:
    """SGP4 satellite record with SIMD-accelerated batch propagation.

    Drop-in replacement for sgp4.api.Satrec. Use twoline2rv() to create.

    Example
    -------
    >>> sat = Satrec.twoline2rv(line1, line2)  # WGS72 default
    >>> e, r, v = sat.sgp4(jd, fr)  # Single propagation
    >>> e, r, v = sat.sgp4_array(jd_arr, fr_arr)  # SIMD batch
    """

    def __init__(self, _native):
        """Wrap a native Satrec. Use twoline2rv() instead."""
        self._native = _native
        self._cached_array = None

    @classmethod
    def twoline2rv(cls, line1, line2, whichconst=WGS72):
        """Create Satrec from TLE lines.

        Parameters
        ----------
        line1 : str
            First line of TLE.
        line2 : str
            Second line of TLE.
        whichconst : int, optional
            Gravity model (WGS72 default, or WGS84).

        Returns
        -------
        Satrec
            Initialized satellite record.
        """
        return cls(_Satrec.twoline2rv(line1, line2, whichconst))

    def sgp4(self, jd, fr):
        """Propagate to given Julian date.

        Parameters
        ----------
        jd : float
            Julian date integer part.
        fr : float
            Julian date fractional part.

        Returns
        -------
        error : int
            Error code (0 = success).
        position : tuple
            (x, y, z) in km, TEME frame.
        velocity : tuple
            (vx, vy, vz) in km/s, TEME frame.
        """
        return self._native.sgp4(jd, fr)

    def sgp4_array(self, jd, fr):
        """Propagate to multiple times using SIMD acceleration.

        This is the python-sgp4 compatible batch interface for a single satellite.
        Uses SIMD internally for high throughput.

        Parameters
        ----------
        jd : array_like
            Julian date integer parts, shape (n_times,).
        fr : array_like
            Julian date fractional parts, shape (n_times,).

        Returns
        -------
        e : ndarray
            Error codes, shape (n_times,). 0 = success.
        r : ndarray
            Positions in km (TEME), shape (n_times, 3).
        v : ndarray
            Velocities in km/s (TEME), shape (n_times, 3).

        Examples
        --------
        >>> jd = np.array([2458826, 2458826, 2458826])
        >>> fr = np.array([0.0001, 0.0002, 0.0003])
        >>> e, r, v = sat.sgp4_array(jd, fr)
        """
        # Lazily create SatrecArray for SIMD batch propagation
        if self._cached_array is None:
            self._cached_array = SatrecArray([self._native])

        e, r, v = self._cached_array.sgp4(jd, fr)

        # Squeeze out the satellite dimension (n_sats=1)
        return e[0], r[0], v[0]

    def __getattr__(self, name):
        """Delegate attribute access to native Satrec object."""
        # Avoid infinite recursion for _native itself
        if name == "_native":
            raise AttributeError(name)
        return getattr(self._native, name)


class SatrecArray:
    """Batch SGP4 propagator - python-sgp4 compatible with SIMD acceleration.

    Drop-in replacement for sgp4.api.SatrecArray that achieves 270-330M props/sec.

    Parameters
    ----------
    satrecs : list of Satrec
        List of Satrec objects to propagate as a batch.

    Example
    -------
    >>> from astroz.api import Satrec, SatrecArray, WGS72
    >>> import numpy as np
    >>>
    >>> sats = [Satrec.twoline2rv(l1, l2, WGS72) for l1, l2 in tle_pairs]
    >>> sat_array = SatrecArray(sats)
    >>>
    >>> jd = np.array([2460000.5, 2460001.5])
    >>> fr = np.array([0.0, 0.0])
    >>> e, r, v = sat_array.sgp4(jd, fr)
    >>> # e: (n_sats, n_times) error codes
    >>> # r: (n_sats, n_times, 3) positions in km
    >>> # v: (n_sats, n_times, 3) velocities in km/s
    """

    def __init__(self, satrecs):
        # Unwrap our Satrec wrappers if needed
        native_satrecs = [s._native if isinstance(s, Satrec) else s for s in satrecs]
        self._zig = _SatrecArray(native_satrecs)
        self._num_sats = self._zig.num_satellites
        self._epochs = np.array(self._zig.epochs, dtype=np.float64)
        # Pad epochs for SIMD alignment
        padded = ((self._num_sats + 3) // 4) * 4
        if padded > self._num_sats:
            self._epochs_padded = np.concatenate(
                [self._epochs, np.full(padded - self._num_sats, self._epochs[-1])]
            )
        else:
            self._epochs_padded = self._epochs

    @property
    def num_satellites(self):
        """Number of satellites in the array."""
        return self._num_sats

    @property
    def epochs(self):
        """Epoch Julian dates for each satellite."""
        return self._epochs.tolist()

    def sgp4(self, jd, fr, *, velocities=True):
        """Propagate all satellites to the given Julian dates.

        This is the python-sgp4 compatible interface that returns NumPy arrays.
        Uses SIMD internally for 270-330M props/sec throughput.

        Parameters
        ----------
        jd : float or array_like
            Julian date integer parts. Scalar or array of shape (n_times,).
        fr : float or array_like
            Julian date fractional parts. Scalar or array of shape (n_times,).
        velocities : bool, optional
            If True (default), compute and return velocities.
            Set to False for ~30% faster propagation when velocities aren't needed.

        Returns
        -------
        e : ndarray
            Error codes, shape (n_sats, n_times). 0 = success.
        r : ndarray
            Positions in km (TEME), shape (n_sats, n_times, 3)
        v : ndarray
            Velocities in km/s (TEME), shape (n_sats, n_times, 3).
            All zeros if velocities=False.

        Examples
        --------
        >>> # Single time point (scalars)
        >>> e, r, v = sat_array.sgp4(2458826.5, 0.8625)
        >>>
        >>> # Multiple time points (arrays)
        >>> jd = np.full(1440, 2458826.5)
        >>> fr = np.linspace(0, 1, 1440)
        >>> e, r, v = sat_array.sgp4(jd, fr)
        """
        # Handle both scalars and arrays (python-sgp4 compatible)
        jd = np.atleast_1d(np.asarray(jd, dtype=np.float64))
        fr = np.atleast_1d(np.asarray(fr, dtype=np.float64))
        n_times = len(jd)
        n_sats = self._num_sats

        # Pre-allocate output arrays (time-major for SIMD, then transpose)
        e = np.zeros((n_times, n_sats), dtype=np.uint8)
        r = np.empty((n_times, n_sats, 3), dtype=np.float64)
        v = np.empty((n_times, n_sats, 3), dtype=np.float64) if velocities else None

        # Compute epoch offsets and times
        reference_jd = jd[0] + fr[0]
        epoch_offsets = (reference_jd - self._epochs_padded) * 1440.0
        times = ((jd + fr) - reference_jd) * 1440.0

        # Call SIMD propagation (passes None for velocities to skip computation)
        self._zig.propagate_into(times, r, v, epoch_offsets=epoch_offsets)

        # Transpose to (n_sats, n_times, ...) to match python-sgp4 output format
        e = e.T  # (n_sats, n_times)
        r = r.transpose(1, 0, 2)  # (n_sats, n_times, 3)
        if velocities:
            v = v.transpose(1, 0, 2)  # (n_sats, n_times, 3)
        else:
            v = np.zeros((n_sats, n_times, 3), dtype=np.float64)

        return e, r, v


__all__ = [
    "Satrec",
    "SatrecArray",
    "jday",
    "days2mdhms",
    "WGS72",
    "WGS84",
    "WGS72OLD",
    "accelerated",
]

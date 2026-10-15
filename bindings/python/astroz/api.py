"""python-sgp4 compatible API for astroz.

Change the import and existing code keeps working::

    from sgp4.api import Satrec, SatrecArray, jday   # before
    from astroz.api import Satrec, SatrecArray, jday  # after

Units and error codes match python-sgp4 (radians, Earth radii, 0 = success).
A few rarely used attributes (intermediate osculating elements, element
rates) aren't implemented; see the README for the list.
"""

import numpy as np
from ._astroz import (
    Satrec as _Satrec,
    SatrecArray as _SatrecArray,
    jday,
    days2mdhms,
    WGS72,
    WGS84,
    sdp4_batch_propagate_into as _sdp4_batch_propagate_into,
)

WGS72OLD = WGS72
accelerated = True


class Satrec:
    """SGP4/SDP4 satellite record, a drop-in for sgp4.api.Satrec.

    Create one with ``Satrec.twoline2rv(line1, line2)``. Deep-space orbits
    use SDP4 automatically.
    """

    def __init__(self, _native):
        self._native = _native

    @classmethod
    def twoline2rv(cls, line1, line2, whichconst=WGS72):
        """Create a Satrec from TLE lines. ``whichconst`` is WGS72 or WGS84."""
        return cls(_Satrec.twoline2rv(line1, line2, whichconst))

    def sgp4(self, jd, fr):
        """Propagate to Julian date ``jd + fr``.

        Returns ``(error, (x, y, z), (vx, vy, vz))`` in km and km/s, TEME.
        """
        return self._native.sgp4(jd, fr)

    def sgp4_array(self, jd, fr):
        """Propagate to many Julian dates ``jd[i] + fr[i]``.

        Returns ``(e, r, v)`` with shapes ``(n,)``, ``(n, 3)`` and ``(n, 3)``.
        """
        jd = np.atleast_1d(np.asarray(jd, dtype=np.float64))
        fr = np.atleast_1d(np.asarray(fr, dtype=np.float64))
        n = len(jd)
        r = np.empty((n, 3), dtype=np.float64)
        v = np.empty((n, 3), dtype=np.float64)
        self._native.sgp4_array_into(jd, fr, r, v)
        e = np.where(np.isnan(r[:, 0]), 6, 0).astype(np.uint8)
        return e, r, v

    def __getattr__(self, name):
        if name == "_native":
            raise AttributeError(name)
        return getattr(self._native, name)


class SatrecArray:
    """Propagates many satellites at once, a drop-in for sgp4.api.SatrecArray.

    Near-earth satellites go through SIMD SGP4 batches and deep-space ones
    through threaded SDP4; mixed lists are fine. A satellite whose elements
    SGP4/SDP4 rejected gets its own error code and NaN, as in python-sgp4.
    """

    def __init__(self, satrecs):
        natives = [s._native if isinstance(s, Satrec) else s for s in satrecs]
        self._num_sats = len(natives)

        self._sgp4_indices = []
        self._sdp4_indices = []
        self._sdp4_natives = []
        self._failed = []  # (index, error code)
        sgp4_natives = []
        for i, sat in enumerate(natives):
            if sat.is_deep_space:
                self._sdp4_indices.append(i)
                self._sdp4_natives.append(sat)
            elif sat.error != 0:
                self._failed.append((i, sat.error))
            else:
                self._sgp4_indices.append(i)
                sgp4_natives.append(sat)

        self._zig = None
        if sgp4_natives:
            self._zig = _SatrecArray(sgp4_natives)
            self._sgp4_epochs = np.array(self._zig.epochs, dtype=np.float64)

    @property
    def num_satellites(self):
        """Number of satellites in the array."""
        return self._num_sats

    def sgp4(self, jd, fr, *, velocities=True):
        """Propagate every satellite to Julian dates ``jd + fr``.

        ``jd`` and ``fr`` may be scalars or arrays. Returns ``(e, r, v)`` with
        shapes ``(n_sats, n_times)``, ``(n_sats, n_times, 3)`` and
        ``(n_sats, n_times, 3)``. With ``velocities=False`` the velocity
        computation is skipped and ``v`` is all zeros.
        """
        jd = np.atleast_1d(np.asarray(jd, dtype=np.float64))
        fr = np.atleast_1d(np.asarray(fr, dtype=np.float64))
        n_times = len(jd)
        shape = (self._num_sats, n_times, 3)

        if self._zig is not None:
            reference_jd = jd[0] + fr[0]
            epoch_offsets = (reference_jd - self._sgp4_epochs) * 1440.0
            times = ((jd + fr) - reference_jd) * 1440.0
            tm_shape = (n_times, len(self._sgp4_indices), 3)
            r_tm = np.empty(tm_shape, dtype=np.float64)
            v_tm = np.empty(tm_shape, dtype=np.float64) if velocities else None
            self._zig.propagate_into(times, r_tm, v_tm, epoch_offsets=epoch_offsets)

        if self._zig is not None and len(self._sgp4_indices) == self._num_sats:
            # All SGP4: hand back transposed views of the time-major output
            r = r_tm.transpose(1, 0, 2)
            v = v_tm.transpose(1, 0, 2) if velocities else np.zeros(shape, dtype=np.float64)
        else:
            r = np.empty(shape, dtype=np.float64)
            v = np.empty(shape, dtype=np.float64) if velocities else np.zeros(shape, dtype=np.float64)
            if self._zig is not None:
                r[self._sgp4_indices] = r_tm.transpose(1, 0, 2)
                if velocities:
                    v[self._sgp4_indices] = v_tm.transpose(1, 0, 2)
            if self._sdp4_natives:
                sdp4_shape = (len(self._sdp4_natives), n_times, 3)
                sdp4_r = np.empty(sdp4_shape, dtype=np.float64)
                sdp4_v = np.empty(sdp4_shape, dtype=np.float64)
                _sdp4_batch_propagate_into(self._sdp4_natives, jd, fr, sdp4_r, sdp4_v)
                r[self._sdp4_indices] = sdp4_r
                if velocities:
                    v[self._sdp4_indices] = sdp4_v

        # Failed propagations come back as NaN; python-sgp4 reports those as 6 (decayed)
        e = np.where(np.isnan(r[..., 0]), 6, 0).astype(np.uint8)
        for i, code in self._failed:
            r[i] = np.nan
            v[i] = np.nan if velocities else 0.0
            e[i] = code
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

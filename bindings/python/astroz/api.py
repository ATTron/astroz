"""python-sgp4 compatible API for astroz.

Change the import and existing code keeps working::

    from sgp4.api import Satrec, SatrecArray, jday   # before
    from astroz.api import Satrec, SatrecArray, jday  # after

``sgp4.exporter`` and ``sgp4.omm`` have counterparts in ``astroz.exporter`` and
``astroz.omm``. Units, gravity constants and error codes match python-sgp4. A
few rarely used attributes (osculating elements after propagation, element
rates) aren't implemented; see the README for the list.
"""

import numpy as np
from ._astroz import (
    Satrec as _Satrec,
    SatrecArray as _SatrecArray,
    jday,
    days2mdhms,
    WGS72OLD,
    WGS72,
    WGS84,
    sdp4_batch_propagate_into as _sdp4_batch_propagate_into,
)

accelerated = True


def _fill_errors(natives, jd, fr, r, v, velocities=True):
    """Error codes for a batch result, as python-sgp4 reports them.

    The SIMD paths only return NaN when a propagation fails, so rerun just those
    (satellite, time) entries one at a time: that gives the real code, and for a
    decayed satellite (6) the position python-sgp4 still returns. ``r`` and ``v``
    are (n_sats, n_times, 3) and are updated in place.
    """
    e = np.zeros(r.shape[:2], dtype=np.uint8)
    for i, k in np.argwhere(np.isnan(r[..., 0])):
        code, rk, vk = natives[i].sgp4(jd[k], fr[k])
        e[i, k] = code
        r[i, k] = rk
        if velocities:
            v[i, k] = vk
    return e


def _time_arrays(jd, fr):
    """jd and fr as equal-length, contiguous float64 arrays (either may be a scalar)."""
    jd, fr = np.broadcast_arrays(np.asarray(jd, dtype=np.float64), np.asarray(fr, dtype=np.float64))
    return np.ascontiguousarray(np.atleast_1d(jd)), np.ascontiguousarray(np.atleast_1d(fr))


class Satrec:
    """SGP4/SDP4 satellite record, a drop-in for sgp4.api.Satrec.

    Create one with ``Satrec.twoline2rv(line1, line2)``, or ``Satrec()`` then
    ``sgp4init(...)`` from mean elements. Deep-space orbits use SDP4 automatically.
    """

    def __init__(self, _native=None):
        self._native = _native if _native is not None else _Satrec()

    @classmethod
    def twoline2rv(cls, line1, line2, whichconst=WGS72):
        """Create a Satrec from TLE lines. ``whichconst`` is WGS72OLD, WGS72 or WGS84."""
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
        jd, fr = _time_arrays(jd, fr)
        n = len(jd)
        r = np.empty((n, 3), dtype=np.float64)
        v = np.empty((n, 3), dtype=np.float64)
        self._native.sgp4_array_into(jd, fr, r, v)
        e = _fill_errors([self._native], jd, fr, r[None], v[None])[0]
        return e, r, v

    def __getattr__(self, name):
        if name == "_native":
            raise AttributeError(name)
        return getattr(self._native, name)

    def __setattr__(self, name, value):
        if name == "_native":
            object.__setattr__(self, name, value)
        else:
            setattr(self._native, name, value)


class SatrecArray:
    """Propagates many satellites at once, a drop-in for sgp4.api.SatrecArray.

    Near-earth satellites go through SIMD SGP4 batches and deep-space ones
    through threaded SDP4; mixed lists are fine. Error codes match python-sgp4,
    including for a satellite whose elements SGP4/SDP4 rejected.
    """

    def __init__(self, satrecs):
        natives = [s._native if isinstance(s, Satrec) else s for s in satrecs]
        self._natives = natives
        self._num_sats = len(natives)

        self._sgp4_indices = []
        self._sdp4_indices = []
        self._sdp4_natives = []
        sgp4_natives = []
        for i, sat in enumerate(natives):
            if sat.is_deep_space:
                self._sdp4_indices.append(i)
                self._sdp4_natives.append(sat)
            elif sat.error == 0:  # rows of satellites rejected at init stay NaN
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
        jd, fr = _time_arrays(jd, fr)
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
            # NaN rows stay for satellites rejected at init; _fill_errors codes them
            r = np.full(shape, np.nan)
            v = np.full(shape, np.nan) if velocities else np.zeros(shape, dtype=np.float64)
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

        e = _fill_errors(self._natives, jd, fr, r, v, velocities)
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

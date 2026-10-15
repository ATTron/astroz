"""astroz: satellite propagation, conjunction screening and orbital mechanics.

There are two entry points:

- this module: ``propagate`` and ``screen`` take a ``Constellation`` or
  anything it can load (a CelesTrak group, NORAD IDs, a URL, a file, or raw
  TLE / OMM JSON text)::

      import astroz
      positions = astroz.propagate("starlink", [0, 60, 120])

- ``astroz.api``: a drop-in for ``sgp4.api`` (``Satrec``, ``SatrecArray``,
  ``jday``).

OMM JSON is detected automatically.

Docs: https://attron.github.io/astroz/
"""

import math
from typing import NamedTuple

import numpy as np
from datetime import datetime, timezone
from pathlib import Path
from ._astroz import (
    Tle,
    __version__,
    Satrec as _Satrec,
    Sgp4Constellation as _Sgp4Constellation,
    WGS72 as _WGS72,
    screen_candidates as _screen_candidates,
    sdp4_batch_propagate_into as _sdp4_batch_propagate_into,
    hohmann_transfer,
    bi_elliptic_transfer,
    lambert,
    orbital_velocity,
    orbital_period,
    escape_velocity,
    propagate_numerical,
    EARTH_MU,
    EARTH_R_EQ,
    EARTH_J2,
    SUN_MU,
    MOON_MU,
)

_CELESTRAK_ALIASES = {
    "all": "active",
    "iss": "stations",
    "gps": "gps-ops",
    "glonass": "glo-ops",
}


def _fetch_url(url):
    import urllib.request

    req = urllib.request.Request(url, headers={"User-Agent": "astroz"})
    return urllib.request.urlopen(req, timeout=60).read().decode("utf-8")


def _celestrak_url(group=None, norad_id=None):
    if norad_id is not None:
        ids = ",".join(map(str, norad_id)) if isinstance(norad_id, (list, tuple)) else str(norad_id)
        query = f"CATNR={ids}"
    else:
        query = f"GROUP={_CELESTRAK_ALIASES.get(group.lower(), group)}"
    return f"https://celestrak.org/NORAD/elements/gp.php?{query}&FORMAT=TLE"


def _load_tle_text(source, norad_id=None):
    """Returns (text, format) where format is "tle" or "json"."""
    if norad_id is not None:
        return _fetch_url(_celestrak_url(norad_id=norad_id)), "tle"
    if source is None:
        raise ValueError("Must specify 'source' or 'norad_id'")
    if source.startswith(("http://", "https://")):
        text = _fetch_url(source)
        fmt = "json" if text.lstrip().startswith(("[", "{")) else "tle"
        return text, fmt
    if "1 " in source and "2 " in source:
        return source, "tle"
    if source.lstrip().startswith(("[", "{")):
        return source, "json"
    if Path(source).exists():
        text = Path(source).read_text()
        fmt = "json" if text.lstrip().startswith(("[", "{")) else "tle"
        return text, fmt
    return _fetch_url(_celestrak_url(group=source)), "tle"


def _parse_tle_pairs(tle_text):
    """(line1, line2) pairs from TLE text; name lines and strays are skipped."""
    lines = [l.strip() for l in tle_text.strip().splitlines() if l.strip()]
    pairs = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("1 ") and i + 1 < len(lines) and lines[i + 1].startswith("2 "):
            pairs.append((lines[i], lines[i + 1]))
            i += 2
        else:
            i += 1
    return pairs


def _start_jd(start_time):
    """Julian date of start_time (default: now, UTC)."""
    if start_time is None:
        start_time = datetime.now(timezone.utc)
    return 2440587.5 + (start_time.timestamp() / 86400.0)


def _gmst_rad(jd):
    """Mirror of WorldCoordinateSystem.julianToGmst for an array of Julian dates."""
    d = np.asarray(jd, dtype=np.float64) - 2451545.0
    t = d / 36525.0
    gmst = 280.46061837 + 360.98564736629 * d + 0.000387933 * t * t - t * t * t / 38710000.0
    return np.deg2rad(np.mod(gmst, 360.0))


_EARTH_ROTATION_RATE = 7.2921150e-5  # rad/s, constants.earth.rotationRate
_WGS84_A = 6378.137
_WGS84_E2 = 6.69437999014e-3


def _teme_to_ecef_into(pos, vel, jd, sats):
    """Rotate pos[:, sats] (and vel) from TEME to ECEF in place. Velocity is
    the derivative of the ECEF track (rotation minus omega x r), matching the
    Zig writeOutput path."""
    theta = _gmst_rad(jd)[:, None]
    c, s = np.cos(theta), np.sin(theta)
    p = pos[:, sats]
    x, y = p[..., 0].copy(), p[..., 1].copy()
    p[..., 0] = c * x + s * y
    p[..., 1] = c * y - s * x
    if vel is not None:
        v = vel[:, sats]
        vx, vy = v[..., 0].copy(), v[..., 1].copy()
        v[..., 0] = c * vx + s * vy + _EARTH_ROTATION_RATE * p[..., 1]
        v[..., 1] = c * vy - s * vx - _EARTH_ROTATION_RATE * p[..., 0]


def _ecef_to_geodetic_into(pos, sats):
    """WGS84 ECEF km -> (lat deg, lon deg, alt km) in place; mirrors
    WorldCoordinateSystem.ecefToGeodeticDeg (iterative, 10 rounds)."""
    p = pos[:, sats]
    x, y, z = p[..., 0].copy(), p[..., 1].copy(), p[..., 2].copy()
    lon = np.arctan2(y, x)
    r = np.sqrt(x * x + y * y)
    lat = np.arctan2(z, r * (1.0 - _WGS84_E2))
    for _ in range(10):
        sin_lat = np.sin(lat)
        n = _WGS84_A / np.sqrt(1.0 - _WGS84_E2 * sin_lat * sin_lat)
        lat = np.arctan2(z + _WGS84_E2 * n * sin_lat, r)
    sin_lat = np.sin(lat)
    n = _WGS84_A / np.sqrt(1.0 - _WGS84_E2 * sin_lat * sin_lat)
    p[..., 0] = np.degrees(lat)
    p[..., 1] = np.degrees(lon)
    p[..., 2] = r / np.cos(lat) - n


class Constellation:
    """Satellites parsed once, for repeated ``propagate`` / ``screen`` calls.

    Near-earth satellites run through SIMD SGP4 batches and deep-space ones
    through SDP4; a constellation can mix both.

    Parameters
    ----------
    source : str, optional
        Orbital element data source. Can be one of:

        - CelesTrak group name (e.g., ``"starlink"``, ``"gps"``, ``"stations"``)
        - URL to TLE or OMM JSON data
        - Local file path (TLE or OMM JSON, auto-detected)
        - Raw TLE text string
        - Raw OMM JSON string (single object or array)

    norad_id : int or list of int, optional
        NORAD catalog ID(s) to fetch from CelesTrak. If provided, ``source``
        is ignored.

    Attributes
    ----------
    num_satellites : int
        Number of satellites in the constellation.
    epochs : list of float
        Epoch for each satellite as Julian Date.

    Examples
    --------
    Load from CelesTrak group name:

    >>> c = Constellation("starlink")
    >>> c.num_satellites
    6547

    Load specific satellites by NORAD ID:

    >>> iss = Constellation(norad_id=25544)
    >>> iss.num_satellites
    1

    Load from OMM JSON (supports 6+ digit NORAD IDs):

    >>> import json
    >>> omm = [{"NORAD_CAT_ID": 25544, "EPOCH": "2026-04-15T12:00:00.000000", ...}]
    >>> c = Constellation(json.dumps(omm))

    Propagate multiple times efficiently:

    >>> c = Constellation("gps")
    >>> pos1 = propagate(c, [0, 60, 120])
    >>> pos2 = propagate(c, [180, 240, 300])

    See Also
    --------
    propagate : Propagate satellites to specified times.
    screen : Screen for conjunction events.
    """

    def __init__(self, source=None, *, norad_id=None):
        text, fmt = _load_tle_text(source, norad_id)
        if fmt == "json":
            sats = _Satrec.from_omm(text, _WGS72)
        else:
            sats = [_Satrec.twoline2rv(l1, l2, _WGS72) for l1, l2 in _parse_tle_pairs(text)]

        # Output columns: SGP4 satellites, then SDP4, then any whose elements
        # SGP4/SDP4 rejected (e.g. already decayed), which always come out NaN.
        # The batch is built from the same Satrecs, so it can't classify differently.
        sgp4 = [s for s in sats if not s.is_deep_space and s.error == 0]
        self._sdp4_satrecs = [s for s in sats if s.is_deep_space]
        self._failed_satrecs = [s for s in sats if not s.is_deep_space and s.error != 0]
        self._satrecs = sgp4 + self._sdp4_satrecs + self._failed_satrecs  # by output column
        self._n_sgp4 = len(sgp4)
        self._total_sats = len(sats)
        self._zig = _Sgp4Constellation(sgp4, gravity_model=_WGS72) if sgp4 else None

    @property
    def num_satellites(self):
        """Number of satellites in the constellation."""
        return self._total_sats

    @property
    def epochs(self):
        """TLE epoch for each satellite as Julian Date."""
        epochs = []
        if self._zig is not None:
            epochs.extend(self._zig.epochs)
        for sat in self._sdp4_satrecs + self._failed_satrecs:
            epochs.append(sat.jdsatepoch + sat.jdsatepochF)
        return epochs


def propagate(
    source, times, *, start_time=None, output="ecef", velocities=False, norad_id=None
):
    """Propagate one or more satellites to the given times.

    Parameters
    ----------
    source : Constellation, str, or None
        Orbital element data source. Can be:

        - ``Constellation`` object (most efficient for repeated calls)
        - CelesTrak group name (e.g., ``"starlink"``, ``"gps"``)
        - URL to TLE or OMM JSON data
        - Local file path (TLE or OMM JSON, auto-detected)
        - Raw TLE text string or OMM JSON string
        - ``None`` if using ``norad_id``

    times : array-like
        Times at which to compute positions, in minutes from ``start_time``.
        Will be converted to a contiguous float64 array.
    start_time : datetime, optional
        Reference time for the ``times`` array. Defaults to current UTC time.
    output : {"ecef", "teme", "geodetic"}, default "ecef"
        Output coordinate frame:

        - ``"ecef"``: Earth-Centered Earth-Fixed (x, y, z in km)
        - ``"teme"``: True Equator Mean Equinox (x, y, z in km)
        - ``"geodetic"``: (latitude in deg, longitude in deg, altitude in km)

    velocities : bool, default False
        If True, also return velocity vectors.
    norad_id : int or list of int, optional
        NORAD catalog ID(s) to fetch from CelesTrak. If provided, ``source``
        is ignored.

    Returns
    -------
    positions : ndarray
        Array of shape ``(n_times, n_satellites, 3)`` containing positions.
        Units depend on ``output`` parameter.
    velocities : ndarray, optional
        Only returned if ``velocities=True``. Array of shape
        ``(n_times, n_satellites, 3)`` containing velocities in km/s.

    Examples
    --------
    Basic propagation with CelesTrak group:

    >>> positions = propagate("starlink", [0, 60, 120])
    >>> positions.shape
    (3, 6547, 3)  # 3 times, 6547 satellites, xyz

    Get positions and velocities:

    >>> pos, vel = propagate("gps", [0, 30, 60], velocities=True)
    >>> vel.shape
    (3, 31, 3)

    Output in geodetic coordinates (lat, lon, alt):

    >>> geo = propagate("stations", [0], output="geodetic")
    >>> lat, lon, alt = geo[0, 0]  # First time, first satellite

    Propagate specific satellite by NORAD ID:

    >>> iss_pos = propagate(None, range(1440), norad_id=25544)

    Reuse a Constellation to avoid re-parsing:

    >>> c = Constellation("starlink")
    >>> pos1 = propagate(c, [0, 60])
    >>> pos2 = propagate(c, [120, 180])

    See Also
    --------
    Constellation : Pre-parsed TLE data for efficient reuse.
    screen : Screen for conjunction events.
    """
    if isinstance(source, Constellation):
        const = source
    else:
        const = Constellation(source, norad_id=norad_id)

    times = np.ascontiguousarray(times, dtype=np.float64)
    n_sats, n_times = const.num_satellites, len(times)
    shape = (n_times, n_sats, 3)
    start = _start_jd(start_time)

    pos = np.empty(shape, dtype=np.float64)
    vel = np.empty(shape, dtype=np.float64) if velocities else None
    n_sgp4 = const._n_sgp4

    # SGP4 fills columns [0, n_sgp4); output_stride=n_sats leaves room for SDP4
    if const._zig is not None and n_sgp4 > 0:
        sgp4_offsets = (start - np.array(const._zig.epochs, dtype=np.float64)) * 1440.0
        const._zig.propagate_into(
            times,
            pos,
            vel,
            epoch_offsets=sgp4_offsets,
            satellite_mask=None,
            output=output,
            reference_jd=start,
            time_major=True,
            output_stride=n_sats,
        )

    # SDP4 fills columns [n_sgp4, n_sats). It only emits TEME, so ECEF and
    # geodetic conversion happen here.
    if const._sdp4_satrecs:
        vel_out = vel if velocities else np.empty(shape, dtype=np.float64)
        _sdp4_batch_propagate_into(
            const._sdp4_satrecs,
            np.full(n_times, start),
            times / 1440.0,
            pos,
            vel_out,
            output_stride=n_sats,
            sat_offset=n_sgp4,
        )
        if output != "teme":
            sdp4 = slice(n_sgp4, n_sats)
            _teme_to_ecef_into(pos, vel, start + times / 1440.0, sdp4)
            if output == "geodetic":
                _ecef_to_geodetic_into(pos, sdp4)

    if const._failed_satrecs:
        failed = slice(n_sats - len(const._failed_satrecs), n_sats)
        pos[:, failed] = np.nan
        if velocities:
            vel[:, failed] = np.nan

    return (pos, vel) if velocities else pos


_SCREEN_STATE_BUDGET = 256 * 1024 * 1024  # bytes of positions + velocities per chunk
_MU_WGS72 = 398600.8
_TCA_TOL_SEC = 1e-4
_TCA_MAX_ITER = 60


class Conjunctions(NamedTuple):
    """Close approaches found by ``screen``, one entry per encounter.

    ``i`` and ``j`` are satellite indices (``i < j``) in the constellation's
    output order, ``time`` is the time of closest approach in minutes from
    ``start_time``, ``distance`` the miss distance in km and ``speed`` the
    relative speed in km/s at that moment.
    """

    i: np.ndarray
    j: np.ndarray
    time: np.ndarray
    distance: np.ndarray
    speed: np.ndarray


def screen(
    source,
    duration,
    *,
    threshold=5.0,
    step=1.0,
    target=None,
    start_time=None,
    norad_id=None,
):
    """Find every close approach within ``threshold`` km over a time window.

    Satellites are propagated on a grid ``step`` minutes apart and every pair is
    checked with a bound that accounts for the motion between grid points, so an
    encounter that falls between samples is still found. Each candidate is then
    refined to the exact time of closest approach. A coarser ``step`` is faster
    but sends more candidates to refinement; it never misses encounters.

    Parameters
    ----------
    source : Constellation, str, or None
        Anything ``Constellation`` accepts, or ``None`` with ``norad_id``.
    duration : float
        Length of the window in minutes, starting at ``start_time``.
    threshold : float, default 5.0
        Miss distance in km.
    step : float, default 1.0
        Grid spacing in minutes.
    target : int or sequence of int, optional
        Only report encounters involving these satellite indices.
    start_time : datetime, optional
        Start of the window. Defaults to now (UTC).
    norad_id : int or list of int, optional
        NORAD catalog ID(s) to fetch from CelesTrak instead of ``source``.

    Returns
    -------
    Conjunctions
        Arrays ``i``, ``j``, ``time``, ``distance`` and ``speed``, sorted by time.

    Examples
    --------
    >>> c = Constellation("stations")
    >>> events = screen(c, 1440, threshold=10.0)
    >>> for i, j, t, d in zip(events.i, events.j, events.time, events.distance):
    ...     print(f"{i} and {j} pass {d:.2f} km apart at t+{t:.1f} min")

    Notes
    -----
    The search assumes at most one range minimum per grid interval, which only
    fails for co-moving pairs (docked or formation-flying objects closing at
    well under 10 m/s).
    """
    const = source if isinstance(source, Constellation) else Constellation(source, norad_id=norad_id)
    n = const.num_satellites
    dt = step * 60.0
    n_steps = max(int(math.ceil(duration * 60.0 / dt)), 1)
    primary = None
    if target is not None:
        primary = np.zeros(n, dtype=bool)
        primary[np.atleast_1d(target)] = True

    # Propagate in chunks that share their end step, so no interval is skipped;
    # duplicates from the shared step are removed below.
    chunk = int(np.clip(_SCREEN_STATE_BUDGET // max(n * 48, 1), 2, 480))
    found = []
    start = 0
    while start < n_steps:
        stop = min(start + chunk, n_steps)
        times = np.arange(start, stop + 1) * step
        pos, vel = propagate(const, times, start_time=start_time, output="teme", velocities=True)
        cand = np.frombuffer(_screen_candidates(pos, vel, n, dt, threshold, primary), dtype=np.uint32)
        cand = cand.reshape(-1, 3).astype(np.int64)
        cand[:, 2] += start
        found.append(cand)
        start = stop

    cand = np.unique(np.concatenate(found), axis=0)
    i, j, t_sec, distance, speed, converged = _refine_tca(const._satrecs, _start_jd(start_time), cand, dt)
    # The grid is rounded up to whole steps, so drop anything past the window
    keep = np.flatnonzero(converged & (distance <= threshold) & (t_sec <= duration * 60.0))
    keep = keep[np.argsort(t_sec[keep], kind="stable")]
    return Conjunctions(i[keep], j[keep], t_sec[keep] / 60.0, distance[keep], speed[keep])


def _states_at(satrecs, jd0, sat_idx, t_sec):
    """States of satellite ``sat_idx[k]`` at ``jd0 + t_sec[k]``; NaN on failure."""
    pos = np.empty((len(sat_idx), 3))
    vel = np.empty((len(sat_idx), 3))
    jd_base = math.floor(jd0)
    fr_base = jd0 - jd_base
    order = np.argsort(sat_idx, kind="stable")
    for grp in np.split(order, np.flatnonzero(np.diff(sat_idx[order])) + 1):
        r = np.empty((len(grp), 3))
        v = np.empty((len(grp), 3))
        jd = np.full(len(grp), float(jd_base))
        fr = np.ascontiguousarray(fr_base + t_sec[grp] / 86400.0)
        satrecs[sat_idx[grp[0]]].sgp4_array_into(jd, fr, r, v)
        pos[grp] = r
        vel[grp] = v
    return pos, vel


def _relative(satrecs, jd0, i, j, t):
    ri, vi = _states_at(satrecs, jd0, i, t)
    rj, vj = _states_at(satrecs, jd0, j, t)
    rho, rho_dot = ri - rj, vi - vj
    acc = -_MU_WGS72 * (ri / np.linalg.norm(ri, axis=1, keepdims=True) ** 3
                        - rj / np.linalg.norm(rj, axis=1, keepdims=True) ** 3)
    f = np.einsum("ij,ij->i", rho, rho_dot)
    fp = np.einsum("ij,ij->i", rho_dot, rho_dot) + np.einsum("ij,ij->i", rho, acc)
    return f, fp


def _refine_tca(satrecs, jd0, cand, dt):
    """Root of f = rho . rho_dot in each candidate interval, where f(a) < 0 <= f(b).

    Safeguarded Newton on all candidates at once (f' with a two-body rho_ddot),
    bisecting whenever a step leaves the bracket or stalls, to ``_TCA_TOL_SEC``.
    """
    i, j = cand[:, 0], cand[:, 1]
    if len(cand) == 0:
        empty = np.empty(0)
        return i, j, empty, empty, empty, np.empty(0, dtype=bool)
    a = cand[:, 2] * dt
    b = a + dt
    t = a + 0.5 * dt
    width = b - a
    active = np.ones(len(cand), dtype=bool)
    for _ in range(_TCA_MAX_ITER):
        idx = np.flatnonzero(active)
        if len(idx) == 0:
            break
        f, fp = _relative(satrecs, jd0, i[idx], j[idx], t[idx])
        neg = f < 0
        a[idx] = np.where(neg, t[idx], a[idx])
        b[idx] = np.where(neg, b[idx], t[idx])
        with np.errstate(divide="ignore", invalid="ignore"):
            t_new = t[idx] - f / fp
        # Bisect when Newton leaves the bracket or stalls (flat f on slow pairs), so
        # the bracket at least halves every other iteration.
        stalled = (b[idx] - a[idx]) > 0.5 * width[idx]
        width[idx] = b[idx] - a[idx]
        bad = ~np.isfinite(t_new) | (t_new <= a[idx]) | (t_new >= b[idx]) | stalled
        t_new = np.where(bad, 0.5 * (a[idx] + b[idx]), t_new)
        step = np.abs(t_new - t[idx])
        t[idx] = t_new
        active[idx] = (step > _TCA_TOL_SEC) & (b[idx] - a[idx] > _TCA_TOL_SEC)

    ri, vi = _states_at(satrecs, jd0, i, t)
    rj, vj = _states_at(satrecs, jd0, j, t)
    return i, j, t, np.linalg.norm(ri - rj, axis=1), np.linalg.norm(vi - vj, axis=1), ~active


__all__ = [
    "__version__",
    "Tle",
    "Constellation",
    "propagate",
    "screen",
    "Conjunctions",
    "hohmann_transfer",
    "bi_elliptic_transfer",
    "lambert",
    "orbital_velocity",
    "orbital_period",
    "escape_velocity",
    "propagate_numerical",
    "EARTH_MU",
    "EARTH_R_EQ",
    "EARTH_J2",
    "SUN_MU",
    "MOON_MU",
]

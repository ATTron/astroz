"""astroz - High-performance astrodynamics library."""

import numpy as np
from datetime import datetime, timezone
from pathlib import Path
from ._astroz import (
    Tle,
    Sgp4 as _Sgp4,
    Sgp4Constellation as _Sgp4Constellation,
    WGS84,
    min_distances as _min_distances,
    coarse_screen as _coarse_screen,
)

# Aliases for CelesTrak group names (only where they differ)
_CELESTRAK_ALIASES = {
    "all": "active",
    "iss": "stations",
    "gps": "gps-ops",
    "glonass": "glo-ops",
}


class Sgp4(_Sgp4):
    """SGP4 orbit propagator for a single satellite."""

    def propagate_batch(self, times):
        """Batch propagation returning (positions, velocities) arrays."""
        times = np.asarray(times, dtype=np.float64)
        pos = np.empty((len(times), 3), dtype=np.float64)
        vel = np.empty((len(times), 3), dtype=np.float64)
        self.propagate_into(times, pos, vel)
        return pos, vel


def _fetch_url(url):
    """Fetch URL content."""
    import urllib.request

    req = urllib.request.Request(url, headers={"User-Agent": "astroz"})
    return urllib.request.urlopen(req, timeout=60).read().decode("utf-8")


def _celestrak_url(group=None, norad_id=None):
    """Build CelesTrak URL."""
    if norad_id is not None:
        ids = (
            ",".join(str(i) for i in norad_id)
            if isinstance(norad_id, (list, tuple))
            else str(norad_id)
        )
        return f"https://celestrak.org/NORAD/elements/gp.php?CATNR={ids}&FORMAT=tle"
    group_name = _CELESTRAK_ALIASES.get(group.lower(), group)
    return f"https://celestrak.org/NORAD/elements/gp.php?GROUP={group_name}&FORMAT=tle"


def _parse_tle_metadata(tle_text):
    """Parse TLE text into list of metadata dicts."""
    lines = tle_text.strip().replace("\r", "").split("\n")
    metadata = []
    i = 0
    while i < len(lines) - 2:
        name = lines[i].strip()
        line1 = lines[i + 1].strip()
        line2 = lines[i + 2].strip()
        if line1.startswith("1 ") and line2.startswith("2 "):
            mean_motion = float(line2[52:63])
            period = 1440.0 / mean_motion
            if period <= 225:  # Near-earth only
                metadata.append(
                    {
                        "name": name,
                        "norad_id": int(line1[2:7]),
                        "inclination": float(line2[8:16]),
                        "eccentricity": float("0." + line2[26:33].strip()),
                        "period": period,
                        "mean_motion": mean_motion,
                    }
                )
            i += 3
        else:
            i += 1
    return metadata


class Constellation:
    """SGP4 constellation propagator with automatic buffer reuse.

    Example:
        constellation = Constellation("starlink")
        positions = constellation.propagate(np.arange(1440))
    """

    def __init__(self, source=None, *, norad_id=None, with_metadata=False):
        """Load constellation from TLE source.

        Args:
            source: File path, URL, TLE string, or CelesTrak group name
            norad_id: NORAD catalog ID or list of IDs
            with_metadata: If True, parse satellite metadata
        """
        if norad_id is not None:
            tle_text = _fetch_url(_celestrak_url(norad_id=norad_id))
        elif source is None:
            raise ValueError("Must specify 'source' or 'norad_id'")
        elif source.startswith(("http://", "https://")):
            tle_text = _fetch_url(source)
        elif Path(source).exists():
            tle_text = Path(source).read_text()
        elif "1 " in source and "2 " in source:
            tle_text = source
        else:
            tle_text = _fetch_url(_celestrak_url(group=source))

        self._zig = _Sgp4Constellation.from_tle_text(tle_text)
        self._metadata = _parse_tle_metadata(tle_text) if with_metadata else None
        self._pos_buffer = None
        self._buffer_shape = None

    @property
    def num_satellites(self):
        return self._zig.num_satellites

    @property
    def metadata(self):
        return self._metadata

    @property
    def epochs(self):
        return self._zig.epochs

    @property
    def batch_size(self):
        return self._zig.batch_size

    def propagate(
        self, times, *, start_time=None, output="ecef", velocities=False, out=None
    ):
        """Propagate constellation and return positions.

        Args:
            times: minutes from start_time (e.g., np.arange(1440) for 1 day)
            start_time: datetime (defaults to now UTC)
            output: "ecef" (default), "teme", or "geodetic"
            velocities: if True, also return velocities
            out: optional pre-allocated output array

        Returns:
            positions (num_times, num_sats, 3), or (positions, velocities) tuple
        """
        times = np.ascontiguousarray(times, dtype=np.float64)
        n_sats, n_times = self._zig.num_satellites, len(times)
        shape = (n_times, n_sats, 3)

        if start_time is None:
            start_time = datetime.now(timezone.utc)

        start_jd = 2440587.5 + (start_time.timestamp() / 86400.0)
        epochs = np.array(self._zig.epochs, dtype=np.float64)
        epoch_offsets = (start_jd - epochs) * 1440.0

        # Use provided buffer, cached buffer, or allocate new
        if out is not None:
            pos = out
        elif self._buffer_shape == shape:
            pos = self._pos_buffer
        else:
            self._pos_buffer = np.empty(shape, dtype=np.float64)
            self._buffer_shape = shape
            pos = self._pos_buffer

        vel = np.empty(shape, dtype=np.float64) if velocities else None

        self._zig.propagate_into(
            times,
            pos,
            vel,
            epoch_offsets=epoch_offsets,
            satellite_mask=None,
            output=output,
            reference_jd=start_jd,
            time_major=True,
        )
        return (pos, vel) if velocities else pos


def min_distances(positions, pairs, velocities=None):
    """Find minimum distance between satellite pairs across all time steps."""
    positions = np.ascontiguousarray(positions, dtype=np.float64)
    pairs = np.ascontiguousarray(pairs, dtype=np.uint32)
    num_sats = positions.shape[0]

    if velocities is not None:
        velocities = np.ascontiguousarray(velocities, dtype=np.float64)
        result = _min_distances(positions, pairs, velocities, num_sats)
        return (
            np.array(result[0]),
            np.array(result[1], dtype=np.uint32),
            np.array(result[2]),
        )

    result = _min_distances(positions, pairs, None, num_sats)
    return np.array(result[0]), np.array(result[1], dtype=np.uint32)


def coarse_screen(positions, threshold=10.0):
    """Find satellite pairs within threshold distance at any time step."""
    positions = np.ascontiguousarray(positions, dtype=np.float64)
    return _coarse_screen(positions, positions.shape[0], threshold, None)


# Backward compatibility aliases
def load_constellation(source=None, *, norad_id=None, with_metadata=False):
    """Load constellation. Returns Constellation object (or tuple with metadata)."""
    c = Constellation(source, norad_id=norad_id, with_metadata=with_metadata)
    if with_metadata:
        return c, c.metadata
    return c


__version__ = "0.4.4"

__all__ = [
    "__version__",
    "Tle",
    "Sgp4",
    "WGS84",
    "Constellation",
    "load_constellation",
    "min_distances",
    "coarse_screen",
]

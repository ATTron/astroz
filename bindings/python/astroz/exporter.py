"""Export a Satrec as TLE lines or an OMM record, a drop-in for sgp4.exporter."""

from datetime import datetime, timedelta
from math import pi

_deg2rad = pi / 180.0
_xpdotp = 1440.0 / (2.0 * pi)  # rad/min -> rev/day


def compute_checksum(line):
    """TLE checksum: digits count as themselves, minus signs as 1, everything else 0."""
    return sum(int(c) if c.isdigit() else c == "-" for c in line[:68]) % 10


def export_tle(satrec):
    """The two TLE lines for ``satrec``, formatted as python-sgp4 writes them."""
    line1 = "".join((
        "1 ",
        satrec.satnum_str,
        (satrec.classification.strip() or "U"),
        " ",
        f"{satrec.intldesg:8} ",
        f"{satrec.epochyr % 100:02d}{satrec.epochdays:012.8f} ",
        f"{satrec.ndot * (_xpdotp * 1440.0): 8.8f}".replace("0", "", 1) + " ",
        _abbreviate_rate(satrec.nddot * _xpdotp * 20736000.0, "-0"),
        _abbreviate_rate(satrec.bstar * 10.0, "+0"),
        f"{satrec.ephtype} {satrec.elnum:4}",
    ))
    line1 += str(compute_checksum(line1))

    for name, value, top in (("Inclination", satrec.inclo, pi), ("RAAN", satrec.nodeo, 2 * pi),
                             ("Argument of Perigee", satrec.argpo, 2 * pi), ("Mean Anomaly", satrec.mo, 2 * pi)):
        if not 0 <= value <= top:
            raise ValueError(f"{name} must be between 0 and {'pi' if top == pi else '2 pi'}, got {value!r}")

    line2 = "".join((
        "2 ",
        satrec.satnum_str,
        f" {satrec.inclo / _deg2rad:8.4f} ",
        f"{satrec.nodeo / _deg2rad:8.4f} ",
        f"{satrec.ecco:8.7f}".replace("0.", "") + " ",
        f"{satrec.argpo / _deg2rad:8.4f} ",
        f"{satrec.mo / _deg2rad:8.4f} ",
        f"{satrec.no_kozai * _xpdotp:11.8f}",
        str(satrec.revnum).rjust(5),
    ))
    line2 += str(compute_checksum(line2))
    return line1, line2


def export_omm(satrec, object_name):
    """An OMM record (dict of CCSDS field names) for ``satrec``."""
    launch_year = int(satrec.intldesg[:2])
    launch_year += 1900 + (launch_year < 57) * 100
    epoch = datetime(2000, 1, 1, 12) + timedelta(days=(satrec.jdsatepoch - 2451545.0) + satrec.jdsatepochF)
    return {
        "OBJECT_NAME": object_name,
        "OBJECT_ID": f"{launch_year}-{satrec.intldesg[2:]}",
        "CENTER_NAME": "EARTH",
        "REF_FRAME": "TEME",
        "TIME_SYSTEM": "UTC",
        "MEAN_ELEMENT_THEORY": "SGP4",
        "EPOCH": epoch.strftime("%Y-%m-%dT%H:%M:%S.%f"),
        "MEAN_MOTION": satrec.no_kozai * _xpdotp,
        "ECCENTRICITY": satrec.ecco,
        "INCLINATION": satrec.inclo / _deg2rad,
        "RA_OF_ASC_NODE": satrec.nodeo / _deg2rad,
        "ARG_OF_PERICENTER": satrec.argpo / _deg2rad,
        "MEAN_ANOMALY": satrec.mo / _deg2rad,
        "EPHEMERIS_TYPE": satrec.ephtype,
        "CLASSIFICATION_TYPE": satrec.classification,
        "NORAD_CAT_ID": satrec.satnum,
        "ELEMENT_SET_NO": satrec.elnum,
        "REV_AT_EPOCH": satrec.revnum,
        "BSTAR": satrec.bstar,
        "MEAN_MOTION_DOT": satrec.ndot * (_xpdotp * 1440.0),
        "MEAN_MOTION_DDOT": satrec.nddot * (_xpdotp * 1440.0 * 1440),
    }


def _abbreviate_rate(value, zero_exponent):
    """TLE implied-decimal exponent field, e.g. 0.00027310 -> ' 27310-3 '."""
    return (f"{value: 4.4e} ".replace(".", "").replace("e+00", zero_exponent)
            .replace("e-0", "-").replace("e+0", "+"))

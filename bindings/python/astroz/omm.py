"""Read CCSDS OMM records from CSV or XML into a Satrec, a drop-in for sgp4.omm.

For OMM JSON, ``astroz.Constellation`` and ``_astroz.Satrec.from_omm`` read the
text directly.
"""

import csv
import xml.etree.ElementTree as ET
from datetime import datetime
from math import pi

from .api import WGS72

_epoch0 = datetime(1949, 12, 31)
_to_radians = pi / 180.0
_ndot_units = 1036800.0 / pi  # rev/day^2 -> rad/min^2
_nddot_units = 2985984000.0 / 2.0 / pi  # rev/day^3 -> rad/min^3


def parse_csv(file):
    """One dict of OMM fields per CSV row."""
    return csv.DictReader(file)


def parse_xml(file):
    """One dict of OMM fields per <segment> of an OMM XML file."""
    root = ET.parse(file).getroot()
    for segment in root.findall(".//segment"):
        data = segment.find("data")
        fields = {}
        for element in segment.find("metadata"), data.find("meanElements"), data.find("tleParameters"):
            fields.update((field.tag, field.text) for field in element)
        yield fields


def initialize(sat, fields, gravconst=WGS72):
    """Set up ``sat`` (a ``Satrec()``) from one OMM record, as python-sgp4 does."""
    sat.classification = fields["CLASSIFICATION_TYPE"]
    sat.intldesg = (fields["OBJECT_ID"] or "")[2:].replace("-", "")
    sat.ephtype = int(fields["EPHEMERIS_TYPE"])
    sat.elnum = int(fields["ELEMENT_SET_NO"])
    sat.revnum = int(fields["REV_AT_EPOCH"])

    epoch_text = fields["EPOCH"]
    try:
        epoch_datetime = datetime.strptime(epoch_text, "%Y-%m-%dT%H:%M:%S.%f")
    except ValueError:  # no fractional seconds, or a trailing Z
        epoch_datetime = datetime.fromisoformat(epoch_text.rstrip("Z"))
    epoch = (epoch_datetime - _epoch0).total_seconds() / 86400.0

    sat.sgp4init(
        gravconst, "i", int(fields["NORAD_CAT_ID"]), epoch,
        float(fields["BSTAR"]),
        float(fields["MEAN_MOTION_DOT"]) / _ndot_units,
        float(fields["MEAN_MOTION_DDOT"]) / _nddot_units,
        float(fields["ECCENTRICITY"]),
        float(fields["ARG_OF_PERICENTER"]) * _to_radians,
        float(fields["INCLINATION"]) * _to_radians,
        float(fields["MEAN_ANOMALY"]) * _to_radians,
        float(fields["MEAN_MOTION"]) / 720.0 * pi,
        float(fields["RA_OF_ASC_NODE"]) * _to_radians,
    )

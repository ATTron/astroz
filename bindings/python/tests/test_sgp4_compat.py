"""astroz.api, astroz.exporter and astroz.omm against python-sgp4 on the same inputs."""

import csv
import io

import numpy as np
import pytest

sgp4_api = pytest.importorskip("sgp4.api")
from sgp4 import exporter as ref_exporter  # noqa: E402
from sgp4 import omm as ref_omm  # noqa: E402
from sgp4.io import compute_checksum  # noqa: E402

from astroz import exporter, omm  # noqa: E402
from astroz.api import WGS72, WGS72OLD, WGS84, Satrec, SatrecArray  # noqa: E402


def _checksummed(line):
    return line + str(compute_checksum(line))


TLES = {
    "ISS": ("1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995",
            "2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123"),
    "GPS (SDP4)": ("1 20413U 90005A   24186.00000000  .00000012  00000+0  10000-3 0  9992",
                   "2 20413  55.4408  61.4858 0112981 129.5765 231.5553  2.00561730104446"),
    "Molniya (SDP4, resonant)": ("1 08195U 75081A   06176.33215444  .00000099  00000-0  11873-3 0   813",
                                 "2 08195  64.1586 279.0717 6877146 264.7651  20.2257  2.00491383225656"),
    "Alpha-5 100404": ("1 A0404U 26192A   26263.38590896  .00031477  00000+0  57128-3 0  9992",
                       "2 A0404  53.1597 307.1840 0001009  61.6382 298.4729 15.49331404  5885"),
    # Heavy drag: decays within days, so it exercises the error codes
    "decaying": (_checksummed("1 99001U 26001A   26272.50000000  .01000000  00000+0  50000-1 0  999"),
                 _checksummed("2 99001  51.6000 100.0000 0010000  90.0000 270.0000 16.10000000    10")),
}
ATTRS = ("satnum satnum_str classification intldesg epochyr epochdays jdsatepoch jdsatepochF ndot nddot "
         "bstar inclo nodeo ecco argpo mo no_kozai elnum revnum ephtype method error a alta altp").split()
DAYS = np.array([0, 0.5, 1, 3, 7.0])


def _times(sat):
    return np.full(len(DAYS), sat.jdsatepoch), sat.jdsatepochF + DAYS


def _same_states(e, r, ref_e, ref_r, tol=1e-6):
    assert (e == ref_e).all()
    assert (np.isnan(r) == np.isnan(ref_r)).all()
    ok = ~np.isnan(ref_r)
    assert np.abs(r[ok] - ref_r[ok]).max(initial=0) < tol  # km


def test_gravity_constants_match():
    assert (WGS72OLD, WGS72, WGS84) == (sgp4_api.WGS72OLD, sgp4_api.WGS72, sgp4_api.WGS84)
    with pytest.raises(ValueError):
        Satrec.twoline2rv(*TLES["ISS"], 3)


@pytest.mark.parametrize("name", TLES)
@pytest.mark.parametrize("gravity", [WGS72OLD, WGS72, WGS84])
def test_twoline2rv_matches(name, gravity):
    sat, ref = Satrec.twoline2rv(*TLES[name], gravity), sgp4_api.Satrec.twoline2rv(*TLES[name], gravity)
    for attr in ATTRS:
        got, want = getattr(sat, attr), getattr(ref, attr)
        assert got == pytest.approx(want, rel=1e-12, abs=1e-15) if isinstance(want, float) else got == want, attr
    _same_states(*sat.sgp4_array(*_times(ref))[:2], *ref.sgp4_array(*_times(ref))[:2])
    assert exporter.export_tle(sat) == ref_exporter.export_tle(ref)


def test_scalar_sgp4_reports_decay_like_python_sgp4():
    sat, ref = Satrec.twoline2rv(*TLES["decaying"]), sgp4_api.Satrec.twoline2rv(*TLES["decaying"])
    codes = []
    for d in DAYS:
        e, r, _ = sat.sgp4(ref.jdsatepoch, ref.jdsatepochF + d)
        ref_e, ref_r, _ = ref.sgp4(ref.jdsatepoch, ref.jdsatepochF + d)
        _same_states(np.array([e]), np.array([r]), np.array([ref_e]), np.array([ref_r]))
        codes.append(e)
    assert 6 in codes or 1 in codes  # the case really does fail


def test_satrec_array_matches():
    sats = [Satrec.twoline2rv(*t) for t in TLES.values()]
    refs = [sgp4_api.Satrec.twoline2rv(*t) for t in TLES.values()]
    jd, fr = _times(refs[0])
    e, r, _ = SatrecArray(sats).sgp4(jd, fr)
    ref_e, ref_r, _ = sgp4_api.SatrecArray(refs).sgp4(jd, fr)
    # The SDP4 batch path's SIMD trig puts deep-space orbits metres, not millimetres, from the scalar result
    _same_states(e, r, ref_e, ref_r, tol=0.01)


# astroz also takes the catalog number as a TLE field string; python-sgp4 only an int
@pytest.mark.parametrize("satnum,ref_satnum", [(25544, 25544), (100000, 100000), ("A0404", 100404)])
def test_sgp4init_matches(satnum, ref_satnum):
    ref = sgp4_api.Satrec.twoline2rv(*TLES["ISS"])
    elements = (ref.jdsatepoch - 2433281.5 + ref.jdsatepochF, ref.bstar, ref.ndot, ref.nddot, ref.ecco,
                ref.argpo, ref.inclo, ref.mo, ref.no_kozai, ref.nodeo)
    sat, want = Satrec(), sgp4_api.Satrec()
    sat.sgp4init(WGS72, "i", satnum, *elements)
    want.sgp4init(WGS72, "i", ref_satnum, *elements)
    for attr in ("satnum", "satnum_str", "epochyr", "epochdays", "jdsatepoch", "jdsatepochF", "method", "error"):
        assert getattr(sat, attr) == pytest.approx(getattr(want, attr)), attr
    _same_states(*sat.sgp4_array(*_times(want))[:2], *want.sgp4_array(*_times(want))[:2])
    assert exporter.export_tle(sat) == ref_exporter.export_tle(want)


def test_metadata_set_before_sgp4init_is_kept():
    sat = Satrec()
    sat.classification, sat.intldesg, sat.elnum, sat.revnum, sat.ephtype = "C", "26001A", 42, 7, 0
    ref = sgp4_api.Satrec.twoline2rv(*TLES["ISS"])
    sat.sgp4init(WGS72, "i", 99999, 27000.5, ref.bstar, ref.ndot, ref.nddot, ref.ecco, ref.argpo,
                 ref.inclo, ref.mo, ref.no_kozai, ref.nodeo)
    assert (sat.classification, sat.intldesg, sat.elnum, sat.revnum) == ("C", "26001A", 42, 7)
    with pytest.raises(ValueError):
        sat.sgp4init(WGS72, "a", 99999, 27000.5, 0, 0, 0, 0.001, 0, 0.9, 0, 0.06, 0)


def test_omm_csv_matches():
    refs = [sgp4_api.Satrec.twoline2rv(*t) for t in TLES.values()]
    records = [ref_exporter.export_omm(r, name) for r, name in zip(refs, TLES)]
    text = io.StringIO()
    writer = csv.DictWriter(text, fieldnames=list(records[0]))
    writer.writeheader()
    writer.writerows(records)

    for got_fields, want_fields in zip(omm.parse_csv(io.StringIO(text.getvalue())),
                                       ref_omm.parse_csv(io.StringIO(text.getvalue()))):
        sat, want = Satrec(), sgp4_api.Satrec()
        omm.initialize(sat, got_fields)
        ref_omm.initialize(want, want_fields)
        assert (sat.satnum, sat.intldesg, sat.classification) == (want.satnum, want.intldesg, want.classification)
        _same_states(*sat.sgp4_array(*_times(want))[:2], *want.sgp4_array(*_times(want))[:2])

#!/usr/bin/env python3
"""Render a day of the active satellite catalog as an interactive Cesium page."""

import json
import numpy as np
import time as time_module
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone

COLORS = {
    "ISS": [255, 255, 255],
    "Starlink": [100, 180, 255],
    "OneWeb": [255, 160, 100],
    "Planet": [130, 200, 130],
    "Spire": [240, 220, 140],
    "Iridium": [200, 140, 220],
    "SSO": [255, 130, 130],
    "LEO53": [140, 220, 200],
    "Other": [180, 180, 190],
}


def classify_satellite(name, inc):
    name_upper = name.upper()
    if "ISS" in name_upper and "ZARYA" in name_upper:
        return "ISS"
    if "STARLINK" in name_upper:
        return "Starlink"
    if "ONEWEB" in name_upper:
        return "OneWeb"
    if "PLANET" in name_upper or "FLOCK" in name_upper:
        return "Planet"
    if "SPIRE" in name_upper or "LEMUR" in name_upper:
        return "Spire"
    if "IRIDIUM" in name_upper:
        return "Iridium"
    if 96 <= inc <= 99:
        return "SSO"
    if 50 <= inc <= 55:
        return "LEO53"
    return "Other"


def near_earth_tles(data):
    """Near-earth satellites from 3-line TLE text: their TLE text and display metadata by NORAD id."""
    lines = data.strip().split("\n")
    tle_lines, meta = [], {}
    for name, line1, line2 in zip(lines, lines[1:], lines[2:]):
        if not (line1.startswith("1 ") and line2.startswith("2 ")):
            continue
        inc, mean_motion = float(line2[8:16]), float(line2[52:63])
        period = 1440.0 / mean_motion
        if period > 225:  # near-earth only
            continue
        tle_lines += [line1, line2]
        satnum = int(line1[2:7])  # CelesTrak TLEs never use Alpha-5
        meta[satnum] = {
            "name": name.strip(),
            "constellation": classify_satellite(name, inc),
            "data": {
                "norad_id": satnum,
                "inclination": inc,
                "eccentricity": float("0." + line2[26:33].strip()),
                "period": period,
                "raan": float(line2[17:25].strip()),
                "mean_motion": mean_motion,
            },
        }
    return "\n".join(tle_lines), meta


def generate_legend_html(counts):
    items = [
        f'<label class="legend-item toggle iss-item" data-constellation="ISS"><input type="checkbox" checked><div class="legend-dot iss-dot">🛰</div>ISS<span class="legend-count">{counts.get("ISS", 0)}</span></label>'
    ]
    for name, rgb in COLORS.items():
        if name == "ISS":
            continue
        label = "LEO 53°" if name == "LEO53" else name
        color = ",".join(map(str, rgb))
        items.append(
            f'<label class="legend-item toggle" data-constellation="{name}"><input type="checkbox" checked><div class="legend-dot" style="background:rgb({color})"></div>{label}<span class="legend-count">{counts.get(name, 0):,}</span></label>'
        )
    return "\n            ".join(items)


def load_active_catalog():
    import urllib.request

    url = "https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=tle"
    req = urllib.request.Request(url, headers={"User-Agent": "astroz"})
    return urllib.request.urlopen(req, timeout=60).read().decode("utf-8")


def main():
    import astroz

    tle_text, meta = near_earth_tles(load_active_catalog())
    const = astroz.Constellation(tle_text)
    start_time = datetime.now(timezone.utc).replace(second=0, microsecond=0)
    num_times = 1440  # 1 day at 1-minute resolution
    times = np.arange(num_times, dtype=np.float64)

    astroz.propagate(const, times[:10], start_time=start_time)  # warm up threads
    t0 = time_module.perf_counter()
    positions = astroz.propagate(const, times, start_time=start_time, output="ecef")
    prop_time = time_module.perf_counter() - t0

    # Cesium points are Earth-fixed, in meters. Columns follow the
    # constellation's order, so look the metadata up by catalog number, and
    # drop any satellite SGP4 rejected (its column is NaN).
    keep = np.isfinite(positions).all(axis=(0, 2))
    satellites = [meta[n] for n, ok in zip(const.satnums, keep) if ok]
    positions = positions[:, keep] * 1000.0
    num_sats = len(satellites)

    total_props = num_sats * num_times
    prop_throughput = total_props / prop_time / 1e6
    print(f"{num_sats:,} satellites, start {start_time.isoformat()}")
    print(f"Propagated {total_props:,} positions to ECEF in {prop_time:.2f}s ({prop_throughput:.1f}M/s)")

    counts = Counter(s["constellation"] for s in satellites)

    import gzip
    import io

    positions_flat = positions.reshape(-1).astype(np.float32)

    gz_buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=gz_buffer, mode="wb", compresslevel=9) as gz:
        gz.write(positions_flat.tobytes())
    gz_data = gz_buffer.getvalue()

    chunk_size = 50 * 1024 * 1024
    chunk_names = []
    for i, start in enumerate(range(0, len(gz_data), chunk_size)):
        suffix = f"a{'abcdefghijklmnopqrstuvwxyz'[i]}"
        chunk_names.append(suffix)
        chunk_path = Path(f"cesium_positions.bin.gz.part_{suffix}")
        chunk_path.write_bytes(gz_data[start : start + chunk_size])
    print(f"Positions: {len(gz_data) / 1e6:.1f} MB gzipped in {len(chunk_names)} chunk(s)")

    template_data = {
        "{{CHUNK_NAMES}}": json.dumps(chunk_names),
        "{{NUM_SATS}}": f"{num_sats:,}",
        "{{NUM_SATS_RAW}}": str(num_sats),
        "{{NUM_FRAMES}}": str(num_times),
        "{{NUM_FRAMES_MINUS_1}}": str(num_times - 1),
        "{{START_FRAME}}": "0",
        "{{THROUGHPUT}}": f"{prop_throughput:.1f}",
        "{{TOTAL_PROPS}}": f"{total_props:,}",
        "{{PROP_TIME_MS}}": f"{prop_time * 1000:.0f}",
        "{{REFERENCE_EPOCH}}": start_time.isoformat(),
        "{{SAT_COLORS}}": json.dumps(
            [COLORS.get(s["constellation"], [150, 150, 150]) for s in satellites],
            separators=(",", ":"),
        ),
        "{{SAT_CONSTELLATIONS}}": json.dumps(
            [s["constellation"] for s in satellites], separators=(",", ":")
        ),
        "{{SAT_NAMES}}": json.dumps(
            [s["name"] for s in satellites], separators=(",", ":")
        ),
        "{{SAT_DATA}}": json.dumps(
            [s["data"] for s in satellites], separators=(",", ":")
        ),
        "{{LEGEND_ITEMS}}": generate_legend_html(counts),
        "{{TLE_FETCH_DATE}}": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
    }

    template = (Path(__file__).parent / "cesium_template.html").read_text()
    html = template
    for key, value in template_data.items():
        html = html.replace(key, value)

    Path("cesium_fast.html").write_text(html)
    print(f"Wrote cesium_fast.html ({len(html) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()

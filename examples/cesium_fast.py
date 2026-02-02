#!/usr/bin/env python3
"""
Generate an interactive Cesium visualization of the satellite catalog.
Uses astroz for high-performance SGP4 propagation.
"""

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


def classify_satellite(name, inc, period):
    """Classify satellite into constellation category."""
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


def parse_tles_and_metadata(data):
    """Parse TLE data into (line1, line2) pairs and satellite metadata."""
    lines = data.strip().split("\n")
    tle_pairs = []
    satellites = []
    i = 0
    while i < len(lines) - 2:
        name, line1, line2 = (
            lines[i].strip(),
            lines[i + 1].strip(),
            lines[i + 2].strip(),
        )
        if line1.startswith("1 ") and line2.startswith("2 "):
            inc = float(line2[8:16])
            mean_motion = float(line2[52:63])
            period = 1440.0 / mean_motion
            if period <= 225:  # Near-earth only (SGP4 limitation)
                tle_pairs.append((line1, line2))
                satellites.append(
                    {
                        "name": name,
                        "constellation": classify_satellite(name, inc, period),
                        "data": {
                            "norad_id": line1[2:7].strip(),
                            "inclination": inc,
                            "eccentricity": float("0." + line2[26:33].strip()),
                            "period": period,
                            "raan": float(line2[17:25].strip()),
                            "mean_motion": mean_motion,
                        },
                    }
                )
            i += 3
        else:
            i += 1
    return tle_pairs, satellites


def generate_legend_html(counts):
    """Generate legend HTML items."""
    items = [
        f'<label class="legend-item toggle iss-item" data-constellation="ISS"><input type="checkbox" checked><div class="legend-dot iss-dot">🛰</div>ISS<span class="legend-count">{counts.get("ISS", 0)}</span></label>'
    ]
    for name, color in [
        ("Starlink", "100,180,255"),
        ("OneWeb", "255,160,100"),
        ("Planet", "130,200,130"),
        ("Spire", "240,220,140"),
        ("Iridium", "200,140,220"),
        ("SSO", "255,130,130"),
        ("LEO53", "140,220,200"),
        ("Other", "180,180,190"),
    ]:
        label = "LEO 53°" if name == "LEO53" else name
        items.append(
            f'<label class="legend-item toggle" data-constellation="{name}"><input type="checkbox" checked><div class="legend-dot" style="background:rgb({color})"></div>{label}<span class="legend-count">{counts.get(name, 0):,}</span></label>'
        )
    return "\n            ".join(items)


def load_tle_text(source):
    """Load TLE text from CelesTrak."""
    import urllib.request

    aliases = {
        "all": "active",
        "iss": "stations",
        "gps": "gps-ops",
        "glonass": "glo-ops",
    }
    group = aliases.get(source.lower(), source)
    url = f"https://celestrak.org/NORAD/elements/gp.php?GROUP={group}&FORMAT=tle"
    req = urllib.request.Request(url, headers={"User-Agent": "astroz"})
    return urllib.request.urlopen(req, timeout=60).read().decode("utf-8")


def main():
    from astroz.api import Satrec, SatrecArray, jday, WGS72

    print("=" * 60)
    print("  Cesium Satellite Visualization (python-sgp4 compatible API)")
    print("=" * 60)

    # Load all active satellites with metadata
    print("\n[1/4] Loading catalog...", end=" ", flush=True)
    tle_text = load_tle_text("all")
    tle_pairs, satellites = parse_tles_and_metadata(tle_text)

    # Create Satrec objects (python-sgp4 compatible)
    satrecs = []
    valid_satellites = []
    for i, (line1, line2) in enumerate(tle_pairs):
        try:
            sat = Satrec.twoline2rv(line1, line2, WGS72)
            if sat.error == 0:
                satrecs.append(sat)
                valid_satellites.append(satellites[i])
        except Exception:
            pass

    satellites = valid_satellites
    sat_array = SatrecArray(satrecs)
    num_sats = len(satrecs)
    print(f"done - {num_sats:,} satellites")

    # Choose a start time: current UTC time rounded down to the minute
    now = datetime.now(timezone.utc)
    start_time = now.replace(second=0, microsecond=0)
    print(f"       Start time: {start_time.isoformat()}")

    # --- Propagation ---
    print("[2/4] Propagating...", end=" ", flush=True)

    num_times = 1440  # 1 day at 1-minute resolution
    t0_total = time_module.perf_counter()

    # Convert start_time to Julian Date arrays (python-sgp4 compatible)
    base_jd, base_fr = jday(
        start_time.year,
        start_time.month,
        start_time.day,
        start_time.hour,
        start_time.minute,
        float(start_time.second),
    )
    # Create time arrays: jd stays constant, fr increments by minutes
    jd = np.full(num_times, base_jd, dtype=np.float64)
    fr = base_fr + np.arange(num_times, dtype=np.float64) / 1440.0  # minutes to days

    # Warmup
    for _ in range(3):
        sat_array.sgp4(jd[:10], fr[:10], velocities=False)

    t0_prop = time_module.perf_counter()
    errors, positions, _ = sat_array.sgp4(jd, fr, velocities=False)
    prop_time = time_module.perf_counter() - t0_prop

    # positions shape from sgp4: (num_sats, num_times, 3) -> transpose to (num_times, num_sats, 3)
    # and convert km → meters
    t0_conv = time_module.perf_counter()
    positions = np.transpose(positions, (1, 0, 2)) * 1000.0
    conv_time = time_module.perf_counter() - t0_conv

    total_time = time_module.perf_counter() - t0_total
    total_props = num_sats * num_times
    prop_throughput = total_props / prop_time / 1e6
    effective_throughput = total_props / total_time / 1e6

    print("done")
    print(
        f"       SGP4 + ECEF:      {prop_time:.2f}s ({prop_throughput:.1f}M props/sec)"
    )
    print(f"       km → m:           {conv_time:.3f}s")
    print(
        f"       Total:            {total_time:.2f}s ({effective_throughput:.1f}M props/sec effective)"
    )

    # Prepare template data
    counts = Counter(s["constellation"] for s in satellites)

    # Save positions as gzipped binary
    import gzip
    import io

    print("[3/4] Binary data...", end=" ", flush=True)
    positions_flat = positions.reshape(-1).astype(np.float32)

    gz_buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=gz_buffer, mode="wb", compresslevel=9) as gz:
        gz.write(positions_flat.tobytes())
    gz_data = gz_buffer.getvalue()
    gz_size_mb = len(gz_data) / 1e6
    print(f"{gz_size_mb:.1f} MB gzipped...", end=" ", flush=True)

    # Split into 50MB chunks
    chunk_size = 50 * 1024 * 1024
    chunk_names = []
    for i, start in enumerate(range(0, len(gz_data), chunk_size)):
        suffix = f"a{'abcdefghijklmnopqrstuvwxyz'[i]}"
        chunk_names.append(suffix)
        chunk_path = Path(f"cesium_positions.bin.gz.part_{suffix}")
        chunk_path.write_bytes(gz_data[start : start + chunk_size])
    print(f"{len(chunk_names)} chunks")

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
        "{{CONV_TIME_MS}}": f"{conv_time * 1000:.0f}",
        "{{TOTAL_TIME_MS}}": f"{total_time * 1000:.0f}",
        "{{EFFECTIVE_THROUGHPUT}}": f"{effective_throughput:.1f}",
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

    # Generate HTML
    print("[4/4] Generating HTML...", end=" ", flush=True)
    template = (Path(__file__).parent / "cesium_template.html").read_text()
    html = template
    for key, value in template_data.items():
        html = html.replace(key, value)

    Path("cesium_fast.html").write_text(html)
    print("done")
    print(f"\nCreated: cesium_fast.html ({len(html) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Generate an interactive Cesium visualization of the satellite catalog.
Uses astroz for high-performance SGP4 propagation.
"""

import json
import numpy as np
import urllib.request
import time as time_module
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone, timedelta

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


def download_tles():
    """Download TLEs from CelesTrak with caching."""
    cache_file = Path(__file__).parent / "tle_cache.txt"
    if (
        cache_file.exists()
        and (time_module.time() - cache_file.stat().st_mtime) < 86400
    ):
        print("(cached)", end=" ")
        return cache_file.read_text()

    print("(downloading)", end=" ")
    req = urllib.request.Request(
        "https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=tle",
        headers={"User-Agent": "Mozilla/5.0 (astroz satellite-viz)"},
    )
    data = urllib.request.urlopen(req, timeout=60).read().decode("utf-8")
    cache_file.write_text(data)
    return data


def parse_tle_epoch(line1):
    """Extract epoch from TLE line 1."""
    epoch_str = line1[18:32].strip()
    year = (
        2000 + int(epoch_str[:2])
        if int(epoch_str[:2]) < 57
        else 1957 + int(epoch_str[:2]) - 57
    )
    return datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(
        days=float(epoch_str[2:]) - 1
    )


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


def parse_tles(data):
    """Parse TLE data into satellite list."""
    lines = data.strip().split("\n")
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
                satellites.append(
                    {
                        "name": name,
                        "line1": line1,
                        "line2": line2,
                        "constellation": classify_satellite(name, inc, period),
                        "epoch": parse_tle_epoch(line1),
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
    return satellites


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


def teme_to_ecef_at_time(positions_teme, utc_time):
    """Convert TEME coordinates to ECEF for a specific UTC time.

    Args:
        positions_teme: Array of shape (num_sats, 3) in meters
        utc_time: datetime object for the conversion

    Returns:
        positions_ecef: Array of same shape in ECEF meters
    """
    # Julian date
    jd = 2440587.5 + (utc_time.timestamp() / 86400.0)

    # GMST in radians (simplified formula)
    t = (jd - 2451545.0) / 36525.0  # Julian centuries from J2000
    gmst_deg = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * t**2
    gmst = np.radians(gmst_deg % 360)

    # Rotation matrix (TEME to ECEF is rotation around Z by -GMST)
    cos_g, sin_g = np.cos(gmst), np.sin(gmst)

    x_teme = positions_teme[:, 0]
    y_teme = positions_teme[:, 1]
    z_teme = positions_teme[:, 2]

    # Rotate around Z axis
    positions_ecef = np.empty_like(positions_teme)
    positions_ecef[:, 0] = x_teme * cos_g + y_teme * sin_g
    positions_ecef[:, 1] = -x_teme * sin_g + y_teme * cos_g
    positions_ecef[:, 2] = z_teme

    return positions_ecef


def main():
    from astroz import Tle, Sgp4

    print("=" * 60)
    print("  Cesium Satellite Visualization")
    print("=" * 60)

    # Download and parse
    print("\n[1/4] Loading catalog...", end=" ", flush=True)
    satellites = parse_tles(download_tles())
    print(f"done - {len(satellites):,} satellites")

    # Choose a start time: use the current UTC time rounded down to the minute
    now = datetime.now(timezone.utc)
    start_time = now.replace(second=0, microsecond=0)
    print(f"       Start time: {start_time.isoformat()}")

    # Propagate each satellite for 1 day at 1-minute intervals
    # Each satellite needs its own tsince calculated from its TLE epoch
    print("[2/4] Propagating...", end=" ", flush=True)

    num_times = 1440  # 1 day
    num_sats = len(satellites)

    # Pre-compute UTC times and their GMST values for TEME→ECEF conversion
    utc_times = [start_time + timedelta(minutes=t) for t in range(num_times)]
    gmst_values = np.empty(num_times, dtype=np.float64)
    for t_idx, utc_time in enumerate(utc_times):
        jd = 2440587.5 + (utc_time.timestamp() / 86400.0)
        t = (jd - 2451545.0) / 36525.0
        gmst_deg = (
            280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * t**2
        )
        gmst_values[t_idx] = np.radians(gmst_deg % 360)

    # Propagate all times for each satellite, then convert to ECEF
    # This minimizes Python↔Zig crossings: one propagate_into call per satellite
    positions = np.empty((num_times, num_sats, 3), dtype=np.float64)
    positions_teme = np.empty((num_times, num_sats, 3), dtype=np.float64)
    pos_buffer = np.empty((num_times, 3), dtype=np.float64)
    vel_buffer = np.empty((num_times, 3), dtype=np.float64)

    # Phase 1: Pure SGP4 propagation
    t0_prop = time_module.perf_counter()
    for sat_idx, s in enumerate(satellites):
        tle = Tle(f"{s['line1']}\n{s['line2']}")
        sgp4 = Sgp4(tle)
        tle_epoch = s["epoch"]

        # Compute tsince for all time steps relative to this satellite's TLE epoch
        epoch_offset = (start_time - tle_epoch).total_seconds() / 60.0
        times = np.arange(num_times, dtype=np.float64) + epoch_offset

        # Single call to propagate all times (fast, happens in Zig)
        try:
            sgp4.propagate_into(times, pos_buffer, vel_buffer)
            positions_teme[:, sat_idx, :] = pos_buffer
        except:
            positions_teme[:, sat_idx, :] = 0

    prop_time = time_module.perf_counter() - t0_prop
    total_props = num_sats * num_times
    prop_throughput = total_props / prop_time / 1e6

    # Phase 2: TEME→ECEF coordinate conversion
    t0_conv = time_module.perf_counter()
    cos_g = np.cos(gmst_values)[:, np.newaxis]
    sin_g = np.sin(gmst_values)[:, np.newaxis]
    x_teme = positions_teme[:, :, 0]
    y_teme = positions_teme[:, :, 1]
    z_teme = positions_teme[:, :, 2]

    positions[:, :, 0] = (x_teme * cos_g + y_teme * sin_g) * 1000  # km to m
    positions[:, :, 1] = (-x_teme * sin_g + y_teme * cos_g) * 1000
    positions[:, :, 2] = z_teme * 1000
    conv_time = time_module.perf_counter() - t0_conv

    total_time = prop_time + conv_time
    effective_throughput = total_props / total_time / 1e6
    print(f"done")
    print(f"       SGP4 propagation: {prop_time:.2f}s ({prop_throughput:.1f}M props/sec)")
    print(f"       TEME→ECEF:        {conv_time:.2f}s")
    print(f"       Total:            {total_time:.2f}s ({effective_throughput:.1f}M props/sec effective)")

    # Prepare template data
    counts = Counter(s["constellation"] for s in satellites)

    # Save positions as gzipped binary, split into chunks for GitHub (<100MB limit)
    import gzip
    import io

    print("[3/4] Binary data...", end=" ", flush=True)
    positions_flat = positions.reshape(-1).astype(np.float32)

    # Gzip in memory
    gz_buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=gz_buffer, mode="wb", compresslevel=9) as gz:
        gz.write(positions_flat.tobytes())
    gz_data = gz_buffer.getvalue()
    gz_size_mb = len(gz_data) / 1e6
    print(f"{gz_size_mb:.1f} MB gzipped...", end=" ", flush=True)

    # Split into 50MB chunks (matching `split -b 50M` naming: aa, ab, ac, ...)
    chunk_size = 50 * 1024 * 1024  # 50MB
    chunk_names = []
    for i, start in enumerate(range(0, len(gz_data), chunk_size)):
        suffix = f"a{'abcdefghijklmnopqrstuvwxyz'[i]}"  # aa, ab, ac, ...
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

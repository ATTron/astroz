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
    'ISS': [255, 255, 255], 'Starlink': [100, 180, 255], 'OneWeb': [255, 160, 100],
    'Planet': [130, 200, 130], 'Spire': [240, 220, 140], 'Iridium': [200, 140, 220],
    'SSO': [255, 130, 130], 'LEO53': [140, 220, 200], 'Other': [180, 180, 190],
}


def download_tles():
    """Download TLEs from CelesTrak with caching."""
    cache_file = Path(__file__).parent / 'tle_cache.txt'
    if cache_file.exists() and (time_module.time() - cache_file.stat().st_mtime) < 86400:
        print("(cached)", end=" ")
        return cache_file.read_text()

    print("(downloading)", end=" ")
    req = urllib.request.Request(
        'https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=tle',
        headers={'User-Agent': 'Mozilla/5.0 (astroz satellite-viz)'}
    )
    data = urllib.request.urlopen(req, timeout=60).read().decode('utf-8')
    cache_file.write_text(data)
    return data


def parse_tle_epoch(line1):
    """Extract epoch from TLE line 1."""
    epoch_str = line1[18:32].strip()
    year = 2000 + int(epoch_str[:2]) if int(epoch_str[:2]) < 57 else 1957 + int(epoch_str[:2]) - 57
    return datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(days=float(epoch_str[2:]) - 1)


def classify_satellite(name, inc, period):
    """Classify satellite into constellation category."""
    name_upper = name.upper()
    if 'ISS' in name_upper and 'ZARYA' in name_upper: return 'ISS'
    if 'STARLINK' in name_upper: return 'Starlink'
    if 'ONEWEB' in name_upper: return 'OneWeb'
    if 'PLANET' in name_upper or 'FLOCK' in name_upper: return 'Planet'
    if 'SPIRE' in name_upper or 'LEMUR' in name_upper: return 'Spire'
    if 'IRIDIUM' in name_upper: return 'Iridium'
    if 96 <= inc <= 99: return 'SSO'
    if 50 <= inc <= 55: return 'LEO53'
    return 'Other'


def parse_tles(data):
    """Parse TLE data into satellite list."""
    lines = data.strip().split('\n')
    satellites = []
    i = 0
    while i < len(lines) - 2:
        name, line1, line2 = lines[i].strip(), lines[i + 1].strip(), lines[i + 2].strip()
        if line1.startswith('1 ') and line2.startswith('2 '):
            inc = float(line2[8:16])
            mean_motion = float(line2[52:63])
            period = 1440.0 / mean_motion
            if period <= 225:  # Near-earth only (SGP4 limitation)
                satellites.append({
                    'name': name, 'line1': line1, 'line2': line2,
                    'constellation': classify_satellite(name, inc, period),
                    'epoch': parse_tle_epoch(line1),
                    'data': {
                        'norad_id': line1[2:7].strip(),
                        'inclination': inc,
                        'eccentricity': float('0.' + line2[26:33].strip()),
                        'period': period,
                        'raan': float(line2[17:25].strip()),
                        'mean_motion': mean_motion,
                    }
                })
            i += 3
        else:
            i += 1
    return satellites


def generate_legend_html(counts):
    """Generate legend HTML items."""
    items = [f'<label class="legend-item toggle iss-item" data-constellation="ISS"><input type="checkbox" checked><div class="legend-dot iss-dot">🛰</div>ISS<span class="legend-count">{counts.get("ISS", 0)}</span></label>']
    for name, color in [('Starlink', '100,180,255'), ('OneWeb', '255,160,100'), ('Planet', '130,200,130'),
                        ('Spire', '240,220,140'), ('Iridium', '200,140,220'), ('SSO', '255,130,130'),
                        ('LEO53', '140,220,200'), ('Other', '180,180,190')]:
        label = 'LEO 53°' if name == 'LEO53' else name
        items.append(f'<label class="legend-item toggle" data-constellation="{name}"><input type="checkbox" checked><div class="legend-dot" style="background:rgb({color})"></div>{label}<span class="legend-count">{counts.get(name, 0):,}</span></label>')
    return '\n            '.join(items)


def main():
    from astroz import Tle, Sgp4Constellation

    print("=" * 60)
    print("  Cesium Satellite Visualization")
    print("=" * 60)

    # Download and parse
    print("\n[1/3] Loading catalog...", end=" ", flush=True)
    satellites = parse_tles(download_tles())
    print(f"done - {len(satellites):,} satellites")

    # Propagate
    print("[2/3] Propagating...", end=" ", flush=True)
    t0 = time_module.perf_counter()

    tle_objects = [Tle(f"{s['line1']}\n{s['line2']}") for s in satellites]
    constellation = Sgp4Constellation(tle_objects)

    num_times = 120
    times = np.arange(num_times, dtype=np.float64)
    out = np.empty((num_times * constellation.num_batches * 4 * 6,), dtype=np.float64)
    constellation.propagate_into(times, out)

    prop_time = time_module.perf_counter() - t0
    throughput = len(satellites) * num_times / prop_time / 1e6
    print(f"done - {throughput:.1f}M props/sec")

    # Reshape positions
    num_sats = len(satellites)
    reshaped = out.reshape(constellation.num_batches, num_times, 4, 6)
    positions = reshaped.transpose(1, 0, 2, 3).reshape(num_times, -1, 6)[:, :num_sats, :3] * 1000

    # Prepare template data
    epochs = sorted(s['epoch'] for s in satellites)
    ref_epoch = epochs[len(epochs) // 2]
    counts = Counter(s['constellation'] for s in satellites)

    template_data = {
        '{{NUM_SATS}}': f'{num_sats:,}',
        '{{NUM_SATS_RAW}}': str(num_sats),
        '{{NUM_FRAMES}}': str(num_times),
        '{{START_FRAME}}': '10',
        '{{THROUGHPUT}}': f'{throughput:.1f}',
        '{{TOTAL_PROPS}}': f'{num_sats * num_times:,}',
        '{{PROP_TIME_MS}}': f'{prop_time * 1000:.0f}',
        '{{REFERENCE_EPOCH}}': ref_epoch.isoformat(),
        '{{FRAMES}}': json.dumps([positions[t].flatten().tolist() for t in range(num_times)]),
        '{{SAT_COLORS}}': json.dumps([COLORS.get(s['constellation'], [150, 150, 150]) for s in satellites]),
        '{{SAT_CONSTELLATIONS}}': json.dumps([s['constellation'] for s in satellites]),
        '{{SAT_NAMES}}': json.dumps([s['name'] for s in satellites]),
        '{{SAT_DATA}}': json.dumps([s['data'] for s in satellites]),
        '{{LEGEND_ITEMS}}': generate_legend_html(counts),
    }

    # Generate HTML
    print("[3/3] Generating HTML...", end=" ", flush=True)
    template = (Path(__file__).parent / 'cesium_template.html').read_text()
    html = template
    for key, value in template_data.items():
        html = html.replace(key, value)

    Path('cesium_fast.html').write_text(html)
    print("done")
    print(f"\nCreated: cesium_fast.html ({len(html) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()

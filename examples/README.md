# astroz Examples

## Cesium Satellite Visualization

Interactive 3D visualization of the entire near-earth satellite catalog (~13,000 satellites).

**Files:**
- `cesium_fast.py` - Main script (downloads TLEs, propagates, generates HTML)
- `cesium_template.html` - HTML/JS template for the visualization

### Running

```bash
uv run --python 3.12 examples/cesium_fast.py
```

This generates `cesium_fast.html` in the project root. Open it in your browser.

### Features

- Multithreaded SGP4 propagation (~300M propagations/sec on 16 threads with AVX512)
- Color coded constellations (Starlink, OneWeb, Planet, Spire, Iridium, etc.)
- Toggle constellation visibility
- Search satellites by name (auto-labels when ≤10 results)
- Click for orbital details, double-click to track
- Real time mode showing actual current UTC positions

### Keyboard Shortcuts

| Key | Action |
|-----|--------|
| Space | Pause/resume |
| R | Toggle real-time mode |
| Esc | Stop tracking / reset view |

## Other Examples

- `python_sgp4.py` - Basic SGP4 propagation using python-sgp4 compatible API (`astroz.api`)
- `conjunction_screening.py` - Conjunction detection example
- `sgp4_propagation.zig` - Zig SGP4 example
- `orbit_maneuvers.zig` - Hohmann transfer and orbital maneuvers

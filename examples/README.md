# astroz Examples

The full list of examples is in the [top-level README](../README.md#examples). This page covers the Cesium visualization, which needs a bit more setup.

## Cesium Satellite Visualization

A 3D view of the near-earth satellite catalog (about 13,000 satellites).

- `cesium_fast.py` downloads TLEs, propagates them and writes the HTML page
- `cesium_template.html` is the HTML/JS template it fills in

### Running

```bash
uv run --python 3.12 examples/cesium_fast.py
```

This writes `cesium_fast.html` to the current directory. Open it in a browser.

### What it does

- Propagates the catalog with the multithreaded SGP4 path (about 300M propagations/sec on 16 threads with AVX512)
- Colors satellites by constellation (Starlink, OneWeb, Planet, Spire, Iridium, ...) and lets you toggle each one
- Searches by name and labels the matches when there are 10 or fewer
- Shows orbital details on click and tracks a satellite on double-click
- Has a real-time mode that shows current UTC positions

### Keyboard shortcuts

| Key | Action |
|-----|--------|
| Space | Pause/resume |
| R | Toggle real-time mode |
| Esc | Stop tracking / reset view |

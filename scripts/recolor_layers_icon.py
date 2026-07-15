"""Create an orange map-layers icon that is pixel-identical to the cyan one
except for layer colors (cyan/teal family -> orange family).
Navy background and non-layer pixels are left unchanged.
"""
from __future__ import annotations

import colorsys
from pathlib import Path

from PIL import Image

ASSETS = Path(
    r"C:\Users\pjohnson\.cursor\projects\c-Users-pjohnson-Dropbox-Documents-Python-SAT-Planner\assets"
)
SRC = ASSETS / "map_layers_icon.png"
OUT_ORANGE = ASSETS / "map_layers_icon_active.png"
# Also write a matched pair into the project media folder if it exists
PROJECT_ASSETS = Path(__file__).resolve().parents[1] / "media"


def is_layer_pixel(r: int, g: int, b: int, a: int) -> bool:
    """True for cyan/teal layer pixels; False for navy bg / gray shadows."""
    if a < 8:
        return False
    h, s, v = colorsys.rgb_to_hsv(r / 255.0, g / 255.0, b / 255.0)
    # Cyan/teal/aqua hue band roughly 160–200° -> 0.44–0.56 in [0,1]
    # Allow a bit wider for cyan highlights
    in_cyan_hue = 0.40 <= h <= 0.58
    # Layers are relatively saturated vs navy / soft gray shadows
    if in_cyan_hue and s >= 0.18 and v >= 0.25:
        return True
    # Brighter near-cyan highlights (slightly desaturated)
    if in_cyan_hue and s >= 0.10 and v >= 0.55 and (g > r + 15) and (g >= b - 40 or b > r):
        return True
    return False


def cyan_hsv_to_orange(h: float, s: float, v: float) -> tuple[float, float, float]:
    """Map cyan/teal hue band onto an orange band while preserving s/v structure.

    Cyan ~0.50, teal darker ~0.48–0.52.
    Map linearly onto orange: bright top ~0.08 (orange), bottom ~0.04 (deeper orange/amber).
    """
    # Normalize position within cyan band
    h_clamped = min(max(h, 0.40), 0.58)
    t = (h_clamped - 0.40) / (0.58 - 0.40)  # 0..1 across cyan band
    # Invert lightly so brighter cyan-ish maps to brighter orange-ish;
    # then blend with value so dark layers stay deeper orange.
    # Target orange hues: deep ~0.03, bright ~0.09
    orange_h = 0.035 + (1.0 - t) * 0.055  # ~0.035..0.09

    # Nudge: higher value (brighter layers) -> warmer/brighter orange hue
    orange_h = 0.045 + v * 0.05  # ~0.045 (dark) to ~0.095 (bright)
    # Slightly boost saturation so orange reads clearly against navy
    orange_s = min(1.0, s * 1.05 + 0.05)
    # Keep value almost identical so lighting/shadows match
    orange_v = v
    return orange_h, orange_s, orange_v


def recolor(src: Path, dest: Path) -> None:
    im = Image.open(src).convert("RGBA")
    pixels = im.load()
    w, h = im.size
    converted = 0
    for y in range(h):
        for x in range(w):
            r, g, b, a = pixels[x, y]
            if not is_layer_pixel(r, g, b, a):
                continue
            hv, s, v = colorsys.rgb_to_hsv(r / 255.0, g / 255.0, b / 255.0)
            oh, os_, ov = cyan_hsv_to_orange(hv, s, v)
            nr, ng, nb = colorsys.hsv_to_rgb(oh, os_, ov)
            pixels[x, y] = (int(nr * 255), int(ng * 255), int(nb * 255), a)
            converted += 1
    dest.parent.mkdir(parents=True, exist_ok=True)
    im.save(dest)
    print(f"Wrote {dest} ({w}x{h}), converted {converted} layer pixels")


def main() -> None:
    print(f"Source: {SRC} exists={SRC.exists()}")
    if not SRC.exists():
        raise SystemExit(f"Missing source icon: {SRC}")
    src_im = Image.open(SRC)
    print(f"Source size={src_im.size} mode={src_im.mode}")
    recolor(SRC, OUT_ORANGE)
    # Keep cyan master unchanged; confirm orange matches dimensions
    cyan = Image.open(SRC)
    orange = Image.open(OUT_ORANGE)
    assert cyan.size == orange.size, (cyan.size, orange.size)
    print(f"Pair dimensions match: {cyan.size}")
    # Optional project copy
    if PROJECT_ASSETS.exists():
        (PROJECT_ASSETS / "map_layers_icon.png").write_bytes(SRC.read_bytes())
        recolor(SRC, PROJECT_ASSETS / "map_layers_icon_active.png")


if __name__ == "__main__":
    main()

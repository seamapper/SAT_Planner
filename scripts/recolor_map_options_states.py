"""Build map options icon pair from a bold cyan master:
  - off  = light grey accents
  - on   = orange accents
"""
from __future__ import annotations

import colorsys
from pathlib import Path

from PIL import Image

ASSETS = Path(
    r"C:\Users\pjohnson\.cursor\projects\c-Users-pjohnson-Dropbox-Documents-Python-SAT-Planner\assets"
)
MASTER = ASSETS / "map_options_icon.png"
OUT_OFF = ASSETS / "map_options_icon.png"
OUT_ON = ASSETS / "map_options_icon_active.png"
STAGING = ASSETS / "_map_options_master_cyan.png"
PROJECT_MEDIA = Path(
    r"C:\Users\pjohnson\Dropbox\Documents\Python\SAT_Planner\media"
)


def is_accent_pixel(r: int, g: int, b: int, a: int) -> bool:
    if a < 8:
        return False
    h, s, v = colorsys.rgb_to_hsv(r / 255.0, g / 255.0, b / 255.0)
    in_cyan_hue = 0.40 <= h <= 0.58
    if in_cyan_hue and s >= 0.18 and v >= 0.25:
        return True
    if in_cyan_hue and s >= 0.10 and v >= 0.55 and (g > r + 15) and (g >= b - 40 or b > r):
        return True
    if s >= 0.08 and v >= 0.70 and 0.38 <= h <= 0.60 and g >= r and (b >= r - 10):
        return True
    return False


def cyan_to_orange(h: float, s: float, v: float) -> tuple[float, float, float]:
    orange_h = 0.045 + v * 0.05
    orange_s = min(1.0, s * 1.05 + 0.05)
    return orange_h, orange_s, v


def cyan_to_light_grey(h: float, s: float, v: float) -> tuple[float, float, float]:
    grey_h = 0.0
    grey_s = 0.0
    grey_v = min(1.0, 0.62 + v * 0.32)
    return grey_h, grey_s, grey_v


def recolor(src: Image.Image, mapper) -> Image.Image:
    im = src.convert("RGBA").copy()
    pixels = im.load()
    w, h = im.size
    converted = 0
    for y in range(h):
        for x in range(w):
            r, g, b, a = pixels[x, y]
            if not is_accent_pixel(r, g, b, a):
                continue
            hv, s, v = colorsys.rgb_to_hsv(r / 255.0, g / 255.0, b / 255.0)
            nh, ns, nv = mapper(hv, s, v)
            nr, ng, nb = colorsys.hsv_to_rgb(nh, ns, nv)
            pixels[x, y] = (int(nr * 255), int(ng * 255), int(nb * 255), a)
            converted += 1
    print(f"  converted {converted} accent pixels")
    return im


def main() -> None:
    if not MASTER.exists():
        raise SystemExit(f"Missing master: {MASTER}")

    # Preserve cyan master before overwriting off icon
    STAGING.write_bytes(MASTER.read_bytes())
    master = Image.open(STAGING)
    print(f"Master size={master.size} mode={master.mode}")

    print("Building OFF (light grey)...")
    off = recolor(master, cyan_to_light_grey)
    off.save(OUT_OFF)
    print(f"Wrote {OUT_OFF}")

    print("Building ON (orange)...")
    on = recolor(master, cyan_to_orange)
    on.save(OUT_ON)
    print(f"Wrote {OUT_ON}")

    assert Image.open(OUT_OFF).size == Image.open(OUT_ON).size
    print(f"Pair dimensions match: {Image.open(OUT_OFF).size}")

    if PROJECT_MEDIA.exists():
        (PROJECT_MEDIA / "map_options_icon.png").write_bytes(OUT_OFF.read_bytes())
        (PROJECT_MEDIA / "map_options_icon_active.png").write_bytes(OUT_ON.read_bytes())
        print("Copied to media/")


if __name__ == "__main__":
    main()

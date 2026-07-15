"""Create an orange map-options icon that is pixel-identical to the cyan
one except for cyan/teal accents remapped to orange.
"""
from __future__ import annotations

import colorsys
from pathlib import Path

from PIL import Image

ASSETS = Path(
    r"C:\Users\pjohnson\.cursor\projects\c-Users-pjohnson-Dropbox-Documents-Python-SAT-Planner\assets"
)
SRC = ASSETS / "map_options_icon.png"
OUT_ORANGE = ASSETS / "map_options_icon_active.png"
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
    return False


def cyan_hsv_to_orange(h: float, s: float, v: float) -> tuple[float, float, float]:
    orange_h = 0.045 + v * 0.05
    orange_s = min(1.0, s * 1.05 + 0.05)
    return orange_h, orange_s, v


def recolor(src: Path, dest: Path) -> None:
    im = Image.open(src).convert("RGBA")
    pixels = im.load()
    w, h = im.size
    converted = 0
    for y in range(h):
        for x in range(w):
            r, g, b, a = pixels[x, y]
            if not is_accent_pixel(r, g, b, a):
                continue
            hv, s, v = colorsys.rgb_to_hsv(r / 255.0, g / 255.0, b / 255.0)
            oh, os_, ov = cyan_hsv_to_orange(hv, s, v)
            nr, ng, nb = colorsys.hsv_to_rgb(oh, os_, ov)
            pixels[x, y] = (int(nr * 255), int(ng * 255), int(nb * 255), a)
            converted += 1
    dest.parent.mkdir(parents=True, exist_ok=True)
    im.save(dest)
    print(f"Wrote {dest} ({w}x{h}), converted {converted} accent pixels")


def main() -> None:
    print(f"Source: {SRC} exists={SRC.exists()}")
    if not SRC.exists():
        raise SystemExit(f"Missing source icon: {SRC}")
    src_im = Image.open(SRC)
    print(f"Source size={src_im.size} mode={src_im.mode}")
    recolor(SRC, OUT_ORANGE)
    cyan = Image.open(SRC)
    orange = Image.open(OUT_ORANGE)
    assert cyan.size == orange.size, (cyan.size, orange.size)
    print(f"Pair dimensions match: {cyan.size}")
    if PROJECT_MEDIA.exists():
        (PROJECT_MEDIA / "map_options_icon.png").write_bytes(SRC.read_bytes())
        recolor(SRC, PROJECT_MEDIA / "map_options_icon_active.png")


if __name__ == "__main__":
    main()

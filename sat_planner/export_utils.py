"""
Shared survey export writers. All CSV/TXT formats use the same row conventions:
- DDD/DMM/DMS CSV rows: (line_num, line_name, point_label, lat, lon) with lat/lon in decimal degrees.
  Writers convert to deg/min/sec for DMM/DMS as needed.
- DDD/DMM/DMS TXT rows: (point_label, lat, lon) in decimal degrees.
- SIS asciiplan: lines = list of (line_name, points) with points = list of (lat, lon).
- Hypack LNW: lines = list of (line_name, points) with points = list of (lat, lon); written in UTM.
"""
import csv
import datetime
import glob
import math
import os
from xml.sax.saxutils import escape

try:
    from sat_planner.constants import GEOSPATIAL_LIBS_AVAILABLE, pyproj
except ImportError:
    GEOSPATIAL_LIBS_AVAILABLE = False
    pyproj = None


def remove_export_file(path):
    """Remove an export target (and shapefile sidecars) so the next write overwrites cleanly."""
    if not path:
        return
    targets = {path}
    base, ext = os.path.splitext(path)
    if ext.lower() == ".shp":
        # Only remove known shapefile sidecars; do not wipe unrelated files
        # that share the same basename (e.g. .geojson, .gpx, .asciiplan).
        for shp_ext in (".shp", ".shx", ".dbf", ".prj", ".cpg", ".qix", ".sbn", ".sbx", ".xml"):
            targets.add(base + shp_ext)
    for target in targets:
        try:
            if os.path.isfile(target):
                os.remove(target)
        except OSError:
            pass


def _deg_min(d):
    deg = int(d)
    min_val = (abs(d) - abs(deg)) * 60.0
    return deg, min_val


def _deg_min_sec(d):
    deg = int(d)
    total_mins = (abs(d) - abs(deg)) * 60.0
    mins = int(total_mins)
    secs = (total_mins - mins) * 60.0
    return deg, mins, secs


def write_ddd_csv(path, rows, *, newline='', encoding='utf-8'):
    """Write *_DDD.csv. rows: list of (line_num, line_name, point_label, lat, lon)."""
    remove_export_file(path)
    with open(path, 'w', newline=newline, encoding=encoding) as f:
        w = csv.writer(f)
        w.writerow(['Line Number', 'Line Name', 'Point Label', 'Latitude', 'Longitude'])
        for row in rows:
            line_num, line_name, point_label, lat, lon = row
            w.writerow([line_num, line_name, point_label, lat, lon])


def write_ddd_txt(path, rows, *, encoding='utf-8'):
    """Write *_DDD.txt. rows: list of (line_num, line_name, point_label, lat, lon). Writes label lat lon per line."""
    remove_export_file(path)
    with open(path, 'w', encoding=encoding) as f:
        for row in rows:
            _, _, point_label, lat, lon = row
            f.write(f"{point_label} {lat:.6f} {lon:.6f}\n")


def write_dmm_csv(path, rows, *, newline='', encoding='utf-8'):
    """Write *_DMM.csv. rows: list of (line_num, line_name, point_label, lat, lon). Converts to deg/min."""
    remove_export_file(path)
    with open(path, 'w', newline=newline, encoding=encoding) as f:
        w = csv.writer(f)
        w.writerow(['Line Number', 'Line Name', 'Point Label', 'Latitude (Deg)', 'Latitude (Min)', 'Longitude (Deg)', 'Longitude (Min)'])
        for row in rows:
            line_num, line_name, point_label, lat, lon = row
            lat_d, lat_m = _deg_min(lat)
            lon_d, lon_m = _deg_min(lon)
            w.writerow([line_num, line_name, point_label, lat_d, lat_m, lon_d, lon_m])


def write_dmm_txt(path, rows, *, encoding='utf-8'):
    """Write *_DMM.txt. rows: list of (line_num, line_name, point_label, lat, lon). One line: label lat_deg lat_min lon_deg lon_min."""
    remove_export_file(path)
    with open(path, 'w', encoding=encoding) as f:
        for row in rows:
            _, _, point_label, lat, lon = row
            lat_d, lat_m = _deg_min(lat)
            lon_d, lon_m = _deg_min(lon)
            f.write(f"{point_label} {lat_d} {lat_m} {lon_d} {lon_m}\n")


def write_dms_csv(path, rows, *, newline='', encoding='utf-8'):
    """Write *_DMS.csv. rows: list of (line_num, line_name, point_label, lat, lon). Converts to deg/min/sec."""
    remove_export_file(path)
    with open(path, 'w', newline=newline, encoding=encoding) as f:
        w = csv.writer(f)
        w.writerow(['Line Number', 'Line Name', 'Point Label', 'Latitude (Deg)', 'Latitude (Min)', 'Latitude (Sec)', 'Longitude (Deg)', 'Longitude (Min)', 'Longitude (Sec)'])
        for row in rows:
            line_num, line_name, point_label, lat, lon = row
            lat_d, lat_m, lat_s = _deg_min_sec(lat)
            lon_d, lon_m, lon_s = _deg_min_sec(lon)
            w.writerow([line_num, line_name, point_label, lat_d, lat_m, lat_s, lon_d, lon_m, lon_s])


def write_dms_txt(path, rows, *, encoding='utf-8'):
    """Write *_DMS.txt. rows: list of (line_num, line_name, point_label, lat, lon). One line: label lat_deg lat_min lat_sec lon_deg lon_min lon_sec."""
    remove_export_file(path)
    with open(path, 'w', encoding=encoding) as f:
        for row in rows:
            _, _, point_label, lat, lon = row
            lat_d, lat_m, lat_s = _deg_min_sec(lat)
            lon_d, lon_m, lon_s = _deg_min_sec(lon)
            f.write(f"{point_label} {lat_d} {lat_m} {lat_s} {lon_d} {lon_m} {lon_s}\n")


def write_profile_csv(path, distances, elevations, slopes, *, newline="", encoding="utf-8"):
    """Write crossline/pitch-style profile CSV: Distance (m), Elevation (m), Slope (deg).

    ``distances``, ``elevations``, ``slopes`` are equal-length sequences (e.g. numpy arrays).
    Non-finite elevation/slope cells are written as empty fields.
    """
    remove_export_file(path)
    with open(path, "w", newline=newline, encoding=encoding) as f:
        w = csv.writer(f)
        w.writerow(["Distance", "Elevation", "Slope"])

        def _cell(x):
            try:
                v = float(x)
            except (TypeError, ValueError):
                return ""
            if math.isnan(v) or math.isinf(v):
                return ""
            return v

        for d, e, s in zip(distances, elevations, slopes):
            w.writerow([_cell(d), _cell(e), _cell(s)])


def write_line_plan_profile_csv(
    path, distances, elevations, slopes, waypoint_labels, *, newline="", encoding="utf-8"
):
    """Write line-plan profile CSV: Distance (m), Elevation (m), Slope (deg), Waypoint.

    ``waypoint_labels`` is the same length as the numeric arrays; use \"\" for rows
    with no waypoint. Non-finite elevation/slope cells are written as empty fields.
    """
    remove_export_file(path)
    with open(path, "w", newline=newline, encoding=encoding) as f:
        w = csv.writer(f)
        w.writerow(["Distance", "Elevation", "Slope", "Waypoint"])

        def _cell(x):
            try:
                v = float(x)
            except (TypeError, ValueError):
                return ""
            if math.isnan(v) or math.isinf(v):
                return ""
            return v

        for d, e, s, wp in zip(distances, elevations, slopes, waypoint_labels):
            wp_cell = wp if (wp is not None and str(wp).strip() != "") else ""
            w.writerow([_cell(d), _cell(e), _cell(s), wp_cell])


def sanitize_export_basename(name):
    """Make a string safe for use as a Windows filename stem (no path)."""
    bad = '<>:"/\\|?*'
    s = str(name).strip()
    for c in bad:
        s = s.replace(c, "_")
    return s.strip().strip(".")


def write_gpx_per_test_files(export_dir, export_base, tests, *, creator="SAT Planner"):
    """Write one GPX file per test: ``{export_base}_{suffix}.gpx`` with a single ``<trk>``.

    ``tests`` is an iterable of ``(filename_suffix, gpx_track_name, points)`` where
    ``points`` is a list of ``(lat, lon)``. ``gpx_track_name`` is the ``<name>`` inside the GPX.

    Returns a list of basenames (filenames) written successfully.
    """
    base = sanitize_export_basename(export_base)
    written = []
    for suffix, trk_name, pts in tests:
        if not pts:
            continue
        stem = sanitize_export_basename(suffix)
        if not stem:
            continue
        path = os.path.join(export_dir, f"{base}_{stem}.gpx")
        if write_gpx(path, [(trk_name, pts)], creator=creator):
            written.append(os.path.basename(path))
    return written


def write_asciiplan(path, lines, *, encoding='utf-8'):
    """Write SIS .asciiplan. lines: list of (line_name, points) where points is list of (lat, lon) in decimal degrees."""
    remove_export_file(path)
    ts = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    with open(path, 'w', encoding=encoding) as f:
        f.write("DEG\n\n0 0 0 0\n")
        for line_index, (line_name, pts) in enumerate(lines):
            coords = " ".join(f"{lat:.6f} {lon:.6f}" for lat, lon in pts)
            f.write(f'_LINE {line_name} {line_index} {ts} 0 {coords} "\n')


def write_gpx(path, lines, *, creator="SAT Planner", encoding='utf-8'):
    """Write GPX 1.1 with one track per line.
    lines: list of (line_name, points) where points is list of (lat, lon).
    Returns True if written, False otherwise.
    """
    remove_export_file(path)
    ts = datetime.datetime.utcnow().replace(microsecond=0).isoformat() + "Z"
    try:
        with open(path, "w", encoding=encoding) as f:
            f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            f.write(
                '<gpx version="1.1" creator="{creator}" xmlns="http://www.topografix.com/GPX/1/1" '
                'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
                'xsi:schemaLocation="http://www.topografix.com/GPX/1/1 '
                'http://www.topografix.com/GPX/1/1/gpx.xsd">\n'.format(creator=escape(str(creator)))
            )
            f.write(f"  <metadata><time>{ts}</time></metadata>\n")
            for line_name, pts in lines:
                if not pts:
                    continue
                safe_name = escape(str(line_name))
                f.write("  <trk>\n")
                f.write(f"    <name>{safe_name}</name>\n")
                f.write("    <trkseg>\n")
                for lat, lon in pts:
                    f.write(f'      <trkpt lat="{lat:.8f}" lon="{lon:.8f}"></trkpt>\n')
                f.write("    </trkseg>\n")
                f.write("  </trk>\n")
            f.write("</gpx>\n")
    except OSError:
        return False
    return True


def compute_utm_zone_from_points(points):
    """From a list of (lat, lon) in WGS84, return (utm_zone, hemisphere) for LNW/filename use.
    points: list of (lat, lon). Returns (zone 1-60, 'North' or 'South')."""
    if not points:
        return (18, "North")
    lats = [p[0] for p in points]
    lons = [p[1] for p in points]
    lon_centroid = sum(lons) / len(lons)
    lat_centroid = sum(lats) / len(lats)
    zone = int((lon_centroid + 180) / 6) + 1
    zone = max(1, min(60, zone))
    hemisphere = "North" if lat_centroid >= 0 else "South"
    return (zone, hemisphere)


def write_lnw(path, lines, *, encoding='utf-8'):
    """Write HYPACK LNW format (LIN/PTS/UTM) to path.
    lines: list of (line_name, [(lat, lon), ...]) with at least 2 points per line.
    Coordinates are converted to UTM from centroid of all points. Returns True if written, False otherwise."""
    if not GEOSPATIAL_LIBS_AVAILABLE or pyproj is None:
        return False
    all_points = []
    for _name, pts in lines:
        all_points.extend(pts)
    if not all_points:
        return False
    zone, hemisphere = compute_utm_zone_from_points(all_points)
    utm_crs = f"EPSG:326{zone:02d}" if hemisphere == "North" else f"EPSG:327{zone:02d}"
    try:
        transformer = pyproj.Transformer.from_crs("EPSG:4326", utm_crs, always_xy=True)
    except Exception:
        return False
    remove_export_file(path)
    try:
        with open(path, "w", encoding=encoding) as f:
            f.write("LNW 1.0\n")
            for name, pts in lines:
                if len(pts) < 2:
                    continue
                utm_pts = []
                for lat, lon in pts:
                    try:
                        easting, northing = transformer.transform(lon, lat)
                        utm_pts.append((easting, northing))
                    except Exception:
                        return False
                f.write(f"LIN {len(utm_pts)}\n")
                for easting, northing in utm_pts:
                    f.write(f"PTS {easting:.4f} {northing:.4f}\n")
                if name:
                    f.write(f"LNN {name}\n")
                f.write("EOL\n")
    except OSError:
        return False
    return True


EMAIL_LOW_PNG_MAX_DIMENSION_PX = 1280
MAP_EXPORT_WIDTH_IN = 11.0
MAP_EXPORT_HEIGHT_IN = 8.5
MAP_EXPORT_HIGH_DPI = 300
MAP_EXPORT_LOW_DPI = 100


def low_resolution_png_path(path):
    """Return path for email-friendly copy with ``_low`` before the extension."""
    root, ext = os.path.splitext(path)
    if not ext:
        ext = ".png"
    if root.endswith("_low"):
        return path
    return f"{root}_low{ext}"


def _scaled_png(image, max_width, max_height):
    """Resize an image so it fits inside a box. It may grow or shrink."""
    from PIL import Image

    max_width = max(1, int(max_width))
    max_height = max(1, int(max_height))
    scale = min(max_width / float(image.width), max_height / float(image.height))
    size = (max(1, int(round(image.width * scale))), max(1, int(round(image.height * scale))))
    if size == image.size:
        return image
    return image.resize(size, Image.Resampling.LANCZOS)


def _map_pixels_within_limit(map_w, map_h, extra_w, extra_h, limit):
    """Map area after the combined image is fit to ``limit`` on its long side."""
    longest = max(map_w + extra_w, map_h + extra_h, 1)
    scale = min(1.0, float(limit) / float(longest))
    return (map_w * scale) * (map_h * scale)


def composite_colorbar_png(map_png, colorbar_png, horizontal_png=None, *, max_dimension=EMAIL_LOW_PNG_MAX_DIMENSION_PX):
    """Add a colorbar beside or under a map image, in the margin. Returns PNG bytes.

    The side that leaves more of the map inside ``max_dimension`` is used.
    A wide map gains the scale underneath. A tall map gains it on the right.
    """
    import io

    if not map_png or not colorbar_png:
        return map_png
    try:
        from PIL import Image
    except ImportError:
        return map_png
    try:
        with Image.open(io.BytesIO(map_png)) as base_image, Image.open(io.BytesIO(colorbar_png)) as bar_image:
            base = base_image.convert("RGBA")
            vertical = bar_image.convert("RGBA")
            margin = max(8, int(max(base.width, base.height) * 0.012))
            vertical = _scaled_png(vertical, base.width * 0.22, base.height * 0.78)
            right_extra = margin + vertical.width
            horizontal = None
            bottom_extra = 0
            if horizontal_png:
                with Image.open(io.BytesIO(horizontal_png)) as horizontal_image:
                    horizontal = _scaled_png(
                        horizontal_image.convert("RGBA"),
                        base.width * 0.88,
                        base.height * 0.22,
                    )
                    bottom_extra = margin + horizontal.height
            use_bottom = (
                horizontal is not None
                and _map_pixels_within_limit(base.width, base.height, 0, bottom_extra, max_dimension)
                > _map_pixels_within_limit(base.width, base.height, right_extra, 0, max_dimension)
            )
            if use_bottom:
                canvas = Image.new("RGBA", (base.width, base.height + bottom_extra), (255, 255, 255, 255))
                canvas.paste(base, (0, 0))
                x = max(0, (base.width - horizontal.width) // 2)
                canvas.alpha_composite(horizontal, (x, base.height + margin))
            else:
                canvas = Image.new("RGBA", (base.width + right_extra, base.height), (255, 255, 255, 255))
                canvas.paste(base, (0, 0))
                y = max(0, (base.height - vertical.height) // 2)
                canvas.alpha_composite(vertical, (base.width + margin, y))
            buffer = io.BytesIO()
            canvas.convert("RGB").save(buffer, format="PNG")
            return buffer.getvalue()
    except Exception:
        return map_png


def composite_legend_png(map_png, legend_png):
    """Place a legend at the lower left of a map image. Returns PNG bytes."""
    import io

    if not map_png or not legend_png:
        return map_png
    try:
        from PIL import Image
    except ImportError:
        return map_png
    try:
        with Image.open(io.BytesIO(map_png)) as base_image, Image.open(io.BytesIO(legend_png)) as legend_image:
            base = base_image.convert("RGBA")
            legend = legend_image.convert("RGBA")
            max_height = max(1, int(base.height * 0.36))
            max_width = max(1, int(base.width * 0.42))
            scale = min(max_height / float(legend.height), max_width / float(legend.width), 1.0)
            if scale < 1.0:
                legend = legend.resize(
                    (max(1, int(legend.width * scale)), max(1, int(legend.height * scale))),
                    Image.Resampling.LANCZOS,
                )
            margin = max(8, int(base.width * 0.015))
            x = margin
            y = max(0, base.height - legend.height - margin)
            base.alpha_composite(legend, (x, y))
            buffer = io.BytesIO()
            base.convert("RGB").save(buffer, format="PNG")
            return buffer.getvalue()
    except Exception:
        return map_png


def map_export_page_px(dpi):
    """Pixel size of the letter-landscape map page at ``dpi``."""
    return (
        max(1, int(round(MAP_EXPORT_WIDTH_IN * float(dpi)))),
        max(1, int(round(MAP_EXPORT_HEIGHT_IN * float(dpi)))),
    )


def write_png_fitted_to_page(
    path,
    data,
    *,
    width_in=MAP_EXPORT_WIDTH_IN,
    height_in=MAP_EXPORT_HEIGHT_IN,
    dpi=MAP_EXPORT_LOW_DPI,
):
    """Fit image bytes on a white page and write the email-sized ``*_low.png``."""
    import io

    low_path = low_resolution_png_path(path)
    remove_export_file(low_path)
    if not data:
        return None
    page_w = max(1, int(round(float(width_in) * float(dpi))))
    page_h = max(1, int(round(float(height_in) * float(dpi))))
    try:
        from PIL import Image

        with Image.open(io.BytesIO(data)) as image:
            fitted = _scaled_png(image.convert("RGB"), page_w, page_h)
            page = Image.new("RGB", (page_w, page_h), (255, 255, 255))
            page.paste(fitted, ((page_w - fitted.width) // 2, (page_h - fitted.height) // 2))
            page.save(low_path, format="PNG", optimize=True, dpi=(float(dpi), float(dpi)))
        return low_path
    except Exception:
        return None


def write_low_resolution_png_bytes(path, data, *, low_max_dimension_px=EMAIL_LOW_PNG_MAX_DIMENSION_PX):
    """Write map-capture bytes to the email-sized ``*_low.png`` beside ``path``."""
    import io

    low_path = low_resolution_png_path(path)
    remove_export_file(low_path)
    if not data:
        return None
    try:
        from PIL import Image

        with Image.open(io.BytesIO(data)) as im:
            im = im.convert("RGB")
            width, height = im.size
            longest = max(width, height)
            if longest > low_max_dimension_px:
                scale = low_max_dimension_px / float(longest)
                im = im.resize(
                    (max(1, int(width * scale)), max(1, int(height * scale))),
                    Image.Resampling.LANCZOS,
                )
            im.save(low_path, format="PNG", optimize=True)
        return low_path
    except Exception:
        try:
            with open(low_path, "wb") as handle:
                handle.write(data)
            return low_path
        except OSError:
            return None


def _write_low_resolution_png_from_file(path, *, low_max_dimension_px=EMAIL_LOW_PNG_MAX_DIMENSION_PX):
    """Resize an existing PNG to ``*_low.png``; returns low path or None on failure."""
    low_path = low_resolution_png_path(path)
    try:
        from PIL import Image
    except ImportError:
        return None
    try:
        with Image.open(path) as im:
            im = im.convert("RGB")
            w, h = im.size
            longest = max(w, h)
            if longest > low_max_dimension_px:
                scale = low_max_dimension_px / float(longest)
                new_size = (max(1, int(w * scale)), max(1, int(h * scale)))
                im = im.resize(new_size, Image.Resampling.LANCZOS)
            im.save(low_path, format="PNG", optimize=True)
        return low_path
    except Exception:
        return None


def _figure_save_kwargs(*, dpi=300, bbox_inches="tight", facecolor="white"):
    save_kwargs = {"dpi": dpi}
    if bbox_inches is not None:
        save_kwargs["bbox_inches"] = bbox_inches
    if facecolor is not None:
        save_kwargs["facecolor"] = facecolor
    return save_kwargs


def save_export_png(
    fig,
    path,
    *,
    save_high=True,
    save_low=True,
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
    low_max_dimension_px=EMAIL_LOW_PNG_MAX_DIMENSION_PX,
    low_dpi=96,
    low_match_figure=False,
):
    """Save figure PNG(s) per high/low toggles. Returns ``(high_path_or_none, low_path_or_none)``."""
    if not save_high and not save_low:
        return None, None
    if save_high:
        remove_export_file(path)
    if save_low:
        remove_export_file(low_resolution_png_path(path))
    save_kwargs = _figure_save_kwargs(dpi=dpi, bbox_inches=bbox_inches, facecolor=facecolor)
    high_path = path if save_high else None
    low_path = low_resolution_png_path(path) if save_low else None
    if save_high:
        fig.savefig(path, **save_kwargs)
    if save_low:
        if low_match_figure or not (save_high and os.path.isfile(path)):
            low_kwargs = dict(save_kwargs)
            low_kwargs["dpi"] = low_dpi
            fig.savefig(low_path, **low_kwargs)
        else:
            written = _write_low_resolution_png_from_file(path, low_max_dimension_px=low_max_dimension_px)
            if written is None:
                low_kwargs = dict(save_kwargs)
                low_kwargs["dpi"] = low_dpi
                fig.savefig(low_path, **low_kwargs)
    return high_path, low_path


def save_figure_png_with_low_res(
    fig,
    path,
    *,
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
    low_max_dimension_px=EMAIL_LOW_PNG_MAX_DIMENSION_PX,
    low_dpi=96,
):
    """Save figure at full resolution and create an email-friendly ``*_low.png`` copy."""
    _, low_path = save_export_png(
        fig,
        path,
        save_high=True,
        save_low=True,
        dpi=dpi,
        bbox_inches=bbox_inches,
        facecolor=facecolor,
        low_max_dimension_px=low_max_dimension_px,
        low_dpi=low_dpi,
    )
    return low_path


def save_existing_png_with_low_res(path, *, low_max_dimension_px=EMAIL_LOW_PNG_MAX_DIMENSION_PX):
    """Create ``*_low.png`` from an existing PNG (e.g. after a custom savefig)."""
    return _write_low_resolution_png_from_file(path, low_max_dimension_px=low_max_dimension_px)


def png_export_basenames(path, *, include_high=True, include_low=True):
    """Basenames for exported PNGs when the files exist."""
    names = []
    if not path:
        return names
    if include_high and os.path.isfile(path):
        names.append(os.path.basename(path))
    if include_low:
        low_path = low_resolution_png_path(path)
        if os.path.isfile(low_path):
            names.append(os.path.basename(low_path))
    return names

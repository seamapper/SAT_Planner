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
import math
import os
from xml.sax.saxutils import escape

try:
    from sat_planner.constants import GEOSPATIAL_LIBS_AVAILABLE, pyproj
except ImportError:
    GEOSPATIAL_LIBS_AVAILABLE = False
    pyproj = None


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
    with open(path, 'w', newline=newline, encoding=encoding) as f:
        w = csv.writer(f)
        w.writerow(['Line Number', 'Line Name', 'Point Label', 'Latitude', 'Longitude'])
        for row in rows:
            line_num, line_name, point_label, lat, lon = row
            w.writerow([line_num, line_name, point_label, lat, lon])


def write_ddd_txt(path, rows, *, encoding='utf-8'):
    """Write *_DDD.txt. rows: list of (line_num, line_name, point_label, lat, lon). Writes label lat lon per line."""
    with open(path, 'w', encoding=encoding) as f:
        for row in rows:
            _, _, point_label, lat, lon = row
            f.write(f"{point_label} {lat:.6f} {lon:.6f}\n")


def write_dmm_csv(path, rows, *, newline='', encoding='utf-8'):
    """Write *_DMM.csv. rows: list of (line_num, line_name, point_label, lat, lon). Converts to deg/min."""
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
    with open(path, 'w', encoding=encoding) as f:
        for row in rows:
            _, _, point_label, lat, lon = row
            lat_d, lat_m = _deg_min(lat)
            lon_d, lon_m = _deg_min(lon)
            f.write(f"{point_label} {lat_d} {lat_m} {lon_d} {lon_m}\n")


def write_dms_csv(path, rows, *, newline='', encoding='utf-8'):
    """Write *_DMS.csv. rows: list of (line_num, line_name, point_label, lat, lon). Converts to deg/min/sec."""
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
                f.write(f"LNN {name}\n")
                f.write("EOL\n")
    except OSError:
        return False
    return True

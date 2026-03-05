"""
Shared survey export writers. All CSV/TXT formats use the same row conventions:
- DDD/DMM/DMS CSV rows: (line_num, line_name, point_label, lat, lon) with lat/lon in decimal degrees.
  Writers convert to deg/min/sec for DMM/DMS as needed.
- DDD/DMM/DMS TXT rows: (point_label, lat, lon) in decimal degrees.
"""
import csv
import os


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

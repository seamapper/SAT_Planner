"""
Survey file parsers: DDD, DMS, DMM, LNW (line and polyline).
Shared by Calibration, Reference, and Line Planning import.
Expects host to provide _show_message().
"""
import os
import re

from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel,
                             QPushButton, QSpinBox, QComboBox)

from sat_planner.constants import GEOSPATIAL_LIBS_AVAILABLE, fiona, pyproj


# Match the common "UTM<zone><N|S>" tag in filenames, e.g. Survey_UTM18N.lnw,
# survey-utm5s_2024.lnw, foo_UTM01N_bar.lnw. The hemisphere letter must not be
# immediately followed by another letter so "UTM18North" or "UTM18Norway" still
# match (they end with N + non-letter). 1-2 digit zone, case-insensitive.
_UTM_FILENAME_PATTERN = re.compile(r"UTM(\d{1,2})([NS])(?![A-Za-z])", re.IGNORECASE)


def parse_utm_zone_from_filename(file_path):
    """Extract UTM zone and hemisphere from a filename if it follows the
    common ``*_UTM<zone><N|S>*`` convention.

    Args:
        file_path: Full path or bare filename to inspect.

    Returns:
        (zone:int, hemisphere:str) tuple where hemisphere is "North" or "South",
        or (None, None) if no valid match is found. Zone must be 1..60.
    """
    if not file_path:
        return (None, None)
    basename = os.path.basename(str(file_path))
    match = _UTM_FILENAME_PATTERN.search(basename)
    if not match:
        return (None, None)
    try:
        zone = int(match.group(1))
    except (TypeError, ValueError):
        return (None, None)
    if not (1 <= zone <= 60):
        return (None, None)
    hemisphere = "North" if match.group(2).upper() == "N" else "South"
    return (zone, hemisphere)


class UTMZoneDialog(QDialog):
    """Dialog for getting UTM zone and hemisphere from user."""

    def __init__(self, parent, default_zone=18, default_hemisphere="North", detected_from_filename=False):
        super().__init__(parent)
        self.setWindowTitle("UTM Zone Information")
        self.setMinimumWidth(300)
        self.setModal(True)

        self.utm_zone = None
        self.hemisphere = None

        try:
            default_zone_int = int(default_zone) if default_zone is not None else 18
        except (TypeError, ValueError):
            default_zone_int = 18
        if not (1 <= default_zone_int <= 60):
            default_zone_int = 18
        default_hemi = str(default_hemisphere or "North").strip().capitalize()
        if default_hemi not in ("North", "South"):
            default_hemi = "North"

        layout = QVBoxLayout(self)
        layout.setSpacing(10)
        layout.setContentsMargins(15, 15, 15, 15)

        if detected_from_filename:
            instructions_text = (
                f"Detected UTM zone {default_zone_int} {default_hemi} from the file name. "
                "Confirm or change before importing the LNW file:"
            )
        else:
            instructions_text = "Please enter the UTM zone information for the LNW file:"
        instructions = QLabel(instructions_text)
        instructions.setWordWrap(True)
        layout.addWidget(instructions)

        zone_layout = QHBoxLayout()
        zone_label = QLabel("UTM Zone:")
        zone_label.setMinimumWidth(100)
        zone_layout.addWidget(zone_label)

        self.zone_spinbox = QSpinBox()
        self.zone_spinbox.setMinimum(1)
        self.zone_spinbox.setMaximum(60)
        self.zone_spinbox.setValue(default_zone_int)
        zone_layout.addWidget(self.zone_spinbox)
        zone_layout.addStretch()
        layout.addLayout(zone_layout)

        hemisphere_layout = QHBoxLayout()
        hemisphere_label = QLabel("Hemisphere:")
        hemisphere_label.setMinimumWidth(100)
        hemisphere_layout.addWidget(hemisphere_label)

        self.hemisphere_combo = QComboBox()
        self.hemisphere_combo.addItems(["North", "South"])
        self.hemisphere_combo.setCurrentText(default_hemi)
        hemisphere_layout.addWidget(self.hemisphere_combo)
        hemisphere_layout.addStretch()
        layout.addLayout(hemisphere_layout)

        button_layout = QHBoxLayout()
        button_layout.addStretch()

        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        button_layout.addWidget(cancel_btn)

        ok_btn = QPushButton("OK")
        ok_btn.clicked.connect(self._on_ok_clicked)
        button_layout.addWidget(ok_btn)

        layout.addLayout(button_layout)

    def _on_ok_clicked(self):
        self.utm_zone = self.zone_spinbox.value()
        self.hemisphere = self.hemisphere_combo.currentText()
        self.accept()

    def get_utm_info(self):
        return (self.utm_zone, self.hemisphere)

    @classmethod
    def for_file(cls, parent, file_path):
        """Build a UTMZoneDialog with zone/hemisphere pre-filled from the
        filename when possible (e.g. ``Survey_UTM18N.lnw`` -> 18, North).
        Falls back to (18, North) when the pattern is absent."""
        zone, hemi = parse_utm_zone_from_filename(file_path)
        if zone is None:
            return cls(parent)
        return cls(parent, default_zone=zone, default_hemisphere=hemi, detected_from_filename=True)


class SurveyParsersMixin:
    """
    Parsers for survey file formats: DDD, DMS, DMM, LNW.
    Line parsers return list of lines, each [(lat1,lon1), (lat2,lon2)].
    Polyline parsers return list of (lat, lon) in order.
    """

    def _dms_to_decimal_degrees(self, degrees, minutes, seconds, is_latitude=True):
        try:
            sign = -1 if degrees < 0 else 1
            abs_degrees = abs(degrees)
            decimal = abs_degrees + (minutes / 60.0) + (seconds / 3600.0)
            return sign * decimal
        except (ValueError, TypeError, ZeroDivisionError):
            raise ValueError(f"Invalid DMS values: {degrees}° {minutes}' {seconds}\"")

    def _dmm_to_decimal_degrees(self, degrees, decimal_minutes):
        try:
            sign = -1 if degrees < 0 else 1
            abs_degrees = abs(degrees)
            decimal = abs_degrees + (decimal_minutes / 60.0)
            return sign * decimal
        except (ValueError, TypeError, ZeroDivisionError):
            raise ValueError(f"Invalid DMM values: {degrees}° {decimal_minutes}'")

    def _parse_dms_txt_file(self, file_path):
        """Parse *_DMS.txt: even rows, 7 cols. Pairs form lines. Returns list of lines or None."""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            lines = [line.strip() for line in content.replace('\r\n', '\n').replace('\r', '\n').split('\n') if line.strip()]
            if len(lines) < 2:
                self._show_message("error", "Import Error",
                                 f"Expected at least 2 data rows in *_DMS.txt file, found {len(lines)} rows.")
                return None
            if len(lines) % 2 != 0:
                self._show_message("error", "Import Error",
                                 f"Expected an even number of data rows (pairs form lines). Found {len(lines)} rows.")
                return None
            points = []
            for i, line in enumerate(lines):
                parts = line.split()
                if len(parts) < 7:
                    self._show_message("error", "Import Error",
                                     f"Row {i+1} does not have 7 columns (Point, deg_lat, min_lat, sec_lat, deg_lon, min_lon, sec_lon). Found {len(parts)} columns.")
                    return None
                try:
                    deg_lat = float(parts[1])
                    min_lat = float(parts[2])
                    sec_lat = float(parts[3])
                    deg_lon = float(parts[4])
                    min_lon = float(parts[5])
                    sec_lon = float(parts[6])
                    lat = self._dms_to_decimal_degrees(deg_lat, min_lat, sec_lat, is_latitude=True)
                    lon = self._dms_to_decimal_degrees(deg_lon, min_lon, sec_lon, is_latitude=False)
                    points.append((lat, lon))
                except (ValueError, IndexError) as e:
                    self._show_message("error", "Import Error",
                                     f"Row {i+1} has invalid DMS values: {str(e)}")
                    return None
            imported_lines = []
            for i in range(0, len(points), 2):
                if i + 1 < len(points):
                    imported_lines.append([points[i], points[i+1]])
                else:
                    break
            return imported_lines
        except FileNotFoundError:
            self._show_message("error", "Import Error", f"File not found: {file_path}")
            return None
        except Exception as e:
            self._show_message("error", "Import Error", f"Failed to parse *_DMS.txt file: {e}")
            import traceback
            traceback.print_exc()
            return None

    def _parse_dmm_txt_file(self, file_path):
        """Parse *_DMM.txt: even rows, 5 cols. Pairs form lines. Returns list of lines or None."""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            lines = [line.strip() for line in content.replace('\r\n', '\n').replace('\r', '\n').split('\n') if line.strip()]
            if len(lines) < 2:
                self._show_message("error", "Import Error",
                                 f"Expected at least 2 data rows in *_DMM.txt file, found {len(lines)} rows.")
                return None
            if len(lines) % 2 != 0:
                self._show_message("error", "Import Error",
                                 f"Expected an even number of data rows (pairs form lines). Found {len(lines)} rows.")
                return None
            points = []
            for i, line in enumerate(lines):
                parts = line.split()
                if len(parts) < 5:
                    self._show_message("error", "Import Error",
                                     f"Row {i+1} does not have 5 columns (Point, deg_lat, min_lat, deg_lon, min_lon). Found {len(parts)} columns.")
                    return None
                try:
                    deg_lat = float(parts[1])
                    min_lat = float(parts[2])
                    deg_lon = float(parts[3])
                    min_lon = float(parts[4])
                    lat = self._dmm_to_decimal_degrees(deg_lat, min_lat)
                    lon = self._dmm_to_decimal_degrees(deg_lon, min_lon)
                    points.append((lat, lon))
                except (ValueError, IndexError) as e:
                    self._show_message("error", "Import Error",
                                     f"Row {i+1} has invalid DMM values: {str(e)}")
                    return None
            imported_lines = []
            for i in range(0, len(points), 2):
                if i + 1 < len(points):
                    imported_lines.append([points[i], points[i+1]])
                else:
                    break
            return imported_lines
        except FileNotFoundError:
            self._show_message("error", "Import Error", f"File not found: {file_path}")
            return None
        except Exception as e:
            self._show_message("error", "Import Error", f"Failed to parse *_DMM.txt file: {e}")
            import traceback
            traceback.print_exc()
            return None

    def _parse_ddd_txt_file(self, file_path):
        """Parse *_DDD.txt: even rows, 3 cols. Pairs form lines. Returns list of lines or None."""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            lines = [line.strip() for line in content.replace('\r\n', '\n').replace('\r', '\n').split('\n') if line.strip()]
            if len(lines) < 2:
                self._show_message("error", "Import Error",
                                 f"Expected at least 2 data rows in *_DDD.txt file, found {len(lines)} rows.")
                return None
            if len(lines) % 2 != 0:
                self._show_message("error", "Import Error",
                                 f"Expected an even number of data rows (pairs form lines). Found {len(lines)} rows.")
                return None
            points = []
            for i, line in enumerate(lines):
                parts = line.split()
                if len(parts) < 3:
                    self._show_message("error", "Import Error",
                                     f"Row {i+1} does not have 3 columns (Point, Latitude, Longitude). Found {len(parts)} columns.")
                    return None
                try:
                    lat = float(parts[1])
                    lon = float(parts[2])
                    points.append((lat, lon))
                except (ValueError, IndexError) as e:
                    self._show_message("error", "Import Error",
                                     f"Row {i+1} has invalid latitude/longitude values: {parts[1] if len(parts) > 1 else 'N/A'}, {parts[2] if len(parts) > 2 else 'N/A'}")
                    return None
            imported_lines = []
            for i in range(0, len(points), 2):
                if i + 1 < len(points):
                    imported_lines.append([points[i], points[i+1]])
                else:
                    break
            return imported_lines
        except FileNotFoundError:
            self._show_message("error", "Import Error", f"File not found: {file_path}")
            return None
        except Exception as e:
            self._show_message("error", "Import Error", f"Failed to parse *_DDD.txt file: {e}")
            import traceback
            traceback.print_exc()
            return None

    def _parse_ddd_txt_file_as_polyline(self, file_path):
        """Parse *_DDD.txt as single polyline: all rows in order, 3 cols. Returns list of (lat, lon) or None. Requires >= 2 points."""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            lines = [line.strip() for line in content.replace('\r\n', '\n').replace('\r', '\n').split('\n') if line.strip()]
            if len(lines) < 2:
                self._show_message("error", "Import Error",
                                 f"Line plan requires at least 2 points in *_DDD.txt file, found {len(lines)} rows.")
                return None
            points = []
            for i, line in enumerate(lines):
                parts = line.split()
                if len(parts) < 3:
                    self._show_message("error", "Import Error",
                                     f"Row {i+1} does not have 3 columns (Point, Latitude, Longitude). Found {len(parts)} columns.")
                    return None
                try:
                    lat = float(parts[1])
                    lon = float(parts[2])
                    points.append((lat, lon))
                except (ValueError, IndexError):
                    self._show_message("error", "Import Error",
                                     f"Row {i+1} has invalid latitude/longitude values.")
                    return None
            return points
        except FileNotFoundError:
            self._show_message("error", "Import Error", f"File not found: {file_path}")
            return None
        except Exception as e:
            self._show_message("error", "Import Error", f"Failed to parse *_DDD.txt file: {e}")
            return None

    def _parse_dms_txt_file_as_polyline(self, file_path):
        """Parse *_DMS.txt as single polyline: all rows in order, 7 cols. Returns list of (lat, lon) or None. Requires >= 2 points."""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            lines = [line.strip() for line in content.replace('\r\n', '\n').replace('\r', '\n').split('\n') if line.strip()]
            if len(lines) < 2:
                self._show_message("error", "Import Error",
                                 f"Line plan requires at least 2 points in *_DMS.txt file, found {len(lines)} rows.")
                return None
            points = []
            for i, line in enumerate(lines):
                parts = line.split()
                if len(parts) < 7:
                    self._show_message("error", "Import Error",
                                     f"Row {i+1} does not have 7 columns (Point, deg_lat, min_lat, sec_lat, deg_lon, min_lon, sec_lon). Found {len(parts)} columns.")
                    return None
                try:
                    deg_lat, min_lat, sec_lat = float(parts[1]), float(parts[2]), float(parts[3])
                    deg_lon, min_lon, sec_lon = float(parts[4]), float(parts[5]), float(parts[6])
                    lat = self._dms_to_decimal_degrees(deg_lat, min_lat, sec_lat, is_latitude=True)
                    lon = self._dms_to_decimal_degrees(deg_lon, min_lon, sec_lon, is_latitude=False)
                    points.append((lat, lon))
                except (ValueError, IndexError):
                    self._show_message("error", "Import Error", f"Row {i+1} has invalid DMS values.")
                    return None
            return points
        except FileNotFoundError:
            self._show_message("error", "Import Error", f"File not found: {file_path}")
            return None
        except Exception as e:
            self._show_message("error", "Import Error", f"Failed to parse *_DMS.txt file: {e}")
            return None

    def _parse_dmm_txt_file_as_polyline(self, file_path):
        """Parse *_DMM.txt as single polyline: all rows in order, 5 cols. Returns list of (lat, lon) or None. Requires >= 2 points."""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            lines = [line.strip() for line in content.replace('\r\n', '\n').replace('\r', '\n').split('\n') if line.strip()]
            if len(lines) < 2:
                self._show_message("error", "Import Error",
                                 f"Line plan requires at least 2 points in *_DMM.txt file, found {len(lines)} rows.")
                return None
            points = []
            for i, line in enumerate(lines):
                parts = line.split()
                if len(parts) < 5:
                    self._show_message("error", "Import Error",
                                     f"Row {i+1} does not have 5 columns (Point, deg_lat, min_lat, deg_lon, min_lon). Found {len(parts)} columns.")
                    return None
                try:
                    deg_lat, min_lat = float(parts[1]), float(parts[2])
                    deg_lon, min_lon = float(parts[3]), float(parts[4])
                    lat = self._dmm_to_decimal_degrees(deg_lat, min_lat)
                    lon = self._dmm_to_decimal_degrees(deg_lon, min_lon)
                    points.append((lat, lon))
                except (ValueError, IndexError):
                    self._show_message("error", "Import Error", f"Row {i+1} has invalid DMM values.")
                    return None
            return points
        except FileNotFoundError:
            self._show_message("error", "Import Error", f"File not found: {file_path}")
            return None
        except Exception as e:
            self._show_message("error", "Import Error", f"Failed to parse *_DMM.txt file: {e}")
            return None

    def _parse_lnw_file(self, file_path, utm_zone, hemisphere):
        """Parse LNW file. Returns list of lines, each [(lat1,lon1), (lat2,lon2)], or None. Uses first/last point per LIN block."""
        try:
            if not GEOSPATIAL_LIBS_AVAILABLE or pyproj is None:
                self._show_message("error", "Import Error",
                                 "Geospatial libraries not available. Cannot import LNW files.")
                return None
            if hemisphere.lower() == 'north' or hemisphere.lower() == 'n':
                utm_crs = f"EPSG:326{utm_zone:02d}"
            else:
                utm_crs = f"EPSG:327{utm_zone:02d}"
            try:
                transformer = pyproj.Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True)
            except Exception as e:
                self._show_message("error", "Import Error",
                                 f"Failed to create UTM transformer for zone {utm_zone} {hemisphere}: {e}")
                return None
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            lines = [line.strip() for line in content.replace('\r\n', '\n').replace('\r', '\n').split('\n') if line.strip()]
            if not lines:
                self._show_message("error", "Import Error", "LNW file is empty.")
                return None
            if lines[0].upper().startswith('LNW'):
                lines = lines[1:]
            imported_lines = []
            i = 0
            while i < len(lines):
                line = lines[i].upper().strip()
                if line.startswith('LIN'):
                    parts = line.split()
                    if len(parts) < 2:
                        i += 1
                        continue
                    try:
                        num_points = int(parts[1])
                    except ValueError:
                        i += 1
                        continue
                    points_utm = []
                    i += 1
                    points_read = 0
                    while i < len(lines) and points_read < num_points:
                        pts_line = lines[i].upper().strip()
                        if pts_line.startswith('PTS'):
                            parts = pts_line.split()
                            if len(parts) >= 3:
                                try:
                                    easting = float(parts[1])
                                    northing = float(parts[2])
                                    points_utm.append((easting, northing))
                                    points_read += 1
                                except (ValueError, IndexError):
                                    pass
                        elif pts_line.startswith('LNN') or pts_line.startswith('EOL'):
                            break
                        i += 1
                    if len(points_utm) >= 2:
                        points_geo = []
                        for easting, northing in points_utm:
                            try:
                                lon, lat = transformer.transform(easting, northing)
                                points_geo.append((lat, lon))
                            except Exception as e:
                                self._show_message("error", "Import Error",
                                                 f"Failed to convert UTM coordinates: {e}")
                                return None
                        imported_lines.append([points_geo[0], points_geo[-1]])
                    while i < len(lines):
                        if lines[i].upper().strip().startswith('EOL'):
                            i += 1
                            break
                        elif lines[i].upper().strip().startswith('LIN'):
                            break
                        i += 1
                else:
                    i += 1
            if len(imported_lines) == 0:
                self._show_message("error", "Import Error",
                                 "No valid lines found in LNW file.")
                return None
            return imported_lines
        except FileNotFoundError:
            self._show_message("error", "Import Error", f"File not found: {file_path}")
            return None
        except Exception as e:
            self._show_message("error", "Import Error", f"Failed to parse LNW file: {e}")
            import traceback
            traceback.print_exc()
            return None

    def _parse_lnw_file_as_polyline(self, file_path, utm_zone, hemisphere):
        """Parse first LIN block of LNW as single polyline. Returns list of (lat, lon) or None. Requires >= 2 points."""
        try:
            if not GEOSPATIAL_LIBS_AVAILABLE or pyproj is None:
                self._show_message("error", "Import Error",
                                 "Geospatial libraries not available. Cannot import LNW files.")
                return None
            if hemisphere.lower() == 'north' or hemisphere.lower() == 'n':
                utm_crs = f"EPSG:326{utm_zone:02d}"
            else:
                utm_crs = f"EPSG:327{utm_zone:02d}"
            try:
                transformer = pyproj.Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True)
            except Exception as e:
                self._show_message("error", "Import Error",
                                 f"Failed to create UTM transformer for zone {utm_zone} {hemisphere}: {e}")
                return None
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            lines = [line.strip() for line in content.replace('\r\n', '\n').replace('\r', '\n').split('\n') if line.strip()]
            if not lines:
                self._show_message("error", "Import Error", "LNW file is empty.")
                return None
            if lines[0].upper().startswith('LNW'):
                lines = lines[1:]
            i = 0
            while i < len(lines):
                line = lines[i].upper().strip()
                if line.startswith('LIN'):
                    parts = line.split()
                    if len(parts) < 2:
                        i += 1
                        continue
                    try:
                        num_points = int(parts[1])
                    except ValueError:
                        i += 1
                        continue
                    points_utm = []
                    i += 1
                    points_read = 0
                    while i < len(lines) and points_read < num_points:
                        pts_line = lines[i].upper().strip()
                        if pts_line.startswith('PTS') and len(pts_line.split()) >= 3:
                            try:
                                p = pts_line.split()
                                easting, northing = float(p[1]), float(p[2])
                                points_utm.append((easting, northing))
                                points_read += 1
                            except (ValueError, IndexError):
                                pass
                        elif pts_line.startswith('LNN') or pts_line.startswith('EOL'):
                            break
                        i += 1
                    if len(points_utm) < 2:
                        self._show_message("error", "Import Error",
                                         "Line plan requires at least 2 points in the first line block of the LNW file.")
                        return None
                    points_geo = []
                    for easting, northing in points_utm:
                        try:
                            lon, lat = transformer.transform(easting, northing)
                            points_geo.append((lat, lon))
                        except Exception as e:
                            self._show_message("error", "Import Error", f"Failed to convert UTM coordinates: {e}")
                            return None
                    return points_geo
                i += 1
            self._show_message("error", "Import Error", "No valid line block found in LNW file.")
            return None
        except FileNotFoundError:
            self._show_message("error", "Import Error", f"File not found: {file_path}")
            return None
        except Exception as e:
            self._show_message("error", "Import Error", f"Failed to parse LNW file: {e}")
            return None

    # ------------------------------------------------------------------
    # Shapefile / GeoPackage parsing
    # ------------------------------------------------------------------
    # Both formats are read through Fiona, so a single low-level reader
    # handles them. Each feature's geometry is reprojected to EPSG:4326
    # (WGS84 lat/lon) so callers can hand the result straight into the
    # plan-specific UI without further CRS work.
    #
    # Three public helpers are exposed:
    #   _parse_vector_file_as_line_list(path)
    #       -> [[(lat1, lon1), (lat2, lon2)], ...] using the first/last
    #          point of every LineString in file order. Used by
    #          Calibration, Accuracy, and Performance (each of which
    #          then routes the lines through its existing line-
    #          assignment dialog).
    #
    #   _parse_vector_file_as_polyline(path)
    #       -> [(lat, lon), ...] from the first LineString in the file.
    #          Used by Line Planning.
    #
    #   _parse_vector_file_as_polyline_or_polygon(path)
    #       -> [(lat, lon), ...] from the first LineString, or, if none
    #          exists, from the outer ring of the first Polygon. Used
    #          by Backscatter so users can drop in either a centerline
    #          or an area polygon.
    #
    # All three return None on hard failure (after showing a message),
    # or an empty list when the file is structurally fine but contains
    # no usable geometry. The caller's existing "must have >= 2 points"
    # warnings then fire normally.

    @staticmethod
    def _shapefile_supported_extensions():
        """File extensions handled by the shapefile/GPKG parsers."""
        return (".shp", ".gpkg")

    def _open_vector_features(self, file_path):
        """Open a shapefile/GeoPackage and yield (geom_in_wgs84_dict,) tuples.

        Returns a list of geometry dicts in EPSG:4326 (WGS84) ordered as they
        appear in the source. Returns ``None`` on hard failure (with a user
        message already shown). For GeoPackage files containing multiple
        layers the first layer is read; that matches Fiona's default and
        keeps the API simple. If a user needs another layer they can drop
        it into the first slot or open it in GIS software first.
        """
        if not GEOSPATIAL_LIBS_AVAILABLE or fiona is None:
            self._show_message(
                "error",
                "Import Error",
                "Geospatial libraries (fiona) not available. Cannot import shapefile/GeoPackage.",
            )
            return None
        try:
            from fiona.transform import transform_geom
        except Exception as e:
            self._show_message("error", "Import Error", f"Failed to import fiona.transform: {e}")
            return None

        ext = os.path.splitext(file_path)[1].lower()
        if ext == ".shp":
            # Quick sanity-check the bundle. A bare .shp without .shx/.dbf is
            # unusable and Fiona's error is opaque; this gives a clearer one.
            base = os.path.splitext(file_path)[0]
            missing = [s for s in (".shx", ".dbf") if not os.path.exists(base + s)]
            if missing:
                self._show_message(
                    "error",
                    "Import Error",
                    "Shapefile is missing required sidecar file(s): "
                    + ", ".join(missing)
                    + ".\nA shapefile is a bundle of .shp + .shx + .dbf (and ideally .prj). "
                    "Make sure all parts are in the same folder.",
                )
                return None

        try:
            geoms_out = []
            with fiona.open(file_path, "r") as source:
                source_crs = source.crs_wkt or source.crs
                if not source_crs:
                    self._show_message(
                        "warning",
                        "Import Warning",
                        "No CRS metadata found in the file (missing .prj?). "
                        "Coordinates will be assumed to already be WGS84 lon/lat.",
                    )
                for feature in source:
                    geom = None
                    if hasattr(feature, "get"):
                        try:
                            geom = feature.get("geometry")
                        except Exception:
                            geom = None
                    if geom is None:
                        try:
                            geom = feature["geometry"]
                        except Exception:
                            geom = getattr(feature, "geometry", None)
                    if geom is None:
                        continue
                    if hasattr(geom, "__geo_interface__"):
                        geom = geom.__geo_interface__
                    elif not isinstance(geom, dict):
                        try:
                            geom = dict(geom)
                        except Exception:
                            continue
                    if source_crs:
                        try:
                            geom = transform_geom(source_crs, "EPSG:4326", geom)
                        except Exception as e:
                            print(f"Warning: failed to reproject feature: {e}")
                            continue
                    # Fiona >= 1.9 returns ``fiona.Geometry`` objects (not plain
                    # dicts) from both feature access and ``transform_geom``.
                    # Normalize to a GeoJSON-style dict so downstream consumers
                    # can use ``geom["type"]`` / ``geom["coordinates"]`` without
                    # caring which Fiona version produced the value.
                    if not isinstance(geom, dict):
                        if hasattr(geom, "__geo_interface__"):
                            geom = geom.__geo_interface__
                        else:
                            try:
                                geom = dict(geom)
                            except Exception:
                                continue
                    geoms_out.append(geom)
            return geoms_out
        except FileNotFoundError:
            self._show_message("error", "Import Error", f"File not found: {file_path}")
            return None
        except Exception as e:
            self._show_message("error", "Import Error", f"Failed to read shapefile/GeoPackage: {e}")
            return None

    @staticmethod
    def _normalize_geom_type(geom_type):
        """Normalize a GeoJSON / Fiona geometry-type string for matching.

        Strips trailing Z / M / ZM / 25D markers and removes whitespace so
        all of these compare equal to ``LINESTRING``:
        ``LineString``, ``LineStringZ``, ``LineString Z``, ``LineString25D``,
        ``LineString 25D``, ``3D LineString``.
        """
        t = (geom_type or "").upper().replace(" ", "")
        if t.startswith("3D"):
            t = t[2:]
        if t.endswith("ZM"):
            t = t[:-2]
        elif t.endswith("25D"):
            t = t[:-3]
        elif t.endswith("Z") or t.endswith("M"):
            t = t[:-1]
        return t

    # Kept under the old name for any external callers that may still rely
    # on it; new code should call ``_normalize_geom_type``.
    _strip_dimensional_suffix = _normalize_geom_type

    @staticmethod
    def _as_geom_dict(geom):
        """Coerce a Fiona Geometry / shapely object / dict-like into a plain dict."""
        if geom is None:
            return None
        if isinstance(geom, dict):
            return geom
        if hasattr(geom, "__geo_interface__"):
            try:
                return geom.__geo_interface__
            except Exception:
                return None
        try:
            return dict(geom)
        except Exception:
            return None

    @classmethod
    def _iter_linestrings_from_geom(cls, geom):
        """Yield each contained LineString coord list (in lon/lat order).
        Handles LineString, MultiLineString, and GeometryCollection, regardless
        of whether the input is a plain dict or a Fiona/shapely geometry object.
        """
        geom = cls._as_geom_dict(geom)
        if not isinstance(geom, dict):
            return
        gt = cls._normalize_geom_type(geom.get("type"))
        coords = geom.get("coordinates")
        if gt == "LINESTRING":
            if isinstance(coords, list) and len(coords) >= 2:
                yield coords
        elif gt == "MULTILINESTRING":
            if isinstance(coords, list):
                for line_coords in coords:
                    if isinstance(line_coords, list) and len(line_coords) >= 2:
                        yield line_coords
        elif gt == "GEOMETRYCOLLECTION":
            for sub in geom.get("geometries", []) or []:
                yield from cls._iter_linestrings_from_geom(sub)

    @classmethod
    def _iter_polygon_outer_rings(cls, geom):
        """Yield outer-ring coord lists for Polygon / MultiPolygon (lon/lat order)."""
        geom = cls._as_geom_dict(geom)
        if not isinstance(geom, dict):
            return
        gt = cls._normalize_geom_type(geom.get("type"))
        coords = geom.get("coordinates")
        if gt == "POLYGON":
            if isinstance(coords, list) and len(coords) >= 1 and isinstance(coords[0], list) and len(coords[0]) >= 3:
                yield coords[0]
        elif gt == "MULTIPOLYGON":
            if isinstance(coords, list):
                for poly_coords in coords:
                    if (
                        isinstance(poly_coords, list)
                        and len(poly_coords) >= 1
                        and isinstance(poly_coords[0], list)
                        and len(poly_coords[0]) >= 3
                    ):
                        yield poly_coords[0]
        elif gt == "GEOMETRYCOLLECTION":
            for sub in geom.get("geometries", []) or []:
                yield from cls._iter_polygon_outer_rings(sub)

    @staticmethod
    def _coords_lonlat_to_latlon_pairs(coords):
        """Convert a list of (lon, lat[, z]) tuples to [(lat, lon), ...]."""
        out = []
        for c in coords or []:
            if not isinstance(c, (list, tuple)) or len(c) < 2:
                continue
            try:
                lon = float(c[0])
                lat = float(c[1])
            except (TypeError, ValueError):
                continue
            out.append((lat, lon))
        return out

    @classmethod
    def _summarize_geom_types(cls, geoms):
        """Return a comma-separated summary of geometry types in ``geoms``,
        suitable for embedding in user-visible error messages."""
        types_seen = {}
        for g in geoms or []:
            g = cls._as_geom_dict(g)
            if isinstance(g, dict):
                t = str(g.get("type") or "(unknown)")
            else:
                t = type(g).__name__
            types_seen[t] = types_seen.get(t, 0) + 1
        if not types_seen:
            return "no features"
        return ", ".join(f"{n}x {t}" for t, n in types_seen.items())

    def _parse_vector_file_as_line_list(self, file_path):
        """Parse a shapefile/GeoPackage as a list of survey-line endpoint pairs.

        Returns ``[[(lat1, lon1), (lat2, lon2)], ...]`` using the first and
        last vertex of every LineString in the source (matching how LNW/
        DDD/DMM/DMS imports already shape line data). Returns ``None`` on
        hard failure and ``[]`` when the file has no LineString features.
        """
        geoms = self._open_vector_features(file_path)
        if geoms is None:
            return None
        imported_lines = []
        for geom in geoms:
            for line_coords in self._iter_linestrings_from_geom(geom):
                pts = self._coords_lonlat_to_latlon_pairs(line_coords)
                if len(pts) >= 2:
                    imported_lines.append([pts[0], pts[-1]])
        if not imported_lines:
            self._show_message(
                "error",
                "Import Error",
                "No LineString geometry found in the selected shapefile/GeoPackage. "
                f"Found: {self._summarize_geom_types(geoms)}.",
            )
            return None
        return imported_lines

    def _parse_vector_file_as_polyline(self, file_path):
        """Parse the first LineString of a shapefile/GeoPackage as a polyline.

        Returns ``[(lat, lon), ...]`` (preserving every vertex) for the
        first LineString found, or ``None`` on hard failure. Returns ``[]``
        when no LineStrings are present.
        """
        geoms = self._open_vector_features(file_path)
        if geoms is None:
            return None
        for geom in geoms:
            for line_coords in self._iter_linestrings_from_geom(geom):
                pts = self._coords_lonlat_to_latlon_pairs(line_coords)
                if len(pts) >= 2:
                    return pts
        self._show_message(
            "error",
            "Import Error",
            "No LineString geometry found in the selected shapefile/GeoPackage. "
            f"Found: {self._summarize_geom_types(geoms)}.",
        )
        return None

    def _parse_vector_file_as_polyline_or_polygon(self, file_path):
        """Like ``_parse_vector_file_as_polyline``, but if no LineStrings are
        present falls back to the outer ring of the first Polygon /
        MultiPolygon. Used by Backscatter, which can accept either a
        centerline-style polyline or an area polygon as its source geometry.
        """
        geoms = self._open_vector_features(file_path)
        if geoms is None:
            return None
        for geom in geoms:
            for line_coords in self._iter_linestrings_from_geom(geom):
                pts = self._coords_lonlat_to_latlon_pairs(line_coords)
                if len(pts) >= 2:
                    return pts
        for geom in geoms:
            for ring_coords in self._iter_polygon_outer_rings(geom):
                pts = self._coords_lonlat_to_latlon_pairs(ring_coords)
                if len(pts) >= 2:
                    return pts
        self._show_message(
            "error",
            "Import Error",
            "No LineString or Polygon geometry found in the selected shapefile/GeoPackage. "
            f"Found: {self._summarize_geom_types(geoms)}.",
        )
        return None

    def _compute_utm_zone_from_points(self, points):
        """From a list of (lat, lon) in WGS84, return (utm_zone, hemisphere) for LNW export.
        points: list of (lat, lon). Returns (zone 1-60, 'North' or 'South')."""
        if not points or pyproj is None:
            return (18, "North")
        lats = [p[0] for p in points]
        lons = [p[1] for p in points]
        lon_centroid = sum(lons) / len(lons)
        lat_centroid = sum(lats) / len(lats)
        zone = int((lon_centroid + 180) / 6) + 1
        zone = max(1, min(60, zone))
        hemisphere = "North" if lat_centroid >= 0 else "South"
        return (zone, hemisphere)

    def _write_lnw_file(self, file_path, lines_with_names):
        """Write proper HYPACK LNW format (LIN/PTS/UTM) to file_path.
        lines_with_names: list of (line_name, [(lat, lon), ...]) with at least 2 points per line.
        Coordinates are converted to UTM from centroid of all points."""
        if not GEOSPATIAL_LIBS_AVAILABLE or pyproj is None:
            return False
        all_points = []
        for _name, pts in lines_with_names:
            all_points.extend(pts)
        if not all_points:
            return False
        zone, hemisphere = self._compute_utm_zone_from_points(all_points)
        utm_crs = f"EPSG:326{zone:02d}" if hemisphere == "North" else f"EPSG:327{zone:02d}"
        try:
            transformer = pyproj.Transformer.from_crs("EPSG:4326", utm_crs, always_xy=True)
        except Exception:
            return False
        with open(file_path, "w", encoding="utf-8") as f:
            f.write("LNW 1.0\n")
            for name, pts in lines_with_names:
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
        return True

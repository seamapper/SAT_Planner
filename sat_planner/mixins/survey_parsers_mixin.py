"""
Survey file parsers: DDD, DMS, DMM, LNW (line and polyline).
Shared by Calibration, Reference, and Line Planning import.
Expects host to provide _show_message().
"""
import os

from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel,
                             QPushButton, QSpinBox, QComboBox)

from sat_planner.constants import GEOSPATIAL_LIBS_AVAILABLE, pyproj


class UTMZoneDialog(QDialog):
    """Dialog for getting UTM zone and hemisphere from user."""

    def __init__(self, parent):
        super().__init__(parent)
        self.setWindowTitle("UTM Zone Information")
        self.setMinimumWidth(300)
        self.setModal(True)

        self.utm_zone = None
        self.hemisphere = None

        layout = QVBoxLayout(self)
        layout.setSpacing(10)
        layout.setContentsMargins(15, 15, 15, 15)

        instructions = QLabel("Please enter the UTM zone information for the LNW file:")
        instructions.setWordWrap(True)
        layout.addWidget(instructions)

        zone_layout = QHBoxLayout()
        zone_label = QLabel("UTM Zone:")
        zone_label.setMinimumWidth(100)
        zone_layout.addWidget(zone_label)

        self.zone_spinbox = QSpinBox()
        self.zone_spinbox.setMinimum(1)
        self.zone_spinbox.setMaximum(60)
        self.zone_spinbox.setValue(18)
        zone_layout.addWidget(self.zone_spinbox)
        zone_layout.addStretch()
        layout.addLayout(zone_layout)

        hemisphere_layout = QHBoxLayout()
        hemisphere_label = QLabel("Hemisphere:")
        hemisphere_label.setMinimumWidth(100)
        hemisphere_layout.addWidget(hemisphere_label)

        self.hemisphere_combo = QComboBox()
        self.hemisphere_combo.addItems(["North", "South"])
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
            if len(imported_lines) != 4:
                self._show_message("warning", "Import Warning",
                                 f"LNW file contains {len(imported_lines)} lines. Expected 4 lines for calibration survey.")
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

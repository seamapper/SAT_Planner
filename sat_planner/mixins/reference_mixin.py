"""
Reference tab: survey lines validation, export/import, info text, reset, statistics.
"""
import os
import csv
import json

from PyQt6.QtWidgets import (QFileDialog, QDialog, QVBoxLayout, QHBoxLayout,
                             QLabel, QComboBox, QPushButton, QWidget)
from PyQt6.QtGui import QTextCursor

from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas

import numpy as np
from sat_planner.constants import GEOSPATIAL_LIBS_AVAILABLE, pyproj, LineString, fiona
from sat_planner import decimal_degrees_to_ddm
from sat_planner.utils_ui import show_statistics_dialog

try:
    from .survey_parsers_mixin import UTMZoneDialog
except ImportError:
    UTMZoneDialog = None


def _line_heading_degrees(geod, line):
    """Return forward azimuth (0-360) from first point to second point of a 2-point line."""
    if not line or len(line) < 2:
        return None
    (lat1, lon1), (lat2, lon2) = line[0], line[1]
    try:
        fwd_az, _, _ = geod.inv(lon1, lat1, lon2, lat2)
        return fwd_az % 360.0
    except Exception:
        return None


def _suggest_crossline_index(imported_lines):
    """Suggest which line index is the crossline (orientation most different from the median).
    Uses orientation 0-180 so opposite directions (e.g. 0° vs 180°) are treated as the same
    line direction; the crossline is usually perpendicular to the reference lines."""
    if not imported_lines or not pyproj:
        return 0
    geod = pyproj.Geod(ellps="WGS84")
    headings = []
    for line in imported_lines:
        h = _line_heading_degrees(geod, line)
        headings.append(h if h is not None else 0.0)
    if len(headings) < 2:
        return 0
    # Normalize to orientation 0-180 (so 0° and 180° count as same direction)
    headings_arr = np.array(headings)
    orientation = headings_arr % 180.0
    median_orient = np.median(orientation)
    # Angular difference in 0-180 space
    diffs = np.abs(orientation - median_orient)
    diffs = np.minimum(diffs, 180.0 - diffs)
    return int(np.argmax(diffs))


class ReferenceLineAssignmentDialog(QDialog):
    """Dialog for assigning imported lines to crossline vs reference lines."""

    def __init__(self, parent, imported_lines):
        super().__init__(parent)
        self.setWindowTitle("Assign Reference Survey Lines")
        self.setMinimumSize(900, 700)
        self.setModal(True)
        self.imported_lines = imported_lines
        self.assignments = {}

        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(10)

        instructions = QLabel("Assign one line as Crossline and the rest as Reference Lines:")
        instructions.setStyleSheet("font-weight: bold;")
        main_layout.addWidget(instructions)

        self.figure = Figure(figsize=(8, 6))
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self.figure)
        main_layout.addWidget(self.canvas)

        self._plot_lines()

        assignment_widget = QWidget()
        assignment_layout = QVBoxLayout(assignment_widget)
        assignment_layout.setContentsMargins(5, 5, 5, 5)
        assignment_label = QLabel("Line Assignments:")
        assignment_label.setStyleSheet("font-weight: bold;")
        assignment_layout.addWidget(assignment_label)

        suggested_cross = _suggest_crossline_index(imported_lines)
        line_colors = ['red', 'blue', 'green', 'purple', 'orange', 'cyan', 'magenta', 'yellow']
        num_lines = len(imported_lines)
        ref_options = ["Reference Line " + str(i + 1) for i in range(num_lines)]
        self.comboboxes = []
        for i in range(num_lines):
            line_widget = QWidget()
            line_layout = QHBoxLayout(line_widget)
            line_layout.setContentsMargins(0, 0, 0, 0)
            line_label = QLabel(f"Line {i+1}:")
            line_label.setMinimumWidth(60)
            line_layout.addWidget(line_label)
            combo = QComboBox()
            combo.addItem("")
            combo.addItem("Crossline")
            for opt in ref_options:
                combo.addItem(opt)
            if i == suggested_cross:
                combo.setCurrentIndex(1)  # Crossline
            self.comboboxes.append(combo)
            line_layout.addWidget(combo)
            color_label = QLabel("●")
            color_label.setStyleSheet(f"color: {line_colors[i % len(line_colors)]}; font-size: 16pt;")
            color_label.setMinimumWidth(30)
            line_layout.addWidget(color_label)
            assignment_layout.addWidget(line_widget)

        main_layout.addWidget(assignment_widget)

        button_layout = QHBoxLayout()
        button_layout.addStretch()
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        button_layout.addWidget(cancel_btn)
        ok_btn = QPushButton("OK")
        ok_btn.clicked.connect(self._on_ok_clicked)
        button_layout.addWidget(ok_btn)
        main_layout.addLayout(button_layout)

    def _plot_lines(self):
        self.ax.clear()
        line_colors = ['red', 'blue', 'green', 'purple', 'orange', 'cyan', 'magenta', 'yellow']
        all_lats, all_lons = [], []
        for i, line in enumerate(self.imported_lines):
            if len(line) == 2:
                (lat1, lon1), (lat2, lon2) = line
                self.ax.plot([lon1, lon2], [lat1, lat2], color=line_colors[i % len(line_colors)],
                           linewidth=2, label=f'Line {i+1}', zorder=10)
                mid_lat = (lat1 + lat2) / 2
                mid_lon = (lon1 + lon2) / 2
                self.ax.text(mid_lon, mid_lat, f'Line {i+1}', fontsize=10, ha='center', va='center',
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.7), zorder=11)
                all_lats.extend([lat1, lat2])
                all_lons.extend([lon1, lon2])
        if all_lats and all_lons:
            lat_range = max(all_lats) - min(all_lats)
            lon_range = max(all_lons) - min(all_lons)
            pad_lat = max(lat_range * 0.1, 0.01)
            pad_lon = max(lon_range * 0.1, 0.01)
            self.ax.set_xlim(min(all_lons) - pad_lon, max(all_lons) + pad_lon)
            self.ax.set_ylim(min(all_lats) - pad_lat, max(all_lats) + pad_lat)
        self.ax.set_xlabel('Longitude')
        self.ax.set_ylabel('Latitude')
        self.ax.grid(True, alpha=0.3)
        self.ax.legend(loc='upper right')
        self.figure.tight_layout()
        self.canvas.draw()

    def _on_ok_clicked(self):
        from PyQt6.QtWidgets import QMessageBox
        crossline_idx = None
        for i, combo in enumerate(self.comboboxes):
            idx = combo.currentIndex()
            if idx == 0:
                QMessageBox.warning(self, "Assignment Error", "Please assign every line.")
                return
            if idx == 1:  # Crossline
                if crossline_idx is not None:
                    QMessageBox.warning(self, "Assignment Error", "Exactly one line must be Crossline.")
                    return
                crossline_idx = i
                self.assignments[i] = 'crossline'
            else:
                # idx 2 = Reference Line 1, idx 3 = Reference Line 2, ...
                ref_num = idx - 1
                self.assignments[i] = f'ref_{ref_num}'
        if crossline_idx is None:
            QMessageBox.warning(self, "Assignment Error", "Exactly one line must be Crossline.")
            return
        self.accept()

    def get_assignments(self):
        return self.assignments.copy()


class ReferenceMixin:
    """Mixin for reference planning: _validate_inputs, _export_survey_data, _import_survey_files,
    set_ref_info_text, _reset_reference_tab, _calculate_total_survey_time,
    _calculate_reference_survey_statistics, _show_reference_planning_info."""

    def _validate_inputs(self):
        """Validates all input fields and returns a tuple (is_valid, values_dict)."""
        values = {}
        try:
            # Central Latitude
            central_lat_str = self.central_lat_entry.text().strip()
            if not central_lat_str:
                self._show_message("error","Input Error", "Central Latitude cannot be empty.")
                return False, {}
            values['central_lat'] = float(central_lat_str)
            if not (-90 <= values['central_lat'] <= 90):
                self._show_message("error","Input Error", "Central Latitude must be between -90 and 90.")
                return False, {}

            # Central Longitude
            central_lon_str = self.central_lon_entry.text().strip()
            if not central_lon_str:
                self._show_message("error","Input Error", "Central Longitude cannot be empty.")
                return False, {}
            values['central_lon'] = float(central_lon_str)
            if not (-180 <= values['central_lon'] <= 180):
                self._show_message("error","Input Error", "Central Longitude must be between -180 and 180.")
                return False, {}

            # Line Length
            line_length_str = self.line_length_entry.text().strip()
            if not line_length_str:
                self._show_message("error","Input Error", "Line Length cannot be empty.")
                return False, {}
            values['line_length'] = float(line_length_str)
            if not (values['line_length'] > 0):
                self._show_message("error","Input Error", "Line Length must be positive.")
                return False, {}

            # Heading
            heading_str = self.heading_entry.text().strip()
            if not heading_str:
                self._show_message("error","Input Error", "Heading cannot be empty.")
                return False, {}
            values['heading'] = float(heading_str)
            if not (0 <= values['heading'] <= 360):
                self._show_message("error","Input Error", "Heading must be between 0 and 360 degrees.")
                return False, {}

            # Distance Between Lines
            dist_between_lines_str = self.dist_between_lines_entry.text().strip()
            if not dist_between_lines_str:
                self._show_message("error","Input Error", "Distance Between Lines cannot be empty.")
                return False, {}
            values['dist_between_lines'] = float(dist_between_lines_str)
            if not (values['dist_between_lines'] > 0):
                self._show_message("error","Input Error", "Distance Between Lines must be positive.")
                return False, {}

            # Number of Lines
            num_lines_str = self.num_lines_entry.text().strip()
            if not num_lines_str:
                self._show_message("error","Input Error", "Number of Lines cannot be empty.")
                return False, {}
            values['num_lines'] = int(num_lines_str)
            if not (values['num_lines'] > 0):
                self._show_message("error","Input Error", "Number of Lines must be a positive integer.")
                return False, {}

            # Crossline Lead-in/out
            bisect_lead_str = self.bisect_lead_entry.text().strip()
            if not bisect_lead_str:
                self._show_message("error","Input Error", "Crossline Lead-in/out cannot be empty.")
                return False, {}
            values['bisect_lead'] = float(bisect_lead_str)
            if not (values['bisect_lead'] >= 0):
                self._show_message("error","Input Error", "Crossline Lead-in/out cannot be negative.")
                return False, {}

            # Survey Speed
            survey_speed_str = self.survey_speed_entry.text().strip()
            if not survey_speed_str:
                self._show_message("error","Input Error", "Survey Speed cannot be empty.")
                return False, {}
            values['survey_speed'] = float(survey_speed_str)
            if not (values['survey_speed'] > 0):
                self._show_message("error","Input Error", "Survey Speed must be positive.")
                return False, {}

            # Export Name
            values['export_name'] = self.export_name_entry.text().strip()
            if not values['export_name']:
                self._show_message("error","Input Error", "Export Name cannot be empty.")
                return False, {}

            values['offset_direction'] = self.offset_direction_var
            values['line_length_multiplier'] = self.line_length_multiplier  # Get slider value
            values['dist_between_lines_multiplier'] = self.dist_between_lines_multiplier  # Get new slider value

        except ValueError as e:
            self._show_message("error","Input Error", f"Invalid numeric input: {str(e)}")
            return False, {}
        except Exception as e:
            self._show_message("error","Error", f"An unexpected error occurred during input validation: {str(e)}")
            return False, {}

        return True, values

    def _export_survey_data(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot export data.")
            return

        is_valid, values = self._validate_inputs()
        if not is_valid:
            return

        export_name = values['export_name']

        if not self.survey_lines_data and not self.cross_line_data:
            self._show_message("warning","No Data", "No survey lines to export. Generate them first.")
            return

        export_dir = QFileDialog.getExistingDirectory(self, "Select Export Directory", self.last_export_dir)
        if not export_dir:
            return
        self.last_export_dir = export_dir
        self._save_last_export_dir()

        try:
            # --- Export to CSV ---
            # CSV files for decimal degrees always use _DD suffix
            csv_file_path = os.path.join(export_dir, f"{export_name}_DD.csv")
            with open(csv_file_path, 'w', newline='') as csvfile:
                csv_writer = csv.writer(csvfile)
                csv_writer.writerow(['Line Number', 'Point Label', 'Latitude', 'Longitude'])
                # Main survey lines
                for i, line in enumerate(self.survey_lines_data):
                    if i % 2 == 0:  # Even lines - normal order
                        start, end = line[0], line[1]
                    else:  # Odd lines - flipped order
                        start, end = line[1], line[0]
                    start_label = f'L{i+1}S'
                    end_label = f'L{i+1}E'
                    csv_writer.writerow([i + 1, start_label, start[0], start[1]])
                    csv_writer.writerow([i + 1, end_label, end[0], end[1]])
                # Crossline
                if self.cross_line_data:
                    csv_writer.writerow([0, 'CLS', self.cross_line_data[0][0], self.cross_line_data[0][1]])
                    csv_writer.writerow([0, 'CLE', self.cross_line_data[1][0], self.cross_line_data[1][1]])

            # --- Export to DDM format (Decimal Minutes) ---
            ddm_file_path = os.path.join(export_dir, f"{export_name}_DM.csv")
            with open(ddm_file_path, 'w', newline='', encoding='utf-8') as ddmfile:
                ddm_writer = csv.writer(ddmfile)
                ddm_writer.writerow(['Line Number', 'Point Label', 'Latitude with Decimal Minutes', 'Longitude with Decimal Minutes'])
                # Main survey lines
                for i, line in enumerate(self.survey_lines_data):
                    if i % 2 == 0:  # Even lines - normal order
                        start, end = line[0], line[1]
                    else:  # Odd lines - flipped order
                        start, end = line[1], line[0]
                    start_label = f'L{i+1}S'
                    end_label = f'L{i+1}E'
                    start_lat_ddm = decimal_degrees_to_ddm(start[0], is_latitude=True)
                    start_lon_ddm = decimal_degrees_to_ddm(start[1], is_latitude=False)
                    end_lat_ddm = decimal_degrees_to_ddm(end[0], is_latitude=True)
                    end_lon_ddm = decimal_degrees_to_ddm(end[1], is_latitude=False)
                    ddm_writer.writerow([i + 1, start_label, start_lat_ddm, start_lon_ddm])
                    ddm_writer.writerow([i + 1, end_label, end_lat_ddm, end_lon_ddm])
                # Crossline
                if self.cross_line_data:
                    cls_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][0], is_latitude=True)
                    cls_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][1], is_latitude=False)
                    cle_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][0], is_latitude=True)
                    cle_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][1], is_latitude=False)
                    ddm_writer.writerow([0, 'CLS', cls_lat_ddm, cls_lon_ddm])
                    ddm_writer.writerow([0, 'CLE', cle_lat_ddm, cle_lon_ddm])

            # --- Export to DDM text format (Decimal Minutes) ---
            ddm_txt_file_path = os.path.join(export_dir, f"{export_name}_DM.txt")
            with open(ddm_txt_file_path, 'w', encoding='utf-8') as ddm_txt_file:
                # Main survey lines
                for i, line in enumerate(self.survey_lines_data):
                    if i % 2 == 0:  # Even lines - normal order
                        start, end = line[0], line[1]
                    else:  # Odd lines - flipped order
                        start, end = line[1], line[0]
                    start_label = f'L{i+1}S'
                    end_label = f'L{i+1}E'
                    start_lat_ddm = decimal_degrees_to_ddm(start[0], is_latitude=True)
                    start_lon_ddm = decimal_degrees_to_ddm(start[1], is_latitude=False)
                    end_lat_ddm = decimal_degrees_to_ddm(end[0], is_latitude=True)
                    end_lon_ddm = decimal_degrees_to_ddm(end[1], is_latitude=False)
                    ddm_txt_file.write(f"{start_label}, {start_lat_ddm}, {start_lon_ddm}\n")
                    ddm_txt_file.write(f"{end_label}, {end_lat_ddm}, {end_lon_ddm}\n")
                # Crossline
                if self.cross_line_data:
                    cls_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][0], is_latitude=True)
                    cls_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][1], is_latitude=False)
                    cle_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][0], is_latitude=True)
                    cle_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][1], is_latitude=False)
                    ddm_txt_file.write(f"CLS, {cls_lat_ddm}, {cls_lon_ddm}\n")
                    ddm_txt_file.write(f"CLE, {cle_lat_ddm}, {cle_lon_ddm}\n")

            # --- Export Input Parameters to TXT ---
            params_file_path = os.path.join(export_dir, f"{export_name}_params_DD.txt")
            with open(params_file_path, 'w') as f:
                f.write("Survey Plan Input Parameters:\n\n")
                for key, value in values.items():
                    f.write(f"{key.replace('_', ' ').title()}: {value}\n")

            # --- Export to ESRI Shapefile (.shp) ---
            if LineString is not None and fiona is not None:
                schema = {
                    'geometry': 'LineString',
                    'properties': {'line_num': 'int'},
                }
                crs_epsg = 'EPSG:4326'  # WGS 84

                features = []
                # Add main survey lines
                for i, line_coords in enumerate(self.survey_lines_data):
                    # Shapely expects (lon, lat) order
                    shapely_line = LineString([(p[1], p[0]) for p in line_coords])
                    features.append({
                        'geometry': shapely_line.__geo_interface__,
                        'properties': {'line_num': i + 1},
                    })

                # Add bisecting line
                if self.cross_line_data:
                    shapely_cross_line = LineString([(p[1], p[0]) for p in self.cross_line_data])
                    features.append({
                        'geometry': shapely_cross_line.__geo_interface__,
                        'properties': {'line_num': 0},  # Using 0 for crossline
                    })

                shapefile_path = os.path.join(export_dir, f"{export_name}.shp")
                with fiona.open(shapefile_path, 'w', driver='ESRI Shapefile', crs=crs_epsg, schema=schema) as collection:
                    collection.writerecords(features)
                self.set_ref_info_text(
                    f"Data exported successfully to:\n"
                    f"- {os.path.basename(csv_file_path)}\n"
                    f"- {os.path.basename(params_file_path)}\n"
                    f"- {os.path.basename(shapefile_path)} (and associated files)\n"
                    f"in directory: {export_dir}", append=False)
            else:
                self.set_ref_info_text(
                    f"Data exported successfully to:\n"
                    f"- {os.path.basename(csv_file_path)}\n"
                    f"- {os.path.basename(params_file_path)}\n"
                    f"(Shapefile export requires fiona/shapely)\n"
                    f"in directory: {export_dir}", append=False)

        except Exception as e:
            self._show_message("error","Export Error", f"Failed to export data: {e}")

    def _apply_ref_assignments(self, assignments, imported_lines):
        """Set cross_line_data and survey_lines_data from assignment dialog result."""
        self.cross_line_data = []
        self.survey_lines_data = []
        ref_lines = []
        for line_idx, role in assignments.items():
            if role == 'crossline':
                self.cross_line_data = imported_lines[line_idx]
            else:
                ref_lines.append((int(role.split('_')[1]), imported_lines[line_idx]))
        ref_lines.sort(key=lambda x: x[0])
        self.survey_lines_data = [line for _, line in ref_lines]

    def _download_and_load_gmrt_after_ref_import(self):
        """Compute extent from reference + crossline points (with buffer), then use GMRT mixin to download and load."""
        all_points = []
        for line in self.survey_lines_data:
            all_points.extend(line)
        if self.cross_line_data:
            all_points.extend(self.cross_line_data)
        if not all_points:
            self._show_message("warning", "GMRT Download", "No reference line points to compute extent.")
            return
        lats = [p[0] for p in all_points]
        lons = [p[1] for p in all_points]
        min_lat, max_lat = min(lats), max(lats)
        min_lon, max_lon = min(lons), max(lons)
        mid_lat = (min_lat + max_lat) / 2.0
        mid_lon = (min_lon + max_lon) / 2.0
        buffer_deg = 0.5
        if hasattr(self, "ref_gmrt_buffer_spin"):
            try:
                buffer_deg = float(self.ref_gmrt_buffer_spin.value())
            except (ValueError, TypeError):
                pass
        west = mid_lon - buffer_deg
        east = mid_lon + buffer_deg
        south = mid_lat - buffer_deg
        north = mid_lat + buffer_deg
        self._download_gmrt_and_load(
            west, east, south, north,
            resolution=100,
            layer="topo",
            default_filename_prefix="GMRT_Bathy",
            log_func=lambda msg, append=True: self.set_ref_info_text(msg, append=append),
            default_directory=getattr(self, "last_ref_import_dir", None),
        )

    def _compute_crossline_lead_from_import(self):
        """Compute symmetric Crossline Lead-in/out (m) from imported ref lines and crossline.
        Intersects crossline with first and last reference lines, then uses min of the two
        end-distances so the regenerated crossline does not extend past the import.
        Returns the value in meters, or None if not computable."""
        if not pyproj or not LineString:
            return None
        if len(self.survey_lines_data) < 2 or not self.cross_line_data or len(self.cross_line_data) != 2:
            return None
        try:
            geod = pyproj.Geod(ellps="WGS84")
            all_pts = []
            for line in self.survey_lines_data:
                all_pts.extend(line)
            all_pts.extend(self.cross_line_data)
            if hasattr(self, "_compute_utm_zone_from_points"):
                zone, hem = self._compute_utm_zone_from_points(all_pts)
            else:
                lats = [p[0] for p in all_pts]
                lons = [p[1] for p in all_pts]
                zone = max(1, min(60, int((sum(lons) / len(lons) + 180) / 6) + 1))
                hem = "North" if (sum(lats) / len(lats)) >= 0 else "South"
            utm_crs = f"EPSG:326{zone:02d}" if hem == "North" else f"EPSG:327{zone:02d}"
            to_utm = pyproj.Transformer.from_crs("EPSG:4326", utm_crs, always_xy=True)
            from_utm = pyproj.Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True)

            def to_utm_xy(lat, lon):
                e, n = to_utm.transform(lon, lat)
                return (e, n)

            def from_utm_to_lonlat(e, n):
                lon, lat = from_utm.transform(e, n)
                return (lat, lon)

            cross_a, cross_b = self.cross_line_data[0], self.cross_line_data[1]
            cross_utm = LineString([to_utm_xy(cross_a[0], cross_a[1]), to_utm_xy(cross_b[0], cross_b[1])])
            first_line = self.survey_lines_data[0]
            first_utm = LineString([to_utm_xy(first_line[0][0], first_line[0][1]), to_utm_xy(first_line[1][0], first_line[1][1])])
            last_line = self.survey_lines_data[-1]
            last_utm = LineString([to_utm_xy(last_line[0][0], last_line[0][1]), to_utm_xy(last_line[1][0], last_line[1][1])])

            inter_first = cross_utm.intersection(first_utm)
            inter_last = cross_utm.intersection(last_utm)
            if inter_first.is_empty or inter_last.is_empty:
                return None

            def one_point(geom):
                if geom.geom_type == "Point":
                    return (geom.x, geom.y)
                if geom.geom_type in ("MultiPoint", "LineString") and len(geom.coords) > 0:
                    return (geom.coords[0][0], geom.coords[0][1])
                return None

            pt_first_utm = one_point(inter_first)
            pt_last_utm = one_point(inter_last)
            if pt_first_utm is None or pt_last_utm is None:
                return None

            pt_first = from_utm_to_lonlat(pt_first_utm[0], pt_first_utm[1])
            pt_last = from_utm_to_lonlat(pt_last_utm[0], pt_last_utm[1])

            _, _, d_first_a = geod.inv(pt_first[1], pt_first[0], cross_a[1], cross_a[0])
            _, _, d_first_b = geod.inv(pt_first[1], pt_first[0], cross_b[1], cross_b[0])
            lead_in = min(d_first_a, d_first_b)

            _, _, d_last_a = geod.inv(pt_last[1], pt_last[0], cross_a[1], cross_a[0])
            _, _, d_last_b = geod.inv(pt_last[1], pt_last[0], cross_b[1], cross_b[0])
            lead_out = min(d_last_a, d_last_b)

            return min(lead_in, lead_out)
        except Exception:
            return None

    def _ref_import_post_import(self, file_path):
        """Shared post-import: validate, load params, populate fields, plot, message."""
        if not self.survey_lines_data and not self.cross_line_data:
            self._show_message("warning", "Import Warning", "No valid survey lines found in the selected file.")
            return
        # Set Crossline Lead-in/out from imported geometry (symmetric) so regenerated crossline matches import
        bisect_lead = self._compute_crossline_lead_from_import()
        if bisect_lead is not None and hasattr(self, "bisect_lead_entry"):
            self.bisect_lead_entry.blockSignals(True)
            self.bisect_lead_entry.setText(f"{bisect_lead:.1f}")
            self.bisect_lead_entry.blockSignals(False)
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        dir_name = os.path.dirname(file_path)
        metadata_path = os.path.join(dir_name, f"{base_name}_params.json")
        params = None
        if os.path.exists(metadata_path):
            try:
                with open(metadata_path, 'r', encoding='utf-8') as f:
                    params = json.load(f)
            except Exception as e:
                print(f"Warning: Could not load metadata file: {e}")
        try:
            if params:
                if params.get('central_lat') is not None:
                    self.central_lat_entry.clear()
                    self.central_lat_entry.setText(str(params['central_lat']))
                    if hasattr(self, 'pick_center_btn'):
                        self.pick_center_btn.setStyleSheet("")
                if params.get('central_lon') is not None:
                    self.central_lon_entry.clear()
                    self.central_lon_entry.setText(str(params['central_lon']))
                if params.get('line_length') is not None:
                    self.line_length_entry.setText(str(params['line_length']))
                if params.get('heading') is not None:
                    self.heading_entry.setText(str(params['heading']))
                if params.get('dist_between_lines') is not None:
                    self.dist_between_lines_entry.setText(str(params['dist_between_lines']))
                if params.get('num_lines') is not None:
                    self.num_lines_entry.setText(str(params['num_lines']))
                if params.get('bisect_lead') is not None:
                    self.bisect_lead_entry.setText(str(params['bisect_lead']))
                if params.get('survey_speed') is not None:
                    self.survey_speed_entry.setText(str(params['survey_speed']))
                if params.get('crossline_passes') is not None:
                    self.crossline_passes_entry.setText(str(params['crossline_passes']))
                if params.get('export_name'):
                    self.export_name_entry.clear()
                    self.export_name_entry.setText(params['export_name'])
                if params.get('offset_direction') and hasattr(self, 'offset_direction_combo'):
                    idx = self.offset_direction_combo.findText(params['offset_direction'])
                    if idx >= 0:
                        self.offset_direction_combo.setCurrentIndex(idx)
                    self.offset_direction_var = params['offset_direction']
                if params.get('line_length_multiplier') is not None and hasattr(self, '_update_multiplier_label_len'):
                    self.line_length_multiplier.set(params['line_length_multiplier'])
                    self._update_multiplier_label_len(params['line_length_multiplier'])
                if params.get('dist_between_lines_multiplier') is not None and hasattr(self, '_update_multiplier_label_dist'):
                    self.dist_between_lines_multiplier.set(params['dist_between_lines_multiplier'])
                    self._update_multiplier_label_dist(params['dist_between_lines_multiplier'])
            else:
                if pyproj is not None:
                    geod = pyproj.Geod(ellps="WGS84")
                    all_points = []
                    for line in self.survey_lines_data:
                        all_points.extend(line)
                    if self.cross_line_data:
                        all_points.extend(self.cross_line_data)
                    if all_points:
                        all_lats = [p[0] for p in all_points]
                        all_lons = [p[1] for p in all_points]
                        central_lat = (min(all_lats) + max(all_lats)) / 2.0
                        central_lon = (min(all_lons) + max(all_lons)) / 2.0
                        self.central_lat_entry.clear()
                        self.central_lat_entry.setText(f"{central_lat:.6f}")
                        if hasattr(self, 'pick_center_btn'):
                            self.pick_center_btn.setStyleSheet("")
                        self.central_lon_entry.clear()
                        self.central_lon_entry.setText(f"{central_lon:.6f}")
                    if len(self.survey_lines_data) > 0:
                        first_line = self.survey_lines_data[0]
                        try:
                            lat1, lon1 = first_line[0]
                            lat2, lon2 = first_line[1]
                            fwd_az, back_az, dist = geod.inv(lon1, lat1, lon2, lat2)
                            heading = fwd_az % 360
                            self.heading_entry.clear()
                            self.heading_entry.setText(f"{heading:.1f}")
                            self.line_length_entry.clear()
                            self.line_length_entry.setText(f"{dist:.1f}")
                        except Exception:
                            pass
                    if self.survey_lines_data:
                        self.num_lines_entry.clear()
                        self.num_lines_entry.setText(str(len(self.survey_lines_data)))
                    if len(self.survey_lines_data) > 1:
                        try:
                            line1_mid = ((self.survey_lines_data[0][0][0] + self.survey_lines_data[0][1][0]) / 2,
                                       (self.survey_lines_data[0][0][1] + self.survey_lines_data[0][1][1]) / 2)
                            line2_mid = ((self.survey_lines_data[1][0][0] + self.survey_lines_data[1][1][0]) / 2,
                                       (self.survey_lines_data[1][0][1] + self.survey_lines_data[1][1][1]) / 2)
                            _, _, dist = geod.inv(line1_mid[1], line1_mid[0], line2_mid[1], line2_mid[0])
                            self.dist_between_lines_entry.clear()
                            self.dist_between_lines_entry.setText(f"{dist:.1f}")
                        except Exception:
                            pass
                if not self.export_name_entry.text().strip():
                    self.export_name_entry.clear()
                    self.export_name_entry.setText(base_name)
        except Exception as e:
            print(f"Warning: Error populating parameter fields: {e}")
        self._plot_survey_plan(preserve_view_limits=True)
        imported_items = []
        if self.survey_lines_data:
            imported_items.append(f"{len(self.survey_lines_data)} survey line(s)")
        if self.cross_line_data:
            imported_items.append("crossline")
        if imported_items:
            msg = f"Successfully imported: {', '.join(imported_items)}"
            if params:
                msg += " (parameters loaded from metadata)"
            else:
                msg += " (parameters calculated from lines)"
            self.set_ref_info_text(msg)
            if getattr(self, "ref_download_gmrt_checkbox", None) and self.ref_download_gmrt_checkbox.isChecked():
                self._download_and_load_gmrt_after_ref_import()
        else:
            self._show_message("warning", "Import Warning", "No valid survey lines found in the selected file.")

    def _import_survey_files(self):
        """Import reference planning survey lines from CSV, GeoJSON, or DDD/DMS/DMM/LNW text files."""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Survey File to Import",
            self.last_ref_import_dir,
            "Known Survey Files (*_DMS.txt *_DMM.txt *_DDD.txt *_DD.csv *.csv *.geojson *.json *.lnw);;"
            "Hypack LNW files (*.lnw);;Degrees Minutes Seconds (*_DMS.txt);;Degrees Decimal Minutes (*_DMM.txt);;"
            "Decimal Degrees (*_DDD.txt);;Decimal Degree CSV (*_DD.csv);;CSV (*.csv);;GeoJSON (*.geojson);;JSON (*.json)"
        )
        if not file_path:
            return
        import_dir = os.path.dirname(file_path)
        if import_dir and os.path.isdir(import_dir):
            self.last_ref_import_dir = import_dir
            self._save_last_ref_import_dir()

        try:
            self.survey_lines_data = []
            self.cross_line_data = []
            file_ext = os.path.splitext(file_path)[1].lower()
            file_basename = os.path.basename(file_path)
            file_processed = False

            # DDD / DMS / DMM / LNW: use calibration parsers + reference assignment dialog
            if file_ext == '.lnw':
                if UTMZoneDialog is None:
                    self._show_message("error", "Import Error", "UTM zone dialog not available.")
                    return
                utm_dialog = UTMZoneDialog(self)
                if utm_dialog.exec() != QDialog.DialogCode.Accepted:
                    return
                utm_zone, hemisphere = utm_dialog.get_utm_info()
                imported_lines = self._parse_lnw_file(file_path, utm_zone, hemisphere)
                if imported_lines is None:
                    return
                dialog = ReferenceLineAssignmentDialog(self, imported_lines)
                if dialog.exec() != QDialog.DialogCode.Accepted:
                    return
                self._apply_ref_assignments(dialog.get_assignments(), imported_lines)
                file_processed = True

            if not file_processed and file_basename.lower().endswith('_dms.txt'):
                imported_lines = self._parse_dms_txt_file(file_path)
                if imported_lines is None:
                    return
                dialog = ReferenceLineAssignmentDialog(self, imported_lines)
                if dialog.exec() != QDialog.DialogCode.Accepted:
                    return
                self._apply_ref_assignments(dialog.get_assignments(), imported_lines)
                file_processed = True

            if not file_processed and file_basename.lower().endswith('_dmm.txt'):
                imported_lines = self._parse_dmm_txt_file(file_path)
                if imported_lines is None:
                    return
                dialog = ReferenceLineAssignmentDialog(self, imported_lines)
                if dialog.exec() != QDialog.DialogCode.Accepted:
                    return
                self._apply_ref_assignments(dialog.get_assignments(), imported_lines)
                file_processed = True

            if not file_processed and file_basename.lower().endswith('_ddd.txt'):
                imported_lines = self._parse_ddd_txt_file(file_path)
                if imported_lines is None:
                    return
                dialog = ReferenceLineAssignmentDialog(self, imported_lines)
                if dialog.exec() != QDialog.DialogCode.Accepted:
                    return
                self._apply_ref_assignments(dialog.get_assignments(), imported_lines)
                file_processed = True

            if file_processed:
                # Run same post-import logic as CSV/GeoJSON path
                self._ref_import_post_import(file_path)
                return

            if file_ext == '.csv':
                # Import from CSV
                with open(file_path, 'r', encoding='utf-8') as csvfile:
                    csv_reader = csv.DictReader(csvfile)

                    # Dictionary to collect points for each line
                    lines_data = {}
                    crossline_points = {}

                    for row in csv_reader:
                        try:
                            line_num = int(row.get('Line Number', -1))
                            point_label = row.get('Point Label', '').strip()
                            lat = float(row.get('Latitude', 0))
                            lon = float(row.get('Longitude', 0))
                        except (ValueError, TypeError):
                            continue

                        # Handle crossline (line_num == 0)
                        if line_num == 0:
                            if point_label == 'CLS':
                                crossline_points['start'] = (lat, lon)
                            elif point_label == 'CLE':
                                crossline_points['end'] = (lat, lon)
                        else:
                            # Main survey line
                            if line_num not in lines_data:
                                lines_data[line_num] = {}
                            lines_data[line_num][point_label] = (lat, lon)

                    # Process main survey lines
                    for line_num in sorted(lines_data.keys()):
                        points = lines_data[line_num]
                        # Look for L{n}S and L{n}E labels
                        start_label = f'L{line_num}S'
                        end_label = f'L{line_num}E'
                        if start_label in points and end_label in points:
                            self.survey_lines_data.append([points[start_label], points[end_label]])

                    # Process crossline
                    if 'start' in crossline_points and 'end' in crossline_points:
                        self.cross_line_data = [crossline_points['start'], crossline_points['end']]

            elif file_ext in ['.geojson', '.json']:
                # Import from GeoJSON
                with open(file_path, 'r', encoding='utf-8') as f:
                    geojson_data = json.load(f)

                # Ensure geojson_data is a dictionary
                if not isinstance(geojson_data, dict):
                    self._show_message("error","Import Error", "Invalid GeoJSON format: root must be an object")
                    return

                # Handle FeatureCollection
                if geojson_data.get('type') == 'FeatureCollection':
                    features = geojson_data.get('features', [])
                    if not isinstance(features, list):
                        features = []
                elif geojson_data.get('type') == 'Feature':
                    features = [geojson_data]
                else:
                    features = []

                # Sort features by line_num to maintain order
                features_with_num = []
                for feature in features:
                    if not isinstance(feature, dict) or feature.get('type') != 'Feature':
                        continue
                    geometry = feature.get('geometry', {})
                    if not isinstance(geometry, dict) or geometry.get('type') != 'LineString':
                        continue
                    properties = feature.get('properties', {})
                    if not isinstance(properties, dict):
                        properties = {}
                    line_num = properties.get('line_num', 0)
                    features_with_num.append((line_num, feature))

                # Sort by line_num
                features_with_num.sort(key=lambda x: x[0])

                for line_num, feature in features_with_num:
                    if not isinstance(feature, dict):
                        continue
                    geometry = feature.get('geometry', {})
                    if not isinstance(geometry, dict):
                        continue
                    coordinates = geometry.get('coordinates', [])
                    if not isinstance(coordinates, list) or len(coordinates) < 2:
                        continue

                    # GeoJSON uses [lon, lat] format, we need [lat, lon]
                    point1 = (coordinates[0][1], coordinates[0][0])
                    point2 = (coordinates[-1][1], coordinates[-1][0])

                    if line_num == 0:
                        # Crossline
                        self.cross_line_data = [point1, point2]
                    else:
                        # Main survey line
                        self.survey_lines_data.append([point1, point2])

            else:
                self._show_message("error","Import Error", f"Unsupported file format: {file_ext}")
                return

            if not self.survey_lines_data and not self.cross_line_data:
                self._show_message("warning","Import Warning", "No valid survey lines found in the selected file.")
                return

            # Try to load parameters metadata JSON file
            base_name = os.path.splitext(os.path.basename(file_path))[0]
            dir_name = os.path.dirname(file_path)
            metadata_path = os.path.join(dir_name, f"{base_name}_params.json")
            params = None

            if os.path.exists(metadata_path):
                try:
                    with open(metadata_path, 'r', encoding='utf-8') as f:
                        params = json.load(f)
                except Exception as e:
                    print(f"Warning: Could not load metadata file: {e}")

            # Populate parameter fields from metadata or calculate from lines
            try:
                if params:
                    # Use metadata if available
                    if params.get('central_lat') is not None:
                        self.central_lat_entry.clear()
                        self.central_lat_entry.setText(str(params['central_lat']))
                        if hasattr(self, 'pick_center_btn'):
                            self.pick_center_btn.setStyleSheet("")

                    if params.get('central_lon') is not None:
                        self.central_lon_entry.clear()
                        self.central_lon_entry.setText(str(params['central_lon']))

                    if params.get('line_length') is not None:
                        self.line_length_entry.setText(str(params['line_length']))

                    if params.get('heading') is not None:
                        self.heading_entry.setText(str(params['heading']))

                    if params.get('dist_between_lines') is not None:
                        self.dist_between_lines_entry.setText(str(params['dist_between_lines']))

                    if params.get('num_lines') is not None:
                        self.num_lines_entry.setText(str(params['num_lines']))

                    if params.get('bisect_lead') is not None:
                        self.bisect_lead_entry.setText(str(params['bisect_lead']))

                    if params.get('survey_speed') is not None:
                        self.survey_speed_entry.setText(str(params['survey_speed']))

                    if params.get('crossline_passes') is not None:
                        self.crossline_passes_entry.setText(str(params['crossline_passes']))

                    if params.get('export_name'):
                        self.export_name_entry.clear()
                        self.export_name_entry.setText(params['export_name'])

                    if params.get('offset_direction'):
                        direction = params['offset_direction']
                        if hasattr(self, 'offset_direction_combo'):
                            idx = self.offset_direction_combo.findText(direction)
                            if idx >= 0:
                                self.offset_direction_combo.setCurrentIndex(idx)
                        self.offset_direction_var = direction

                    if params.get('line_length_multiplier') is not None:
                        self.line_length_multiplier.set(params['line_length_multiplier'])
                        self._update_multiplier_label_len(params['line_length_multiplier'])

                    if params.get('dist_between_lines_multiplier') is not None:
                        self.dist_between_lines_multiplier.set(params['dist_between_lines_multiplier'])
                        self._update_multiplier_label_dist(params['dist_between_lines_multiplier'])
                else:
                    # Calculate parameters from imported lines
                    if pyproj is not None:
                        geod = pyproj.Geod(ellps="WGS84")

                        # Calculate central lat/lon from all points
                        all_points = []
                        for line in self.survey_lines_data:
                            all_points.extend(line)
                        if self.cross_line_data:
                            all_points.extend(self.cross_line_data)

                        if all_points:
                            all_lats = [p[0] for p in all_points]
                            all_lons = [p[1] for p in all_points]
                            central_lat = (min(all_lats) + max(all_lats)) / 2.0
                            central_lon = (min(all_lons) + max(all_lons)) / 2.0

                            self.central_lat_entry.clear()
                            self.central_lat_entry.setText(f"{central_lat:.6f}")
                            if hasattr(self, 'pick_center_btn'):
                                self.pick_center_btn.setStyleSheet("")

                            self.central_lon_entry.clear()
                            self.central_lon_entry.setText(f"{central_lon:.6f}")

                        # Calculate heading from first line
                        if len(self.survey_lines_data) > 0:
                            first_line = self.survey_lines_data[0]
                            try:
                                lat1, lon1 = first_line[0]
                                lat2, lon2 = first_line[1]
                                fwd_az, back_az, dist = geod.inv(lon1, lat1, lon2, lat2)
                                heading = fwd_az % 360

                                self.heading_entry.clear()
                                self.heading_entry.setText(f"{heading:.1f}")

                                line_length = dist
                                self.line_length_entry.clear()
                                self.line_length_entry.setText(f"{line_length:.1f}")
                            except Exception:
                                pass

                        # Count number of lines
                        if self.survey_lines_data:
                            self.num_lines_entry.clear()
                            self.num_lines_entry.setText(str(len(self.survey_lines_data)))

                        # Calculate distance between lines (approximate)
                        if len(self.survey_lines_data) > 1:
                            try:
                                line1_mid = ((self.survey_lines_data[0][0][0] + self.survey_lines_data[0][1][0]) / 2,
                                           (self.survey_lines_data[0][0][1] + self.survey_lines_data[0][1][1]) / 2)
                                line2_mid = ((self.survey_lines_data[1][0][0] + self.survey_lines_data[1][1][0]) / 2,
                                           (self.survey_lines_data[1][0][1] + self.survey_lines_data[1][1][1]) / 2)
                                _, _, dist = geod.inv(line1_mid[1], line1_mid[0], line2_mid[1], line2_mid[0])
                                self.dist_between_lines_entry.clear()
                                self.dist_between_lines_entry.setText(f"{dist:.1f}")
                            except Exception:
                                pass

                    # Update export name if not set
                    if not self.export_name_entry.text().strip():
                        self.export_name_entry.clear()
                        self.export_name_entry.setText(base_name)
            except Exception as e:
                print(f"Warning: Error populating parameter fields: {e}")

            # Update plot
            self._plot_survey_plan(preserve_view_limits=True)

            # Show success message
            imported_items = []
            if self.survey_lines_data:
                imported_items.append(f"{len(self.survey_lines_data)} survey line(s)")
            if self.cross_line_data:
                imported_items.append("crossline")

            if imported_items:
                msg = f"Successfully imported: {', '.join(imported_items)}"
                if params:
                    msg += " (parameters loaded from metadata)"
                else:
                    msg += " (parameters calculated from lines)"
                self.set_ref_info_text(msg)
                if getattr(self, "ref_download_gmrt_checkbox", None) and self.ref_download_gmrt_checkbox.isChecked():
                    self._download_and_load_gmrt_after_ref_import()
            else:
                self._show_message("warning","Import Warning", "No valid survey lines found in the selected file.")

        except Exception as e:
            self._show_message("error","Import Error", f"Failed to import survey files: {e}")
            import traceback
            traceback.print_exc()

    def set_ref_info_text(self, message, append=False):
        """Add a message to the Activity Log with [Reference] prefix. Maintains up to 200 lines of history."""
        if not hasattr(self, 'activity_log_text'):
            return
        prefixed_message = f"[Reference] {message}"
        self.activity_log_text.setReadOnly(False)
        if append:
            # Prepend new message (newest at top) - only when explicitly requested
            current_text = self.activity_log_text.toPlainText().rstrip('\n')
            new_text = prefixed_message + "\n" + current_text if current_text else prefixed_message + "\n"
            # Split into lines and limit to 200
            lines = new_text.splitlines()
            if len(lines) > 200:
                lines = lines[:200]
            self.activity_log_text.setPlainText("\n".join(lines) + "\n")
            # Move cursor to top
            cursor = self.activity_log_text.textCursor()
            cursor.movePosition(QTextCursor.MoveOperation.Start)
            self.activity_log_text.setTextCursor(cursor)
        else:
            # Append to end (default behavior - newest at bottom)
            self.activity_log_text.append(prefixed_message)
            # Maintain up to 200 lines by removing oldest from top
            num_lines = self.activity_log_text.document().blockCount()
            if num_lines > 200:
                cursor = self.activity_log_text.textCursor()
                cursor.movePosition(QTextCursor.MoveOperation.Start)
                cursor.movePosition(QTextCursor.MoveOperation.Down, QTextCursor.MoveMode.MoveAnchor, num_lines - 200)
                cursor.movePosition(QTextCursor.MoveOperation.Start, QTextCursor.MoveMode.KeepAnchor)
                cursor.removeSelectedText()
            # Move cursor to end and scroll to show newest message
            cursor = self.activity_log_text.textCursor()
            cursor.movePosition(QTextCursor.MoveOperation.End)
            self.activity_log_text.setTextCursor(cursor)
        self.activity_log_text.ensureCursorVisible()
        self.activity_log_text.setReadOnly(True)

    def _reset_reference_tab(self):
        """Reset Reference tab data and entry fields to defaults."""
        # Clear line data
        self.survey_lines_data = []
        self.cross_line_data = []
        self.central_point_coords = (None, None)

        # Reset entry fields to defaults
        if hasattr(self, 'central_lat_entry'):
            self.central_lat_entry.clear()
        if hasattr(self, 'central_lon_entry'):
            self.central_lon_entry.clear()
        if hasattr(self, 'line_length_entry'):
            self.line_length_entry.clear()
        if hasattr(self, 'heading_entry'):
            self.heading_entry.setText("0")
        if hasattr(self, 'dist_between_lines_entry'):
            self.dist_between_lines_entry.clear()
        if hasattr(self, 'num_lines_entry'):
            self.num_lines_entry.setText("5")
        if hasattr(self, 'bisect_lead_entry'):
            self.bisect_lead_entry.setText("100")
        if hasattr(self, 'survey_speed_entry'):
            self.survey_speed_entry.setText("8")
        if hasattr(self, 'crossline_passes_entry'):
            self.crossline_passes_entry.setText("2")
        if hasattr(self, 'export_name_entry'):
            self.export_name_entry.clear()
            try:
                dist = int(float(self.dist_between_lines_entry.text())) if self.dist_between_lines_entry.text() else 0
                heading = int(float(self.heading_entry.text())) if self.heading_entry.text() else 0
                export_name = f"Reference_{dist}m_{heading}deg"
            except Exception:
                export_name = "Reference_0m_0deg"
            self.export_name_entry.setText(export_name)
        if hasattr(self, 'offset_direction_var'):
            self.offset_direction_var = "North"
            if hasattr(self, 'offset_direction_combo'):
                self.offset_direction_combo.setCurrentIndex(0)
        if hasattr(self, 'line_length_multiplier'):
            if hasattr(self.line_length_multiplier, 'set'):
                self.line_length_multiplier.set(8.0)
            else:
                self.line_length_multiplier = 8.0
            if hasattr(self, 'multiplier_slider_len'):
                self.multiplier_slider_len.setValue(int(8.0 * 10))
            if hasattr(self, 'multiplier_label_len'):
                self.multiplier_label_len.setText("8.0")
        if hasattr(self, 'dist_between_lines_multiplier'):
            if hasattr(self.dist_between_lines_multiplier, 'set'):
                self.dist_between_lines_multiplier.set(1.0)
            else:
                self.dist_between_lines_multiplier = 1.0
            if hasattr(self, 'multiplier_slider_dist'):
                self.multiplier_slider_dist.setValue(int(1.0 * 10))
            if hasattr(self, 'multiplier_label_dist'):
                self.multiplier_label_dist.setText("1.0")

        # Redraw plot (without clearing GeoTIFF)
        self._plot_survey_plan(preserve_view_limits=True)

    def _calculate_total_survey_time(self):
        """Calculate total survey time including travel between lines and crossline."""
        try:
            geod = pyproj.Geod(ellps="WGS84")
            speed_knots = float(self.survey_speed_entry.text()) if self.survey_speed_entry.text() else 8.0
            speed_m_per_h = speed_knots * 1852
            turn_time_min = float(self.ref_turn_time_entry.text()) if hasattr(self, 'ref_turn_time_entry') and self.ref_turn_time_entry.text() else 5.0

            # Initialize distance calculations
            main_lines_total_distance = 0
            travel_between_lines_total_distance = 0
            travel_to_crossline_distance = 0
            crossline_single_pass_distance = 0

            # Calculate main lines survey time and distance (all main lines)
            main_lines_survey_time_minutes = 0
            for line in self.survey_lines_data:
                lat1, lon1 = line[0]
                lat2, lon2 = line[1]
                _, _, line_length = geod.inv(lon1, lat1, lon2, lat2)
                main_lines_total_distance += line_length
                line_time_hours = line_length / speed_m_per_h if speed_m_per_h > 0 else 0
                main_lines_survey_time_minutes += line_time_hours * 60

            # Calculate travel time and distance between main lines (zigzag pattern)
            travel_between_lines_minutes = 0
            if len(self.survey_lines_data) > 1:
                for i in range(1, len(self.survey_lines_data)):
                    if (i-1) % 2 == 0:
                        end_lat, end_lon = self.survey_lines_data[i-1][1]
                    else:
                        end_lat, end_lon = self.survey_lines_data[i-1][0]

                    if i % 2 == 0:
                        start_lat, start_lon = self.survey_lines_data[i][0]
                    else:
                        start_lat, start_lon = self.survey_lines_data[i][1]

                    _, _, travel_distance = geod.inv(end_lon, end_lat, start_lon, start_lat)
                    travel_between_lines_total_distance += travel_distance
                    travel_time_hours = travel_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                    # Travel time + turn time between lines
                    travel_between_lines_minutes += travel_time_hours * 60 + turn_time_min

            # Calculate travel time and distance from last main line to crossline
            travel_to_crossline_minutes = 0
            if self.cross_line_data and self.survey_lines_data:
                last_line_index = len(self.survey_lines_data) - 1
                if last_line_index % 2 == 0:
                    last_end_lat, last_end_lon = self.survey_lines_data[last_line_index][1]
                else:
                    last_end_lat, last_end_lon = self.survey_lines_data[last_line_index][0]

                crossline_start_lat, crossline_start_lon = self.cross_line_data[0]
                crossline_end_lat, crossline_end_lon = self.cross_line_data[1]

                _, _, dist_to_crossline_start = geod.inv(last_end_lon, last_end_lat, crossline_start_lon, crossline_start_lat)
                _, _, dist_to_crossline_end = geod.inv(last_end_lon, last_end_lat, crossline_end_lon, crossline_end_lat)

                if dist_to_crossline_start <= dist_to_crossline_end:
                    travel_distance = dist_to_crossline_start
                else:
                    travel_distance = dist_to_crossline_end

                travel_to_crossline_distance = travel_distance
                travel_time_hours = travel_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                # Travel time + turn time to crossline
                travel_to_crossline_minutes = travel_time_hours * 60 + turn_time_min

            # Calculate crossline survey time and distance
            crossline_survey_minutes = 0
            num_passes = 1
            if self.cross_line_data:
                lat1, lon1 = self.cross_line_data[0]
                lat2, lon2 = self.cross_line_data[1]
                _, _, crossline_length = geod.inv(lon1, lat1, lon2, lat2)
                crossline_single_pass_distance = crossline_length
                crossline_time_hours = crossline_length / speed_m_per_h if speed_m_per_h > 0 else 0
                crossline_survey_minutes = crossline_time_hours * 60

                try:
                    num_passes = int(self.crossline_passes_entry.text()) if self.crossline_passes_entry.text() else 2
                    crossline_survey_minutes *= num_passes
                    # Add turn times between crossline passes (num_passes - 1 turns)
                    if num_passes > 1:
                        crossline_survey_minutes += turn_time_min * (num_passes - 1)
                except (ValueError, AttributeError):
                    num_passes = 2
                    crossline_survey_minutes *= num_passes
                    if num_passes > 1:
                        crossline_survey_minutes += turn_time_min * (num_passes - 1)

            # Calculate turn times
            num_turns_between_main_lines = len(self.survey_lines_data) - 1 if len(self.survey_lines_data) > 1 else 0
            num_turns_to_crossline = 1 if (self.cross_line_data and self.survey_lines_data) else 0
            num_turns_between_crossline_passes = (num_passes - 1) if (self.cross_line_data and num_passes > 1) else 0
            total_turn_time_min = turn_time_min * (num_turns_between_main_lines + num_turns_to_crossline + num_turns_between_crossline_passes)

            # Calculate total time and distance
            total_minutes = (main_lines_survey_time_minutes +
                           travel_between_lines_minutes +
                           travel_to_crossline_minutes +
                           crossline_survey_minutes)
            total_hours = total_minutes / 60

            total_distance = (main_lines_total_distance +
                            travel_between_lines_total_distance +
                            travel_to_crossline_distance +
                            (crossline_single_pass_distance * num_passes))

            return {
                'total_minutes': total_minutes,
                'total_hours': total_hours,
                'main_lines_minutes': main_lines_survey_time_minutes,
                'crossline_minutes': crossline_survey_minutes,
                'travel_minutes': travel_between_lines_minutes,
                'travel_to_crossline_minutes': travel_to_crossline_minutes,
                'total_turn_time_min': total_turn_time_min,
                'num_turns_between_main_lines': num_turns_between_main_lines,
                'num_turns_to_crossline': num_turns_to_crossline,
                'num_turns_between_crossline_passes': num_turns_between_crossline_passes,
                'turn_time_per_turn_min': turn_time_min,
                'total_distance_m': total_distance,
                'total_distance_km': total_distance / 1000.0,
                'total_distance_nm': total_distance / 1852.0,
                'main_lines_distance_m': main_lines_total_distance,
                'main_lines_distance_km': main_lines_total_distance / 1000.0,
                'main_lines_distance_nm': main_lines_total_distance / 1852.0,
                'travel_between_lines_distance_m': travel_between_lines_total_distance,
                'travel_between_lines_distance_km': travel_between_lines_total_distance / 1000.0,
                'travel_between_lines_distance_nm': travel_between_lines_total_distance / 1852.0,
                'travel_to_crossline_distance_m': travel_to_crossline_distance,
                'travel_to_crossline_distance_km': travel_to_crossline_distance / 1000.0,
                'travel_to_crossline_distance_nm': travel_to_crossline_distance / 1852.0,
                'crossline_single_pass_distance_m': crossline_single_pass_distance,
                'crossline_single_pass_distance_km': crossline_single_pass_distance / 1000.0,
                'crossline_single_pass_distance_nm': crossline_single_pass_distance / 1852.0,
                'crossline_total_distance_m': crossline_single_pass_distance * num_passes,
                'crossline_total_distance_km': (crossline_single_pass_distance * num_passes) / 1000.0,
                'crossline_total_distance_nm': (crossline_single_pass_distance * num_passes) / 1852.0,
                'num_crossline_passes': num_passes
            }
        except Exception as e:
            print(f"Error calculating total survey time: {e}")
            return {
                'total_minutes': 0,
                'total_hours': 0,
                'main_lines_minutes': 0,
                'crossline_minutes': 0,
                'travel_minutes': 0,
                'travel_to_crossline_minutes': 0,
                'total_turn_time_min': 0,
                'num_turns_between_main_lines': 0,
                'num_turns_to_crossline': 0,
                'num_turns_between_crossline_passes': 0,
                'turn_time_per_turn_min': 5.0,
                'total_distance_m': 0,
                'total_distance_km': 0,
                'total_distance_nm': 0,
                'main_lines_distance_m': 0,
                'main_lines_distance_km': 0,
                'main_lines_distance_nm': 0,
                'travel_between_lines_distance_m': 0,
                'travel_between_lines_distance_km': 0,
                'travel_between_lines_distance_nm': 0,
                'travel_to_crossline_distance_m': 0,
                'travel_to_crossline_distance_km': 0,
                'travel_to_crossline_distance_nm': 0,
                'crossline_single_pass_distance_m': 0,
                'crossline_single_pass_distance_km': 0,
                'crossline_single_pass_distance_nm': 0,
                'crossline_total_distance_m': 0,
                'crossline_total_distance_km': 0,
                'crossline_total_distance_nm': 0,
                'num_crossline_passes': 1
            }

    def _calculate_reference_survey_statistics(self):
        """Calculate comprehensive reference planning survey statistics using the existing calculation function."""
        try:
            stats = self._calculate_total_survey_time()
            if not stats:
                return None

            stats['num_main_lines'] = len(self.survey_lines_data)

            return stats

        except Exception as e:
            print(f"Error calculating reference survey statistics: {e}")
            return None

    def _show_reference_planning_info(self):
        """Display comprehensive reference planning survey statistics in a custom dialog with copy functionality."""
        stats = self._calculate_reference_survey_statistics()
        if not stats:
            self._show_message("warning","Statistics Error", "Unable to calculate reference planning statistics. Please ensure survey lines are generated.")
            return

        # Format the statistics for display
        stats_text = "COMPREHENSIVE REFERENCE PLANNING STATISTICS\n"
        stats_text += "=" * 50 + "\n\n"

        # Survey parameters
        try:
            speed_knots = float(self.survey_speed_entry.text()) if self.survey_speed_entry.text() else 8.0
            stats_text += f"Survey Speed: {speed_knots} knots\n"
        except Exception:
            stats_text += "Survey Speed: 8.0 knots (default)\n"
        try:
            turn_time_min = float(self.ref_turn_time_entry.text()) if hasattr(self, 'ref_turn_time_entry') and self.ref_turn_time_entry.text() else 5.0
            stats_text += f"Turn Time (per turn): {turn_time_min} min\n\n"
        except Exception:
            turn_time_min = 5.0
            stats_text += f"Turn Time (per turn): {turn_time_min} min (default)\n\n"

        # Main survey lines
        stats_text += "MAIN SURVEY LINES\n"
        stats_text += "-" * 20 + "\n"
        stats_text += f"Number of Main Lines: {stats['num_main_lines']}\n"

        # Add heading information for the first main line
        if self.survey_lines_data and len(self.survey_lines_data) > 0:
            try:
                if pyproj is not None:
                    geod = pyproj.Geod(ellps="WGS84")
                    first_line = self.survey_lines_data[0]
                    lat1, lon1 = first_line[0]
                    lat2, lon2 = first_line[1]
                    fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)

                    heading = fwd_az % 360
                    reciprocal_heading = back_az % 360

                    stats_text += f"Heading: {heading:.1f}°\n"
                    stats_text += f"Reciprocal Heading: {reciprocal_heading:.1f}°\n"
                else:
                    stats_text += f"Heading: pyproj not available\n"
            except Exception:
                stats_text += f"Heading: Unable to calculate\n"

        if stats['num_main_lines'] > 0:
            single_line_length_m = stats['main_lines_distance_m'] / stats['num_main_lines']
            single_line_length_km = single_line_length_m / 1000
            single_line_length_nm = single_line_length_m / 1852
            stats_text += f"Length of Main Line: {single_line_length_m:.1f} m ({single_line_length_km:.3f} km, {single_line_length_nm:.3f} nm)\n"
        stats_text += f"Total Distance: {stats['main_lines_distance_m']:.1f} m ({stats['main_lines_distance_km']:.3f} km, {stats['main_lines_distance_nm']:.3f} nm)\n"
        stats_text += f"Survey Time: {stats['main_lines_minutes']:.1f} min\n\n"

        # Travel between lines
        if stats['travel_between_lines_distance_m'] > 0:
            stats_text += "TRAVEL BETWEEN LINES\n"
            stats_text += "-" * 25 + "\n"
            stats_text += f"Total Travel Distance: {stats['travel_between_lines_distance_m']:.1f} m ({stats['travel_between_lines_distance_km']:.3f} km, {stats['travel_between_lines_distance_nm']:.3f} nm)\n"
            stats_text += f"Travel Time: {stats['travel_minutes']:.1f} min\n\n"

        # Crossline information
        if stats['crossline_single_pass_distance_m'] > 0:
            stats_text += "CROSSLINE\n"
            stats_text += "-" * 10 + "\n"
            stats_text += f"Crossline Passes: {stats['num_crossline_passes']}\n"

            if self.cross_line_data and len(self.cross_line_data) >= 2:
                try:
                    if pyproj is not None:
                        geod = pyproj.Geod(ellps="WGS84")
                        lat1, lon1 = self.cross_line_data[0]
                        lat2, lon2 = self.cross_line_data[1]
                        fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)

                        crossline_heading = fwd_az % 360
                        crossline_reciprocal_heading = back_az % 360

                        stats_text += f"Crossline Heading: {crossline_heading:.1f}°\n"
                        stats_text += f"Crossline Reciprocal Heading: {crossline_reciprocal_heading:.1f}°\n"
                    else:
                        stats_text += f"Crossline Heading: pyproj not available\n"
                except Exception:
                    stats_text += f"Crossline Heading: Unable to calculate\n"

            stats_text += f"Crossline (per pass): {stats['crossline_single_pass_distance_m']:.1f} m ({stats['crossline_single_pass_distance_km']:.3f} km, {stats['crossline_single_pass_distance_nm']:.3f} nm)\n"
            stats_text += f"Crossline Total Distance: {stats['crossline_total_distance_m']:.1f} m ({stats['crossline_total_distance_km']:.3f} km, {stats['crossline_total_distance_nm']:.3f} nm)\n"
            stats_text += f"Crossline Time: {stats['crossline_minutes']:.1f} min\n"
            if stats['travel_to_crossline_distance_m'] > 0:
                stats_text += f"Travel to Crossline: {stats['travel_to_crossline_distance_m']:.1f} m ({stats['travel_to_crossline_distance_km']:.3f} km, {stats['travel_to_crossline_distance_nm']:.3f} nm)\n"
            stats_text += f"Travel to Crossline Time: {stats['travel_to_crossline_minutes']:.1f} min\n"
            stats_text += "\n"

        # Survey Time Breakdown
        stats_text += "SURVEY TIME BREAKDOWN\n"
        stats_text += "-" * 25 + "\n"
        stats_text += f"Main Lines Survey: {stats['main_lines_minutes']:.1f} min\n"
        stats_text += f"Crossline Survey: {stats['crossline_minutes']:.1f} min\n"
        stats_text += f"Travel Between Lines: {stats['travel_minutes']:.1f} min\n"
        if stats['travel_to_crossline_minutes'] > 0:
            stats_text += f"Travel to Crossline: {stats['travel_to_crossline_minutes']:.1f} min\n"
        stats_text += f"Total Turn Time: {stats['total_turn_time_min']:.1f} min ({stats['num_turns_between_main_lines'] + stats['num_turns_to_crossline'] + stats['num_turns_between_crossline_passes']} turns × {stats['turn_time_per_turn_min']:.1f} min)\n\n"
        
        # Calculate total survey time (main lines + crossline, without turn times)
        # Crossline minutes includes turn times between passes, so subtract them
        crossline_survey_time_only = stats['crossline_minutes']
        if stats['num_turns_between_crossline_passes'] > 0:
            crossline_survey_time_only -= stats['turn_time_per_turn_min'] * stats['num_turns_between_crossline_passes']
        total_survey_time_min = stats['main_lines_minutes'] + crossline_survey_time_only
        
        # Calculate total transit time (travel times without turn times)
        # Travel minutes include turn times, so subtract them
        travel_between_lines_time_only = stats['travel_minutes']
        if stats['num_turns_between_main_lines'] > 0:
            travel_between_lines_time_only -= stats['turn_time_per_turn_min'] * stats['num_turns_between_main_lines']
        
        travel_to_crossline_time_only = stats['travel_to_crossline_minutes']
        if stats['num_turns_to_crossline'] > 0:
            travel_to_crossline_time_only -= stats['turn_time_per_turn_min'] * stats['num_turns_to_crossline']
        
        total_transit_time_min = travel_between_lines_time_only + travel_to_crossline_time_only
        
        # Total summary
        stats_text += "TOTAL SURVEY SUMMARY\n"
        stats_text += "-" * 25 + "\n"
        stats_text += f"Total Distance: {stats['total_distance_m']:.1f} m ({stats['total_distance_km']:.3f} km, {stats['total_distance_nm']:.3f} nm)\n"
        stats_text += f"Total Survey Time: {total_survey_time_min:.1f} min ({total_survey_time_min / 60.0:.2f} hr)\n"
        stats_text += f"Total Transit Time: {total_transit_time_min:.1f} min ({total_transit_time_min / 60.0:.2f} hr)\n"
        stats_text += f"Total Turn Time: {stats['total_turn_time_min']:.1f} min ({stats['num_turns_between_main_lines'] + stats['num_turns_to_crossline'] + stats['num_turns_between_crossline_passes']} turns × {stats['turn_time_per_turn_min']:.1f} min)\n"
        stats_text += f"Total Time: {stats['total_minutes']:.1f} min ({stats['total_hours']:.2f} hr)\n\n"

        # Survey pattern explanation
        stats_text += "SURVEY PATTERN\n"
        stats_text += "-" * 15 + "\n"
        stats_text += "1. Run main survey lines (zigzag pattern)\n"
        if stats['crossline_single_pass_distance_m'] > 0:
            stats_text += "2. Travel from last main line to crossline start\n"
            stats_text += f"3. Run crossline ({stats['num_crossline_passes']} passes)\n"
        stats_text += "\n"

        # Waypoints with labels
        stats_text += "WAYPOINTS\n"
        stats_text += "-" * 10 + "\n"
        if self.survey_lines_data:
            for i, line in enumerate(self.survey_lines_data):
                if i % 2 == 0:
                    start, end = line[0], line[1]
                else:
                    start, end = line[1], line[0]
                start_label = f'L{i+1}S'
                end_label = f'L{i+1}E'
                start_lat_ddm = decimal_degrees_to_ddm(start[0], is_latitude=True)
                start_lon_ddm = decimal_degrees_to_ddm(start[1], is_latitude=False)
                end_lat_ddm = decimal_degrees_to_ddm(end[0], is_latitude=True)
                end_lon_ddm = decimal_degrees_to_ddm(end[1], is_latitude=False)
                stats_text += f"{start_label}: {start_lat_ddm}, {start_lon_ddm} ({start[0]:.6f}, {start[1]:.6f})\n"
                stats_text += f"{end_label}: {end_lat_ddm}, {end_lon_ddm} ({end[0]:.6f}, {end[1]:.6f})\n"
        if self.cross_line_data:
            cls_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][0], is_latitude=True)
            cls_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][1], is_latitude=False)
            cle_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][0], is_latitude=True)
            cle_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][1], is_latitude=False)
            stats_text += f"CLS: {cls_lat_ddm}, {cls_lon_ddm} ({self.cross_line_data[0][0]:.6f}, {self.cross_line_data[0][1]:.6f})\n"
            stats_text += f"CLE: {cle_lat_ddm}, {cle_lon_ddm} ({self.cross_line_data[1][0]:.6f}, {self.cross_line_data[1][1]:.6f})\n"

        # Create custom dialog window with copy functionality
        show_statistics_dialog(self, "Reference Planning Statistics", stats_text)

    def _update_multiplier_label_len(self, val):
        """Updates the label next to the line length multiplier slider."""
        self.multiplier_label_len.setText(f"{float(val):.1f}")

    def _update_multiplier_label_dist(self, val):
        """Updates the label next to the distance between lines multiplier slider."""
        self.multiplier_label_dist.setText(f"{float(val):.1f}")

    def _on_parameter_changed(self):
        """Handle parameter changes by restarting the auto-regenerate timer."""
        self.auto_regenerate_timer.stop()
        self.auto_regenerate_timer.start(800)  # 800ms debounce delay

    def _auto_regenerate_survey_plan(self):
        """Auto-regenerate survey plan when timer expires, if survey plan already exists."""
        if hasattr(self, 'survey_lines_data') and len(self.survey_lines_data) > 0:
            try:
                is_valid, values = self._validate_inputs()
                if is_valid:
                    # Do not auto-zoom when parameters change; preserve current view
                    self._generate_and_plot(show_success_dialog=False, auto_zoom=False)
            except Exception as e:
                print(f"Auto-regenerate failed: {e}")

    def _on_line_length_or_speed_change(self, *args):
        pass

    def _update_export_name(self):
        try:
            dist = int(float(self.dist_between_lines_entry.text()))
            heading = int(float(self.heading_entry.text()))
            export_name = f"Reference_{dist}m_{heading}deg"
            self.export_name_entry.clear()
            self.export_name_entry.setText(export_name)
        except Exception:
            pass

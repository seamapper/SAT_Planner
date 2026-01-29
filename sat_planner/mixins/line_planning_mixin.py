"""
Line planning tab: draw/edit line, profile, export/import, set_line_info_text,
reset, statistics, handle plot click and motion for line planning.
"""
import os
import csv
import json
import datetime

import numpy as np
from PyQt6.QtWidgets import QFileDialog
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QTextCursor

from sat_planner.constants import GEOSPATIAL_LIBS_AVAILABLE, pyproj, LineString, fiona
from sat_planner import decimal_degrees_to_ddm

try:
    from shapely.geometry import mapping as _shapely_mapping
except ImportError:
    _shapely_mapping = None


class LinePlanningMixin:
    """Mixin for line planning: set_line_info_text, _toggle_line_planning_mode,
    _clear_line_planning, _handle_line_planning_plot_click, _on_line_planning_motion,
    _export_drawn_line, _import_drawn_line, _update_line_planning_button_states,
    _toggle_edit_line_planning_mode, _on_line_planning_handle_*, _reset_line_planning_tab,
    _draw_line_planning_profile, _calculate_line_planning_statistics, _show_line_information,
    _clear_line_planning_profile."""

    def set_line_info_text(self, message, append=False):
        """Add a message to the Activity Log with [Line] prefix. Maintains up to 200 lines of history."""
        if not hasattr(self, 'activity_log_text'):
            return
        prefixed_message = f"[Line] {message}"
        self.activity_log_text.setReadOnly(False)
        if append:
            current_text = self.activity_log_text.toPlainText().rstrip('\n')
            new_text = prefixed_message + "\n" + current_text if current_text else prefixed_message + "\n"
            lines = new_text.splitlines()
            if len(lines) > 200:
                lines = lines[:200]
            self.activity_log_text.setPlainText("\n".join(lines) + "\n")
            cursor = self.activity_log_text.textCursor()
            cursor.movePosition(QTextCursor.MoveOperation.Start)
            self.activity_log_text.setTextCursor(cursor)
        else:
            self.activity_log_text.append(prefixed_message)
            num_lines = self.activity_log_text.document().blockCount()
            if num_lines > 200:
                cursor = self.activity_log_text.textCursor()
                cursor.movePosition(QTextCursor.MoveOperation.Start)
                cursor.movePosition(QTextCursor.MoveOperation.Down, QTextCursor.MoveMode.MoveAnchor, num_lines - 200)
                cursor.movePosition(QTextCursor.MoveOperation.Start, QTextCursor.MoveMode.KeepAnchor)
                cursor.removeSelectedText()
            cursor = self.activity_log_text.textCursor()
            cursor.movePosition(QTextCursor.MoveOperation.End)
            self.activity_log_text.setTextCursor(cursor)
        self.activity_log_text.ensureCursorVisible()
        self.activity_log_text.setReadOnly(True)

    def _toggle_line_planning_mode(self):
        if not GEOSPATIAL_LIBS_AVAILABLE or self.geotiff_data_array is None:
            self._show_message("warning", "No GeoTIFF", "Load a GeoTIFF first to draw a line.")
            return
        self.line_planning_mode = not self.line_planning_mode
        if self.line_planning_mode:
            self.line_planning_points = []
            if self.line_planning_line is not None:
                self.line_planning_line.remove()
                self.line_planning_line = None
            self.line_start_draw_btn.setText("Drawing: Left-click to add, Right-click to finish")
            self.canvas_widget.setCursor(Qt.CursorShape.CrossCursor)
            self.set_line_info_text("Left click to start line, left click to add waypoints, and Right click to end line")
        else:
            self.line_start_draw_btn.setText("Start Drawing Line")
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)

    def _clear_line_planning(self):
        self.line_planning_points = []
        if self.line_planning_line is not None:
            self.line_planning_line.remove()
            self.line_planning_line = None
        if hasattr(self, 'line_planning_temp_line') and self.line_planning_temp_line is not None:
            self.line_planning_temp_line.remove()
            self.line_planning_temp_line = None
        if hasattr(self, 'line_planning_info_text') and self.line_planning_info_text is not None:
            self.line_planning_info_text.set_visible(False)
            self.line_planning_info_text = None
        if self.edit_line_planning_mode:
            self.edit_line_planning_mode = False
            self.line_planning_handles = []
            self.dragging_line_planning_handle = None
            if hasattr(self, 'line_planning_edit_line') and self.line_planning_edit_line is not None:
                self.line_planning_edit_line.remove()
                self.line_planning_edit_line = None
            if hasattr(self, 'line_planning_pick_cid'):
                self.canvas.mpl_disconnect(self.line_planning_pick_cid)
            if hasattr(self, 'line_planning_motion_cid'):
                self.canvas.mpl_disconnect(self.line_planning_motion_cid)
            if hasattr(self, 'line_planning_release_cid'):
                self.canvas.mpl_disconnect(self.line_planning_release_cid)
            self.line_edit_btn.setText("Edit Line Planning")
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
        self.canvas.draw_idle()
        self.line_start_draw_btn.setText("Start Drawing Line")
        if hasattr(self, 'geotiff_dataset_original') and self.geotiff_dataset_original is not None:
            if hasattr(self, 'line_start_draw_btn'):
                self.line_start_draw_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")
        else:
            if hasattr(self, 'line_start_draw_btn'):
                self.line_start_draw_btn.setStyleSheet("")
        self.line_planning_mode = False
        self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
        self._clear_line_planning_profile()
        self._update_line_planning_button_states()

    def _clear_line_planning_profile(self):
        """Clear the line planning profile plot."""
        if hasattr(self, 'profile_ax'):
            self.profile_ax.clear()
            self.profile_ax.set_title("Line Planning Elevation Profile", fontsize=8)
            self.profile_ax.set_xlabel("Distance (m)", fontsize=8)
            self.profile_ax.set_ylabel("Elevation (m)", fontsize=8)
            self.profile_ax.tick_params(axis='both', which='major', labelsize=7)
            if hasattr(self, 'profile_canvas'):
                self.profile_canvas.draw_idle()

    def _handle_line_planning_plot_click(self, event):
        """Handle plot click when in line planning tab and mode. Returns True if event was handled."""
        if not (hasattr(self, 'param_notebook') and self.param_notebook.currentIndex() == 2 and self.line_planning_mode):
            return False
        if event.inaxes != self.ax:
            return False
        if event.button == 1:  # Left click: add point
            lat, lon = event.ydata, event.xdata
            if lat is None or lon is None:
                return True
            self.line_planning_points.append((lat, lon))
            lats = [p[0] for p in self.line_planning_points]
            lons = [p[1] for p in self.line_planning_points]
            if self.line_planning_line is None:
                self.line_planning_line, = self.ax.plot(lons, lats, color='orange', linewidth=2, marker='o', label='Line Planning')
            else:
                self.line_planning_line.set_data(lons, lats)
            if hasattr(self, 'line_planning_temp_line') and self.line_planning_temp_line is not None:
                self.line_planning_temp_line.remove()
                self.line_planning_temp_line = None
            self.canvas.draw_idle()
            return True
        if event.button == 3:  # Right click: finish
            self.line_planning_mode = False
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
            if len(self.line_planning_points) >= 2:
                self.line_start_draw_btn.setText("Line Finished. Start Drawing Line")
                if hasattr(self, 'line_start_draw_btn'):
                    self.line_start_draw_btn.setStyleSheet("")
                try:
                    self._draw_line_planning_profile()
                except Exception as e:
                    print(f"Error drawing line planning profile: {e}")
                try:
                    geod = pyproj.Geod(ellps="WGS84")
                    num_segments = len(self.line_planning_points) - 1
                    total_length_m = 0.0
                    for i in range(1, len(self.line_planning_points)):
                        lat1, lon1 = self.line_planning_points[i-1]
                        lat2, lon2 = self.line_planning_points[i]
                        _, _, d = geod.inv(lon1, lat1, lon2, lat2)
                        if not np.isnan(d) and not np.isinf(d):
                            total_length_m += d
                    total_length_km = total_length_m / 1000.0
                    total_length_nm = total_length_m / 1852.0
                    try:
                        speed_knots = float(self.line_survey_speed_entry.text())
                    except Exception:
                        speed_knots = 8.0
                    speed_m_per_h = speed_knots * 1852
                    time_hours = total_length_m / speed_m_per_h if speed_m_per_h > 0 else 0
                    time_minutes = time_hours * 60
                    summary = (
                        f"Line Summary:\n"
                        f"Segments: {num_segments}\n"
                        f"Total Length: {total_length_m:.1f} m\n"
                        f"Total Length: {total_length_km:.3f} km\n"
                        f"Total Length: {total_length_nm:.3f} nautical miles\n"
                        f"Time to Run: {time_minutes:.1f} min\n"
                        f"Time to Run: {time_hours:.2f} hr"
                    )
                    self.set_line_info_text(summary)
                except Exception as e:
                    self.set_line_info_text(f"Error calculating line summary: {e}")
            else:
                self.line_start_draw_btn.setText("Start Drawing Line")
                if len(self.line_planning_points) > 0:
                    self.line_planning_points = []
                    if hasattr(self, 'line_planning_line') and self.line_planning_line is not None:
                        self.line_planning_line.remove()
                        self.line_planning_line = None
            if hasattr(self, 'line_planning_temp_line') and self.line_planning_temp_line is not None:
                self.line_planning_temp_line.remove()
                self.line_planning_temp_line = None
            if hasattr(self, 'line_planning_info_text') and self.line_planning_info_text is not None:
                self.line_planning_info_text.set_visible(False)
                self.line_planning_info_text = None
            if len(self.line_planning_points) >= 2:
                self._plot_survey_plan(preserve_view_limits=True)
            else:
                self.canvas.draw_idle()
            self._update_line_planning_button_states()
            return True
        return False

    def _export_drawn_line(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded. Cannot export line.")
            return
        if not self.line_planning_points or len(self.line_planning_points) < 2:
            self._show_message("warning", "No Line", "Draw a line with at least two points before exporting.")
            return
        export_dir = QFileDialog.getExistingDirectory(self, "Select Export Directory", self.last_used_dir)
        if not export_dir:
            return
        self.last_used_dir = export_dir
        self._save_last_used_dir()
        if hasattr(self, 'line_export_name_entry'):
            export_name = self.line_export_name_entry.text().strip()
            if not export_name:
                export_name = f"Line_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
                self.line_export_name_entry.setText(export_name)
        else:
            export_name = f"LinePlanning_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
        try:
            csv_file_path = os.path.join(export_dir, f"{export_name}_DD.csv")
            with open(csv_file_path, 'w', newline='') as csvfile:
                csv_writer = csv.writer(csvfile)
                csv_writer.writerow(['Point Label', 'Latitude', 'Longitude'])
                for i, (lat, lon) in enumerate(self.line_planning_points):
                    csv_writer.writerow([i + 1, lat, lon])
            ddm_file_path = os.path.join(export_dir, f"{export_name}_DM.csv")
            with open(ddm_file_path, 'w', newline='', encoding='utf-8') as ddmfile:
                ddm_writer = csv.writer(ddmfile)
                ddm_writer.writerow(['Line Number', 'Point Label', 'Latitude with Decimal Minutes', 'Longitude with Decimal Minutes'])
                for i, (lat, lon) in enumerate(self.line_planning_points):
                    lat_ddm = decimal_degrees_to_ddm(lat, is_latitude=True)
                    lon_ddm = decimal_degrees_to_ddm(lon, is_latitude=False)
                    ddm_writer.writerow([1, i + 1, lat_ddm, lon_ddm])
            ddm_txt_file_path = os.path.join(export_dir, f"{export_name}_DM.txt")
            with open(ddm_txt_file_path, 'w', encoding='utf-8') as ddm_txt_file:
                for i, (lat, lon) in enumerate(self.line_planning_points):
                    lat_ddm = decimal_degrees_to_ddm(lat, is_latitude=True)
                    lon_ddm = decimal_degrees_to_ddm(lon, is_latitude=False)
                    ddm_txt_file.write(f"{i + 1}, {lat_ddm}, {lon_ddm}\n")
            shapefile_path = None
            if LineString and fiona and _shapely_mapping:
                schema = {'geometry': 'LineString', 'properties': {'name': 'str'}}
                crs_epsg = 'EPSG:4326'
                shapely_line = LineString([(lon, lat) for lat, lon in self.line_planning_points])
                shapefile_path = os.path.join(export_dir, f"{export_name}.shp")
                with fiona.open(shapefile_path, 'w', driver='ESRI Shapefile', crs=crs_epsg, schema=schema) as collection:
                    collection.writerecords([{'geometry': _shapely_mapping(shapely_line), 'properties': {'name': export_name}}])
            geojson_file_path = os.path.join(export_dir, f"{export_name}.geojson")
            if LineString and _shapely_mapping:
                shapely_line = LineString([(lon, lat) for lat, lon in self.line_planning_points])
                geojson_feature = {"type": "Feature", "geometry": _shapely_mapping(shapely_line), "properties": {"name": export_name, "points": [{"point_num": i + 1, "lat": lat, "lon": lon} for i, (lat, lon) in enumerate(self.line_planning_points)]}}
                with open(geojson_file_path, 'w') as f:
                    json.dump({"type": "FeatureCollection", "features": [geojson_feature]}, f, indent=2)
            lnw_file_path = os.path.join(export_dir, f"{export_name}.lnw")
            with open(lnw_file_path, 'w') as f:
                f.write("LNW 1.0\n")
                for i, (lat, lon) in enumerate(self.line_planning_points):
                    depth = self._get_depth_at_point(lat, lon) if self.geotiff_data_array is not None else 0.0
                    try:
                        speed = float(self.line_survey_speed_entry.text()) if self.line_survey_speed_entry.text() else 8.0
                    except Exception:
                        speed = 8.0
                    f.write(f"WAYPOINT_{i+1:03d},{lat:.6f},{lon:.6f},{abs(depth):.1f},{speed:.1f},50.0,1,1\n")
            sis_file_path = os.path.join(export_dir, f"{export_name}.asciiplan")
            with open(sis_file_path, 'w') as f:
                f.write("SIS ASCII Plan\n")
                for i, (lat, lon) in enumerate(self.line_planning_points):
                    try:
                        speed_knots = float(self.line_survey_speed_entry.text()) if self.line_survey_speed_entry.text() else 8.0
                    except Exception:
                        speed_knots = 8.0
                    depth = self._get_depth_at_point(lat, lon)
                    f.write(f"WAYPOINT_{i+1:03d}, {lat:.6f}, {lon:.6f}, {depth:.1f}, {speed_knots:.1f}, 1, 1\n")
            map_png_path = os.path.join(export_dir, f"{export_name}_map.png")
            self.figure.savefig(map_png_path, dpi=300, bbox_inches='tight', facecolor='white')
            profile_png_path = None
            if hasattr(self, 'profile_fig') and self.profile_fig is not None:
                profile_png_path = os.path.join(export_dir, f"{export_name}_profile.png")
                self.profile_fig.savefig(profile_png_path, dpi=300, bbox_inches='tight', facecolor='white')
            stats = self._calculate_line_planning_statistics()
            stats_file_path = None
            if stats:
                stats_file_path = os.path.join(export_dir, f"{export_name}_statistics.txt")
                with open(stats_file_path, 'w') as f:
                    f.write("LINE PLANNING STATISTICS\n" + "=" * 30 + "\n\n")
                    f.write(f"Number of Points: {stats['num_points']}\n")
                    f.write(f"Total Distance: {stats['total_distance_m']:.1f} m ({stats['total_distance_km']:.3f} km, {stats['total_distance_nm']:.3f} nm)\n")
                    f.write(f"Survey Speed: {stats['speed_knots']:.1f} knots\n")
                    f.write(f"Estimated Time: {stats['total_time_minutes']:.1f} min ({stats['total_time_hours']:.2f} hr)\n\n")
                    if stats['depth_info']:
                        f.write(stats['depth_info'])
                    if len(stats['segment_distances']) > 1:
                        f.write("\nSEGMENT DISTANCE AND HEADINGS\n" + "-" * 30 + "\n")
                        for i in range(len(self.line_planning_points) - 1):
                            lat1, lon1 = self.line_planning_points[i]
                            lat2, lon2 = self.line_planning_points[i + 1]
                            if pyproj is not None:
                                geod = pyproj.Geod(ellps="WGS84")
                                fwd_az, back_az, dist = geod.inv(lon1, lat1, lon2, lat2)
                                f.write(f"Segment {i+1}: {dist:.1f} m ({dist/1000:.3f} km, {dist/1852:.3f} nm) - Heading: {fwd_az % 360:.1f}°\n")
                            else:
                                f.write(f"Segment {i+1}: {stats['segment_distances'][i]:.1f} m - Heading: Unable to calculate\n")
                        f.write("\nWAYPOINTS\n" + "-" * 10 + "\n")
                        for i, (lat, lon) in enumerate(self.line_planning_points):
                            f.write(f"WP{i+1}: {decimal_degrees_to_ddm(lat, True)}, {decimal_degrees_to_ddm(lon, False)} ({lat:.6f}, {lon:.6f})\n")
            success_msg = f"Line exported successfully to:\n- {os.path.basename(csv_file_path)}\n"
            if shapefile_path:
                success_msg += f"- {os.path.basename(shapefile_path)} (and associated files)\n"
            success_msg += f"- {os.path.basename(geojson_file_path)}\n- {os.path.basename(lnw_file_path)}\n- {os.path.basename(sis_file_path)}\n- {os.path.basename(map_png_path)}\n"
            if profile_png_path:
                success_msg += f"- {os.path.basename(profile_png_path)}\n"
            if stats_file_path:
                success_msg += f"- {os.path.basename(stats_file_path)}\n"
            success_msg += f"in directory: {export_dir}"
            self.set_line_info_text(success_msg, append=False)
        except Exception as e:
            self.set_line_info_text(f"Failed to export drawn line: {e}")

    def _import_drawn_line(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select Line Plan File to Import", self.last_line_import_dir,
            "Decimal Degree CSV files (*_DD.csv);;CSV files (*.csv);;GeoJSON files (*.geojson);;JSON files (*.json);;All files (*.*)")
        if not file_path:
            return
        import_dir = os.path.dirname(file_path)
        if import_dir and os.path.isdir(import_dir):
            self.last_line_import_dir = import_dir
            self._save_last_line_import_dir()
        try:
            self.line_planning_points = []
            file_ext = os.path.splitext(file_path)[1].lower()
            if file_ext == '.csv':
                with open(file_path, 'r', encoding='utf-8') as csvfile:
                    csv_reader = csv.DictReader(csvfile)
                    for row in csv_reader:
                        try:
                            lat = float(row.get('Latitude', 0))
                            lon = float(row.get('Longitude', 0))
                            self.line_planning_points.append((lat, lon))
                        except (ValueError, TypeError):
                            continue
            elif file_ext in ['.geojson', '.json']:
                with open(file_path, 'r', encoding='utf-8') as f:
                    geojson_data = json.load(f)
                features = geojson_data.get('features', []) if geojson_data.get('type') == 'FeatureCollection' else [geojson_data] if geojson_data.get('type') == 'Feature' else []
                for feature in features:
                    if feature.get('type') != 'Feature':
                        continue
                    geometry = feature.get('geometry', {})
                    if geometry.get('type') != 'LineString':
                        continue
                    for coord in geometry.get('coordinates', []):
                        if len(coord) >= 2:
                            self.line_planning_points.append((coord[1], coord[0]))
            else:
                self._show_message("error", "Import Error", f"Unsupported file format: {file_ext}")
                return
            if len(self.line_planning_points) < 2:
                self._show_message("warning", "Import Warning", "Line plan must have at least 2 points.")
                self.line_planning_points = []
                return
            if hasattr(self, 'line_export_name_entry'):
                file_basename = os.path.splitext(os.path.basename(file_path))[0]
                for suffix in ['_DD', '_DM', '.geojson', '.shp', '.lnw']:
                    if file_basename.endswith(suffix):
                        file_basename = file_basename[:-len(suffix)]
                self.line_export_name_entry.setText(file_basename)
            self._update_line_planning_button_states()
            self._plot_survey_plan(preserve_view_limits=True)
            self.set_line_info_text(f"Successfully imported line plan with {len(self.line_planning_points)} points.")
            if hasattr(self, 'line_start_draw_btn') and len(self.line_planning_points) >= 2:
                self.line_start_draw_btn.setStyleSheet("")
        except Exception as e:
            self._show_message("error", "Import Error", f"Failed to import line plan: {e}")

    def _on_line_planning_motion(self, event):
        """Handle mouse motion in line planning mode (temp line and info). Used by _on_mouse_motion."""
        if not (hasattr(self, 'line_planning_mode') and self.line_planning_mode):
            return
        if len(self.line_planning_points) < 1:
            return
        lats = [p[0] for p in self.line_planning_points] + [event.ydata]
        lons = [p[1] for p in self.line_planning_points] + [event.xdata]
        if event.ydata is None or event.xdata is None:
            return
        if hasattr(self, 'line_planning_temp_line') and self.line_planning_temp_line is not None:
            self.line_planning_temp_line.set_data(lons, lats)
        else:
            self.line_planning_temp_line, = self.ax.plot(lons, lats, color='orange', linewidth=2, linestyle='--', alpha=0.7)
        if hasattr(self, 'geotiff_data_array') and self.geotiff_data_array is not None:
            geod = pyproj.Geod(ellps="WGS84")
            total_length_m = 0.0
            for i in range(1, len(lats)):
                _, _, seg_length = geod.inv(lons[i-1], lats[i-1], lons[i], lats[i])
                total_length_m += seg_length
            total_length_nm = total_length_m / 1852.0
            try:
                speed_knots = float(self.line_survey_speed_entry.text()) if hasattr(self, 'line_survey_speed_entry') and self.line_survey_speed_entry.text() else 8.0
            except Exception:
                speed_knots = 8.0
            speed_m_per_h = speed_knots * 1852
            time_hours = total_length_m / speed_m_per_h if speed_m_per_h > 0 else 0
            time_minutes = time_hours * 60
            elevation, slope = self._calculate_slope_at_point(event.ydata, event.xdata)
            try:
                depth_str = f"{abs(elevation):.1f} m" if elevation is not None and not np.isnan(elevation) else "-"
                slope_str = f"{slope:.1f}°" if slope is not None and not np.isnan(slope) else "-"
            except Exception:
                depth_str = "-"
                slope_str = "-"
            info_str = f"Depth: {depth_str}\nSlope: {slope_str}\nLine Length: {total_length_m:.1f} m\nLine Length: {total_length_nm:.3f} nm\nSurvey Time: {time_minutes:.1f} min"
            if hasattr(self, 'mouse_hover_info_text') and self.mouse_hover_info_text is not None:
                self.mouse_hover_info_text.set_text(info_str)
                self.mouse_hover_info_text.set_visible(True)
            else:
                self.mouse_hover_info_text = self.ax.text(0.02, 0.98, info_str, transform=self.ax.transAxes, fontsize=9, va='top', ha='left', bbox=dict(boxstyle='round', facecolor='orange', alpha=0.8), zorder=10)
        self.canvas.draw_idle()

    def _update_line_planning_button_states(self):
        has_line_planning = hasattr(self, 'line_planning_points') and len(self.line_planning_points) >= 2
        if hasattr(self, 'line_edit_btn'):
            self.line_edit_btn.setEnabled(has_line_planning)

    def _toggle_edit_line_planning_mode(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded. Cannot edit line planning.")
            return
        if self.geotiff_dataset_original is None:
            self._show_message("warning", "No GeoTIFF", "Load a GeoTIFF first to edit line planning.")
            return
        if len(self.line_planning_points) < 2:
            self._show_message("warning", "No Line Planning", "Draw a line planning route first before editing it.")
            return
        self.edit_line_planning_mode = not self.edit_line_planning_mode
        if self.edit_line_planning_mode:
            self.line_planning_handles = []
            for i, (lat, lon) in enumerate(self.line_planning_points):
                color = 'red' if i == 0 else 'blue' if i == len(self.line_planning_points) - 1 else 'green'
                handle = self.ax.scatter([lon], [lat], color=color, s=100, marker='o', edgecolors='black', linewidth=2, zorder=10, picker=5)
                self.line_planning_handles.append(handle)
            self.line_planning_pick_cid = self.canvas.mpl_connect('pick_event', self._on_line_planning_handle_pick)
            self.line_planning_motion_cid = self.canvas.mpl_connect('motion_notify_event', self._on_line_planning_handle_motion)
            self.line_planning_release_cid = self.canvas.mpl_connect('button_release_event', self._on_line_planning_handle_release)
            self.line_edit_btn.setText("Click to Stop Editing Line Planning")
            self.canvas_widget.setCursor(Qt.CursorShape.SizeAllCursor)
            self.set_line_info_text("You can drag the red (start), green (waypoints), and blue (end) points to edit the line planning route.")
        else:
            if self.line_planning_handles:
                self._plot_survey_plan(preserve_view_limits=True)
                self.line_planning_handles = []
            if hasattr(self, 'line_planning_pick_cid'):
                self.canvas.mpl_disconnect(self.line_planning_pick_cid)
            if hasattr(self, 'line_planning_motion_cid'):
                self.canvas.mpl_disconnect(self.line_planning_motion_cid)
            if hasattr(self, 'line_planning_release_cid'):
                self.canvas.mpl_disconnect(self.line_planning_release_cid)
            self.line_edit_btn.setText("Edit Line Planning")
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
            self.dragging_line_planning_handle = None
            self._plot_survey_plan(preserve_view_limits=True)
            self._draw_line_planning_profile()
            try:
                geod = pyproj.Geod(ellps="WGS84")
                num_segments = len(self.line_planning_points) - 1
                total_length_m = 0.0
                for i in range(1, len(self.line_planning_points)):
                    lat1, lon1 = self.line_planning_points[i-1]
                    lat2, lon2 = self.line_planning_points[i]
                    _, _, d = geod.inv(lon1, lat1, lon2, lat2)
                    total_length_m += d
                try:
                    speed_knots = float(self.line_survey_speed_entry.text())
                except Exception:
                    speed_knots = 8.0
                speed_m_per_h = speed_knots * 1852
                time_hours = total_length_m / speed_m_per_h if speed_m_per_h > 0 else 0
                time_minutes = time_hours * 60
                summary = (f"Updated Line Planning Summary:\nSegments: {num_segments}\nTotal Length: {total_length_m:.1f} m\nTotal Length: {total_length_m/1000:.3f} km\nTotal Length: {total_length_m/1852:.3f} nautical miles\nTime to Run: {time_minutes:.1f} min\nTime to Run: {time_hours:.2f} hr")
                self.set_line_info_text(summary)
            except Exception as e:
                self.set_line_info_text(f"Error calculating updated line planning summary: {e}")

    def _on_line_planning_handle_pick(self, event):
        if not self.edit_line_planning_mode or event.mouseevent.button != 1:
            return
        for i, handle in enumerate(self.line_planning_handles):
            if event.artist == handle:
                self.dragging_line_planning_handle = i
                break

    def _on_line_planning_handle_motion(self, event):
        if not self.edit_line_planning_mode or self.dragging_line_planning_handle is None or event.inaxes != self.ax or event.xdata is None or event.ydata is None:
            return
        handle_idx = self.dragging_line_planning_handle
        if handle_idx < len(self.line_planning_handles) and self.line_planning_handles[handle_idx] is not None:
            self.line_planning_handles[handle_idx].set_offsets([[event.xdata, event.ydata]])
            self.line_planning_points[handle_idx] = (event.ydata, event.xdata)
        if hasattr(self, 'line_planning_edit_line') and self.line_planning_edit_line is not None:
            lats = [p[0] for p in self.line_planning_points]
            lons = [p[1] for p in self.line_planning_points]
            self.line_planning_edit_line.set_data(lons, lats)
        else:
            lats = [p[0] for p in self.line_planning_points]
            lons = [p[1] for p in self.line_planning_points]
            self.line_planning_edit_line, = self.ax.plot(lons, lats, color='orange', linewidth=2, zorder=9)
        self.canvas.draw_idle()

    def _on_line_planning_handle_release(self, event):
        if not self.edit_line_planning_mode:
            return
        if self.dragging_line_planning_handle is not None:
            if hasattr(self, 'line_planning_edit_line') and self.line_planning_edit_line is not None:
                self.line_planning_edit_line.remove()
                self.line_planning_edit_line = None
            self._plot_survey_plan(preserve_view_limits=True)
            self.line_planning_handles = []
            for i, (lat, lon) in enumerate(self.line_planning_points):
                color = 'red' if i == 0 else 'blue' if i == len(self.line_planning_points) - 1 else 'green'
                handle = self.ax.scatter([lon], [lat], color=color, s=100, marker='o', edgecolors='black', linewidth=2, zorder=10, picker=5)
                self.line_planning_handles.append(handle)
            self.dragging_line_planning_handle = None

    def _reset_line_planning_tab(self):
        self.line_planning_points = []
        self.line_planning_mode = False
        self.edit_line_planning_mode = False
        if hasattr(self, 'line_planning_line') and self.line_planning_line is not None:
            self.line_planning_line.remove()
            self.line_planning_line = None
        if hasattr(self, 'line_planning_temp_line') and self.line_planning_temp_line is not None:
            self.line_planning_temp_line.remove()
            self.line_planning_temp_line = None
        if hasattr(self, 'line_planning_edit_line') and self.line_planning_edit_line is not None:
            self.line_planning_edit_line.remove()
            self.line_planning_edit_line = None
        if hasattr(self, 'line_planning_info_text') and self.line_planning_info_text is not None:
            self.line_planning_info_text.set_visible(False)
            self.line_planning_info_text = None
        self.line_planning_handles = []
        self.dragging_line_planning_handle = None
        if hasattr(self, 'line_survey_speed_entry'):
            self.line_survey_speed_entry.clear()
            self.line_survey_speed_entry.setText("8")
        if hasattr(self, 'line_start_draw_btn'):
            self.line_start_draw_btn.setText("Start Drawing Line")
        if hasattr(self, 'line_edit_btn'):
            self.line_edit_btn.setEnabled(False)
        if hasattr(self, 'canvas_widget'):
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
        self._update_line_planning_button_states()
        self._plot_survey_plan(preserve_view_limits=True)

    def _draw_line_planning_profile(self):
        self.profile_ax.clear()
        for ax in self.profile_fig.get_axes():
            if ax != self.profile_ax:
                ax.remove()
        self.profile_ax.set_title("Line Planning Elevation Profile", fontsize=8)
        self.profile_ax.set_xlabel("Distance (m)", fontsize=8)
        self.profile_ax.set_ylabel("Elevation (m)", fontsize=8)
        self.profile_ax.tick_params(axis='both', which='major', labelsize=7)
        slope_ax = None
        if not (hasattr(self, 'line_planning_points') and len(self.line_planning_points) >= 2):
            self.profile_fig.tight_layout(pad=1.0)
            if hasattr(self, 'profile_canvas'):
                self.profile_canvas.draw_idle()
            return
        geod = pyproj.Geod(ellps="WGS84")
        lats, lons, dists = [], [], [0.0]
        total_dist = 0.0
        waypoint_distances = [0.0]
        for i in range(1, len(self.line_planning_points)):
            lat1, lon1 = self.line_planning_points[i-1]
            lat2, lon2 = self.line_planning_points[i]
            seg_lats = np.linspace(lat1, lat2, 50)
            seg_lons = np.linspace(lon1, lon2, 50)
            if i == 1:
                lats.extend(seg_lats)
                lons.extend(seg_lons)
            else:
                lats.extend(seg_lats[1:])
                lons.extend(seg_lons[1:])
            for j in range(1, 50):
                _, _, d = geod.inv(seg_lons[j-1], seg_lats[j-1], seg_lons[j], seg_lats[j])
                if not np.isnan(d) and not np.isinf(d):
                    total_dist += d
                dists.append(total_dist)
            waypoint_distances.append(total_dist)
        lats = np.array(lats)
        lons = np.array(lons)
        dists = np.array(dists)
        elevations, slopes, _ = self._get_profile_data_from_geotiff(lats, lons)
        if elevations is None:
            if self.geotiff_data_array is None or self.geotiff_extent is None:
                return
            left, right, bottom, top = tuple(self.geotiff_extent)
            nrows, ncols = self.geotiff_data_array.shape
            rows = ((top - lats) / (top - bottom) * (nrows - 1)).clip(0, nrows - 1)
            cols = ((lons - left) / (right - left) * (ncols - 1)).clip(0, ncols - 1)
            elevations, slopes = [], []
            for r, c in zip(rows, cols):
                ir, ic = int(round(r)), int(round(c))
                elev = self.geotiff_data_array[ir, ic]
                elevations.append(elev)
                slope = np.nan
                if 0 < ir < nrows-1 and 0 < ic < ncols-1:
                    center_lat_geotiff = (self.geotiff_extent[2] + self.geotiff_extent[3]) / 2
                    m_per_deg_lat, m_per_deg_lon = 111320.0, 111320.0 * np.cos(np.radians(center_lat_geotiff))
                    res_lat_deg = (self.geotiff_extent[3] - self.geotiff_extent[2]) / nrows
                    res_lon_deg = (self.geotiff_extent[1] - self.geotiff_extent[0]) / ncols
                    dx_m, dy_m = res_lon_deg * m_per_deg_lon, res_lat_deg * m_per_deg_lat
                    window = self.geotiff_data_array[ir-1:ir+2, ic-1:ic+2]
                    if window.shape == (3, 3) and not np.all(np.isnan(window)):
                        dz_dy, dz_dx = np.gradient(window, dy_m, dx_m)
                        slope = np.degrees(np.arctan(np.sqrt(dz_dx[1,1]**2 + dz_dy[1,1]**2)))
                slopes.append(slope)
            elevations = np.array(elevations)
            slopes = np.array(slopes)
        self.profile_ax.plot(dists, elevations, color='orange', lw=1, label='Elevation')
        waypoint_elevations = []
        for i, (lat, lon) in enumerate(self.line_planning_points):
            _, idx = min((abs(lat - lats[j]) + abs(lon - lons[j]), j) for j in range(len(lats)))
            waypoint_elevations.append(elevations[idx])
        self.profile_ax.plot(waypoint_distances, waypoint_elevations, 'o', color='orange', markersize=6, label='Waypoints')
        if np.any(~np.isnan(slopes)):
            slope_ax = self.profile_ax.twinx()
            slope_ax.plot(dists, slopes, color='teal', lw=1, linestyle='--', label='Slope (deg)')
            slope_ax.set_ylabel('Slope (deg)', fontsize=8)
            slope_ax.tick_params(axis='y', labelsize=7)
            slope_ax.grid(False)
        if len(dists) > 0 and not (np.isnan(dists[0]) or np.isinf(dists[0]) or np.isnan(dists[-1]) or np.isinf(dists[-1])):
            self.profile_ax.set_xlim(dists[0], dists[-1])
        if np.any(~np.isnan(elevations)):
            self.profile_ax.set_ylim(np.nanmin(elevations), np.nanmax(elevations))
        if slope_ax and np.any(~np.isnan(slopes)):
            slope_ax.set_ylim(0, np.nanmax(slopes[~np.isnan(slopes)]) * 1.1)
        handles, labels = self.profile_ax.get_legend_handles_labels()
        if slope_ax:
            slope_handles, slope_labels = slope_ax.get_legend_handles_labels()
            handles.extend(slope_handles)
            labels.extend(slope_labels)
        if handles:
            self.profile_ax.legend(handles, labels, fontsize=10, loc='upper right')
        self.profile_fig.tight_layout(pad=1.0)
        self.profile_canvas.draw_idle()

    def _calculate_line_planning_statistics(self):
        if not self.line_planning_points or len(self.line_planning_points) < 2:
            return None
        try:
            if pyproj is not None:
                geod = pyproj.Geod(ellps="WGS84")
                total_distance = 0.0
                segment_distances = []
                for i in range(len(self.line_planning_points) - 1):
                    lat1, lon1 = self.line_planning_points[i]
                    lat2, lon2 = self.line_planning_points[i + 1]
                    fwd_az, back_az, dist = geod.inv(lon1, lat1, lon2, lat2)
                    segment_distances.append(dist)
                    total_distance += dist
                if len(self.line_planning_points) >= 2:
                    lat1, lon1 = self.line_planning_points[0]
                    lat2, lon2 = self.line_planning_points[1]
                    fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)
                    heading, reciprocal_heading = fwd_az % 360, back_az % 360
                else:
                    heading, reciprocal_heading = 0, 0
                try:
                    speed_knots = float(self.line_survey_speed_entry.text()) if self.line_survey_speed_entry.text() else 8.0
                except Exception:
                    speed_knots = 8.0
                speed_ms = speed_knots * 0.514444
                total_time_seconds = total_distance / speed_ms
                total_time_minutes = total_time_seconds / 60
                total_time_hours = total_time_minutes / 60
                depth_info = ""
                if self.geotiff_data_array is not None:
                    depths = [self._get_depth_at_point(lat, lon) for lat, lon in self.line_planning_points]
                    if depths:
                        min_d, max_d, avg_d = min(depths), max(depths), sum(depths) / len(depths)
                        depth_info = f"Depth Range: {min_d:.1f} to {max_d:.1f} m (avg: {avg_d:.1f} m)\n"
                return {'num_points': len(self.line_planning_points), 'total_distance_m': total_distance, 'total_distance_km': total_distance / 1000, 'total_distance_nm': total_distance / 1852, 'heading': heading, 'reciprocal_heading': reciprocal_heading, 'speed_knots': speed_knots, 'total_time_minutes': total_time_minutes, 'total_time_hours': total_time_hours, 'depth_info': depth_info, 'segment_distances': segment_distances}
            return None
        except Exception as e:
            print(f"Error calculating line planning statistics: {e}")
            return None

    def _show_line_information(self):
        stats = self._calculate_line_planning_statistics()
        if not stats:
            self._show_message("warning", "No Line Data", "No line has been drawn. Please draw a line first.")
            return
        stats_text = "LINE PLANNING STATISTICS\n" + "=" * 30 + "\n\n"
        stats_text += f"Number of Points: {stats['num_points']}\n"
        stats_text += f"Total Distance: {stats['total_distance_m']:.1f} m ({stats['total_distance_km']:.3f} km, {stats['total_distance_nm']:.3f} nm)\n"
        stats_text += f"Survey Speed: {stats['speed_knots']:.1f} knots\n"
        stats_text += f"Estimated Time: {stats['total_time_minutes']:.1f} min ({stats['total_time_hours']:.2f} hr)\n\n"
        if stats['depth_info']:
            stats_text += stats['depth_info']
        if len(stats['segment_distances']) > 1:
            stats_text += "\nSEGMENT DISTANCE AND HEADINGS\n" + "-" * 30 + "\n"
            for i in range(len(self.line_planning_points) - 1):
                lat1, lon1 = self.line_planning_points[i]
                lat2, lon2 = self.line_planning_points[i + 1]
                if pyproj is not None:
                    geod = pyproj.Geod(ellps="WGS84")
                    fwd_az, back_az, dist = geod.inv(lon1, lat1, lon2, lat2)
                    stats_text += f"Segment {i+1}: {dist:.1f} m ({dist/1000:.3f} km, {dist/1852:.3f} nm) - Heading: {fwd_az % 360:.1f}°\n"
                else:
                    stats_text += f"Segment {i+1}: {stats['segment_distances'][i]:.1f} m - Heading: Unable to calculate\n"
        if self.line_planning_points:
            stats_text += "\nWAYPOINTS\n" + "-" * 10 + "\n"
            for i, (lat, lon) in enumerate(self.line_planning_points):
                stats_text += f"WP{i+1}: {decimal_degrees_to_ddm(lat, True)}, {decimal_degrees_to_ddm(lon, False)} ({lat:.6f}, {lon:.6f})\n"
        self._show_statistics_dialog("Line Planning Statistics", stats_text)

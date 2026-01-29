"""
Calibration tab: pitch/roll/heading lines, pick/edit modes, handle events,
export/import cal survey files, set_cal_info_text, reset, statistics.
"""
import os
import csv
import json
import re
import datetime

from PyQt6.QtWidgets import QFileDialog
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QTextCursor, QColor, QTextCharFormat

import numpy as np
from sat_planner.constants import GEOSPATIAL_LIBS_AVAILABLE, pyproj
from sat_planner import decimal_degrees_to_ddm
from sat_planner.utils_ui import show_statistics_dialog


class CalibrationMixin:
    """Mixin for calibration: set_cal_info_text, toggle pick/edit pitch/roll,
    pitch/roll handle events, _export_cal_survey_files, _import_cal_survey_files,
    _load_last_cal_import_dir, _save_last_cal_import_dir, _clear_calibration_lines,
    _update_pitch_line_button_states, _update_roll_line_button_states,
    _update_cal_line_times, _reset_calibration_tab,
    _calculate_calibration_survey_statistics, _show_calibration_statistics."""

    def set_cal_info_text(self, message, append=False):
        """Add a message to the Activity Log with [Calibration] prefix. Maintains up to 200 lines of history."""
        if not hasattr(self, 'activity_log_text'):
            return
        prefixed_message = f"[Calibration] {message}"
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

    def _toggle_pick_pitch_line_mode(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot pick pitch line.")
            return
        if self.geotiff_dataset_original is None:
            self._show_message("warning","No GeoTIFF", "Load a GeoTIFF first to pick a pitch line.")
            return
        self.pick_pitch_line_mode = not self.pick_pitch_line_mode
        if self.pick_pitch_line_mode:
            self.pitch_line_points = []
            self.pick_pitch_line_btn.setText("Drawing Pitch Line: Click Start Point")
            if hasattr(self, 'calibration_frame'):
                self.calibration_frame.setFocus()
            self.canvas_widget.setCursor(Qt.CursorShape.CrossCursor)
            self.set_cal_info_text("Click the start point, and then the end point of the pitch line on the plot.")
            if hasattr(self, 'activity_log_text'):
                self.activity_log_text.setReadOnly(False)
                cursor = self.activity_log_text.textCursor()
                cursor.movePosition(QTextCursor.MoveOperation.End)
                char_format = QTextCharFormat()
                char_format.setForeground(QColor(255, 165, 0))
                cursor.setCharFormat(char_format)
                cursor.insertText("[Calibration] Left click the starting and ending points of the pitch line on the plot.\n")
                char_format.setForeground(QColor(0, 0, 0))
                cursor.setCharFormat(char_format)
                self.activity_log_text.setTextCursor(cursor)
                self.activity_log_text.setReadOnly(True)
            if hasattr(self, 'pitch_line_info_text') and self.pitch_line_info_text is not None:
                self.pitch_line_info_text.set_visible(False)
                self.pitch_line_info_text = None
                self.canvas.draw_idle()
        else:
            self.pick_pitch_line_btn.setText("Draw a Pitch Line")
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
            if hasattr(self, 'pitch_line_info_text') and self.pitch_line_info_text is not None:
                self.pitch_line_info_text.set_visible(False)
                self.pitch_line_info_text = None
                self.canvas.draw_idle()
        self._update_cal_line_times()

    def _toggle_edit_pitch_line_mode(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot edit pitch line.")
            return
        if self.geotiff_dataset_original is None:
            self._show_message("warning","No GeoTIFF", "Load a GeoTIFF first to edit a pitch line.")
            return
        if len(self.pitch_line_points) < 2:
            self._show_message("warning","No Pitch Line", "Draw a pitch line first before editing it.")
            return

        self.edit_pitch_line_mode = not self.edit_pitch_line_mode
        if self.edit_pitch_line_mode:
            self.heading_lines = []
            self._plot_survey_plan()

            start_lat, start_lon = self.pitch_line_points[0]
            end_lat, end_lon = self.pitch_line_points[1]

            self.pitch_line_start_handle = self.ax.scatter([start_lon], [start_lat],
                                                         color='red', s=100, marker='o',
                                                         edgecolors='black', linewidth=2,
                                                         zorder=10, picker=5)
            self.pitch_line_end_handle = self.ax.scatter([end_lon], [end_lat],
                                                       color='blue', s=100, marker='o',
                                                       edgecolors='black', linewidth=2,
                                                       zorder=10, picker=5)

            self.pitch_line_pick_cid = self.canvas.mpl_connect('pick_event', self._on_pitch_line_handle_pick)
            self.pitch_line_motion_cid = self.canvas.mpl_connect('motion_notify_event', self._on_pitch_line_handle_motion)
            self.pitch_line_release_cid = self.canvas.mpl_connect('button_release_event', self._on_pitch_line_handle_release)

            self.edit_pitch_line_btn.setText("Click to Stop Editing Pitch Line")
            self.canvas_widget.setCursor(Qt.CursorShape.SizeAllCursor)
            self.set_cal_info_text("You can drag the red (start) and blue (end) points to edit the pitch line. Heading lines have been cleared.")
        else:
            if self.pitch_line_start_handle or self.pitch_line_end_handle:
                self._plot_survey_plan()
                self.pitch_line_start_handle = None
                self.pitch_line_end_handle = None

            if hasattr(self, 'pitch_line_pick_cid'):
                self.canvas.mpl_disconnect(self.pitch_line_pick_cid)
            if hasattr(self, 'pitch_line_motion_cid'):
                self.canvas.mpl_disconnect(self.pitch_line_motion_cid)
            if hasattr(self, 'pitch_line_release_cid'):
                self.canvas.mpl_disconnect(self.pitch_line_release_cid)

            self.edit_pitch_line_btn.setText("Edit Pitch Line")
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
            self.dragging_pitch_line_handle = None

            self._plot_survey_plan()
            self._update_cal_line_offset_from_pitch_line()
            self._draw_pitch_line_profile()
            self._update_cal_export_name_from_pitch_line()
            self._update_cal_line_times()

            self.set_cal_info_text("Pitch line editing completed.")

    def _toggle_pick_roll_line_mode(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot pick roll line.")
            return
        if self.geotiff_dataset_original is None:
            self._show_message("warning","No GeoTIFF", "Load a GeoTIFF first to pick a roll line.")
            return
        self.pick_roll_line_mode = not self.pick_roll_line_mode
        if self.pick_roll_line_mode:
            self.roll_line_points = []
            self.pick_roll_line_btn.setText("Drawing Roll Line: Click Start Point")
            if hasattr(self, 'calibration_frame'):
                self.calibration_frame.setFocus()
            self.canvas_widget.setCursor(Qt.CursorShape.CrossCursor)
            self.set_cal_info_text("Click the start point, and then the end point of the roll line on the plot.")
            if hasattr(self, 'activity_log_text'):
                self.activity_log_text.setReadOnly(False)
                cursor = self.activity_log_text.textCursor()
                cursor.movePosition(QTextCursor.MoveOperation.End)
                char_format = QTextCharFormat()
                char_format.setForeground(QColor(255, 165, 0))
                cursor.setCharFormat(char_format)
                cursor.insertText("[Calibration] Left click the starting and ending points of the roll line on the plot.\n")
                char_format.setForeground(QColor(0, 0, 0))
                cursor.setCharFormat(char_format)
                self.activity_log_text.setTextCursor(cursor)
                self.activity_log_text.setReadOnly(True)
        else:
            self.pick_roll_line_btn.setText("Draw a Roll Line")
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
        self._update_cal_line_times()

    def _toggle_edit_roll_line_mode(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot edit roll line.")
            return
        if self.geotiff_dataset_original is None:
            self._show_message("warning","No GeoTIFF", "Load a GeoTIFF first to edit a roll line.")
            return
        if len(self.roll_line_points) < 2:
            self._show_message("warning","No Roll Line", "Draw a roll line first before editing it.")
            return

        self.edit_roll_line_mode = not self.edit_roll_line_mode
        if self.edit_roll_line_mode:
            (lat1, lon1), (lat2, lon2) = self.roll_line_points
            self.roll_line_start_handle = self.ax.scatter([lon1], [lat1], c='red', s=100, zorder=10, picker=True)
            self.roll_line_end_handle = self.ax.scatter([lon2], [lat2], c='blue', s=100, zorder=10, picker=True)

            self.roll_line_pick_cid = self.canvas.mpl_connect('pick_event', self._on_roll_line_handle_pick)
            self.roll_line_motion_cid = self.canvas.mpl_connect('motion_notify_event', self._on_roll_line_handle_motion)
            self.roll_line_release_cid = self.canvas.mpl_connect('button_release_event', self._on_roll_line_handle_release)

            self.edit_roll_line_btn.setText("Click to Stop Editing Roll Line")
            self.canvas_widget.setCursor(Qt.CursorShape.SizeAllCursor)
            self.set_cal_info_text("You can drag the red (start) and blue (end) points to edit the roll line.")
        else:
            if self.roll_line_start_handle or self.roll_line_end_handle:
                self._plot_survey_plan()
                self.roll_line_start_handle = None
                self.roll_line_end_handle = None

            if hasattr(self, 'roll_line_pick_cid'):
                self.canvas.mpl_disconnect(self.roll_line_pick_cid)
            if hasattr(self, 'roll_line_motion_cid'):
                self.canvas.mpl_disconnect(self.roll_line_motion_cid)
            if hasattr(self, 'roll_line_release_cid'):
                self.canvas.mpl_disconnect(self.roll_line_release_cid)

            self.edit_roll_line_btn.setText("Edit Roll Line")
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
            self.dragging_roll_line_handle = None

            self._update_cal_line_times()

            self.set_cal_info_text("Roll line editing completed.")

    def _on_pitch_line_handle_pick(self, event):
        if not self.edit_pitch_line_mode:
            return
        if event.mouseevent.button != 1:
            return
        if event.artist == self.pitch_line_start_handle:
            self.dragging_pitch_line_handle = 'start'
        elif event.artist == self.pitch_line_end_handle:
            self.dragging_pitch_line_handle = 'end'

    def _on_pitch_line_handle_motion(self, event):
        if not self.edit_pitch_line_mode or not self.dragging_pitch_line_handle:
            return
        if event.inaxes != self.ax:
            return
        if event.xdata is None or event.ydata is None:
            return

        if self.dragging_pitch_line_handle == 'start' and self.pitch_line_start_handle is not None:
            self.pitch_line_start_handle.set_offsets([[event.xdata, event.ydata]])
            self.pitch_line_points[0] = (event.ydata, event.xdata)
        elif self.dragging_pitch_line_handle == 'end' and self.pitch_line_end_handle is not None:
            self.pitch_line_end_handle.set_offsets([[event.xdata, event.ydata]])
            self.pitch_line_points[1] = (event.ydata, event.xdata)

        if hasattr(self, 'pitch_line_edit_line') and self.pitch_line_edit_line is not None:
            self.pitch_line_edit_line.set_data(
                [self.pitch_line_points[0][1], self.pitch_line_points[1][1]],
                [self.pitch_line_points[0][0], self.pitch_line_points[1][0]]
            )
        else:
            (lat1, lon1), (lat2, lon2) = self.pitch_line_points
            self.pitch_line_edit_line, = self.ax.plot([lon1, lon2], [lat1, lat2], color='orange', linewidth=2, zorder=9)

        if hasattr(self, 'pitch_line_tooltip') and self.pitch_line_tooltip is not None:
            depth = self._get_depth_at_point(event.ydata, event.xdata)
            slope = self._calculate_slope_at_point(event.ydata, event.xdata)[1]
            if depth is not None and slope is not None:
                tooltip_text = f"Depth: {depth:.1f} m\nSlope: {slope:.1f}°"
                if self.pitch_line_tooltip_text is not None:
                    self.pitch_line_tooltip_text.set_text(tooltip_text)

        self.canvas.draw_idle()

    def _on_pitch_line_handle_release(self, event):
        if not self.edit_pitch_line_mode:
            return

        if self.dragging_pitch_line_handle:
            if hasattr(self, 'pitch_line_edit_line') and self.pitch_line_edit_line is not None:
                self.pitch_line_edit_line.remove()
                self.pitch_line_edit_line = None
            self._plot_survey_plan(preserve_view_limits=True)

            start_lat, start_lon = self.pitch_line_points[0]
            end_lat, end_lon = self.pitch_line_points[1]

            self.pitch_line_start_handle = self.ax.scatter([start_lon], [start_lat],
                                                         color='red', s=100, marker='o',
                                                         edgecolors='black', linewidth=2,
                                                         zorder=10, picker=5)
            self.pitch_line_end_handle = self.ax.scatter([end_lon], [end_lat],
                                                       color='blue', s=100, marker='o',
                                                       edgecolors='black', linewidth=2,
                                                       zorder=10, picker=5)

            self._update_cal_line_offset_from_pitch_line()
            self._draw_pitch_line_profile()
            self._update_cal_export_name_from_pitch_line()
            self._update_cal_line_times()

            try:
                if pyproj is not None:
                    geod = pyproj.Geod(ellps="WGS84")
                    (lat1, lon1), (lat2, lon2) = self.pitch_line_points
                    az12, az21, dist_m = geod.inv(lon1, lat1, lon2, lat2)
                    dist_nm = dist_m / 1852.0
                    speed_knots = float(self.cal_survey_speed_entry.text()) if self.cal_survey_speed_entry.text() else 8.0
                    speed_m_per_h = speed_knots * 1852
                    time_hours = dist_m / speed_m_per_h if speed_m_per_h > 0 else 0
                    time_minutes = time_hours * 60
                    summary = (
                        f"Updated Pitch Line:\n"
                        f"Length: {dist_m:.1f} m ({dist_nm:.3f} nm)\n"
                        f"Azimuth: {az12:.1f}°\n"
                        f"Survey Time: {time_minutes:.1f} min"
                    )
                    self.set_cal_info_text(summary)
            except Exception as e:
                self.set_cal_info_text(f"Error calculating updated pitch line summary: {e}")

            self.dragging_pitch_line_handle = None

    def _on_roll_line_handle_pick(self, event):
        if not self.edit_roll_line_mode:
            return
        if event.mouseevent.button != 1:
            return
        if event.artist == self.roll_line_start_handle:
            self.dragging_roll_line_handle = 'start'
        elif event.artist == self.roll_line_end_handle:
            self.dragging_roll_line_handle = 'end'

    def _on_roll_line_handle_motion(self, event):
        if not self.edit_roll_line_mode or not self.dragging_roll_line_handle:
            return
        if event.inaxes != self.ax:
            return
        if event.xdata is None or event.ydata is None:
            return

        if self.dragging_roll_line_handle == 'start' and self.roll_line_start_handle is not None:
            self.roll_line_start_handle.set_offsets([[event.xdata, event.ydata]])
            self.roll_line_points[0] = (event.ydata, event.xdata)
        elif self.dragging_roll_line_handle == 'end' and self.roll_line_end_handle is not None:
            self.roll_line_end_handle.set_offsets([[event.xdata, event.ydata]])
            self.roll_line_points[1] = (event.ydata, event.xdata)

        if hasattr(self, 'roll_line_edit_line') and self.roll_line_edit_line is not None:
            self.roll_line_edit_line.set_data(
                [self.roll_line_points[0][1], self.roll_line_points[1][1]],
                [self.roll_line_points[0][0], self.roll_line_points[1][0]]
            )
        else:
            (lat1, lon1), (lat2, lon2) = self.roll_line_points
            self.roll_line_edit_line, = self.ax.plot([lon1, lon2], [lat1, lat2], color='purple', linewidth=2, zorder=9)

        self.canvas.draw_idle()

    def _on_roll_line_handle_release(self, event):
        if not self.edit_roll_line_mode:
            return

        if self.dragging_roll_line_handle:
            if hasattr(self, 'roll_line_edit_line') and self.roll_line_edit_line is not None:
                self.roll_line_edit_line.remove()
                self.roll_line_edit_line = None
            self._plot_survey_plan(preserve_view_limits=True)

            start_lat, start_lon = self.roll_line_points[0]
            end_lat, end_lon = self.roll_line_points[1]

            self.roll_line_start_handle = self.ax.scatter([start_lon], [start_lat],
                                                        color='red', s=100, marker='o',
                                                        edgecolors='black', linewidth=2,
                                                        zorder=10, picker=5)
            self.roll_line_end_handle = self.ax.scatter([end_lon], [end_lat],
                                                      color='blue', s=100, marker='o',
                                                      edgecolors='black', linewidth=2,
                                                      zorder=10, picker=5)

            self._update_cal_line_times()

            try:
                if pyproj is not None:
                    geod = pyproj.Geod(ellps="WGS84")
                    (lat1, lon1), (lat2, lon2) = self.roll_line_points
                    az12, az21, dist_m = geod.inv(lon1, lat1, lon2, lat2)
                    dist_nm = dist_m / 1852.0
                    speed_knots = float(self.cal_survey_speed_entry.text()) if self.cal_survey_speed_entry.text() else 8.0
                    speed_m_per_h = speed_knots * 1852
                    time_hours = dist_m / speed_m_per_h if speed_m_per_h > 0 else 0
                    time_minutes = time_hours * 60
                    summary = (
                        f"Updated Roll Line:\n"
                        f"Length: {dist_m:.1f} m ({dist_nm:.3f} nm)\n"
                        f"Azimuth: {az12:.1f}°\n"
                        f"Survey Time: {time_minutes:.1f} min"
                    )
                    self.set_cal_info_text(summary)
            except Exception as e:
                self.set_cal_info_text(f"Error calculating updated roll line summary: {e}")

            self.dragging_roll_line_handle = None

    def _load_last_cal_import_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_cal_import_dir' in config and os.path.isdir(config['last_cal_import_dir']):
                    self.last_cal_import_dir = config['last_cal_import_dir']
        except Exception:
            pass

    def _save_last_cal_import_dir(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
            config['last_cal_import_dir'] = self.last_cal_import_dir
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _clear_calibration_lines(self):
        self.pitch_line_points = []
        self.heading_lines = []
        self.roll_line_points = []

        if hasattr(self, 'pitch_shallowest_depth_label'):
            self.pitch_shallowest_depth_label.setText("-")
        if hasattr(self, 'pitch_max_depth_label'):
            self.pitch_max_depth_label.setText("-")
        if hasattr(self, 'pitch_mean_depth_label'):
            self.pitch_mean_depth_label.setText("-")
        if hasattr(self, 'pitch_median_depth_label'):
            self.pitch_median_depth_label.setText("-")

        self._plot_survey_plan(preserve_view_limits=True)
        self._draw_pitch_line_profile()
        self._update_cal_line_times()
        self._update_cal_export_name_from_pitch_line()
        self._update_pitch_line_button_states()
        self._update_roll_line_button_states()
        self._update_line_planning_button_states()
        if hasattr(self, 'geotiff_dataset_original') and self.geotiff_dataset_original is not None:
            if hasattr(self, 'pick_pitch_line_btn'):
                self.pick_pitch_line_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")
        else:
            if hasattr(self, 'pick_pitch_line_btn'):
                self.pick_pitch_line_btn.setStyleSheet("")
        if hasattr(self, 'add_heading_lines_btn'):
            self.add_heading_lines_btn.setStyleSheet("")
        if hasattr(self, 'pick_roll_line_btn'):
            self.pick_roll_line_btn.setStyleSheet("")

    def _update_pitch_line_button_states(self):
        has_pitch_line = hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2
        if hasattr(self, 'edit_pitch_line_btn'):
            self.edit_pitch_line_btn.setEnabled(has_pitch_line)
        if hasattr(self, 'add_heading_lines_btn'):
            self.add_heading_lines_btn.setEnabled(has_pitch_line)

    def _update_roll_line_button_states(self):
        has_roll_line = hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2
        if hasattr(self, 'edit_roll_line_btn'):
            self.edit_roll_line_btn.setEnabled(has_roll_line)

    def _update_cal_line_times(self):
        try:
            if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2:
                (lat1, lon1), (lat2, lon2) = self.pitch_line_points
                if pyproj is not None:
                    geod = pyproj.Geod(ellps="WGS84")
                    _, _, dist = geod.inv(lon1, lat1, lon2, lat2)
                    speed_knots = float(self.cal_survey_speed_entry.text()) if self.cal_survey_speed_entry.text() else 8.0
                    speed_m_per_h = speed_knots * 1852
                    time_hours = dist / speed_m_per_h if speed_m_per_h > 0 else 0
                    time_minutes = time_hours * 60
        except Exception:
            pass
        try:
            if hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2:
                (lat1, lon1), (lat2, lon2) = self.roll_line_points
                if pyproj is not None:
                    geod = pyproj.Geod(ellps="WGS84")
                    _, _, dist = geod.inv(lon1, lat1, lon2, lat2)
                    speed_knots = float(self.cal_survey_speed_entry.text()) if self.cal_survey_speed_entry.text() else 8.0
                    speed_m_per_h = speed_knots * 1852
                    time_hours = dist / speed_m_per_h if speed_m_per_h > 0 else 0
                    time_minutes = time_hours * 60
        except Exception:
            pass

    def _reset_calibration_tab(self):
        self.pitch_line_points = []
        self.heading_lines = []
        self.roll_line_points = []

        self.pick_pitch_line_mode = False
        self.edit_pitch_line_mode = False
        self.pick_roll_line_mode = False
        self.edit_roll_line_mode = False

        self.pitch_line_start_handle = None
        self.pitch_line_end_handle = None
        self.dragging_pitch_line_handle = None
        self.pitch_line_temp_line = None
        self.pitch_line_annotation = None
        self.pitch_line_tooltip = None
        self.pitch_line_tooltip_text = None
        self.pitch_line_edit_line = None
        self.roll_line_start_handle = None
        self.roll_line_end_handle = None
        self.dragging_roll_line_handle = None
        self.roll_line_edit_line = None

        if hasattr(self, 'cal_survey_speed_entry'):
            self.cal_survey_speed_entry.setText("8")
        if hasattr(self, 'cal_line_offset_entry'):
            self.cal_line_offset_entry.clear()
        if hasattr(self, 'cal_export_name_entry'):
            self.cal_export_name_entry.clear()

        self._update_pitch_line_button_states()
        self._update_roll_line_button_states()

        if hasattr(self, 'geotiff_dataset_original') and self.geotiff_dataset_original is not None:
            if hasattr(self, 'pick_pitch_line_btn'):
                self.pick_pitch_line_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")
        else:
            if hasattr(self, 'pick_pitch_line_btn'):
                self.pick_pitch_line_btn.setStyleSheet("")
        if hasattr(self, 'add_heading_lines_btn'):
            self.add_heading_lines_btn.setStyleSheet("")
        if hasattr(self, 'pick_roll_line_btn'):
            self.pick_roll_line_btn.setStyleSheet("")

        self._plot_survey_plan(preserve_view_limits=True)
        if hasattr(self, '_draw_pitch_line_profile'):
            self._draw_pitch_line_profile()
        if hasattr(self, '_update_cal_line_times'):
            self._update_cal_line_times()

    def _export_cal_survey_files(self):
        try:
            import fiona
            from shapely.geometry import LineString
            from shapely.geometry import mapping
        except ImportError:
            self._show_message("warning","Missing Libraries", "fiona and shapely are required for shapefile export.")
            return
        lines = []
        line_num = 1
        if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2:
            lines.append((line_num, 'Pitch', self.pitch_line_points))
            line_num += 1
        if hasattr(self, 'heading_lines') and len(self.heading_lines) == 2:
            lines.append((line_num, 'Heading1', self.heading_lines[0]))
            line_num += 1
            lines.append((line_num, 'Heading2', self.heading_lines[1]))
            line_num += 1
        if hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2:
            lines.append((line_num, 'Roll', self.roll_line_points))
            line_num += 1
        if not lines:
            self._show_message("warning","No Data", "No calibration lines to export. Define at least one line.")
            return
        export_name = self.cal_export_name_entry.text().strip() or "calibration_survey"
        export_dir = QFileDialog.getExistingDirectory(self, "Select Export Directory", self.last_export_dir)
        if not export_dir:
            return
        self.last_export_dir = export_dir
        self._save_last_export_dir()
        try:
            csv_file_path = os.path.join(export_dir, f"{export_name}_DD.csv")
            with open(csv_file_path, 'w', newline='', encoding='utf-8') as csvfile:
                csv_writer = csv.writer(csvfile)
                csv_writer.writerow(['Line Number', 'Point Label', 'Line Name', 'Latitude', 'Longitude'])
                for num, name, pts in lines:
                    if name.lower().startswith('pitch'):
                        start_label, end_label = 'PLS', 'PLE'
                    elif name.lower().startswith('roll'):
                        start_label, end_label = 'RLS', 'RLE'
                    elif name.lower().startswith('heading'):
                        m = re.search(r'(\d+)', name)
                        n = m.group(1) if m else ''
                        start_label, end_label = f'H{n}S', f'H{n}E'
                    else:
                        start_label, end_label = 'START', 'END'
                    csv_writer.writerow([num, start_label, name, pts[0][0], pts[0][1]])
                    csv_writer.writerow([num, end_label, name, pts[1][0], pts[1][1]])
            ddm_file_path = os.path.join(export_dir, f"{export_name}_DM.csv")
            with open(ddm_file_path, 'w', newline='', encoding='utf-8') as ddmfile:
                ddm_writer = csv.writer(ddmfile)
                ddm_writer.writerow(['Line Number', 'Point Label', 'Latitude with Decimal Minutes', 'Longitude with Decimal Minutes'])
                for num, name, pts in lines:
                    if name.lower().startswith('pitch'):
                        start_label, end_label = 'PLS', 'PLE'
                    elif name.lower().startswith('roll'):
                        start_label, end_label = 'RLS', 'RLE'
                    elif name.lower().startswith('heading'):
                        m = re.search(r'(\d+)', name)
                        n = m.group(1) if m else ''
                        start_label, end_label = f'H{n}S', f'H{n}E'
                    else:
                        start_label, end_label = 'START', 'END'
                    start_lat_ddm = decimal_degrees_to_ddm(pts[0][0], is_latitude=True)
                    start_lon_ddm = decimal_degrees_to_ddm(pts[0][1], is_latitude=False)
                    end_lat_ddm = decimal_degrees_to_ddm(pts[1][0], is_latitude=True)
                    end_lon_ddm = decimal_degrees_to_ddm(pts[1][1], is_latitude=False)
                    ddm_writer.writerow([num, start_label, start_lat_ddm, start_lon_ddm])
                    ddm_writer.writerow([num, end_label, end_lat_ddm, end_lon_ddm])
            schema = {'geometry': 'LineString', 'properties': {'line_num': 'int', 'line_name': 'str'}}
            crs_epsg = 'EPSG:4326'
            features = []
            for num, name, pts in lines:
                shapely_line = LineString([(p[1], p[0]) for p in pts])
                features.append({'geometry': mapping(shapely_line), 'properties': {'line_num': num, 'line_name': name}})
            shapefile_path = os.path.join(export_dir, f"{export_name}.shp")
            with fiona.open(shapefile_path, 'w', driver='ESRI Shapefile', crs=crs_epsg, schema=schema) as collection:
                collection.writerecords(features)
            geojson_file_path = os.path.join(export_dir, f"{export_name}.geojson")
            geojson_features = []
            for num, name, pts in lines:
                geojson_features.append({
                    "type": "Feature",
                    "geometry": {"type": "LineString", "coordinates": [[pts[0][1], pts[0][0]], [pts[1][1], pts[1][0]]]},
                    "properties": {"line_num": num, "line_name": name}
                })
            with open(geojson_file_path, 'w', encoding='utf-8') as f:
                json.dump({"type": "FeatureCollection", "features": geojson_features}, f, indent=2)
            lnw_file_path = os.path.join(export_dir, f"{export_name}.lnw")
            with open(lnw_file_path, 'w') as f:
                f.write("LNW 1.0\n")
                for num, name, pts in lines:
                    try:
                        speed_knots = float(self.cal_survey_speed_entry.text()) if self.cal_survey_speed_entry.text() else 8.0
                    except Exception:
                        speed_knots = 8.0
                    depth1 = self._get_depth_at_point(pts[0][0], pts[0][1])
                    depth2 = self._get_depth_at_point(pts[1][0], pts[1][1])
                    f.write(f"{name}_001, {pts[0][0]:.6f}, {pts[0][1]:.6f}, {depth1:.1f}, {speed_knots:.1f}, 50.0, {num}, {num}\n")
                    f.write(f"{name}_002, {pts[1][0]:.6f}, {pts[1][1]:.6f}, {depth2:.1f}, {speed_knots:.1f}, 50.0, {num}, {num}\n")
            sis_file_path = os.path.join(export_dir, f"{export_name}.asciiplan")
            with open(sis_file_path, 'w') as f:
                f.write("SIS ASCII Plan\n")
                for num, name, pts in lines:
                    try:
                        speed_knots = float(self.cal_survey_speed_entry.text()) if self.cal_survey_speed_entry.text() else 8.0
                    except Exception:
                        speed_knots = 8.0
                    depth1 = self._get_depth_at_point(pts[0][0], pts[0][1])
                    depth2 = self._get_depth_at_point(pts[1][0], pts[1][1])
                    f.write(f"{name}_001, {pts[0][0]:.6f}, {pts[0][1]:.6f}, {depth1:.1f}, {speed_knots:.1f}, {num}, {num}\n")
                    f.write(f"{name}_002, {pts[1][0]:.6f}, {pts[1][1]:.6f}, {depth2:.1f}, {speed_knots:.1f}, {num}, {num}\n")
            txt_file_path = os.path.join(export_dir, f"{export_name}_DD.txt")
            with open(txt_file_path, 'w', encoding='utf-8') as f:
                for num, name, pts in lines:
                    if name.lower().startswith('pitch'):
                        start_label, end_label = 'PLS', 'PLE'
                    elif name.lower().startswith('roll'):
                        start_label, end_label = 'RLS', 'RLE'
                    elif name.lower().startswith('heading'):
                        m = re.search(r'(\d+)', name)
                        n = m.group(1) if m else ''
                        start_label, end_label = f'H{n}S', f'H{n}E'
                    else:
                        start_label, end_label = 'START', 'END'
                    f.write(f"{start_label} {pts[0][0]:.6f} {pts[0][1]:.6f}\n")
                    f.write(f"{end_label} {pts[1][0]:.6f} {pts[1][1]:.6f}\n")
            stats = self._calculate_calibration_survey_statistics()
            stats_file_path = None
            if stats:
                stats_file_path = os.path.join(export_dir, f"{export_name}_stats.txt")
                with open(stats_file_path, 'w', encoding='utf-8') as f:
                    f.write("COMPREHENSIVE CALIBRATION SURVEY STATISTICS\n")
                    f.write("=" * 50 + "\n\n")
                    f.write(f"Calibration Survey: {export_name}\n")
                    f.write(f"Export Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                    f.write(f"Total Distance: {stats['total_distance_m']:.1f} m\n")
                    f.write(f"Total Time: {stats['total_time_min']:.1f} min\n")
            try:
                params = {}
                try:
                    params['line_offset'] = float(self.cal_line_offset_entry.text())
                except Exception:
                    params['line_offset'] = None
                try:
                    params['export_name'] = self.cal_export_name_entry.text().strip()
                except Exception:
                    params['export_name'] = export_name
                try:
                    params['survey_speed'] = float(self.cal_survey_speed_entry.text()) if self.cal_survey_speed_entry.text() else 8.0
                except Exception:
                    params['survey_speed'] = 8.0
                json_metadata_path = os.path.join(export_dir, f"{export_name}_params.json")
                with open(json_metadata_path, 'w', encoding='utf-8') as f:
                    json.dump(params, f, indent=2)
            except Exception as e:
                print(f"Warning: Could not export metadata: {e}")
                json_metadata_path = None
            map_png_path = os.path.join(export_dir, f"{export_name}_map.png")
            self.figure.savefig(map_png_path, dpi=300, bbox_inches='tight', facecolor='white')
            profile_png_path = None
            if hasattr(self, 'profile_fig') and self.profile_fig is not None:
                profile_png_path = os.path.join(export_dir, f"{export_name}_profile.png")
                self.profile_fig.savefig(profile_png_path, dpi=300, bbox_inches='tight', facecolor='white')
            success_msg = f"Survey exported successfully to:\n- {os.path.basename(csv_file_path)}\n- {os.path.basename(shapefile_path)}\n- {os.path.basename(geojson_file_path)}\n- {os.path.basename(lnw_file_path)}\n- {os.path.basename(sis_file_path)}\n"
            if stats_file_path:
                success_msg += f"- {os.path.basename(stats_file_path)}\n"
            if json_metadata_path:
                success_msg += f"- {os.path.basename(json_metadata_path)}\n"
            success_msg += f"- {os.path.basename(map_png_path)}\n"
            if profile_png_path:
                success_msg += f"- {os.path.basename(profile_png_path)}\n"
            success_msg += f"in directory: {export_dir}"
            self.set_cal_info_text(success_msg, append=False)
        except Exception as e:
            self._show_message("error","Export Error", f"Failed to export calibration survey files: {e}")

    def _import_cal_survey_files(self):
        """Import calibration survey lines from CSV or GeoJSON file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Survey File to Import",
            self.last_cal_import_dir,
            "Decimal Degree CSV files (*_DD.csv);;CSV files (*.csv);;GeoJSON files (*.geojson);;JSON files (*.json);;All files (*.*)"
        )
        if not file_path:
            return
        import_dir = os.path.dirname(file_path)
        if import_dir and os.path.isdir(import_dir):
            self.last_cal_import_dir = import_dir
            self._save_last_cal_import_dir()
        try:
            self.pitch_line_points = []
            self.heading_lines = []
            self.roll_line_points = []
            file_ext = os.path.splitext(file_path)[1].lower()
            if file_ext == '.csv':
                with open(file_path, 'r', encoding='utf-8') as csvfile:
                    csv_reader = csv.DictReader(csvfile)
                    lines_data = {}
                    for row in csv_reader:
                        line_name = row.get('Line Name', '').strip()
                        point_label = row.get('Point Label', '').strip()
                        try:
                            lat = float(row.get('Latitude', 0))
                            lon = float(row.get('Longitude', 0))
                        except (ValueError, TypeError):
                            continue
                        if not line_name:
                            continue
                        if line_name not in lines_data:
                            lines_data[line_name] = {}
                        lines_data[line_name][point_label] = (lat, lon)
                    for line_name, points in lines_data.items():
                        line_name_lower = line_name.lower()
                        if line_name_lower.startswith('pitch'):
                            if 'PLS' in points and 'PLE' in points:
                                self.pitch_line_points = [points['PLS'], points['PLE']]
                        elif line_name_lower.startswith('roll'):
                            if 'RLS' in points and 'RLE' in points:
                                self.roll_line_points = [points['RLS'], points['RLE']]
                        elif line_name_lower.startswith('heading'):
                            match = re.search(r'heading\s*(\d+)', line_name_lower)
                            if match:
                                heading_num = int(match.group(1))
                                start_label = f'H{heading_num}S'
                                end_label = f'H{heading_num}E'
                                heading_idx = min(heading_num - 1, 1)
                            else:
                                start_labels = [p for p in points.keys() if p.endswith('S')]
                                end_labels = [p for p in points.keys() if p.endswith('E')]
                                if start_labels and end_labels:
                                    start_label, end_label = start_labels[0], end_labels[0]
                                    heading_idx = 0
                                else:
                                    continue
                            if start_label in points and end_label in points:
                                while len(self.heading_lines) < 2:
                                    self.heading_lines.append([])
                                if heading_idx < 2:
                                    self.heading_lines[heading_idx] = [points[start_label], points[end_label]]
            elif file_ext in ['.geojson', '.json']:
                with open(file_path, 'r', encoding='utf-8') as f:
                    geojson_data = json.load(f)
                if not isinstance(geojson_data, dict):
                    self._show_message("error","Import Error", "Invalid GeoJSON format.")
                    return
                if geojson_data.get('type') == 'FeatureCollection':
                    features = geojson_data.get('features', [])
                elif geojson_data.get('type') == 'Feature':
                    features = [geojson_data]
                else:
                    features = []
                for feature in features:
                    if not isinstance(feature, dict) or feature.get('type') != 'Feature':
                        continue
                    geometry = feature.get('geometry', {})
                    if not isinstance(geometry, dict) or geometry.get('type') != 'LineString':
                        continue
                    coordinates = geometry.get('coordinates', [])
                    if not isinstance(coordinates, list) or len(coordinates) != 2:
                        continue
                    point1 = (coordinates[0][1], coordinates[0][0])
                    point2 = (coordinates[1][1], coordinates[1][0])
                    properties = feature.get('properties', {}) or {}
                    line_name = str(properties.get('line_name', '')).strip().lower()
                    if line_name.startswith('pitch'):
                        self.pitch_line_points = [point1, point2]
                    elif line_name.startswith('roll'):
                        self.roll_line_points = [point1, point2]
                    elif line_name.startswith('heading'):
                        match = re.search(r'heading\s*(\d+)', line_name)
                        heading_idx = min(int(match.group(1)) - 1, 1) if match else 0
                        while len(self.heading_lines) < 2:
                            self.heading_lines.append([])
                        if heading_idx < 2:
                            self.heading_lines[heading_idx] = [point1, point2]
            else:
                self._show_message("error","Import Error", f"Unsupported file format: {file_ext}")
                return
            params = None
            base_path = os.path.splitext(file_path)[0]
            metadata_path = os.path.join(os.path.dirname(file_path), os.path.basename(base_path) + "_params.json")
            if os.path.exists(metadata_path):
                try:
                    with open(metadata_path, 'r', encoding='utf-8') as f:
                        params = json.load(f)
                except Exception:
                    pass
            if params and isinstance(params, dict):
                try:
                    if params.get('line_offset') is not None:
                        self.cal_line_offset_entry.setText(str(params['line_offset']))
                    if params.get('export_name'):
                        self.cal_export_name_entry.setText(params['export_name'])
                    if params.get('survey_speed') is not None:
                        self.cal_survey_speed_entry.setText(str(params['survey_speed']))
                except Exception:
                    pass
            if not params or not isinstance(params, dict) or params.get('line_offset') is None:
                self._update_cal_line_offset_from_pitch_line()
            if not params or not isinstance(params, dict) or not params.get('export_name'):
                self._update_cal_export_name_from_pitch_line()
            self._update_pitch_line_button_states()
            self._update_roll_line_button_states()
            if hasattr(self, 'add_heading_lines_btn'):
                has_heading_lines = any(len(line) == 2 for line in self.heading_lines)
                self.add_heading_lines_btn.setEnabled(has_heading_lines)
            for btn_attr in ('pick_pitch_line_btn', 'add_heading_lines_btn', 'pick_roll_line_btn'):
                if hasattr(self, btn_attr):
                    getattr(self, btn_attr).setStyleSheet("")
            self._update_cal_line_times()
            self._plot_survey_plan(preserve_view_limits=True)
            if hasattr(self, '_draw_current_profile'):
                self._draw_current_profile()
            elif hasattr(self, '_draw_pitch_line_profile'):
                self._draw_pitch_line_profile()
            imported_lines = []
            if len(self.pitch_line_points) == 2:
                imported_lines.append("Pitch line")
            if len(self.roll_line_points) == 2:
                imported_lines.append("Roll line")
            for i, heading_line in enumerate(self.heading_lines):
                if len(heading_line) == 2:
                    imported_lines.append(f"Heading line {i+1}")
            if imported_lines:
                self.set_cal_info_text(f"Successfully imported: {', '.join(imported_lines)}")
            else:
                self._show_message("warning","Import Warning", "No valid calibration lines found in the selected file.")
        except Exception as e:
            self._show_message("error","Import Error", f"Failed to import calibration survey files: {e}")
            import traceback
            traceback.print_exc()

    def _calculate_calibration_survey_statistics(self):
        """Calculate comprehensive calibration survey statistics including pitch, roll, heading lines and travel distances."""
        try:
            if pyproj is None:
                return None
            geod = pyproj.Geod(ellps="WGS84")
            speed_knots = float(self.cal_survey_speed_entry.text()) if self.cal_survey_speed_entry.text() else 8.0
            speed_m_per_h = speed_knots * 1852

            stats = {
                'pitch_line_distance_m': 0.0,
                'pitch_line_time_min': 0.0,
                'roll_line_distance_m': 0.0,
                'roll_line_time_min': 0.0,
                'heading1_distance_m': 0.0,
                'heading1_time_min': 0.0,
                'heading2_distance_m': 0.0,
                'heading2_time_min': 0.0,
                'travel_pitch_to_roll_m': 0.0,
                'travel_pitch_to_roll_min': 0.0,
                'travel_roll_to_heading1_m': 0.0,
                'travel_roll_to_heading1_min': 0.0,
                'travel_heading1_to_heading2_m': 0.0,
                'travel_heading1_to_heading2_min': 0.0,
                'total_distance_m': 0.0,
                'total_time_min': 0.0,
                'total_time_hours': 0.0
            }

            if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2:
                (lat1, lon1), (lat2, lon2) = self.pitch_line_points
                _, _, pitch_distance = geod.inv(lon1, lat1, lon2, lat2)
                stats['pitch_line_distance_m'] = pitch_distance
                pitch_time_hours = pitch_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                stats['pitch_line_time_min'] = pitch_time_hours * 60 * 2

            if hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2:
                (lat1, lon1), (lat2, lon2) = self.roll_line_points
                _, _, roll_distance = geod.inv(lon1, lat1, lon2, lat2)
                stats['roll_line_distance_m'] = roll_distance
                roll_time_hours = roll_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                stats['roll_line_time_min'] = roll_time_hours * 60 * 2

            if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 1:
                (lat1, lon1), (lat2, lon2) = self.heading_lines[0]
                _, _, heading1_distance = geod.inv(lon1, lat1, lon2, lat2)
                stats['heading1_distance_m'] = heading1_distance
                heading1_time_hours = heading1_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                stats['heading1_time_min'] = heading1_time_hours * 60

            if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 2:
                (lat1, lon1), (lat2, lon2) = self.heading_lines[1]
                _, _, heading2_distance = geod.inv(lon1, lat1, lon2, lat2)
                stats['heading2_distance_m'] = heading2_distance
                heading2_time_hours = heading2_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                stats['heading2_time_min'] = heading2_time_hours * 60

            if (hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2 and
                hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2):
                pitch_start_lat, pitch_start_lon = self.pitch_line_points[0]
                roll_start_lat, roll_start_lon = self.roll_line_points[0]
                _, _, travel_distance = geod.inv(pitch_start_lon, pitch_start_lat, roll_start_lon, roll_start_lat)
                stats['travel_pitch_to_roll_m'] = travel_distance
                travel_time_hours = travel_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                stats['travel_pitch_to_roll_min'] = travel_time_hours * 60

            if (hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2 and
                hasattr(self, 'heading_lines') and len(self.heading_lines) >= 1):
                roll_start_lat, roll_start_lon = self.roll_line_points[0]
                heading1_start_lat, heading1_start_lon = self.heading_lines[0][0]
                _, _, travel_distance = geod.inv(roll_start_lon, roll_start_lat, heading1_start_lon, heading1_start_lat)
                stats['travel_roll_to_heading1_m'] = travel_distance
                travel_time_hours = travel_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                stats['travel_roll_to_heading1_min'] = travel_time_hours * 60

            if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 2:
                heading1_end_lat, heading1_end_lon = self.heading_lines[0][1]
                heading2_start_lat, heading2_start_lon = self.heading_lines[1][0]
                _, _, travel_distance = geod.inv(heading1_end_lon, heading1_end_lat, heading2_start_lon, heading2_start_lat)
                stats['travel_heading1_to_heading2_m'] = travel_distance
                travel_time_hours = travel_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                stats['travel_heading1_to_heading2_min'] = travel_time_hours * 60

            stats['total_distance_m'] = (
                stats['pitch_line_distance_m'] * 2 +
                stats['roll_line_distance_m'] * 2 +
                stats['heading1_distance_m'] +
                stats['heading2_distance_m'] +
                stats['travel_pitch_to_roll_m'] +
                stats['travel_roll_to_heading1_m'] +
                stats['travel_heading1_to_heading2_m']
            )
            stats['total_time_min'] = (
                stats['pitch_line_time_min'] +
                stats['roll_line_time_min'] +
                stats['heading1_time_min'] +
                stats['heading2_time_min'] +
                stats['travel_pitch_to_roll_min'] +
                stats['travel_roll_to_heading1_min'] +
                stats['travel_heading1_to_heading2_min']
            )
            stats['total_time_hours'] = stats['total_time_min'] / 60.0
            stats['total_distance_km'] = stats['total_distance_m'] / 1000.0
            stats['total_distance_nm'] = stats['total_distance_m'] / 1852.0
            stats['pitch_line_distance_km'] = stats['pitch_line_distance_m'] / 1000.0
            stats['pitch_line_distance_nm'] = stats['pitch_line_distance_m'] / 1852.0
            stats['roll_line_distance_km'] = stats['roll_line_distance_m'] / 1000.0
            stats['roll_line_distance_nm'] = stats['roll_line_distance_m'] / 1852.0
            stats['heading1_distance_km'] = stats['heading1_distance_m'] / 1000.0
            stats['heading1_distance_nm'] = stats['heading1_distance_m'] / 1852.0
            stats['heading2_distance_km'] = stats['heading2_distance_m'] / 1000.0
            stats['heading2_distance_nm'] = stats['heading2_distance_m'] / 1852.0
            stats['travel_pitch_to_roll_km'] = stats['travel_pitch_to_roll_m'] / 1000.0
            stats['travel_pitch_to_roll_nm'] = stats['travel_pitch_to_roll_m'] / 1852.0
            stats['travel_roll_to_heading1_km'] = stats['travel_roll_to_heading1_m'] / 1000.0
            stats['travel_roll_to_heading1_nm'] = stats['travel_roll_to_heading1_m'] / 1852.0
            stats['travel_heading1_to_heading2_km'] = stats['travel_heading1_to_heading2_m'] / 1000.0
            stats['travel_heading1_to_heading2_nm'] = stats['travel_heading1_to_heading2_m'] / 1852.0

            return stats

        except Exception as e:
            print(f"Error calculating calibration survey statistics: {e}")
            return None

    def _show_calibration_statistics(self):
        """Display comprehensive calibration survey statistics in a custom dialog."""
        stats = self._calculate_calibration_survey_statistics()
        if not stats:
            self._show_message("warning","Statistics Error", "Unable to calculate calibration survey statistics. Please ensure all required lines are drawn.")
            return

        stats_text = "COMPREHENSIVE CALIBRATION SURVEY STATISTICS\n"
        stats_text += "=" * 50 + "\n\n"
        try:
            speed_knots = float(self.cal_survey_speed_entry.text()) if self.cal_survey_speed_entry.text() else 8.0
            stats_text += f"Survey Speed: {speed_knots} knots\n\n"
        except Exception:
            stats_text += "Survey Speed: 8.0 knots (default)\n\n"

        stats_text += "SURVEY LINE DISTANCES AND TIMES\n"
        stats_text += "-" * 35 + "\n"
        if stats['pitch_line_distance_m'] > 0:
            stats_text += f"Pitch Line (2 passes):\n"
            if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2 and pyproj is not None:
                try:
                    geod = pyproj.Geod(ellps="WGS84")
                    (lat1, lon1), (lat2, lon2) = self.pitch_line_points
                    fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)
                    stats_text += f"  Heading: {fwd_az % 360:.1f}°\n"
                    stats_text += f"  Reciprocal Heading: {back_az % 360:.1f}°\n"
                except Exception:
                    pass
            stats_text += f"  Distance (per pass): {stats['pitch_line_distance_m']:.1f} m ({stats['pitch_line_distance_km']:.3f} km, {stats['pitch_line_distance_nm']:.3f} nm)\n"
            stats_text += f"  Time (total): {stats['pitch_line_time_min']:.1f} min\n\n"
        if stats['roll_line_distance_m'] > 0:
            stats_text += f"Roll Line (2 passes):\n"
            if hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2 and pyproj is not None:
                try:
                    geod = pyproj.Geod(ellps="WGS84")
                    (lat1, lon1), (lat2, lon2) = self.roll_line_points
                    fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)
                    stats_text += f"  Heading: {fwd_az % 360:.1f}°\n"
                    stats_text += f"  Reciprocal Heading: {back_az % 360:.1f}°\n"
                except Exception:
                    pass
            stats_text += f"  Distance (per pass): {stats['roll_line_distance_m']:.1f} m\n"
            stats_text += f"  Time (total): {stats['roll_line_time_min']:.1f} min\n\n"
        if stats['heading1_distance_m'] > 0:
            stats_text += f"Heading Line 1: {stats['heading1_distance_m']:.1f} m, {stats['heading1_time_min']:.1f} min\n\n"
        if stats['heading2_distance_m'] > 0:
            stats_text += f"Heading Line 2: {stats['heading2_distance_m']:.1f} m, {stats['heading2_time_min']:.1f} min\n\n"
        stats_text += "TRAVEL DISTANCES AND TIMES\n"
        stats_text += "-" * 30 + "\n"
        if stats['travel_pitch_to_roll_m'] > 0:
            stats_text += f"Pitch → Roll: {stats['travel_pitch_to_roll_m']:.1f} m, {stats['travel_pitch_to_roll_min']:.1f} min\n\n"
        if stats['travel_roll_to_heading1_m'] > 0:
            stats_text += f"Roll → Heading1: {stats['travel_roll_to_heading1_m']:.1f} m\n\n"
        if stats['travel_heading1_to_heading2_m'] > 0:
            stats_text += f"Heading1 → Heading2: {stats['travel_heading1_to_heading2_m']:.1f} m\n\n"
        stats_text += "TOTAL SURVEY SUMMARY\n"
        stats_text += "-" * 25 + "\n"
        stats_text += f"Total Distance: {stats['total_distance_m']:.1f} m ({stats['total_distance_km']:.3f} km, {stats['total_distance_nm']:.3f} nm)\n"
        stats_text += f"Total Time: {stats['total_time_min']:.1f} min ({stats['total_time_hours']:.2f} hr)\n"
        if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2:
            pitch_start = self.pitch_line_points[0]
            pitch_end = self.pitch_line_points[1]
            stats_text += f"\nPitch Start: {decimal_degrees_to_ddm(pitch_start[0], True)}, {decimal_degrees_to_ddm(pitch_start[1], False)}\n"
            stats_text += f"Pitch End: {decimal_degrees_to_ddm(pitch_end[0], True)}, {decimal_degrees_to_ddm(pitch_end[1], False)}\n"
        if hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2:
            roll_start = self.roll_line_points[0]
            roll_end = self.roll_line_points[1]
            stats_text += f"Roll Start: {decimal_degrees_to_ddm(roll_start[0], True)}, {decimal_degrees_to_ddm(roll_start[1], False)}\n"
            stats_text += f"Roll End: {decimal_degrees_to_ddm(roll_end[0], True)}, {decimal_degrees_to_ddm(roll_end[1], False)}\n"

        show_statistics_dialog(self, "Calibration Planning Statistics", stats_text)

    def _update_cal_line_offset_from_pitch_line(self):
        """Calculate heading line offset from median depth along pitch line."""
        if (self.geotiff_data_array is None or self.geotiff_extent is None or
            not hasattr(self, 'pitch_line_points') or len(self.pitch_line_points) != 2):
            self.cal_line_offset_entry.clear()
            return
        (lat1, lon1), (lat2, lon2) = self.pitch_line_points
        lats = np.linspace(lat1, lat2, 100)
        lons = np.linspace(lon1, lon2, 100)
        left, right, bottom, top = tuple(self.geotiff_extent)
        nrows, ncols = self.geotiff_data_array.shape
        rows = ((top - lats) / (top - bottom) * (nrows - 1)).clip(0, nrows - 1)
        cols = ((lons - left) / (right - left) * (ncols - 1)).clip(0, ncols - 1)
        elevations = []
        for r, c in zip(rows, cols):
            ir, ic = int(round(r)), int(round(c))
            elevations.append(self.geotiff_data_array[ir, ic])
        elevations = np.array(elevations)
        if np.any(~np.isnan(elevations)):
            median_val = np.nanmedian(elevations)
            mean_val = np.nanmean(elevations)
            min_val = np.nanmin(elevations)
            max_val = np.nanmax(elevations)
            valid_count = np.sum(~np.isnan(elevations))
            offset_value = abs(median_val)
            self.cal_line_offset_entry.setText(f"{offset_value:.2f}")
            if hasattr(self, 'pitch_shallowest_depth_label'):
                self.pitch_shallowest_depth_label.setText(f"{abs(max_val):.2f}")
            if hasattr(self, 'pitch_max_depth_label'):
                self.pitch_max_depth_label.setText(f"{abs(min_val):.2f}")
            if hasattr(self, 'pitch_mean_depth_label'):
                self.pitch_mean_depth_label.setText(f"{abs(mean_val):.2f}")
            if hasattr(self, 'pitch_median_depth_label'):
                self.pitch_median_depth_label.setText(f"{abs(median_val):.2f}")
            if hasattr(self, 'set_cal_info_text'):
                self.set_cal_info_text(
                    f"Heading Line Offset calculated: {offset_value:.2f} m (median of {valid_count} points along pitch line). "
                    f"Depth range: {abs(max_val):.2f} m (shallowest) to {abs(min_val):.2f} m (deepest). "
                    f"Mean depth: {abs(mean_val):.2f} m, Median depth: {abs(median_val):.2f} m.",
                    append=False
                )
        else:
            self.cal_line_offset_entry.setText("-")
            if hasattr(self, 'pitch_shallowest_depth_label'):
                self.pitch_shallowest_depth_label.setText("-")
            if hasattr(self, 'pitch_max_depth_label'):
                self.pitch_max_depth_label.setText("-")
            if hasattr(self, 'pitch_mean_depth_label'):
                self.pitch_mean_depth_label.setText("-")
            if hasattr(self, 'pitch_median_depth_label'):
                self.pitch_median_depth_label.setText("-")
            if hasattr(self, 'set_cal_info_text'):
                self.set_cal_info_text(
                    "Heading Line Offset: Could not calculate (no valid elevation data along pitch line).",
                    append=False
                )

    def _add_heading_lines_from_pitch_line(self):
        """Add heading lines north and south of pitch line at offset distance."""
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded. Cannot add heading lines.")
            return
        if self.geotiff_dataset_original is None:
            self._show_message("warning", "No GeoTIFF", "Load a GeoTIFF first.")
            return
        if not hasattr(self, 'pitch_line_points') or len(self.pitch_line_points) != 2:
            self._show_message("warning", "No Pitch Line", "Pick a pitch line first.")
            return
        try:
            offset_val = float(self.cal_line_offset_entry.text())
            if offset_val <= 0:
                raise ValueError
        except Exception:
            self._show_message("warning", "Invalid Offset", "Line Offset must be a positive number.")
            return
        (lat1, lon1), (lat2, lon2) = self.pitch_line_points
        geod = pyproj.Geod(ellps="WGS84")
        az12, az21, dist = geod.inv(lon1, lat1, lon2, lat2)
        perp_az_north = (az12 + 90) % 360
        perp_az_south = (az12 - 90) % 360
        n1_lon, n1_lat, _ = geod.fwd(lon1, lat1, perp_az_north, offset_val)
        n2_lon, n2_lat, _ = geod.fwd(lon2, lat2, perp_az_north, offset_val)
        s1_lon, s1_lat, _ = geod.fwd(lon1, lat1, perp_az_south, offset_val)
        s2_lon, s2_lat, _ = geod.fwd(lon2, lat2, perp_az_south, offset_val)
        self.heading_lines = [
            [(n1_lat, n1_lon), (n2_lat, n2_lon)],
            [(s1_lat, s1_lon), (s2_lat, s2_lon)]
        ]
        self._plot_survey_plan(preserve_view_limits=True)
        if hasattr(self, 'add_heading_lines_btn'):
            self.add_heading_lines_btn.setStyleSheet("")
        if hasattr(self, 'set_cal_info_text'):
            self.set_cal_info_text("Heading lines have been added north and south of the pitch line.")

    def _update_cal_export_name_from_pitch_line(self):
        """Update calibration export name from pitch line and offset."""
        if not hasattr(self, 'pitch_line_points') or len(self.pitch_line_points) != 2:
            return
        try:
            offset_val = float(self.cal_line_offset_entry.text())
        except Exception:
            return
        (lat1, lon1), (lat2, lon2) = self.pitch_line_points
        try:
            geod = pyproj.Geod(ellps="WGS84")
            az12, az21, dist = geod.inv(lon1, lat1, lon2, lat2)
            heading = int(round(az12)) % 360
        except Exception:
            heading = 0
        export_name = f"Cal_{int(round(offset_val))}m_{heading}deg"
        self.cal_export_name_entry.setText(export_name)

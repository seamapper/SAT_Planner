"""
Calibration tab: pitch/roll/heading lines, pick/edit modes, handle events,
export/import cal survey files, set_cal_info_text, reset, statistics.
"""
import os
import csv
import json
import re
import datetime

from PyQt6.QtWidgets import (QFileDialog, QDialog, QVBoxLayout, QHBoxLayout, 
                             QLabel, QComboBox, QPushButton, QWidget, QLineEdit, QSpinBox, QCheckBox)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QTextCursor, QColor, QTextCharFormat

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from sat_planner.constants import GEOSPATIAL_LIBS_AVAILABLE, pyproj
from sat_planner import decimal_degrees_to_ddm, export_utils
from sat_planner.utils_ui import show_statistics_dialog

try:
    from .survey_parsers_mixin import UTMZoneDialog
except ImportError:
    UTMZoneDialog = None

# Combo index to assignment type: 0=unset, 1=pitch, 2=roll, 3=heading1, 4=heading2
_CAL_ASSIGNMENT_TO_COMBO = {'pitch': 1, 'roll': 2, 'heading1': 3, 'heading2': 4}


def _suggest_calibration_assignments_from_metadata(file_path, file_basename, num_lines):
    """If the file has point labels (PLS/PLE, RLS/RLE, H1S/H1E, H2S/H2E), return suggested {line_idx: role} or None.
    Only for *_DDD.txt, *_DMM.txt, *_DMS.txt. Requires exactly 4 lines."""
    if num_lines != 4:
        return None
    file_basename_lower = file_basename.lower()
    if not (file_basename_lower.endswith('_ddd.txt') or file_basename_lower.endswith('_dmm.txt') or file_basename_lower.endswith('_dms.txt')):
        return None
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except Exception:
        return None
    lines = [line.strip() for line in content.replace('\r\n', '\n').replace('\r', '\n').split('\n') if line.strip()]
    # Each calibration line = 2 rows. First row's first token is start label (PLS, RLS, H1S, H2S).
    if len(lines) < 8:
        return None
    def label_to_role(first_token):
        t = (first_token or '').upper()
        if t in ('PLS', 'PLE'):
            return 'pitch'
        if t in ('RLS', 'RLE'):
            return 'roll'
        if t in ('H1S', 'H1E'):
            return 'heading1'
        if t in ('H2S', 'H2E'):
            return 'heading2'
        return None
    suggestion = {}
    for line_idx in range(4):
        row_idx = line_idx * 2
        if row_idx >= len(lines):
            return None
        parts = lines[row_idx].split()
        if not parts:
            return None
        role = label_to_role(parts[0])
        if role is None:
            return None
        suggestion[line_idx] = role
    if len(suggestion) != 4 or len(set(suggestion.values())) != 4:
        return None
    return suggestion


def _suggest_calibration_assignments_from_geometry(imported_lines):
    """Suggest calibration line roles from geometry: 3 parallel lines (Pitch + Heading1 + Heading2), 1 Roll.
    Pitch = spatially in the middle between the two heading lines. Heading1 vs Heading2 = file order.
    Returns {line_idx: 'pitch'|'roll'|'heading1'|'heading2'} or None."""
    if not imported_lines or len(imported_lines) != 4 or not GEOSPATIAL_LIBS_AVAILABLE or pyproj is None:
        return None
    try:
        geod = pyproj.Geod(ellps="WGS84")
    except Exception:
        return None
    # Compute for each line: length (m), orientation (0-180), midpoint (lat, lon)
    orientations = []
    lengths = []
    midpoints = []
    for line in imported_lines:
        if len(line) != 2:
            return None
        (lat1, lon1), (lat2, lon2) = line[0], line[1]
        _, _, dist = geod.inv(lon1, lat1, lon2, lat2)
        fwd_az, _, _ = geod.inv(lon1, lat1, lon2, lat2)
        orient = fwd_az % 180.0
        mid = ((lat1 + lat2) / 2.0, (lon1 + lon2) / 2.0)
        orientations.append(orient)
        lengths.append(dist)
        midpoints.append(mid)
    orientations = np.array(orientations)
    # Roll = the line whose orientation differs most from the median (the one not parallel to the other 3)
    median_orient = np.median(orientations)
    diffs = np.abs(orientations - median_orient)
    diffs = np.minimum(diffs, 180.0 - diffs)
    roll_idx = int(np.argmax(diffs))
    parallel_indices = [i for i in range(4) if i != roll_idx]
    if len(parallel_indices) != 3:
        return None
    # Common orientation of the 3 parallel lines (average)
    theta = float(np.mean(orientations[parallel_indices]))
    theta_rad = np.deg2rad(theta)
    # Across direction = perpendicular to line direction (theta + 90)
    across_rad = theta_rad + np.pi / 2.0
    across_lat = np.cos(across_rad)
    across_lon = np.sin(across_rad)
    centroid_lat = np.mean([midpoints[i][0] for i in parallel_indices])
    centroid_lon = np.mean([midpoints[i][1] for i in parallel_indices])
    # Across position: project each midpoint onto the across direction (approx for ordering)
    across_positions = []
    for i in parallel_indices:
        lat, lon = midpoints[i]
        pos = (lat - centroid_lat) * across_lat + (lon - centroid_lon) * across_lon * np.cos(np.deg2rad(centroid_lat))
        across_positions.append((i, pos))
    across_positions.sort(key=lambda x: x[1])
    # Middle (index 1) = Pitch, outer two = Heading1 and Heading2 by file order
    pitch_idx = across_positions[1][0]
    heading_indices = [across_positions[0][0], across_positions[2][0]]
    heading_indices.sort()
    heading1_idx, heading2_idx = heading_indices[0], heading_indices[1]
    return {
        pitch_idx: 'pitch',
        roll_idx: 'roll',
        heading1_idx: 'heading1',
        heading2_idx: 'heading2',
    }


class ReverseLineDirectionDialog(QDialog):
    """Dialog to choose which calibration line(s) to reverse (flip start/end points)."""
    LINE_KEYS = ['pitch', 'roll', 'heading1', 'heading2']
    LINE_LABELS = {'pitch': 'Pitch Line', 'roll': 'Roll Line', 'heading1': 'Heading Line 1', 'heading2': 'Heading Line 2'}

    def __init__(self, parent, available_lines):
        super().__init__(parent)
        self.setWindowTitle("Reverse Line Direction")
        self.setModal(True)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Select which line(s) to reverse (start and end points will be swapped):"))
        self.checkboxes = {}
        for key in self.LINE_KEYS:
            if key not in available_lines:
                continue
            cb = QCheckBox(self.LINE_LABELS[key])
            self.checkboxes[key] = cb
            layout.addWidget(cb)
        if not self.checkboxes:
            layout.addWidget(QLabel("No calibration lines available to reverse."))
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        btn_layout.addWidget(cancel_btn)
        ok_btn = QPushButton("Reverse")
        ok_btn.clicked.connect(self._on_reverse)
        btn_layout.addWidget(ok_btn)
        layout.addLayout(btn_layout)

    def _on_reverse(self):
        from PyQt6.QtWidgets import QMessageBox
        if not any(cb.isChecked() for cb in self.checkboxes.values()):
            QMessageBox.warning(self, "Reverse Line Direction", "Select at least one line to reverse.")
            return
        self.accept()

    def get_selected_lines(self):
        """Return list of line keys that are checked."""
        return [key for key, cb in self.checkboxes.items() if cb.isChecked()]


class LineAssignmentDialog(QDialog):
    """Dialog for assigning imported lines to calibration line types (pitch, roll, heading)."""
    
    def __init__(self, parent, imported_lines, suggested_assignments=None):
        """
        Args:
            parent: Parent widget
            imported_lines: List of 4 lines, each as [(lat1, lon1), (lat2, lon2)]
            suggested_assignments: Optional dict {line_idx: 'pitch'|'roll'|'heading1'|'heading2'} to pre-fill combos
        """
        super().__init__(parent)
        self.setWindowTitle("Assign Calibration Lines")
        self.setMinimumSize(900, 700)
        self.setModal(True)
        
        self.imported_lines = imported_lines
        self.assignments = {}  # Will store {line_idx: assignment_type}
        self._suggested_assignments = suggested_assignments
        
        # Create main layout
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(10)
        
        # Instructions label
        instructions = QLabel("Assign each imported line to a calibration line type:")
        instructions.setStyleSheet("font-weight: bold;")
        main_layout.addWidget(instructions)
        
        # Create map figure and canvas
        self.figure = Figure(figsize=(8, 6))
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self.figure)
        main_layout.addWidget(self.canvas)
        
        # Plot all lines on the map
        self._plot_lines()
        
        # Assignment controls
        assignment_widget = QWidget()
        assignment_layout = QVBoxLayout(assignment_widget)
        assignment_layout.setContentsMargins(5, 5, 5, 5)
        
        assignment_label = QLabel("Line Assignments:")
        assignment_label.setStyleSheet("font-weight: bold;")
        assignment_layout.addWidget(assignment_label)
        
        # Create dropdowns for each line
        self.comboboxes = []
        line_colors = ['red', 'blue', 'green', 'purple', 'orange', 'cyan', 'magenta', 'yellow']
        num_lines = len(imported_lines)
        for i in range(num_lines):
            line_widget = QWidget()
            line_layout = QHBoxLayout(line_widget)
            line_layout.setContentsMargins(0, 0, 0, 0)
            
            line_label = QLabel(f"Line {i+1}:")
            line_label.setMinimumWidth(60)
            line_layout.addWidget(line_label)
            
            combo = QComboBox()
            combo.addItems(["", "Pitch Line", "Roll Line", "Heading Line 1", "Heading Line 2"])
            # Use a closure to properly capture the line index
            def make_handler(line_idx):
                return lambda idx: self._on_assignment_changed(line_idx, idx)
            combo.currentIndexChanged.connect(make_handler(i))
            self.comboboxes.append(combo)
            line_layout.addWidget(combo)
            
            # Color indicator
            color_label = QLabel("●")
            color_style = line_colors[i % len(line_colors)]  # Cycle through colors if more than 4 lines
            color_label.setStyleSheet(f"color: {color_style}; font-size: 16pt;")
            color_label.setMinimumWidth(30)
            line_layout.addWidget(color_label)
            
            assignment_layout.addWidget(line_widget)
        
        # Apply suggestion if provided (pre-fill combos)
        if self._suggested_assignments:
            for line_idx, role in self._suggested_assignments.items():
                if 0 <= line_idx < len(self.comboboxes) and role in _CAL_ASSIGNMENT_TO_COMBO:
                    self.comboboxes[line_idx].setCurrentIndex(_CAL_ASSIGNMENT_TO_COMBO[role])
                    self.assignments[line_idx] = role
        
        main_layout.addWidget(assignment_widget)
        
        # Buttons
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
        """Plot all imported lines on the map with different colors."""
        self.ax.clear()
        
        line_colors = ['red', 'blue', 'green', 'purple', 'orange', 'cyan', 'magenta', 'yellow']
        all_lats = []
        all_lons = []
        
        for i, line in enumerate(self.imported_lines):
            if len(line) == 2:
                (lat1, lon1), (lat2, lon2) = line
                self.ax.plot([lon1, lon2], [lat1, lat2], 
                           color=line_colors[i], linewidth=2, 
                           label=f'Line {i+1}', zorder=10)
                # Add label at midpoint
                mid_lat = (lat1 + lat2) / 2
                mid_lon = (lon1 + lon2) / 2
                self.ax.text(mid_lon, mid_lat, f'Line {i+1}', 
                           fontsize=10, ha='center', va='center',
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.7),
                           zorder=11)
                all_lats.extend([lat1, lat2])
                all_lons.extend([lon1, lon2])
        
        if all_lats and all_lons:
            # Set axis limits with padding
            lat_range = max(all_lats) - min(all_lats)
            lon_range = max(all_lons) - min(all_lons)
            padding_lat = max(lat_range * 0.1, 0.01)
            padding_lon = max(lon_range * 0.1, 0.01)
            
            self.ax.set_xlim(min(all_lons) - padding_lon, max(all_lons) + padding_lon)
            self.ax.set_ylim(min(all_lats) - padding_lat, max(all_lats) + padding_lat)
        
        self.ax.set_xlabel('Longitude')
        self.ax.set_ylabel('Latitude')
        self.ax.grid(True, alpha=0.3)
        self.ax.legend(loc='upper right')
        self.figure.tight_layout()
        self.canvas.draw()
    
    def _on_assignment_changed(self, line_idx, combo_idx):
        """Handle assignment change (optional: could highlight line on map)."""
        # Map combo index to assignment type
        assignment_map = {
            0: None,
            1: 'pitch',
            2: 'roll',
            3: 'heading1',
            4: 'heading2'
        }
        self.assignments[line_idx] = assignment_map.get(combo_idx)
    
    def _on_ok_clicked(self):
        """Validate assignments and accept dialog if valid."""
        # Check that exactly 4 lines are assigned (pitch, roll, heading1, heading2)
        assigned_types = []
        for i, combo in enumerate(self.comboboxes):
            if combo.currentIndex() == 0:  # Empty selection
                continue  # Allow unassigned lines, but we need exactly 4 assigned
            assignment_type = ['', 'pitch', 'roll', 'heading1', 'heading2'][combo.currentIndex()]
            if assignment_type in assigned_types:
                self._show_error(f"Each line type can only be assigned once. '{assignment_type}' is assigned multiple times.")
                return
            assigned_types.append(assignment_type)
            self.assignments[i] = assignment_type
        
        # Validate that we have exactly 4 assignments (pitch, roll, heading1, heading2)
        if len(assigned_types) != 4:
            self._show_error("Please assign exactly 4 lines: Pitch Line, Roll Line, Heading Line 1, and Heading Line 2.")
            return
        
        if 'pitch' not in assigned_types or 'roll' not in assigned_types or 'heading1' not in assigned_types or 'heading2' not in assigned_types:
            self._show_error("Please assign all required line types: Pitch Line, Roll Line, Heading Line 1, and Heading Line 2.")
            return
        
        self.accept()
    
    def _show_error(self, message):
        """Show error message to user."""
        from PyQt6.QtWidgets import QMessageBox
        msg = QMessageBox(self)
        msg.setIcon(QMessageBox.Icon.Warning)
        msg.setWindowTitle("Assignment Error")
        msg.setText(message)
        msg.exec()
    
    def get_assignments(self):
        """Return the assignments as a dict mapping line indices to assignment types."""
        return self.assignments.copy()


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
        self._update_reverse_line_button_states()

    def _update_roll_line_button_states(self):
        has_roll_line = hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2
        if hasattr(self, 'edit_roll_line_btn'):
            self.edit_roll_line_btn.setEnabled(has_roll_line)
        self._update_reverse_line_button_states()

    def _update_reverse_line_button_states(self):
        has_pitch = hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2
        has_roll = hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2
        has_h1 = hasattr(self, 'heading_lines') and len(self.heading_lines) >= 1
        has_h2 = hasattr(self, 'heading_lines') and len(self.heading_lines) >= 2
        if hasattr(self, 'reverse_line_direction_btn'):
            self.reverse_line_direction_btn.setEnabled(has_pitch or has_roll or has_h1 or has_h2)

    def _on_reverse_line_direction_clicked(self):
        available = []
        if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2:
            available.append('pitch')
        if hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2:
            available.append('roll')
        if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 1:
            available.append('heading1')
        if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 2:
            available.append('heading2')
        if not available:
            self._show_message("info", "Reverse Line Direction", "No calibration lines to reverse.")
            return
        dialog = ReverseLineDirectionDialog(self, available)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        selected = dialog.get_selected_lines()
        if not selected:
            return
        for line_type in selected:
            self._reverse_calibration_line(line_type, redraw=False)
        self._plot_survey_plan(preserve_view_limits=True)
        if hasattr(self, '_draw_pitch_line_profile'):
            self._draw_pitch_line_profile()
        if hasattr(self, '_draw_roll_line_profile'):
            self._draw_roll_line_profile()
        if hasattr(self, '_update_cal_line_times'):
            self._update_cal_line_times()
        stats = self._calculate_calibration_survey_statistics()
        labels = [ReverseLineDirectionDialog.LINE_LABELS.get(k, k) for k in selected]
        if stats:
            self.set_cal_info_text(
                f"Reversed: {', '.join(labels)}. Calibration Test Info has been updated (travel and total time may have changed).",
                append=False
            )
        else:
            self.set_cal_info_text(f"Reversed: {', '.join(labels)}.", append=False)

    def _reverse_calibration_line(self, line_type, *, redraw=True):
        """Swap start and end points for the given calibration line. If redraw is True, redraws plot and updates info."""
        if line_type == 'pitch':
            if getattr(self, 'edit_pitch_line_mode', False):
                self._toggle_edit_pitch_line_mode()
            if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2:
                self.pitch_line_points = [self.pitch_line_points[1], self.pitch_line_points[0]]
        elif line_type == 'roll':
            if getattr(self, 'edit_roll_line_mode', False):
                self._toggle_edit_roll_line_mode()
            if hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2:
                self.roll_line_points = [self.roll_line_points[1], self.roll_line_points[0]]
        elif line_type == 'heading1' and hasattr(self, 'heading_lines') and len(self.heading_lines) >= 1:
            self.heading_lines[0] = [self.heading_lines[0][1], self.heading_lines[0][0]]
        elif line_type == 'heading2' and hasattr(self, 'heading_lines') and len(self.heading_lines) >= 2:
            self.heading_lines[1] = [self.heading_lines[1][1], self.heading_lines[1][0]]
        else:
            return
        if redraw:
            self._plot_survey_plan(preserve_view_limits=True)
            if hasattr(self, '_draw_pitch_line_profile'):
                self._draw_pitch_line_profile()
            if hasattr(self, '_draw_roll_line_profile'):
                self._draw_roll_line_profile()
            if hasattr(self, '_update_cal_line_times'):
                self._update_cal_line_times()
            label = ReverseLineDirectionDialog.LINE_LABELS.get(line_type, line_type)
            stats = self._calculate_calibration_survey_statistics()
            if stats:
                self.set_cal_info_text(
                    f"{label} direction reversed. Calibration Test Info has been updated (travel and total time may have changed).",
                    append=False
                )
            else:
                self.set_cal_info_text(f"{label} direction reversed.", append=False)

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
            # --- Build common rows and write DDD/DMM/DMS CSV and TXT via export_utils ---
            cal_rows = []
            for num, name, pts in lines:
                if name.lower().startswith('pitch'):
                    line_name, start_label, end_label = 'Pitch', 'PLS', 'PLE'
                elif name.lower().startswith('roll'):
                    line_name, start_label, end_label = 'Roll', 'RLS', 'RLE'
                elif name.lower().startswith('heading'):
                    m = re.search(r'(\d+)', name)
                    n = m.group(1) if m else '1'
                    line_name, start_label, end_label = f'Heading{n}', f'H{n}S', f'H{n}E'
                else:
                    line_name, start_label, end_label = name, 'START', 'END'
                cal_rows.append((num, line_name, start_label, pts[0][0], pts[0][1]))
                cal_rows.append((num, line_name, end_label, pts[1][0], pts[1][1]))

            csv_file_path = os.path.join(export_dir, f"{export_name}_DDD.csv")
            export_utils.write_ddd_csv(csv_file_path, cal_rows)
            ddm_file_path = os.path.join(export_dir, f"{export_name}_DMM.csv")
            export_utils.write_dmm_csv(ddm_file_path, cal_rows)
            dms_file_path = os.path.join(export_dir, f"{export_name}_DMS.csv")
            export_utils.write_dms_csv(dms_file_path, cal_rows)
            ddm_txt_file_path = os.path.join(export_dir, f"{export_name}_DMM.txt")
            export_utils.write_dmm_txt(ddm_txt_file_path, cal_rows)
            dms_txt_file_path = os.path.join(export_dir, f"{export_name}_DMS.txt")
            export_utils.write_dms_txt(dms_txt_file_path, cal_rows)

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
            lnw_file_path = None
            lnw_lines = [(name, list(pts)) for _num, name, pts in lines]
            if lnw_lines:
                all_pts = [p for _name, pts in lnw_lines for p in pts]
                zone, hem = export_utils.compute_utm_zone_from_points(all_pts)
                utm_suffix = f"_UTM{zone}{'N' if hem == 'North' else 'S'}"
                lnw_file_path = os.path.join(export_dir, f"{export_name}{utm_suffix}.lnw")
                if not export_utils.write_lnw(lnw_file_path, lnw_lines):
                    lnw_file_path = None
            sis_file_path = os.path.join(export_dir, f"{export_name}.asciiplan")
            cal_ascii_lines = [(name, list(pts)) for _num, name, pts in lines]
            export_utils.write_asciiplan(sis_file_path, cal_ascii_lines)
            txt_file_path = os.path.join(export_dir, f"{export_name}_DDD.txt")
            export_utils.write_ddd_txt(txt_file_path, cal_rows)
            stats_file_path = os.path.join(export_dir, f"{export_name}_info.txt")
            stats = self._calculate_calibration_survey_statistics()
            if stats:
                stats_text = self._format_calibration_statistics_text(stats, include_export_date=True, export_name=export_name)
                if stats_text:
                    with open(stats_file_path, 'w', encoding='utf-8') as f:
                        f.write(stats_text)
                else:
                    with open(stats_file_path, 'w', encoding='utf-8') as f:
                        f.write("Calibration survey info.\nNo statistics available for this export.\n")
            else:
                with open(stats_file_path, 'w', encoding='utf-8') as f:
                    f.write("Calibration survey info.\nNo statistics available for this export.\n")
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
            success_msg = f"Survey exported successfully to:\n- {os.path.basename(csv_file_path)}\n- {os.path.basename(shapefile_path)}\n- {os.path.basename(geojson_file_path)}\n"
            if lnw_file_path:
                success_msg += f"- {os.path.basename(lnw_file_path)}\n"
            success_msg += f"- {os.path.basename(sis_file_path)}\n"
            success_msg += f"- {os.path.basename(ddm_txt_file_path)}\n"
            success_msg += f"- {os.path.basename(dms_txt_file_path)}\n"
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
        """Import calibration survey lines from CSV, GeoJSON, LNW, DMS, DMM, or DDD text file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Survey File to Import",
            self.last_cal_import_dir,
            "Known Calibration Files (*_DMS.txt *_DMM.txt *_DDD.txt *_DDD.csv *.csv *.geojson *.json *.lnw);;Hypack LNW files (*.lnw);;Degrees Minutes Seconds Text files (*_DMS.txt);;Degrees Decimal Minutes Text files (*_DMM.txt);;Decimal Degrees Text files (*_DDD.txt);;Decimal Degree CSV files (*_DDD.csv);;CSV files (*.csv);;GeoJSON files (*.geojson);;JSON files (*.json)"
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
            file_basename = os.path.basename(file_path)
            file_processed = False
            
            # Handle *.lnw format (Hypack LNW)
            if file_ext == '.lnw':
                # Show UTM zone dialog
                utm_dialog = UTMZoneDialog(self)
                if utm_dialog.exec() != QDialog.DialogCode.Accepted:
                    return  # User cancelled
                
                utm_zone, hemisphere = utm_dialog.get_utm_info()
                
                imported_lines = self._parse_lnw_file(file_path, utm_zone, hemisphere)
                if imported_lines is None:
                    return  # Error already shown
                suggested = _suggest_calibration_assignments_from_geometry(imported_lines)
                dialog = LineAssignmentDialog(self, imported_lines, suggested)
                if dialog.exec() != QDialog.DialogCode.Accepted:
                    return  # User cancelled
                
                assignments = dialog.get_assignments()
                
                # Map assignments to calibration data structures
                self.pitch_line_points = []
                self.roll_line_points = []
                self.heading_lines = []
                
                for line_idx, assignment_type in assignments.items():
                    if assignment_type == 'pitch':
                        self.pitch_line_points = imported_lines[line_idx]
                    elif assignment_type == 'roll':
                        self.roll_line_points = imported_lines[line_idx]
                    elif assignment_type == 'heading1':
                        if len(self.heading_lines) < 1:
                            self.heading_lines.append([])
                        self.heading_lines[0] = imported_lines[line_idx]
                    elif assignment_type == 'heading2':
                        if len(self.heading_lines) < 2:
                            while len(self.heading_lines) < 2:
                                self.heading_lines.append([])
                        self.heading_lines[1] = imported_lines[line_idx]
                file_processed = True
            
            # Handle *_DMS.txt format (Degrees Minutes Seconds)
            if not file_processed and file_basename.lower().endswith('_dms.txt'):
                imported_lines = self._parse_dms_txt_file(file_path)
                if imported_lines is None:
                    return  # Error already shown
                if len(imported_lines) != 4:
                    self._show_message("error", "Calibration Import",
                                     f"Calibration survey requires exactly 4 lines (8 points). This file has {len(imported_lines)} lines.")
                    return
                file_basename = os.path.basename(file_path)
                suggested = _suggest_calibration_assignments_from_metadata(file_path, file_basename, len(imported_lines))
                if suggested is None:
                    suggested = _suggest_calibration_assignments_from_geometry(imported_lines)
                dialog = LineAssignmentDialog(self, imported_lines, suggested)
                if dialog.exec() != QDialog.DialogCode.Accepted:
                    return  # User cancelled
                
                assignments = dialog.get_assignments()
                
                # Map assignments to calibration data structures
                self.pitch_line_points = []
                self.roll_line_points = []
                self.heading_lines = []
                
                for line_idx, assignment_type in assignments.items():
                    if assignment_type == 'pitch':
                        self.pitch_line_points = imported_lines[line_idx]
                    elif assignment_type == 'roll':
                        self.roll_line_points = imported_lines[line_idx]
                    elif assignment_type == 'heading1':
                        if len(self.heading_lines) < 1:
                            self.heading_lines.append([])
                        self.heading_lines[0] = imported_lines[line_idx]
                    elif assignment_type == 'heading2':
                        if len(self.heading_lines) < 2:
                            while len(self.heading_lines) < 2:
                                self.heading_lines.append([])
                        self.heading_lines[1] = imported_lines[line_idx]
                file_processed = True
            
            # Handle *_DMM.txt format (Degrees and decimal minutes)
            if not file_processed and file_basename.lower().endswith('_dmm.txt'):
                imported_lines = self._parse_dmm_txt_file(file_path)
                if imported_lines is None:
                    return  # Error already shown
                if len(imported_lines) != 4:
                    self._show_message("error", "Calibration Import",
                                     f"Calibration survey requires exactly 4 lines (8 points). This file has {len(imported_lines)} lines.")
                    return
                file_basename = os.path.basename(file_path)
                suggested = _suggest_calibration_assignments_from_metadata(file_path, file_basename, len(imported_lines))
                if suggested is None:
                    suggested = _suggest_calibration_assignments_from_geometry(imported_lines)
                dialog = LineAssignmentDialog(self, imported_lines, suggested)
                if dialog.exec() != QDialog.DialogCode.Accepted:
                    return  # User cancelled
                
                assignments = dialog.get_assignments()
                self.pitch_line_points = []
                self.roll_line_points = []
                self.heading_lines = []
                
                for line_idx, assignment_type in assignments.items():
                    if assignment_type == 'pitch':
                        self.pitch_line_points = imported_lines[line_idx]
                    elif assignment_type == 'roll':
                        self.roll_line_points = imported_lines[line_idx]
                    elif assignment_type == 'heading1':
                        if len(self.heading_lines) < 1:
                            self.heading_lines.append([])
                        self.heading_lines[0] = imported_lines[line_idx]
                    elif assignment_type == 'heading2':
                        if len(self.heading_lines) < 2:
                            while len(self.heading_lines) < 2:
                                self.heading_lines.append([])
                        self.heading_lines[1] = imported_lines[line_idx]
                file_processed = True
            
            # Handle *_DDD.txt format (Decimal Degrees)
            if not file_processed and file_basename.lower().endswith('_ddd.txt'):
                imported_lines = self._parse_ddd_txt_file(file_path)
                if imported_lines is None:
                    return  # Error already shown
                if len(imported_lines) != 4:
                    self._show_message("error", "Calibration Import",
                                     f"Calibration survey requires exactly 4 lines (8 points). This file has {len(imported_lines)} lines.")
                    return
                file_basename = os.path.basename(file_path)
                suggested = _suggest_calibration_assignments_from_metadata(file_path, file_basename, len(imported_lines))
                if suggested is None:
                    suggested = _suggest_calibration_assignments_from_geometry(imported_lines)
                dialog = LineAssignmentDialog(self, imported_lines, suggested)
                if dialog.exec() != QDialog.DialogCode.Accepted:
                    return  # User cancelled
                
                assignments = dialog.get_assignments()
                
                # Map assignments to calibration data structures
                self.pitch_line_points = []
                self.roll_line_points = []
                self.heading_lines = []
                
                for line_idx, assignment_type in assignments.items():
                    if assignment_type == 'pitch':
                        self.pitch_line_points = imported_lines[line_idx]
                    elif assignment_type == 'roll':
                        self.roll_line_points = imported_lines[line_idx]
                    elif assignment_type == 'heading1':
                        if len(self.heading_lines) < 1:
                            self.heading_lines.append([])
                        self.heading_lines[0] = imported_lines[line_idx]
                    elif assignment_type == 'heading2':
                        if len(self.heading_lines) < 2:
                            while len(self.heading_lines) < 2:
                                self.heading_lines.append([])
                        self.heading_lines[1] = imported_lines[line_idx]
                file_processed = True
            elif not file_processed and file_ext == '.csv':
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
                file_processed = True
            elif not file_processed:
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
            # Suggested export name: from params, else from loaded GeoTIFF, else defer if GMRT download, else pitch-line or generic
            if not params or not isinstance(params, dict) or not params.get('export_name'):
                self._set_cal_export_name_after_import()
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
                if getattr(self, 'cal_download_gmrt_checkbox', None) and self.cal_download_gmrt_checkbox.isChecked():
                    self._download_and_load_gmrt_after_import()
            else:
                self._show_message("warning","Import Warning", "No valid calibration lines found in the selected file.")
        except Exception as e:
            self._show_message("error","Import Error", f"Failed to import calibration survey files: {e}")
            import traceback
            traceback.print_exc()

    def _download_and_load_gmrt_after_import(self):
        """Compute extent from calibration lines (1° buffer), then use GMRT mixin to download and load."""
        all_points = []
        if len(self.pitch_line_points) == 2:
            all_points.extend(self.pitch_line_points)
        if len(self.roll_line_points) == 2:
            all_points.extend(self.roll_line_points)
        for line in self.heading_lines:
            if len(line) == 2:
                all_points.extend(line)
        if not all_points:
            self._show_message("warning", "GMRT Download", "No calibration line points to compute extent.")
            return
        lats = [p[0] for p in all_points]
        lons = [p[1] for p in all_points]
        min_lat, max_lat = min(lats), max(lats)
        min_lon, max_lon = min(lons), max(lons)
        mid_lat = (min_lat + max_lat) / 2.0
        mid_lon = (min_lon + max_lon) / 2.0
        buffer_deg = 0.5
        if hasattr(self, 'cal_gmrt_buffer_spin'):
            try:
                buffer_deg = float(self.cal_gmrt_buffer_spin.value())
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
            log_func=lambda msg, append=True: self.set_cal_info_text(msg, append=append),
            default_directory=getattr(self, "last_cal_import_dir", None),
        )

    def _calculate_calibration_survey_statistics(self):
        """Calculate comprehensive calibration survey statistics including pitch, roll, heading lines and travel distances."""
        try:
            if pyproj is None:
                return None
            geod = pyproj.Geod(ellps="WGS84")
            speed_knots = float(self.cal_survey_speed_entry.text()) if self.cal_survey_speed_entry.text() else 8.0
            speed_m_per_h = speed_knots * 1852
            turn_time_min = float(self.cal_turn_time_entry.text()) if hasattr(self, 'cal_turn_time_entry') and self.cal_turn_time_entry.text() else 5.0

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
                # 2 passes + 1 turn time between passes
                stats['pitch_line_time_min'] = pitch_time_hours * 60 * 2 + turn_time_min

            if hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2:
                (lat1, lon1), (lat2, lon2) = self.roll_line_points
                _, _, roll_distance = geod.inv(lon1, lat1, lon2, lat2)
                stats['roll_line_distance_m'] = roll_distance
                roll_time_hours = roll_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                # 2 passes + 1 turn time between passes
                stats['roll_line_time_min'] = roll_time_hours * 60 * 2 + turn_time_min

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
                # Travel time + turn time
                stats['travel_pitch_to_roll_min'] = travel_time_hours * 60 + turn_time_min

            if (hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2 and
                hasattr(self, 'heading_lines') and len(self.heading_lines) >= 1):
                roll_start_lat, roll_start_lon = self.roll_line_points[0]
                heading1_start_lat, heading1_start_lon = self.heading_lines[0][0]
                _, _, travel_distance = geod.inv(roll_start_lon, roll_start_lat, heading1_start_lon, heading1_start_lat)
                stats['travel_roll_to_heading1_m'] = travel_distance
                travel_time_hours = travel_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                # Travel time + turn time
                stats['travel_roll_to_heading1_min'] = travel_time_hours * 60 + turn_time_min

            if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 2:
                heading1_end_lat, heading1_end_lon = self.heading_lines[0][1]
                heading2_start_lat, heading2_start_lon = self.heading_lines[1][0]
                _, _, travel_distance = geod.inv(heading1_end_lon, heading1_end_lat, heading2_start_lon, heading2_start_lat)
                stats['travel_heading1_to_heading2_m'] = travel_distance
                travel_time_hours = travel_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                # Travel time + turn time
                stats['travel_heading1_to_heading2_min'] = travel_time_hours * 60 + turn_time_min

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
            
            # Calculate Total Survey Time (2 pitch passes + 2 roll passes + 2 heading passes, without turn times)
            total_survey_time_min = 0.0
            if stats['pitch_line_time_min'] > 0:
                total_survey_time_min += stats['pitch_line_time_min'] - turn_time_min  # Remove 1 turn time from pitch (between passes)
            if stats['roll_line_time_min'] > 0:
                total_survey_time_min += stats['roll_line_time_min'] - turn_time_min  # Remove 1 turn time from roll (between passes)
            total_survey_time_min += stats['heading1_time_min'] + stats['heading2_time_min']
            stats['total_survey_time_min'] = total_survey_time_min
            stats['total_survey_time_hours'] = total_survey_time_min / 60.0
            
            # Calculate Total Transit Time (transit times without turn times)
            total_transit_time_min = 0.0
            if stats['travel_pitch_to_roll_min'] > 0:
                total_transit_time_min += stats['travel_pitch_to_roll_min'] - turn_time_min
            if stats['travel_roll_to_heading1_min'] > 0:
                total_transit_time_min += stats['travel_roll_to_heading1_min'] - turn_time_min
            if stats['travel_heading1_to_heading2_min'] > 0:
                total_transit_time_min += stats['travel_heading1_to_heading2_min'] - turn_time_min
            stats['total_transit_time_min'] = total_transit_time_min
            stats['total_transit_time_hours'] = total_transit_time_min / 60.0
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

    def _format_calibration_statistics_text(self, stats, include_export_date=False, export_name=None):
        """Format comprehensive calibration survey statistics as text string."""
        if not stats:
            return None
        
        stats_text = "COMPREHENSIVE CALIBRATION SURVEY STATISTICS\n"
        stats_text += "=" * 50 + "\n\n"
        if include_export_date and export_name:
            stats_text += f"Calibration Survey: {export_name}\n"
            stats_text += f"Export Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
        try:
            speed_knots = float(self.cal_survey_speed_entry.text()) if self.cal_survey_speed_entry.text() else 8.0
            stats_text += f"Survey Speed: {speed_knots} knots\n"
        except Exception:
            stats_text += "Survey Speed: 8.0 knots (default)\n"
        try:
            turn_time_min = float(self.cal_turn_time_entry.text()) if hasattr(self, 'cal_turn_time_entry') and self.cal_turn_time_entry.text() else 5.0
            stats_text += f"Turn Time (per turn): {turn_time_min} min\n\n"
        except Exception:
            turn_time_min = 5.0
            stats_text += f"Turn Time (per turn): {turn_time_min} min (default)\n\n"

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
            stats_text += f"Roll → Heading1: {stats['travel_roll_to_heading1_m']:.1f} m, {stats['travel_roll_to_heading1_min']:.1f} min\n\n"
        if stats['travel_heading1_to_heading2_m'] > 0:
            stats_text += f"Heading1 → Heading2: {stats['travel_heading1_to_heading2_m']:.1f} m, {stats['travel_heading1_to_heading2_min']:.1f} min\n\n"
        # Calculate total turn time
        num_turns = 0
        if stats['pitch_line_distance_m'] > 0:
            num_turns += 1  # Turn between pitch passes
        if stats['roll_line_distance_m'] > 0:
            num_turns += 1  # Turn between roll passes
        if stats['travel_pitch_to_roll_m'] > 0:
            num_turns += 1  # Turn between pitch and roll
        if stats['travel_roll_to_heading1_m'] > 0:
            num_turns += 1  # Turn between roll and heading1
        if stats['travel_heading1_to_heading2_m'] > 0:
            num_turns += 1  # Turn between heading1 and heading2
        
        total_turn_time_min = num_turns * turn_time_min
        
        stats_text += "TOTAL SURVEY SUMMARY\n"
        stats_text += "-" * 25 + "\n"
        stats_text += f"Total Distance: {stats['total_distance_m']:.1f} m ({stats['total_distance_km']:.3f} km, {stats['total_distance_nm']:.3f} nm)\n"
        stats_text += f"Total Survey Time: {stats['total_survey_time_min']:.1f} min ({stats['total_survey_time_hours']:.2f} hr)\n"
        stats_text += f"Total Transit Time: {stats['total_transit_time_min']:.1f} min ({stats['total_transit_time_hours']:.2f} hr)\n"
        stats_text += f"Total Turn Time: {total_turn_time_min:.1f} min ({num_turns} turns × {turn_time_min:.1f} min)\n"
        stats_text += f"Total Time: {stats['total_time_min']:.1f} min ({stats['total_time_hours']:.2f} hr)\n"

        # Calibration Waypoints (DMM)
        dmm_heading = "Calibration Waypoints (DMM)"
        stats_text += f"\n{dmm_heading}\n"
        stats_text += "-" * len(dmm_heading) + "\n"
        if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2:
            pitch_start = self.pitch_line_points[0]
            pitch_end = self.pitch_line_points[1]
            stats_text += f"Pitch Start: {decimal_degrees_to_ddm(pitch_start[0], True)}, {decimal_degrees_to_ddm(pitch_start[1], False)}\n"
            stats_text += f"Pitch End: {decimal_degrees_to_ddm(pitch_end[0], True)}, {decimal_degrees_to_ddm(pitch_end[1], False)}\n"
        if hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2:
            roll_start = self.roll_line_points[0]
            roll_end = self.roll_line_points[1]
            stats_text += f"Roll Start: {decimal_degrees_to_ddm(roll_start[0], True)}, {decimal_degrees_to_ddm(roll_start[1], False)}\n"
            stats_text += f"Roll End: {decimal_degrees_to_ddm(roll_end[0], True)}, {decimal_degrees_to_ddm(roll_end[1], False)}\n"
        if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 1:
            h1_start = self.heading_lines[0][0]
            h1_end = self.heading_lines[0][1]
            stats_text += f"Heading1 Start: {decimal_degrees_to_ddm(h1_start[0], True)}, {decimal_degrees_to_ddm(h1_start[1], False)}\n"
            stats_text += f"Heading1 End: {decimal_degrees_to_ddm(h1_end[0], True)}, {decimal_degrees_to_ddm(h1_end[1], False)}\n"
        if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 2:
            h2_start = self.heading_lines[1][0]
            h2_end = self.heading_lines[1][1]
            stats_text += f"Heading2 Start: {decimal_degrees_to_ddm(h2_start[0], True)}, {decimal_degrees_to_ddm(h2_start[1], False)}\n"
            stats_text += f"Heading2 End: {decimal_degrees_to_ddm(h2_end[0], True)}, {decimal_degrees_to_ddm(h2_end[1], False)}\n"

        # Calibration Waypoints (DDD)
        ddd_heading = "Calibration Waypoints (DDD)"
        stats_text += f"\n{ddd_heading}\n"
        stats_text += "-" * len(ddd_heading) + "\n"
        if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2:
            pitch_start = self.pitch_line_points[0]
            pitch_end = self.pitch_line_points[1]
            stats_text += f"Pitch Start: {pitch_start[0]:.6f}, {pitch_start[1]:.6f}\n"
            stats_text += f"Pitch End: {pitch_end[0]:.6f}, {pitch_end[1]:.6f}\n"
        if hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2:
            roll_start = self.roll_line_points[0]
            roll_end = self.roll_line_points[1]
            stats_text += f"Roll Start: {roll_start[0]:.6f}, {roll_start[1]:.6f}\n"
            stats_text += f"Roll End: {roll_end[0]:.6f}, {roll_end[1]:.6f}\n"
        if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 1:
            h1_start = self.heading_lines[0][0]
            h1_end = self.heading_lines[0][1]
            stats_text += f"Heading1 Start: {h1_start[0]:.6f}, {h1_start[1]:.6f}\n"
            stats_text += f"Heading1 End: {h1_end[0]:.6f}, {h1_end[1]:.6f}\n"
        if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 2:
            h2_start = self.heading_lines[1][0]
            h2_end = self.heading_lines[1][1]
            stats_text += f"Heading2 Start: {h2_start[0]:.6f}, {h2_start[1]:.6f}\n"
            stats_text += f"Heading2 End: {h2_end[0]:.6f}, {h2_end[1]:.6f}\n"

        return stats_text

    def _show_calibration_statistics(self):
        """Display comprehensive calibration survey statistics in a custom dialog."""
        stats = self._calculate_calibration_survey_statistics()
        if not stats:
            self._show_message("warning","Statistics Error", "Unable to calculate calibration survey statistics. Please ensure all required lines are drawn.")
            return

        stats_text = self._format_calibration_statistics_text(stats)
        if stats_text:
            show_statistics_dialog(self, "Calibration Survey Info", stats_text)

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
        
        # Get pitch line points (needed for validation and line creation)
        (lat1, lon1), (lat2, lon2) = self.pitch_line_points
        
        # Check if offset is more than 4 times the shallowest depth
        shallowest_depth = None
        # Try to get shallowest depth from label
        if hasattr(self, 'pitch_shallowest_depth_label'):
            try:
                shallowest_depth_text = self.pitch_shallowest_depth_label.text().strip()
                if shallowest_depth_text and shallowest_depth_text != "-":
                    shallowest_depth = float(shallowest_depth_text)
            except (ValueError, AttributeError):
                pass
        
        # If not available from label, recalculate it
        if shallowest_depth is None and self.geotiff_data_array is not None and self.geotiff_extent is not None:
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
                max_val = np.nanmax(elevations)
                shallowest_depth = abs(max_val)
        
        # Warn if offset is more than 2 times shallowest depth
        if shallowest_depth is not None and shallowest_depth > 0:
            if offset_val > 2 * shallowest_depth:
                self._show_message(
                    "warning",
                    "Large Offset Warning",
                    f"Heading Line Offset ({offset_val:.2f} m) is more than 2 times the shallowest depth "
                    f"({shallowest_depth:.2f} m). This may result in heading lines that are too far from "
                    f"the pitch line for effective calibration."
                )
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

    def _set_cal_export_name_after_import(self):
        """Set suggested calibration export name after import using same convention as drawing: Cal_{offset}m_{heading}deg from pitch line."""
        if not hasattr(self, 'cal_export_name_entry'):
            return
        # If Download GMRT is checked, leave empty; GMRT load completion will set name from pitch line (offset + heading) after GeoTIFF loads
        if getattr(self, 'cal_download_gmrt_checkbox', None) and self.cal_download_gmrt_checkbox.isChecked():
            return
        # Use same naming as drawn survey: depth (offset) from pitch line + orientation (heading) from pitch line
        self._update_cal_line_offset_from_pitch_line()  # fill offset from GeoTIFF if loaded
        self._update_cal_export_name_from_pitch_line()
        if not self.cal_export_name_entry.text().strip():
            self.cal_export_name_entry.setText(f"Cal_import_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}")

    def _update_cal_export_name_from_pitch_line(self):
        """Update calibration export name from pitch line and offset (same convention as when drawing: Cal_{offset}m_{heading}deg)."""
        if not hasattr(self, 'pitch_line_points') or len(self.pitch_line_points) != 2:
            return
        try:
            offset_val = float(self.cal_line_offset_entry.text().strip() or "0")
        except (ValueError, TypeError):
            offset_val = 0
        (lat1, lon1), (lat2, lon2) = self.pitch_line_points
        try:
            geod = pyproj.Geod(ellps="WGS84")
            az12, az21, dist = geod.inv(lon1, lat1, lon2, lat2)
            heading = int(round(az12)) % 360
        except Exception:
            heading = 0
        export_name = f"Cal_{int(round(offset_val))}m_{heading}deg"
        self.cal_export_name_entry.setText(export_name)

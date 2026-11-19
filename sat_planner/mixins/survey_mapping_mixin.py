"""
Survey tab: draw a line, move its vertices, then add parallel offset lines.
"""
import csv
import datetime
import json
import math
import os

import numpy as np

from PyQt6.QtWidgets import (
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QFileDialog,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSizePolicy,
    QSpinBox,
    QVBoxLayout,
)
from PyQt6.QtCore import Qt

from sat_planner import decimal_degrees_to_ddm, export_utils
from sat_planner.utils_ui import show_statistics_dialog
from sat_planner.constants import GEOSPATIAL_LIBS_AVAILABLE, LineString, fiona, pyproj
from sat_planner.utils_geo import lat_lon_decimal_from_survey_csv_row

try:
    from shapely.geometry import mapping as _shapely_mapping
except ImportError:
    _shapely_mapping = None

try:
    from .survey_parsers_mixin import UTMZoneDialog
except ImportError:
    UTMZoneDialog = None


class SurveyMappingMixin:
    """A survey line the user can reshape, plus optional parallel lines."""

    def _build_survey_tab(self):
        self.survey_draw_mode = False
        self.survey_edit_mode = False
        self.survey_center_points = []
        self.survey_offset_m = 100.0
        self.survey_offset_count = 1
        self.survey_offset_side = "both"
        self.survey_offsets_applied = False
        self.survey_start_edge = None
        self.survey_start_vertex = None
        self.survey_pick_start_mode = False
        self.survey_mask_points = []
        self.survey_mask_draft = []
        self.survey_mask_draw_mode = False
        self.survey_temp_line = None
        self.survey_vertex_handles = []
        self.survey_dragging_vertex = None
        self.survey_line_artists = []

        layout = QVBoxLayout(self.survey_frame)
        group = QGroupBox("Mapping Survey")
        group_layout = QVBoxLayout(group)

        self.survey_draw_btn = QPushButton("Draw Line")
        self.survey_draw_btn.setToolTip(
            "Left-click to add vertices. Right-click to finish the line."
        )
        self.survey_draw_btn.clicked.connect(self._toggle_survey_draw)
        group_layout.addWidget(self.survey_draw_btn)

        self.survey_move_btn = QPushButton("Move Vertices")
        self.survey_move_btn.setEnabled(False)
        self.survey_move_btn.setToolTip("Drag the vertices until the line is where you want it.")
        self.survey_move_btn.clicked.connect(self._toggle_survey_edit)
        group_layout.addWidget(self.survey_move_btn)

        self.survey_reverse_btn = QPushButton("Reverse Direction")
        self.survey_reverse_btn.setEnabled(False)
        self.survey_reverse_btn.setToolTip("Swap the start and end. Port and starboard swap with them.")
        self.survey_reverse_btn.clicked.connect(self._reverse_survey_centerline)
        group_layout.addWidget(self.survey_reverse_btn)

        distance_row = QHBoxLayout()
        distance_row.addWidget(QLabel("Spacing"))
        self.survey_offset_spin = QDoubleSpinBox()
        self.survey_offset_spin.setRange(1.0, 100000.0)
        self.survey_offset_spin.setDecimals(1)
        self.survey_offset_spin.setSingleStep(10.0)
        self.survey_offset_spin.setValue(self.survey_offset_m)
        self.survey_offset_spin.setSuffix(" m")
        self.survey_offset_spin.setToolTip("Distance between neighboring offset lines")
        self.survey_offset_spin.valueChanged.connect(self._on_survey_offset_changed)
        distance_row.addWidget(self.survey_offset_spin, 1)
        group_layout.addLayout(distance_row)

        count_row = QHBoxLayout()
        count_row.addWidget(QLabel("Lines"))
        self.survey_count_spin = QSpinBox()
        self.survey_count_spin.setRange(1, 200)
        self.survey_count_spin.setValue(self.survey_offset_count)
        self.survey_count_spin.setToolTip("How many offset lines to draw on the chosen side")
        self.survey_count_spin.valueChanged.connect(self._on_survey_count_changed)
        count_row.addWidget(self.survey_count_spin, 1)
        group_layout.addLayout(count_row)

        side_row = QHBoxLayout()
        side_row.addWidget(QLabel("Side"))
        self.survey_side_combo = QComboBox()
        self.survey_side_combo.addItem("Port", "port")
        self.survey_side_combo.addItem("Starboard", "starboard")
        self.survey_side_combo.addItem("Both", "both")
        self.survey_side_combo.setCurrentIndex(2)
        self.survey_side_combo.setToolTip("Port is left of the line direction. Starboard is right.")
        self.survey_side_combo.currentIndexChanged.connect(self._on_survey_side_changed)
        side_row.addWidget(self.survey_side_combo, 1)
        group_layout.addLayout(side_row)

        self.survey_add_offsets_btn = QPushButton("Add Offset Lines")
        self.survey_add_offsets_btn.setEnabled(False)
        self.survey_add_offsets_btn.setToolTip(
            "Draw the offset lines from the line as it is now."
        )
        self.survey_add_offsets_btn.clicked.connect(self._add_survey_offsets)
        group_layout.addWidget(self.survey_add_offsets_btn)

        self.survey_select_beginning_btn = QPushButton("Select Beginning")
        self.survey_select_beginning_btn.setEnabled(False)
        self.survey_select_beginning_btn.setToolTip(
            "Click, then click an end of an outside line. That point becomes line 1, point 1."
        )
        self.survey_select_beginning_btn.clicked.connect(self._survey_select_beginning)
        group_layout.addWidget(self.survey_select_beginning_btn)

        self.survey_clear_offsets_btn = QPushButton("Clear Offset Lines")
        self.survey_clear_offsets_btn.setEnabled(False)
        self.survey_clear_offsets_btn.clicked.connect(self._clear_survey_offsets)
        group_layout.addWidget(self.survey_clear_offsets_btn)

        self.survey_draw_mask_btn = QPushButton("Draw Mask")
        self.survey_draw_mask_btn.setToolTip(
            "Left-click to add corners. Right-click to close the polygon. Lines outside it are clipped."
        )
        self.survey_draw_mask_btn.clicked.connect(self._toggle_survey_mask_draw)
        group_layout.addWidget(self.survey_draw_mask_btn)

        self.survey_remove_mask_btn = QPushButton("Remove Mask")
        self.survey_remove_mask_btn.setEnabled(False)
        self.survey_remove_mask_btn.setToolTip("Remove the mask and show the original lines.")
        self.survey_remove_mask_btn.clicked.connect(self._remove_survey_mask)
        group_layout.addWidget(self.survey_remove_mask_btn)

        hint = QLabel(
            "Draw the line and move its vertices first. "
            "Offset lines are added from that line. "
            "Select Beginning, then click the starting point on the map. "
            "Port is green, to the left. Starboard is red, to the right."
        )
        hint.setWordWrap(True)
        group_layout.addWidget(hint)

        self.survey_clear_btn = QPushButton("Clear Survey")
        self.survey_clear_btn.clicked.connect(self._clear_survey_plan)
        group_layout.addWidget(self.survey_clear_btn)

        self.survey_status_label = QLabel("Draw a line, then move its vertices.")
        self.survey_status_label.setWordWrap(True)
        group_layout.addWidget(self.survey_status_label)

        self.survey_profile_select_combo = QComboBox()
        self.survey_profile_select_combo.addItem("L1")
        self.survey_profile_select_combo.setToolTip(
            "Starts on the line you drew. The other lines are listed by their survey names."
        )
        self.survey_profile_user_choice = None
        self.survey_profile_select_combo.currentIndexChanged.connect(self._on_survey_profile_select_changed)

        layout.addWidget(group)
        layout.addWidget(self._build_survey_info_group())
        layout.addStretch(1)

    def _build_survey_info_group(self):
        """Depth and slope across every survey line, plus survey speed."""
        group = QGroupBox("Survey Info")
        group.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        grid = QGridLayout(group)
        grid.setSpacing(3)
        grid.setContentsMargins(9, 9, 9, 9)
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 2)
        rows = (
            ("Shallowest Depth (m):", "survey_shallowest_depth_label"),
            ("Maximum Depth (m):", "survey_max_depth_label"),
            ("Minimum Slope (deg):", "survey_min_slope_label"),
            ("Maximum Slope (deg):", "survey_max_slope_label"),
        )
        for row, (caption, name) in enumerate(rows):
            grid.addWidget(QLabel(caption), row, 0)
            label = QLabel("-")
            setattr(self, name, label)
            grid.addWidget(label, row, 1)
        grid.addWidget(QLabel("Survey Speed (knots):"), len(rows), 0)
        self.survey_speed_entry = QLineEdit("8")
        grid.addWidget(self.survey_speed_entry, len(rows), 1)
        self.survey_show_info_btn = QPushButton("Show Survey Info")
        self.survey_show_info_btn.clicked.connect(self._show_survey_information)
        grid.addWidget(self.survey_show_info_btn, len(rows) + 1, 0, 1, 2)
        return group

    def _survey_tab_is_current(self):
        notebook = getattr(self, "param_notebook", None)
        index = getattr(self, "survey_tab_index", None)
        return notebook is not None and index is not None and notebook.currentIndex() == index

    def _survey_view_intersects_geotiff(self):
        """True when the current map window is a local view of the loaded grid."""
        extent = getattr(self, "geotiff_extent", None)
        ax = getattr(self, "ax", None)
        if extent is None or ax is None:
            return False
        try:
            west, east = sorted((float(extent[0]), float(extent[1])))
            south, north = sorted((float(extent[2]), float(extent[3])))
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            view_west, view_east = sorted((float(xlim[0]), float(xlim[1])))
            view_south, view_north = sorted((float(ylim[0]), float(ylim[1])))
        except (TypeError, ValueError):
            return False
        if east <= west or north <= south:
            return False
        local = (view_east - view_west) < 300 and (view_north - view_south) < 150
        hits = view_west <= east and view_east >= west and view_south <= north and view_north >= south
        return local and hits

    def _ensure_survey_map_geotiff(self):
        """Draw the loaded planning grid on the map while Survey is open."""
        if not self._survey_tab_is_current():
            return
        self._overview_map_active = False
        if not getattr(self, "show_bathymetry_var", True):
            return
        if getattr(self, "geotiff_data_array", None) is None or getattr(self, "geotiff_extent", None) is None:
            return
        if not self._survey_view_intersects_geotiff():
            if hasattr(self, "_zoom_to_geotiff"):
                self._zoom_to_geotiff()
            return
        has_image = any(
            getattr(self, name, None) is not None
            for name in ("geotiff_hillshade_plot", "geotiff_image_plot", "depth_mode_image_plot")
        )
        if not has_image and hasattr(self, "_plot_survey_plan"):
            self._plot_survey_plan(preserve_view_limits=True)
        elif hasattr(self, "_sync_maplibre_bathymetry"):
            self._sync_maplibre_bathymetry()
        if hasattr(self, "_fit_maplibre_to_axes"):
            self._fit_maplibre_to_axes()

    def _toggle_survey_draw(self):
        if getattr(self, "survey_draw_mode", False):
            self._stop_survey_draw(keep_line=True)
            if hasattr(self, "_plot_survey_plan"):
                self._plot_survey_plan(preserve_view_limits=True)
            return
        self._stop_survey_edit(replot=False)
        if getattr(self, "survey_mask_draw_mode", False):
            self._stop_survey_mask_draw(keep=False)
        self.survey_draw_mode = True
        self.survey_center_points = []
        self.survey_offsets_applied = False
        self._survey_clear_start_end()
        self._stop_survey_pick_start()
        self._survey_remove_temp_line()
        self.survey_draw_btn.setText("Left-click to add, right-click to finish")
        self.survey_draw_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")
        if hasattr(self, "canvas_widget"):
            self.canvas_widget.setCursor(Qt.CursorShape.CrossCursor)
        self._survey_refresh_controls()
        if hasattr(self, "_plot_survey_plan"):
            self._plot_survey_plan(preserve_view_limits=True)

    def _stop_survey_draw(self, keep_line=True):
        self.survey_draw_mode = False
        self._survey_remove_temp_line()
        points = getattr(self, "survey_center_points", None) or []
        if not keep_line or len(points) < 2:
            self.survey_center_points = []
            self.survey_offsets_applied = False
            self._survey_clear_start_end()
            self._stop_survey_pick_start()
        if hasattr(self, "survey_draw_btn"):
            self.survey_draw_btn.setText("Draw Line")
            self.survey_draw_btn.setStyleSheet("")
        if hasattr(self, "canvas_widget") and not getattr(self, "survey_edit_mode", False):
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
        self._survey_refresh_controls()

    def _toggle_survey_edit(self):
        points = getattr(self, "survey_center_points", None) or []
        if len(points) < 2 or getattr(self, "survey_draw_mode", False):
            return
        self.survey_edit_mode = not getattr(self, "survey_edit_mode", False)
        self.survey_dragging_vertex = None
        if self.survey_edit_mode:
            self.survey_move_btn.setText("Stop Moving Vertices")
            self.survey_move_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")
            if hasattr(self, "canvas_widget"):
                self.canvas_widget.setCursor(Qt.CursorShape.SizeAllCursor)
        else:
            self._stop_survey_edit(replot=False)
        if hasattr(self, "_plot_survey_plan"):
            self._plot_survey_plan(preserve_view_limits=True)
        self._survey_refresh_controls()

    def _stop_survey_edit(self, replot=False):
        self.survey_edit_mode = False
        self.survey_dragging_vertex = None
        self.survey_vertex_handles = []
        if hasattr(self, "survey_move_btn"):
            self.survey_move_btn.setText("Move Vertices")
            self.survey_move_btn.setStyleSheet("")
        if hasattr(self, "canvas_widget") and not getattr(self, "survey_draw_mode", False):
            if getattr(self, "survey_pick_start_mode", False):
                self.canvas_widget.setCursor(Qt.CursorShape.CrossCursor)
            else:
                self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
        if replot and hasattr(self, "_plot_survey_plan"):
            self._plot_survey_plan(preserve_view_limits=True)

    def _on_survey_offset_changed(self, value):
        self.survey_offset_m = float(value)
        self._survey_replot_if_offsets()

    def _on_survey_count_changed(self, value):
        self.survey_offset_count = int(value)
        self._survey_replot_if_offsets()

    def _on_survey_side_changed(self, _index=None):
        combo = getattr(self, "survey_side_combo", None)
        side = combo.currentData() if combo is not None else "both"
        self.survey_offset_side = side or "both"
        self._survey_replot_if_offsets()

    def _survey_replot_if_offsets(self):
        if getattr(self, "survey_offsets_applied", False) and hasattr(self, "_plot_survey_plan"):
            self._plot_survey_plan(preserve_view_limits=True)
        self._survey_refresh_controls()

    def _add_survey_offsets(self):
        points = getattr(self, "survey_center_points", None) or []
        if len(points) < 2:
            self.survey_status_label.setText("Draw a line before adding offset lines.")
            return
        self.survey_offsets_applied = True
        self._survey_begin_pick_start()
        if hasattr(self, "_plot_survey_plan"):
            self._plot_survey_plan(preserve_view_limits=True)
        self._survey_refresh_controls()

    def _survey_select_beginning(self):
        """Start the map click that chooses line 1, point 1. Works the first time and again later."""
        points = getattr(self, "survey_center_points", None) or []
        if len(points) < 2 or not getattr(self, "survey_offsets_applied", False):
            self.survey_status_label.setText("Add offset lines before selecting the beginning.")
            return
        self._survey_begin_pick_start()
        if hasattr(self, "_plot_survey_plan"):
            self._plot_survey_plan(preserve_view_limits=True)
        self._survey_refresh_controls()

    def _clear_survey_offsets(self):
        self.survey_offsets_applied = False
        self._survey_clear_start_end()
        self._stop_survey_pick_start()
        if hasattr(self, "_plot_survey_plan"):
            self._plot_survey_plan(preserve_view_limits=True)
        self._survey_refresh_controls()

    def _reverse_survey_centerline(self):
        points = getattr(self, "survey_center_points", None) or []
        if len(points) < 2:
            return
        self.survey_center_points = list(reversed(points))
        self._survey_keep_start_end_after_reverse()
        if hasattr(self, "_plot_survey_plan"):
            self._plot_survey_plan(preserve_view_limits=True)
        self._survey_refresh_controls()

    def _clear_survey_plan(self):
        self.survey_center_points = []
        self.survey_offsets_applied = False
        self.survey_mask_points = []
        self.survey_mask_draft = []
        self._survey_clear_start_end()
        self._stop_survey_pick_start()
        self._stop_survey_mask_draw(keep=False)
        self._stop_survey_edit(replot=False)
        self._stop_survey_draw(keep_line=False)
        self.survey_status_label.setText("Draw a line, then move its vertices.")
        if hasattr(self, "_plot_survey_plan"):
            self._plot_survey_plan(preserve_view_limits=True)

    def _survey_refresh_controls(self):
        points = getattr(self, "survey_center_points", None) or []
        drawing = bool(getattr(self, "survey_draw_mode", False))
        ready = len(points) >= 2 and not drawing
        if hasattr(self, "survey_reverse_btn"):
            self.survey_reverse_btn.setEnabled(ready)
        if hasattr(self, "survey_move_btn"):
            self.survey_move_btn.setEnabled(ready)
        if hasattr(self, "survey_add_offsets_btn"):
            self.survey_add_offsets_btn.setEnabled(ready)
        offsets_ready = ready and bool(getattr(self, "survey_offsets_applied", False))
        if hasattr(self, "survey_select_beginning_btn"):
            self.survey_select_beginning_btn.setEnabled(offsets_ready)
        if hasattr(self, "_set_button_text_emphasis"):
            self._set_button_text_emphasis(
                getattr(self, "survey_select_beginning_btn", None),
                offsets_ready and not self._survey_start_is_chosen(),
            )
        if hasattr(self, "survey_clear_offsets_btn"):
            self.survey_clear_offsets_btn.setEnabled(bool(getattr(self, "survey_offsets_applied", False)))
        if hasattr(self, "survey_remove_mask_btn"):
            self.survey_remove_mask_btn.setEnabled(len(getattr(self, "survey_mask_points", None) or []) >= 3)
        if hasattr(self, "_set_button_text_emphasis"):
            self._set_button_text_emphasis(getattr(self, "survey_show_info_btn", None), ready)
        if hasattr(self, "_update_shared_survey_io_ui") and self._survey_tab_is_current():
            self._update_shared_survey_io_ui()
        if not hasattr(self, "survey_status_label"):
            return
        if getattr(self, "survey_mask_draw_mode", False):
            count = len(getattr(self, "survey_mask_draft", None) or [])
            if count < 3:
                self.survey_status_label.setText(
                    "Click mask corners. Right-click to close the polygon."
                )
            else:
                self.survey_status_label.setText(
                    f"{count} mask corners. Click to add another, or right-click to close."
                )
            return
        if drawing:
            if len(points) == 0:
                self.survey_status_label.setText("Click the first vertex. Right-click to finish.")
            else:
                self.survey_status_label.setText(
                    f"{len(points)} vertices. Click to add another, or right-click to finish."
                )
            return
        if not ready:
            self.survey_status_label.setText("Draw a line, then move its vertices.")
            return
        length_m = self._survey_centerline_meters()
        length = f"Line is {length_m / 1852.0:.2f} nm, {len(points)} vertices. " if length_m else ""
        if not getattr(self, "survey_offsets_applied", False):
            self.survey_status_label.setText(f"{length}Move the vertices, then add offset lines.")
            return
        if not self._survey_start_is_chosen():
            if getattr(self, "survey_edit_mode", False):
                self.survey_status_label.setText(f"{length}Stop moving vertices, then click Select Beginning.")
            elif getattr(self, "survey_pick_start_mode", False):
                self.survey_status_label.setText(
                    f"{length}Click the starting point of the survey. It becomes line 1, point 1."
                )
            else:
                self.survey_status_label.setText(f"{length}Click Select Beginning, then click the starting point.")
            return
        side = getattr(self, "survey_offset_side", "both")
        distance = float(getattr(self, "survey_offset_m", 0) or 0)
        count = int(getattr(self, "survey_offset_count", 1) or 1)
        if side == "port":
            where = "to port"
        elif side == "starboard":
            where = "to starboard"
        else:
            where = "to port and starboard"
        edge, vertex = self._survey_start_choice()
        end_word = "start" if vertex == "start" else "end"
        self.survey_status_label.setText(
            f"{length}{count} line{'s' if count != 1 else ''} at {distance:.1f} m spacing {where}. "
            f"Numbered from the {edge} side at the {end_word} of the drawn line."
        )

    def _handle_survey_plot_click(self, event):
        """Add a vertex while drawing, or choose the survey start after offset lines are added."""
        if self._survey_mask_click(event):
            return True
        if self._survey_pick_start_click(event):
            return True
        if not getattr(self, "survey_draw_mode", False) or not self._survey_tab_is_current():
            return False
        button = int(getattr(event, "button", 1) or 1)
        if button == 3:
            if len(getattr(self, "survey_center_points", None) or []) >= 2:
                self._stop_survey_draw(keep_line=True)
                if hasattr(self, "_plot_survey_plan"):
                    self._plot_survey_plan(preserve_view_limits=True)
            else:
                self.survey_status_label.setText("Add one more vertex, then right-click to finish.")
            return True
        if button != 1:
            return True
        if hasattr(self, "_main_map_click_lonlat"):
            lon, lat = self._main_map_click_lonlat(event)
        else:
            lon, lat = getattr(event, "xdata", None), getattr(event, "ydata", None)
        if lon is None or lat is None:
            return True
        self.survey_center_points.append((float(lat), float(lon)))
        self._survey_remove_temp_line()
        if hasattr(self, "_plot_survey_plan"):
            self._plot_survey_plan(preserve_view_limits=True)
        self._survey_refresh_controls()
        return True

    def _on_survey_draw_motion(self, event):
        """Rubber-band from the last vertex to the pointer."""
        if getattr(self, "survey_mask_draw_mode", False) and self._survey_tab_is_current():
            return self._survey_mask_motion(event)
        if not getattr(self, "survey_draw_mode", False) or not self._survey_tab_is_current():
            return False
        points = getattr(self, "survey_center_points", None) or []
        if not points or getattr(self, "ax", None) is None:
            return False
        lon = getattr(event, "xdata", None)
        lat = getattr(event, "ydata", None)
        if lon is None or lat is None:
            return False
        start_lat, start_lon = points[-1]
        line = getattr(self, "survey_temp_line", None)
        if line is not None and line in self.ax.lines:
            line.set_data([start_lon, lon], [start_lat, lat])
            return True
        self.survey_temp_line, = self.ax.plot(
            [start_lon, lon],
            [start_lat, lat],
            color="#1f4e79",
            linewidth=2,
            linestyle="--",
        )
        return True

    def _survey_remove_temp_line(self):
        line = getattr(self, "survey_temp_line", None)
        if line is not None:
            try:
                line.remove()
            except (NotImplementedError, ValueError, AttributeError):
                pass
        self.survey_temp_line = None

    def _on_survey_vertex_pick(self, event):
        if not getattr(self, "survey_edit_mode", False):
            return
        mouse = getattr(event, "mouseevent", None)
        if mouse is not None and int(getattr(mouse, "button", 1) or 1) != 1:
            return
        for index, handle in enumerate(getattr(self, "survey_vertex_handles", None) or []):
            if event.artist == handle:
                self.survey_dragging_vertex = index
                return

    def _on_survey_vertex_motion(self, event):
        index = getattr(self, "survey_dragging_vertex", None)
        if not getattr(self, "survey_edit_mode", False) or index is None:
            return
        lon = getattr(event, "xdata", None)
        lat = getattr(event, "ydata", None)
        if lon is None or lat is None:
            return
        points = getattr(self, "survey_center_points", None) or []
        if index >= len(points):
            return
        points[index] = (float(lat), float(lon))
        handles = getattr(self, "survey_vertex_handles", None) or []
        if index < len(handles) and handles[index] is not None:
            handles[index].set_offsets([[float(lon), float(lat)]])
        self._survey_update_drawn_lines()

    def _on_survey_vertex_release(self, event):
        if not getattr(self, "survey_edit_mode", False):
            return
        if getattr(self, "survey_dragging_vertex", None) is None:
            return
        self.survey_dragging_vertex = None
        if hasattr(self, "_plot_survey_plan"):
            self._plot_survey_plan(preserve_view_limits=True)
        self._survey_refresh_controls()

    def _survey_update_drawn_lines(self):
        """Move the drawn line with a vertex while it is being dragged."""
        points = getattr(self, "survey_center_points", None) or []
        artists = getattr(self, "survey_line_artists", None) or []
        if len(points) < 2 or not artists:
            return
        artists[0].set_data([point[1] for point in points], [point[0] for point in points])
        if not getattr(self, "survey_offsets_applied", False):
            return
        for artist, (_side, _index, line) in zip(artists[1:], self._survey_offset_records()):
            artist.set_data([point[1] for point in line], [point[0] for point in line])

    def _survey_clear_start_end(self):
        self.survey_start_edge = None
        self.survey_start_vertex = None
        self.survey_profile_user_choice = None

    def _survey_start_choice(self):
        """The end that becomes line 1, point 1. Port and the drawn start are the fallback."""
        edge = getattr(self, "survey_start_edge", None)
        vertex = getattr(self, "survey_start_vertex", None)
        if edge not in ("port", "starboard"):
            edge = "port"
        if vertex not in ("start", "finish"):
            vertex = "start"
        return edge, vertex

    def _survey_start_is_chosen(self):
        edge = getattr(self, "survey_start_edge", None)
        vertex = getattr(self, "survey_start_vertex", None)
        return edge in ("port", "starboard") and vertex in ("start", "finish")

    def _survey_begin_pick_start(self):
        """After the offset lines exist, the next map click chooses line 1, point 1."""
        self._stop_survey_edit(replot=False)
        self.survey_pick_start_mode = True
        self._survey_clear_start_end()
        if hasattr(self, "canvas_widget") and self._survey_tab_is_current():
            self.canvas_widget.setCursor(Qt.CursorShape.CrossCursor)

    def _stop_survey_pick_start(self):
        self.survey_pick_start_mode = False
        if (
            hasattr(self, "canvas_widget")
            and not getattr(self, "survey_draw_mode", False)
            and not getattr(self, "survey_edit_mode", False)
        ):
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)

    def _survey_start_candidates(self):
        """The two ends of each outside line. A click on one of these starts the survey."""
        geometry = self._survey_geometry_lines()
        if len(geometry) < 2:
            return []
        candidates = []
        for edge, index in (("port", 0), ("starboard", -1)):
            _key, line = geometry[index]
            if len(line) < 2:
                continue
            candidates.append((edge, "start", line[0]))
            candidates.append((edge, "finish", line[-1]))
        return candidates

    def _survey_start_from_click(self, lat, lon):
        ax = getattr(self, "ax", None)
        if ax is None:
            return None
        x_span = abs(float(ax.get_xlim()[1]) - float(ax.get_xlim()[0])) or 1.0
        y_span = abs(float(ax.get_ylim()[1]) - float(ax.get_ylim()[0])) or 1.0
        best = None
        best_dist = 0.06 ** 2
        for edge, vertex, (point_lat, point_lon) in self._survey_start_candidates():
            dist = ((lon - point_lon) / x_span) ** 2 + ((lat - point_lat) / y_span) ** 2
            if dist < best_dist:
                best_dist = dist
                best = (edge, vertex)
        return best

    def _survey_pick_start_click(self, event):
        if not getattr(self, "survey_pick_start_mode", False) or not self._survey_tab_is_current():
            return False
        if getattr(self, "survey_edit_mode", False):
            return False
        button = int(getattr(event, "button", 1) or 1)
        if button != 1:
            return True
        if hasattr(self, "_main_map_click_lonlat"):
            lon, lat = self._main_map_click_lonlat(event)
        else:
            lon, lat = getattr(event, "xdata", None), getattr(event, "ydata", None)
        if lon is None or lat is None:
            return True
        choice = self._survey_start_from_click(float(lat), float(lon))
        if choice is None:
            self.survey_status_label.setText("Click an end of an outside line. That point becomes line 1, point 1.")
            return True
        self.survey_start_edge, self.survey_start_vertex = choice
        self._stop_survey_pick_start()
        if hasattr(self, "_plot_survey_plan"):
            self._plot_survey_plan(preserve_view_limits=True)
        self._survey_refresh_controls()
        return True

    def _survey_keep_start_end_after_reverse(self):
        """Reversing the drawn line swaps port with starboard, so keep the same geographic end."""
        edge, vertex = self._survey_start_choice()
        if getattr(self, "survey_start_edge", None) not in ("port", "starboard"):
            return
        self.survey_start_edge = "starboard" if edge == "port" else "port"
        self.survey_start_vertex = "finish" if vertex == "start" else "start"

    def _survey_offset_records(self):
        """Offset lines in drawn-line order: (side, index, points). Index 1 is nearest the drawn line."""
        if not getattr(self, "survey_offsets_applied", False):
            return []
        points = getattr(self, "survey_center_points", None) or []
        geod = self._survey_geod()
        if len(points) < 2 or geod is None:
            return []
        spacing = float(getattr(self, "survey_offset_m", 0) or 0)
        count = int(getattr(self, "survey_offset_count", 1) or 0)
        if spacing <= 0 or count < 1:
            return []
        side = getattr(self, "survey_offset_side", "both")
        records = []
        for index in range(1, count + 1):
            distance = spacing * index
            if side in ("port", "both"):
                records.append((
                    "port",
                    index,
                    self._survey_offset_polyline(geod, points, -90.0, distance),
                ))
            if side in ("starboard", "both"):
                records.append((
                    "starboard",
                    index,
                    self._survey_offset_polyline(geod, points, 90.0, distance),
                ))
        return [record for record in records if len(record[2]) >= 2]

    def _toggle_survey_mask_draw(self):
        if getattr(self, "survey_mask_draw_mode", False):
            self._stop_survey_mask_draw(keep=len(getattr(self, "survey_mask_draft", None) or []) >= 3)
            if hasattr(self, "_plot_survey_plan"):
                self._plot_survey_plan(preserve_view_limits=True)
            return
        if getattr(self, "survey_draw_mode", False):
            self._stop_survey_draw(keep_line=len(getattr(self, "survey_center_points", None) or []) >= 2)
        self._stop_survey_edit(replot=False)
        self._stop_survey_pick_start()
        self.survey_mask_draw_mode = True
        self.survey_mask_draft = []
        self._survey_remove_temp_line()
        self.survey_draw_mask_btn.setText("Left-click to add, right-click to close")
        self.survey_draw_mask_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")
        if hasattr(self, "canvas_widget"):
            self.canvas_widget.setCursor(Qt.CursorShape.CrossCursor)
        self._survey_refresh_controls()
        if hasattr(self, "_plot_survey_plan"):
            self._plot_survey_plan(preserve_view_limits=True)

    def _stop_survey_mask_draw(self, keep=False):
        draft = list(getattr(self, "survey_mask_draft", None) or [])
        self.survey_mask_draw_mode = False
        self.survey_mask_draft = []
        self._survey_remove_temp_line()
        if keep and len(draft) >= 3:
            self.survey_mask_points = draft
            if hasattr(self, "_plot_survey_plan"):
                self._plot_survey_plan(preserve_view_limits=True)
        button = getattr(self, "survey_draw_mask_btn", None)
        if button is not None:
            button.setText("Draw Mask")
            button.setStyleSheet("")
        if hasattr(self, "canvas_widget") and not getattr(self, "survey_edit_mode", False):
            if getattr(self, "survey_pick_start_mode", False) and self._survey_tab_is_current():
                self.canvas_widget.setCursor(Qt.CursorShape.CrossCursor)
            else:
                self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
        self._survey_refresh_controls()

    def _remove_survey_mask(self):
        self.survey_mask_points = []
        self.survey_mask_draft = []
        self._stop_survey_mask_draw(keep=False)
        if hasattr(self, "_plot_survey_plan"):
            self._plot_survey_plan(preserve_view_limits=True)
        self._survey_refresh_controls()

    def _survey_mask_click(self, event):
        if not getattr(self, "survey_mask_draw_mode", False) or not self._survey_tab_is_current():
            return False
        button = int(getattr(event, "button", 1) or 1)
        draft = list(getattr(self, "survey_mask_draft", None) or [])
        if button == 3:
            if len(draft) >= 3:
                self._stop_survey_mask_draw(keep=True)
                if hasattr(self, "_plot_survey_plan"):
                    self._plot_survey_plan(preserve_view_limits=True)
            else:
                self.survey_status_label.setText("Add at least three corners, then right-click to close.")
            return True
        if button != 1:
            return True
        if hasattr(self, "_main_map_click_lonlat"):
            lon, lat = self._main_map_click_lonlat(event)
        else:
            lon, lat = getattr(event, "xdata", None), getattr(event, "ydata", None)
        if lon is None or lat is None:
            return True
        draft.append((float(lat), float(lon)))
        self.survey_mask_draft = draft
        self._survey_remove_temp_line()
        if hasattr(self, "_plot_survey_plan"):
            self._plot_survey_plan(preserve_view_limits=True)
        self._survey_refresh_controls()
        return True

    def _survey_mask_motion(self, event):
        draft = list(getattr(self, "survey_mask_draft", None) or [])
        ax = getattr(self, "ax", None)
        if not draft or ax is None:
            return False
        lon = getattr(event, "xdata", None)
        lat = getattr(event, "ydata", None)
        if lon is None or lat is None:
            return False
        start_lat, start_lon = draft[-1]
        line = getattr(self, "survey_temp_line", None)
        if line is not None and line in ax.lines:
            line.set_data([start_lon, lon], [start_lat, lat])
            return True
        self.survey_temp_line, = ax.plot(
            [start_lon, lon],
            [start_lat, lat],
            color="#6c3483",
            linewidth=1.8,
            linestyle="--",
        )
        return True

    def _survey_mask_ring(self):
        ring = []
        for lat, lon in getattr(self, "survey_mask_points", None) or []:
            try:
                ring.append((float(lat), float(lon)))
            except (TypeError, ValueError):
                continue
        return ring if len(ring) >= 3 else []

    def _survey_mask_pieces(self, points):
        """Portions of a line that lie inside the mask. The drawn line itself is unchanged."""
        line = [(float(lat), float(lon)) for lat, lon in (points or [])]
        if len(line) < 2:
            return []
        ring = self._survey_mask_ring()
        if not ring:
            return [line]
        try:
            from shapely.geometry import LineString as ShapelyLine
            from shapely.geometry import Polygon as ShapelyPolygon
        except ImportError:
            return [line]
        try:
            shape = ShapelyLine([(lon, lat) for lat, lon in line])
            polygon = ShapelyPolygon([(lon, lat) for lat, lon in ring])
            if not polygon.is_valid:
                polygon = polygon.buffer(0)
            if polygon.is_empty:
                return [line]
            hit = shape.intersection(polygon)
        except Exception:
            return [line]
        parts = []
        if hit.is_empty:
            return []
        geoms = [hit] if hit.geom_type == "LineString" else []
        if hit.geom_type == "MultiLineString":
            geoms = list(hit.geoms)
        elif hit.geom_type == "GeometryCollection":
            geoms = [geom for geom in hit.geoms if geom.geom_type == "LineString"]
        for geom in geoms:
            coords = [(float(lat), float(lon)) for lon, lat in geom.coords]
            if len(coords) >= 2:
                parts.append(coords)
        if not parts:
            return []
        try:
            from shapely.geometry import Point as ShapelyPoint
            parts.sort(key=lambda piece: shape.project(ShapelyPoint(piece[0][1], piece[0][0])))
        except Exception:
            pass
        return parts

    def _survey_clipped_line(self, points):
        pieces = self._survey_mask_pieces(points)
        if not pieces:
            return []
        return max(pieces, key=len)

    def _survey_geometry_lines(self):
        """Lines from the port edge to the starboard edge, each in drawn-line vertex order."""
        points = list(getattr(self, "survey_center_points", None) or [])
        if len(points) < 2:
            return []
        records = self._survey_offset_records()
        port = [record for record in records if record[0] == "port"]
        starboard = [record for record in records if record[0] == "starboard"]
        geometry = []
        for side, index, line in reversed(port):
            geometry.append(((side, index), line))
        geometry.append((("source", 0), points))
        for side, index, line in starboard:
            geometry.append(((side, index), line))
        clipped = []
        for key, line in geometry:
            piece = self._survey_clipped_line(line)
            if len(piece) >= 2:
                clipped.append((key, piece))
        return clipped

    def _survey_run_lines(self):
        """(name, points in run order, key). Names follow the sinuous survey order."""
        points = list(getattr(self, "survey_center_points", None) or [])
        if len(points) < 2:
            return []
        if not getattr(self, "survey_offsets_applied", False):
            piece = self._survey_clipped_line(points)
            if len(piece) < 2:
                return []
            return [("L1", piece, ("source", 0))]
        geometry = self._survey_geometry_lines()
        if not geometry:
            if self._survey_mask_ring():
                return []
            return [("L1", points, ("source", 0))]
        if not self._survey_start_is_chosen():
            return [("", list(line), key) for key, line in geometry]
        edge, vertex = self._survey_start_choice()
        if edge == "starboard":
            geometry = list(reversed(geometry))
        forward_first = vertex == "start"
        run = []
        for number, (key, line) in enumerate(geometry, start=1):
            forward = forward_first if number % 2 == 1 else not forward_first
            ordered = list(line) if forward else list(reversed(line))
            run.append((f"L{number}", ordered, key))
        return run

    def _survey_named_lines(self):
        return [(name, line) for name, line, _key in self._survey_run_lines() if name]

    def _survey_extra_mask_lines(self):
        """Pieces left inside the mask after the numbered line keeps the longest one."""
        if len(self._survey_mask_ring()) < 3:
            return []
        names = {key: name for name, _line, key in self._survey_run_lines()}
        sources = [(("source", 0), list(getattr(self, "survey_center_points", None) or []))]
        if getattr(self, "survey_offsets_applied", False):
            sources.extend(
                ((side, index), list(line))
                for side, index, line in self._survey_offset_records()
            )
        extras = []
        for key, line in sources:
            pieces = self._survey_mask_pieces(line)
            if len(pieces) < 2:
                continue
            longest = max(range(len(pieces)), key=lambda index: len(pieces[index]))
            name = names.get(key) or "survey"
            for index, piece in enumerate(pieces):
                if index == longest or len(piece) < 2:
                    continue
                extras.append((name, piece))
        return extras

    def _survey_color_for_key(self, key):
        side = key[0] if key else "source"
        if side == "port":
            return "#1e8449"
        if side == "starboard":
            return "#c0392b"
        return "#1f4e79"

    def _survey_point_label(self, line_name, point_number):
        """Point names are L, the line number, P, and the point number, as in L1P1."""
        return f"{line_name}P{int(point_number)}"

    def _survey_geod(self):
        module = pyproj
        if module is None:
            try:
                import pyproj as module
            except Exception:
                return None
        try:
            return module.Geod(ellps="WGS84")
        except Exception:
            return None

    def _survey_centerline_meters(self):
        points = getattr(self, "survey_center_points", None) or []
        geod = self._survey_geod()
        if len(points) < 2 or geod is None:
            return 0.0
        total = 0.0
        for (lat1, lon1), (lat2, lon2) in zip(points, points[1:]):
            try:
                _azimuth, _back, distance = geod.inv(lon1, lat1, lon2, lat2)
            except Exception:
                continue
            if distance >= 0.5:
                total += float(distance)
        return total

    def _survey_offset_polyline(self, geod, points, turn_deg, distance):
        """Shift every vertex so the copy stays parallel to the line."""
        azimuths = []
        for (lat1, lon1), (lat2, lon2) in zip(points, points[1:]):
            try:
                azimuth, _back, length = geod.inv(lon1, lat1, lon2, lat2)
            except Exception:
                azimuth, length = 0.0, 0.0
            if length < 0.05 and azimuths:
                azimuth = azimuths[-1]
            azimuths.append(azimuth % 360.0)
        if not azimuths:
            return []

        def normal(azimuth):
            return (azimuth + turn_deg) % 360.0

        def vector(bearing):
            radians = math.radians(bearing)
            return math.sin(radians), math.cos(radians)

        shifted = []
        last = len(points) - 1
        for index, (lat, lon) in enumerate(points):
            if index == 0:
                bearing = normal(azimuths[0])
                step = distance
            elif index == last:
                bearing = normal(azimuths[-1])
                step = distance
            else:
                east_in, north_in = vector(normal(azimuths[index - 1]))
                east_out, north_out = vector(normal(azimuths[index]))
                east = east_in + east_out
                north = north_in + north_out
                magnitude = math.hypot(east, north)
                if magnitude < 1e-6:
                    bearing = normal(azimuths[index])
                    step = distance
                else:
                    bearing = math.degrees(math.atan2(east / magnitude, north / magnitude)) % 360.0
                    step = distance / max(magnitude / 2.0, 0.3)
            try:
                lon2, lat2, _back = geod.fwd(lon, lat, bearing, step)
            except Exception:
                lon2, lat2 = lon, lat
            shifted.append((lat2, lon2))
        return shifted

    def _plot_survey_mapping_lines(self):
        """Draw the survey line, and offset lines after they have been added."""
        self.survey_line_artists = []
        self.survey_vertex_handles = []
        points = getattr(self, "survey_center_points", None) or []
        ax = getattr(self, "ax", None)
        if ax is None:
            return
        self._survey_plot_mask(ax)
        if not points:
            return
        lats = [point[0] for point in points]
        lons = [point[1] for point in points]
        if len(points) == 1:
            ax.plot(lons, lats, color="#1f4e79", marker="o", linestyle="None", markersize=6)
            self._survey_annotate_points("L1", points, "#1f4e79")
            return
        names = {key: (name, ordered) for name, ordered, key in self._survey_run_lines()}
        source_name, source_ordered = names.get(("source", 0), ("", list(points)))
        for piece in self._survey_mask_pieces(points):
            plotted = ax.plot(
                [point[1] for point in piece],
                [point[0] for point in piece],
                color="#1f4e79",
                linewidth=2.4,
                linestyle="-",
                marker="o",
                markersize=5,
                label=source_name or "_nolegend_",
            )
            self.survey_line_artists.extend(plotted)
            if source_name and hasattr(self, "_tag_line_label"):
                self._tag_line_label(plotted, source_name)
            source_name = ""
        if names.get(("source", 0), ("",))[0]:
            self._survey_annotate_points(names[("source", 0)][0], source_ordered, "#1f4e79")
        for side, index, line in self._survey_offset_records():
            color = self._survey_color_for_key((side, index))
            name, ordered = names.get((side, index), ("", list(line)))
            for piece in self._survey_mask_pieces(line):
                plotted = ax.plot(
                    [point[1] for point in piece],
                    [point[0] for point in piece],
                    color=color,
                    linewidth=2.0,
                    linestyle="--",
                    label=name or "_nolegend_",
                )
                self.survey_line_artists.extend(plotted)
                if name and hasattr(self, "_tag_line_label"):
                    self._tag_line_label(plotted, name)
                name = ""
            if ordered and names.get((side, index), ("",))[0]:
                self._survey_annotate_points(names[(side, index)][0], ordered, color)
        if getattr(self, "survey_pick_start_mode", False):
            for _edge, _vertex, (lat, lon) in self._survey_start_candidates():
                ax.plot(
                    [lon],
                    [lat],
                    linestyle="None",
                    marker="o",
                    markersize=12,
                    color="#e67e22",
                    markeredgecolor="#e67e22",
                    markerfacecolor="none",
                    markeredgewidth=2,
                    label="_nolegend_",
                )
                ax.annotate(
                    "Start",
                    (lon, lat),
                    xytext=(8, 8),
                    textcoords="offset points",
                    fontsize=8,
                    color="#e67e22",
                    weight="bold",
                )
        if getattr(self, "survey_edit_mode", False):
            self._survey_create_vertex_handles()

    def _survey_create_vertex_handles(self):
        ax = getattr(self, "ax", None)
        points = getattr(self, "survey_center_points", None) or []
        self.survey_vertex_handles = []
        if ax is None or len(points) < 2:
            return
        last = len(points) - 1
        for index, (lat, lon) in enumerate(points):
            if index == 0:
                color = "#c0392b"
            elif index == last:
                color = "#1f4e79"
            else:
                color = "#1e8449"
            handle = ax.scatter(
                [lon],
                [lat],
                color=color,
                s=90,
                marker="o",
                edgecolors="black",
                linewidths=1.5,
                zorder=10,
                picker=5,
            )
            self.survey_vertex_handles.append(handle)

    def _survey_plot_mask(self, ax):
        """Draw the finished mask, and the corners of a new mask while it is being drawn."""
        finished = self._survey_mask_ring()
        draft = [
            (float(lat), float(lon))
            for lat, lon in (getattr(self, "survey_mask_draft", None) or [])
        ]
        if getattr(self, "survey_mask_draw_mode", False):
            self._survey_plot_mask_ring(ax, finished, closed=True)
            self._survey_plot_mask_ring(ax, draft, closed=False)
            return
        self._survey_plot_mask_ring(ax, finished, closed=True)

    def _survey_plot_mask_ring(self, ax, ring, closed):
        if len(ring) == 1 and not closed:
            ax.plot(
                [ring[0][1]],
                [ring[0][0]],
                linestyle="None",
                marker="o",
                color="#6c3483",
                markersize=6,
                zorder=3,
                label="_nolegend_",
            )
            return
        if len(ring) < 2:
            return
        lons = [point[1] for point in ring]
        lats = [point[0] for point in ring]
        if closed:
            lons.append(ring[0][1])
            lats.append(ring[0][0])
        ax.plot(
            lons,
            lats,
            color="#6c3483",
            linewidth=2.0,
            linestyle="-" if closed else "--",
            marker="o",
            markersize=5,
            zorder=3,
            label="_nolegend_",
        )

    def _survey_annotate_points(self, line_name, line, color):
        """Label each vertex as the line name plus its point number."""
        ax = getattr(self, "ax", None)
        if ax is None:
            return
        for index, (lat, lon) in enumerate(line, start=1):
            ax.annotate(
                self._survey_point_label(line_name, index),
                (lon, lat),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=8,
                color=color,
                weight="bold",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7),
            )

    def _survey_profile_lines(self):
        """(name, points) in sinuous run order, after a starting point has been clicked."""
        if getattr(self, "survey_offsets_applied", False) and not self._survey_start_is_chosen():
            points = list(getattr(self, "survey_center_points", None) or [])
            return [("Drawn line", points)] if len(points) >= 2 else []
        return self._survey_named_lines()

    def _survey_profile_default_name(self):
        """The profile list opens on the line the user drew."""
        if getattr(self, "survey_offsets_applied", False) and not self._survey_start_is_chosen():
            return "Drawn line"
        for name, _points, key in self._survey_run_lines():
            if key == ("source", 0) and name:
                return name
        return "L1"

    def _on_survey_profile_select_changed(self, *_args):
        if not self._survey_tab_is_current():
            return
        combo = getattr(self, "survey_profile_select_combo", None)
        if combo is not None and combo.count():
            self.survey_profile_user_choice = combo.currentText().strip()
        self._draw_survey_selected_profile()

    def _sync_survey_profile_select_options(self):
        """Keep a chosen line, and otherwise show the line that was drawn."""
        combo = getattr(self, "survey_profile_select_combo", None)
        if combo is None:
            return
        lines = self._survey_profile_lines()
        options = [name for name, _points in lines] or ["L1"]
        choice = getattr(self, "survey_profile_user_choice", None)
        if choice not in options:
            choice = self._survey_profile_default_name()
        combo.blockSignals(True)
        combo.clear()
        combo.addItems(options)
        if choice in options:
            combo.setCurrentText(choice)
        else:
            combo.setCurrentIndex(0)
        combo.blockSignals(False)

    def _draw_survey_selected_profile(self):
        """Draw the main line profile, or the offset line chosen in the list."""
        if not hasattr(self, "profile_ax") or not hasattr(self, "profile_canvas"):
            return
        self._sync_survey_profile_select_options()
        if getattr(self, "survey_offsets_applied", False) and not self._survey_start_is_chosen():
            points = list(getattr(self, "survey_center_points", None) or [])
            name, key = "Drawn line", ("source", 0)
        else:
            combo = getattr(self, "survey_profile_select_combo", None)
            selection = combo.currentText().strip() if combo is not None and combo.count() else self._survey_profile_default_name()
            chosen = None
            for run_name, run_points, run_key in self._survey_run_lines():
                if run_name and run_name == selection:
                    chosen = (run_name, run_points, run_key)
                    break
            if chosen is None:
                default_name = self._survey_profile_default_name()
                for run_name, run_points, run_key in self._survey_run_lines():
                    if run_name == default_name:
                        chosen = (run_name, run_points, run_key)
                        break
            if chosen is None:
                lines = [(run_name, run_points, run_key) for run_name, run_points, run_key in self._survey_run_lines() if run_name]
                chosen = lines[0] if lines else ("L1", [], ("source", 0))
            name, points, key = chosen
        color = self._survey_color_for_key(key)
        self.profile_ax.clear()
        self.profile_ax.set_axis_on()
        for axis in self.profile_fig.get_axes():
            if axis is not self.profile_ax:
                axis.remove()
        self.profile_ax.set_title(f"{name} Elevation Profile", fontsize=8)
        self.profile_ax.set_xlabel("Distance (m)", fontsize=8)
        self.profile_ax.set_ylabel("Elevation (m)", fontsize=8)
        self.profile_ax.tick_params(axis="both", which="major", labelsize=7)
        profile = self._line_plan_profile_arrays(points) if hasattr(self, "_line_plan_profile_arrays") else None
        if profile is None:
            self._update_survey_info_labels(points)
            self._show_profile_planning_placeholder()
            return
        dists, elevations, slopes, waypoint_distances = profile
        self.profile_ax.plot(dists, elevations, color=color, lw=1, label="Elevation")
        waypoint_elevations = []
        for distance in waypoint_distances:
            index = int(np.argmin(np.abs(dists - distance)))
            waypoint_elevations.append(elevations[index])
        self.profile_ax.plot(
            waypoint_distances,
            waypoint_elevations,
            "o",
            color=color,
            markersize=6,
            label="Vertices",
        )
        for index, (distance, elevation) in enumerate(zip(waypoint_distances, waypoint_elevations), start=1):
            if not str(name).startswith("L"):
                break
            self.profile_ax.annotate(
                self._survey_point_label(name, index),
                (distance, elevation),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=7,
                color=color,
            )
        slope_ax = None
        if self._show_slope_on_profile() and np.any(~np.isnan(slopes)):
            slope_ax = self.profile_ax.twinx()
            slope_ax.plot(dists, slopes, color="teal", lw=1, linestyle="--", label="Slope (deg)")
            slope_ax.set_ylabel("Slope (deg)", fontsize=8)
            slope_ax.tick_params(axis="y", labelsize=7)
            slope_ax.grid(False)
            slope_ax.set_ylim(0, np.nanmax(slopes[~np.isnan(slopes)]) * 1.1)
        if len(dists) > 1 and np.isfinite(dists[0]) and np.isfinite(dists[-1]):
            self.profile_ax.set_xlim(dists[0], dists[-1])
        if np.any(~np.isnan(elevations)):
            self.profile_ax.set_ylim(np.nanmin(elevations), np.nanmax(elevations))
        handles, labels = self.profile_ax.get_legend_handles_labels()
        if slope_ax is not None:
            slope_handles, slope_labels = slope_ax.get_legend_handles_labels()
            handles.extend(slope_handles)
            labels.extend(slope_labels)
        if handles:
            self.profile_ax.legend(handles, labels, fontsize=10, loc="upper right")
        self._remember_profile_track(
            dists,
            self._last_profile_sample_lats,
            self._last_profile_sample_lons,
            elevations,
            slopes,
        )
        self.profile_fig.tight_layout(pad=1.0)
        self.profile_canvas.draw_idle()
        self._update_survey_info_labels(points)

    def _zoom_to_survey_plan(self):
        points = list(getattr(self, "survey_center_points", None) or [])
        for _side, _index, line in self._survey_offset_records():
            points.extend(line)
        if len(points) < 2 or not hasattr(self, "_apply_map_zoom_limits_and_reload_geotiff"):
            return
        lats = [point[0] for point in points]
        lons = [point[1] for point in points]
        pad_lat = max((max(lats) - min(lats)) * 0.15, 0.002)
        pad_lon = max((max(lons) - min(lons)) * 0.15, 0.002)
        self._apply_map_zoom_limits_and_reload_geotiff(
            (min(lons) - pad_lon, max(lons) + pad_lon),
            (min(lats) - pad_lat, max(lats) + pad_lat),
        )

    def _survey_export_name(self):
        entry = getattr(self, "survey_export_name_entry", None)
        name = entry.text().strip() if entry is not None else ""
        if name:
            return name
        name = f"Survey_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
        if entry is not None:
            entry.setText(name)
        return name

    def _survey_export_lines(self):
        """Numbered lines in sinuous run order."""
        return self._survey_named_lines()

    def _survey_export_rows(self, lines):
        rows = []
        for line_number, (name, line) in enumerate(lines, start=1):
            for point_number, (lat, lon) in enumerate(line, start=1):
                rows.append((line_number, name, self._survey_point_label(name, point_number), lat, lon))
        return rows

    def _export_mapping_survey(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded. Cannot export the survey.")
            return
        if getattr(self, "survey_offsets_applied", False) and not self._survey_start_is_chosen():
            self._show_message("warning", "Choose a Start", "Click the starting point of the survey before exporting.")
            return
        lines = self._survey_export_lines()
        if not lines:
            self._show_message("warning", "No Line", "Draw a line with at least two points before exporting.")
            return
        export_dir = self._select_export_directory(getattr(self, "last_used_dir", None))
        if not export_dir:
            return
        self.last_used_dir = export_dir
        if hasattr(self, "_save_last_used_dir"):
            self._save_last_used_dir()
        export_name = self._survey_export_name()
        try:
            exported_geotiff_path = (
                self._maybe_export_survey_geotiff(export_dir, export_name)
                if hasattr(self, "_maybe_export_survey_geotiff")
                else None
            )
            params_geotiff_path = (
                self._resolve_export_params_geotiff_path(exported_geotiff_path, export_dir=export_dir)
                if hasattr(self, "_resolve_export_params_geotiff_path")
                else getattr(self, "current_geotiff_path", None)
            )
            enabled = self._export_type_enabled if hasattr(self, "_export_type_enabled") else (lambda _key: True)
            export_shapefile = enabled("esri_shapefile")
            export_gpkg = enabled("gpkg")
            export_sis = enabled("sis_asciiplan")
            export_gpx = enabled("gpx")
            export_text_csv = enabled("text_csv")
            export_text_txt = enabled("text_txt")
            export_hypack = enabled("hypack_lnw")
            export_kml = enabled("kml")
            export_map_png = self._export_map_png_enabled() if hasattr(self, "_export_map_png_enabled") else True
            rows = self._survey_export_rows(lines)
            written = []

            def note(path):
                if path and os.path.exists(path):
                    written.append(os.path.basename(path))

            csv_path = os.path.join(export_dir, f"{export_name}_DDD.csv")
            dmm_path = os.path.join(export_dir, f"{export_name}_DMM.csv")
            dms_path = os.path.join(export_dir, f"{export_name}_DMS.csv")
            txt_path = os.path.join(export_dir, f"{export_name}_DDD.txt")
            dmm_txt_path = os.path.join(export_dir, f"{export_name}_DMM.txt")
            dms_txt_path = os.path.join(export_dir, f"{export_name}_DMS.txt")
            if export_text_csv:
                export_utils.write_ddd_csv(csv_path, rows)
                export_utils.write_dmm_csv(dmm_path, rows)
                export_utils.write_dms_csv(dms_path, rows)
                note(csv_path)
                note(dmm_path)
                note(dms_path)
            if export_text_txt:
                export_utils.write_ddd_txt(txt_path, rows)
                export_utils.write_dmm_txt(dmm_txt_path, rows)
                export_utils.write_dms_txt(dms_txt_path, rows)
                note(txt_path)
                note(dmm_txt_path)
                note(dms_txt_path)

            collection_props = self._survey_export_properties(params_geotiff_path, export_dir)
            features = []
            shapely_features = []
            for name, line in lines:
                props = {
                    "name": name,
                    "role": "survey",
                }
                props.update(collection_props)
                if LineString and _shapely_mapping:
                    geometry = _shapely_mapping(LineString([(lon, lat) for lat, lon in line]))
                    features.append({"type": "Feature", "geometry": geometry, "properties": props})
                    shapely_features.append({"geometry": geometry, "properties": {"name": name}})
            for name, line in self._survey_extra_mask_lines():
                props = {"name": name, "role": "survey"}
                props.update(collection_props)
                if LineString and _shapely_mapping:
                    geometry = _shapely_mapping(LineString([(lon, lat) for lat, lon in line]))
                    features.append({"type": "Feature", "geometry": geometry, "properties": props})
                    shapely_features.append({"geometry": geometry, "properties": {"name": name}})
            mask_ring = self._survey_mask_ring()
            if mask_ring:
                closed = [[float(lon), float(lat)] for lat, lon in mask_ring]
                if closed[0] != closed[-1]:
                    closed.append(closed[0])
                features.append({
                    "type": "Feature",
                    "geometry": {"type": "Polygon", "coordinates": [closed]},
                    "properties": {"name": "mask", "role": "mask"},
                })
            geojson_path = os.path.join(export_dir, f"{export_name}.geojson")
            export_utils.remove_export_file(geojson_path)
            with open(geojson_path, "w", encoding="utf-8") as handle:
                json.dump(
                    {"type": "FeatureCollection", "properties": collection_props, "features": features},
                    handle,
                    indent=2,
                )
            note(geojson_path)

            shapefile_path = os.path.join(export_dir, f"{export_name}.shp")
            if shapely_features and (export_shapefile or export_gpkg) and fiona and _shapely_mapping:
                schema = {"geometry": "LineString", "properties": {"name": "str"}}
                if export_shapefile:
                    export_utils.remove_export_file(shapefile_path)
                    with fiona.open(shapefile_path, "w", driver="ESRI Shapefile", crs="EPSG:4326", schema=schema) as collection:
                        collection.writerecords(shapely_features)
                    note(shapefile_path)
                if hasattr(self, "_write_gpkg_if_enabled"):
                    self._write_gpkg_if_enabled(shapefile_path, schema, shapely_features, crs="EPSG:4326")
                    gpkg_path = os.path.splitext(shapefile_path)[0] + ".gpkg"
                    note(gpkg_path)
            self._survey_export_mask(export_dir, export_name, export_shapefile, export_gpkg, note)

            named_lines = [(name, list(line)) for name, line in lines]
            if export_hypack:
                zone, hem = export_utils.compute_utm_zone_from_points(lines[0][1])
                suffix = f"_UTM{zone}{'N' if hem == 'North' else 'S'}"
                lnw_path = os.path.join(export_dir, f"{export_name}{suffix}.lnw")
                if export_utils.write_lnw(lnw_path, named_lines):
                    note(lnw_path)
            if export_sis:
                sis_path = os.path.join(export_dir, f"{export_name}.asciiplan")
                export_utils.write_asciiplan(sis_path, named_lines)
                note(sis_path)
            if export_gpx:
                gpx_path = os.path.join(export_dir, f"{export_name}.gpx")
                if export_utils.write_gpx(gpx_path, named_lines, creator="SAT Planner Survey"):
                    note(gpx_path)
            if export_kml:
                kml_path = os.path.join(export_dir, f"{export_name}.kml")
                if export_utils.write_kml(kml_path, named_lines, document_name=export_name):
                    note(kml_path)
            if export_map_png and hasattr(self, "_save_export_map_png"):
                if hasattr(self, "_hide_map_hover_tooltip_for_export"):
                    self._hide_map_hover_tooltip_for_export()
                map_path = os.path.join(export_dir, f"{export_name}_map.png")
                self._save_export_map_png(map_path, dpi=300, bbox_inches="tight", facecolor="white")
                if hasattr(self, "_map_png_export_basenames"):
                    written.extend(self._map_png_export_basenames(map_path))
                else:
                    note(map_path)
            info_path = os.path.join(export_dir, f"{export_name}_info.txt")
            self._survey_write_info_file(info_path, lines)
            note(info_path)
            params_path = os.path.join(export_dir, f"{export_name}_params.json")
            export_utils.remove_export_file(params_path)
            with open(params_path, "w", encoding="utf-8") as handle:
                json.dump(collection_props, handle, indent=2)
            note(params_path)
            if exported_geotiff_path and os.path.isfile(exported_geotiff_path):
                note(exported_geotiff_path)
            self.survey_status_label.setText(f"Exported {export_name}.")
            self._show_message(
                "info",
                "Export Survey",
                "Survey exported to:\n" + "\n".join(f"- {name}" for name in written),
            )
        except Exception as exc:
            self._show_message("error", "Export Error", f"Failed to export the survey.\n{exc}")

    def _survey_write_info_file(self, path, lines):
        """Write the full survey: every line in run order, transits, and all points."""
        export_utils.remove_export_file(path)
        parts = ["SURVEY STATISTICS\n" + "=" * 30 + "\n"]
        if getattr(self, "survey_offsets_applied", False):
            side = getattr(self, "survey_offset_side", "both") or "both"
            spacing = float(getattr(self, "survey_offset_m", 0) or 0)
            count = int(getattr(self, "survey_offset_count", 1) or 1)
            parts.append(f"Offset spacing: {spacing:.1f} m\n")
            parts.append(f"Offset lines per side: {count}\n")
            parts.append(f"Offset side: {side}\n\n")
        parts.append(self._survey_full_info_text(lines))
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("".join(parts))

    def _survey_speed_knots(self):
        entry = getattr(self, "survey_speed_entry", None)
        try:
            speed = float(entry.text()) if entry is not None and entry.text().strip() else 8.0
        except (TypeError, ValueError):
            speed = 8.0
        return speed if speed > 0 else 8.0

    def _survey_selected_profile_line(self):
        """The line named in the profile list, or the line that was drawn."""
        lines = self._survey_profile_lines()
        combo = getattr(self, "survey_profile_select_combo", None)
        selection = combo.currentText().strip() if combo is not None and combo.count() else ""
        for name, points in lines:
            if name == selection:
                return name, points
        default_name = self._survey_profile_default_name()
        for name, points in lines:
            if name == default_name:
                return name, points
        if lines:
            return lines[0]
        return "L1", []

    def _survey_info_lines(self):
        """Every survey line, including the one that was drawn."""
        lines = []
        for _key, points in self._survey_geometry_lines():
            if len(points) >= 2:
                lines.append(list(points))
        if lines:
            return lines
        if self._survey_mask_ring():
            return []
        points = list(getattr(self, "survey_center_points", None) or [])
        return [points] if len(points) >= 2 else []

    def _update_survey_info_labels(self, points=None):
        """Shallowest depth, maximum depth, and slope range across every survey line."""
        del points
        shallow = getattr(self, "survey_shallowest_depth_label", None)
        maximum = getattr(self, "survey_max_depth_label", None)
        slope_min = getattr(self, "survey_min_slope_label", None)
        slope_max = getattr(self, "survey_max_slope_label", None)
        if shallow is None:
            return
        elevations_all = []
        slopes_all = []
        if hasattr(self, "_line_plan_profile_arrays"):
            for line in self._survey_info_lines():
                profile = self._line_plan_profile_arrays(line)
                if profile is None:
                    continue
                _dists, elevations, slopes, _waypoints = profile
                elevations_all.append(np.asarray(elevations, dtype=float))
                slopes_all.append(np.asarray(slopes, dtype=float))
        stats = None
        if elevations_all:
            elevations = np.concatenate(elevations_all)
            slopes = np.concatenate(slopes_all)
            if np.any(~np.isnan(elevations)):
                finite = np.isfinite(slopes) & ~np.isnan(slopes)
                stats = {
                    "shallowest": abs(float(np.nanmax(elevations))),
                    "maximum": abs(float(np.nanmin(elevations))),
                    "slope_min_deg": float(np.min(slopes[finite])) if np.any(finite) else None,
                    "slope_max_deg": float(np.max(slopes[finite])) if np.any(finite) else None,
                }
        shallow.setText(f"{stats['shallowest']:.2f}" if stats else "-")
        if maximum is not None:
            maximum.setText(f"{stats['maximum']:.2f}" if stats else "-")
        if slope_min is not None:
            slope_min.setText(
                f"{stats['slope_min_deg']:.2f}" if stats and stats["slope_min_deg"] is not None else "-"
            )
        if slope_max is not None:
            slope_max.setText(
                f"{stats['slope_max_deg']:.2f}" if stats and stats["slope_max_deg"] is not None else "-"
            )

    def _show_survey_information(self):
        if getattr(self, "survey_offsets_applied", False) and not self._survey_start_is_chosen():
            self._show_message("warning", "Choose a Start", "Click the starting point of the survey before showing survey info.")
            return
        lines = self._survey_export_lines()
        if not lines:
            self._show_message("warning", "No Line Data", "Draw a line before showing survey info.")
            return
        show_statistics_dialog(self, "Survey Info", self._survey_full_info_text(lines))

    def _survey_pair_metrics(self, start, end):
        """Geodesic distance in meters and heading in degrees from start to end."""
        geod = self._survey_geod()
        if geod is None or start is None or end is None:
            return 0.0, 0.0
        try:
            azimuth, _back, distance = geod.inv(start[1], start[0], end[1], end[0])
        except Exception:
            return 0.0, 0.0
        if distance != distance:
            return 0.0, 0.0
        return float(distance), float(azimuth) % 360.0

    def _survey_minutes(self, meters, speed_knots):
        speed_m_per_s = float(speed_knots) * 0.514444
        if speed_m_per_s <= 0:
            return 0.0
        return (float(meters) / speed_m_per_s) / 60.0

    def _survey_format_distance(self, meters):
        return f"{meters:.1f} m ({meters / 1000.0:.3f} km, {meters / 1852.0:.3f} nm)"

    def _survey_format_time(self, minutes):
        return f"{minutes:.1f} min ({minutes / 60.0:.2f} hr)"

    def _survey_full_info_text(self, lines):
        """Full survey: line times, transits between lines, total time, and every point in order."""
        speed_knots = self._survey_speed_knots()
        named = [(name, list(points)) for name, points in lines if len(points) >= 2]
        parts = [
            "SURVEY INFO\n" + "=" * 30 + "\n",
            f"Number of Lines: {len(named)}\n",
            f"Survey Speed: {speed_knots:.1f} knots\n",
        ]
        line_rows = []
        line_distance = 0.0
        for name, points in named:
            segments = []
            distance = 0.0
            for start, end in zip(points, points[1:]):
                length, heading = self._survey_pair_metrics(start, end)
                segments.append((length, heading))
                distance += length
            line_distance += distance
            line_rows.append((name, points, segments, distance, self._survey_minutes(distance, speed_knots)))
        line_minutes = self._survey_minutes(line_distance, speed_knots)
        parts.append("\nLINE RUNNING\n" + "-" * 30 + "\n")
        for name, _points, _segments, distance, minutes in line_rows:
            parts.append(f"{name}: {self._survey_format_distance(distance)} — {self._survey_format_time(minutes)}\n")
        parts.append(f"Line running time: {self._survey_format_time(line_minutes)}\n")
        transits = []
        transit_distance = 0.0
        for (name, points, _segments, _distance, _minutes), (next_name, next_points, *_rest) in zip(line_rows, line_rows[1:]):
            length, heading = self._survey_pair_metrics(points[-1], next_points[0])
            minutes = self._survey_minutes(length, speed_knots)
            transits.append((name, next_name, length, heading, minutes))
            transit_distance += length
        transit_minutes = self._survey_minutes(transit_distance, speed_knots)
        parts.append("\nTRANSITS BETWEEN LINES\n" + "-" * 30 + "\n")
        if transits:
            for name, next_name, length, heading, minutes in transits:
                parts.append(
                    f"{name} to {next_name}: {self._survey_format_distance(length)}"
                    f" — Heading {heading:.1f}° — {self._survey_format_time(minutes)}\n"
                )
        else:
            parts.append("No transits. The survey is one line.\n")
        parts.append(f"Transit time: {self._survey_format_time(transit_minutes)}\n")
        total_minutes = line_minutes + transit_minutes
        parts.append("\nTOTAL SURVEY TIME\n" + "-" * 30 + "\n")
        parts.append(f"Line running: {self._survey_format_time(line_minutes)}\n")
        parts.append(f"Transits: {self._survey_format_time(transit_minutes)}\n")
        parts.append(f"Total: {self._survey_format_time(total_minutes)}\n")
        depth_lines = [line for _name, line in named]
        depth = self._survey_lines_depth_info(depth_lines)
        if depth:
            parts.append("\nDEPTH AND SLOPE (ALL LINES)\n" + "-" * 30 + "\n")
            parts.append(depth if depth.endswith("\n") else depth + "\n")
        for name, _points, segments, _distance, _minutes in line_rows:
            if not segments:
                continue
            parts.append(f"\n{name} SEGMENTS\n" + "-" * 30 + "\n")
            for index, (distance, heading) in enumerate(segments, start=1):
                parts.append(
                    f"Segment {index}: {self._survey_format_distance(distance)} — Heading: {heading:.1f}°\n"
                )
        dmm_heading = "Points in Survey Order (DMM)"
        parts.append(f"\n{dmm_heading}\n" + "-" * len(dmm_heading) + "\n")
        for name, points, _segments, _distance, _minutes in line_rows:
            for index, (lat, lon) in enumerate(points, start=1):
                label = self._survey_point_label(name, index)
                parts.append(
                    f"{label}: {decimal_degrees_to_ddm(lat, True)}, {decimal_degrees_to_ddm(lon, False)}\n"
                )
        ddd_heading = "Points in Survey Order (DDD)"
        parts.append(f"\n{ddd_heading}\n" + "-" * len(ddd_heading) + "\n")
        for name, points, _segments, _distance, _minutes in line_rows:
            for index, (lat, lon) in enumerate(points, start=1):
                label = self._survey_point_label(name, index)
                parts.append(f"{label}: {lat:.6f}, {lon:.6f}\n")
        return "".join(parts)

    def _survey_lines_depth_info(self, lines):
        """Shallowest depth, maximum depth, and slope range across every line."""
        if not hasattr(self, "_line_plan_profile_arrays"):
            return ""
        elevations_all = []
        slopes_all = []
        for points in lines:
            profile = self._line_plan_profile_arrays(points)
            if profile is None:
                continue
            _dists, elevations, slopes, _waypoints = profile
            elevations_all.append(np.asarray(elevations, dtype=float))
            slopes_all.append(np.asarray(slopes, dtype=float))
        if not elevations_all:
            return ""
        elevations = np.concatenate(elevations_all)
        slopes = np.concatenate(slopes_all)
        if not np.any(~np.isnan(elevations)):
            return ""
        finite = np.isfinite(slopes) & ~np.isnan(slopes)
        stats = {
            "shallowest": abs(float(np.nanmax(elevations))),
            "maximum": abs(float(np.nanmin(elevations))),
            "slope_min_deg": float(np.min(slopes[finite])) if np.any(finite) else None,
            "slope_max_deg": float(np.max(slopes[finite])) if np.any(finite) else None,
            "valid_count": int(np.sum(~np.isnan(elevations))),
        }
        if hasattr(self, "_format_line_survey_info_text"):
            return self._format_line_survey_info_text(stats)
        return ""

    def _survey_line_info_block(self, points, speed_knots, line_name="L1"):
        """Distance, time, depth, headings, and vertices for one survey line."""
        geod = self._survey_geod()
        segment_distances = []
        headings = []
        total = 0.0
        if geod is not None:
            for (lat1, lon1), (lat2, lon2) in zip(points, points[1:]):
                try:
                    azimuth, _back, distance = geod.inv(lon1, lat1, lon2, lat2)
                except Exception:
                    continue
                segment_distances.append(float(distance))
                headings.append(float(azimuth) % 360.0)
                total += float(distance)
        speed_m_per_s = speed_knots * 0.514444
        minutes = (total / speed_m_per_s) / 60.0 if speed_m_per_s > 0 else 0.0
        lines = [
            f"Number of Points: {len(points)}\n",
            f"Total Distance: {total:.1f} m ({total / 1000.0:.3f} km, {total / 1852.0:.3f} nm)\n",
            f"Survey Speed: {speed_knots:.1f} knots\n",
            f"Estimated Time: {minutes:.1f} min ({minutes / 60.0:.2f} hr)\n\n",
        ]
        depth = self._survey_line_depth_info(points)
        if depth:
            lines.append(depth)
            if not depth.endswith("\n"):
                lines.append("\n")
        if segment_distances:
            lines.append("\nSEGMENT DISTANCE AND HEADINGS\n" + "-" * 30 + "\n")
            for index, (distance, heading) in enumerate(zip(segment_distances, headings), start=1):
                lines.append(
                    f"Segment {index}: {distance:.1f} m ({distance / 1000.0:.3f} km, {distance / 1852.0:.3f} nm)"
                    f" - Heading: {heading:.1f}°\n"
                )
        dmm_heading = "Line Waypoints (DMM)"
        lines.append(f"\n{dmm_heading}\n" + "-" * len(dmm_heading) + "\n")
        for index, (lat, lon) in enumerate(points, start=1):
            label = self._survey_point_label(line_name, index)
            lines.append(
                f"{label}: {decimal_degrees_to_ddm(lat, True)}, {decimal_degrees_to_ddm(lon, False)}\n"
            )
        ddd_heading = "Line Waypoints (DDD)"
        lines.append(f"\n{ddd_heading}\n" + "-" * len(ddd_heading) + "\n")
        for index, (lat, lon) in enumerate(points, start=1):
            label = self._survey_point_label(line_name, index)
            lines.append(f"{label}: {lat:.6f}, {lon:.6f}\n")
        return "".join(lines)

    def _survey_line_depth_info(self, points):
        if not hasattr(self, "_line_plan_profile_arrays"):
            return ""
        profile = self._line_plan_profile_arrays(points)
        if profile is None:
            return ""
        _dists, elevations, slopes, _waypoints = profile
        elevations = np.asarray(elevations, dtype=float)
        slopes = np.asarray(slopes, dtype=float)
        if not np.any(~np.isnan(elevations)):
            return ""
        finite = np.isfinite(slopes) & ~np.isnan(slopes)
        stats = {
            "shallowest": abs(float(np.nanmax(elevations))),
            "maximum": abs(float(np.nanmin(elevations))),
            "slope_min_deg": float(np.min(slopes[finite])) if np.any(finite) else None,
            "slope_max_deg": float(np.max(slopes[finite])) if np.any(finite) else None,
            "valid_count": int(np.sum(~np.isnan(elevations))),
        }
        if hasattr(self, "_format_line_survey_info_text"):
            return self._format_line_survey_info_text(stats)
        return ""

    def _survey_export_properties(self, geotiff_path, export_dir):
        properties = {
            "plan_type": "survey",
            "offset_m": float(getattr(self, "survey_offset_m", 100.0) or 100.0),
            "offset_count": int(getattr(self, "survey_offset_count", 1) or 1),
            "offset_side": getattr(self, "survey_offset_side", "both") or "both",
            "offsets_applied": bool(getattr(self, "survey_offsets_applied", False)),
            "start_edge": getattr(self, "survey_start_edge", None),
            "start_vertex": getattr(self, "survey_start_vertex", None),
            "centerline_latlon": [
                [float(lat), float(lon)]
                for lat, lon in (getattr(self, "survey_center_points", None) or [])
            ],
            "mask_latlon": [
                [float(lat), float(lon)]
                for lat, lon in self._survey_mask_ring()
            ],
            "geotiff_path": geotiff_path,
            "geotiff_nan_value": float(getattr(self, "geotiff_nan_value", -99999.0)),
        }
        if hasattr(self, "_stamp_topography_export_path"):
            self._stamp_topography_export_path(properties, export_dir)
        if hasattr(self, "_add_geotiff_viz_params_to_params"):
            self._add_geotiff_viz_params_to_params(properties)
        return properties

    def _import_mapping_survey(self):
        if hasattr(self, "_gmrt_is_downloading") and self._gmrt_is_downloading():
            self._gmrt_cancel_active_download()
            return
        start_dir = getattr(self, "last_survey_import_dir", None) or getattr(self, "last_used_dir", None)
        file_path, _selected = QFileDialog.getOpenFileName(
            self,
            "Select Survey File to Import",
            start_dir or "",
            "Known Survey Files (*_DMS.txt *_DMM.txt *_DDD.txt *.lnw *_DDD.csv *_DMM.csv *_DMS.csv *.csv *.geojson *.json *.gpx *.shp *.gpkg);;"
            "Decimal Degrees (*_DDD.txt);;Degrees Minutes Seconds (*_DMS.txt);;Degrees Decimal Minutes (*_DMM.txt);;"
            "Hypack LNW (*.lnw);;Shapefile (*.shp);;GeoPackage (*.gpkg);;Decimal Degree CSV (*_DDD.csv);;DMM CSV (*_DMM.csv);;DMS CSV (*_DMS.csv);;CSV (*.csv);;GeoJSON (*.geojson);;JSON (*.json);;GPX (*.gpx);;All files (*.*)",
        )
        if not file_path:
            return
        import_dir = os.path.dirname(file_path)
        if import_dir and os.path.isdir(import_dir):
            self.last_survey_import_dir = import_dir
        try:
            points, settings = self._survey_points_from_file(file_path)
            if points is None:
                return
            if len(points) < 2:
                self._show_message("warning", "Import Warning", "The survey line needs at least two points.")
                return
            self._stop_survey_edit(replot=False)
            self._stop_survey_draw(keep_line=False)
            saved_line = self._survey_centerline_from_settings(settings)
            if saved_line:
                points = saved_line
            self.survey_center_points = points
            self.survey_offsets_applied = False
            self.survey_mask_points = []
            self.survey_mask_draft = []
            self._survey_clear_start_end()
            self._stop_survey_pick_start()
            self._survey_apply_saved_settings(settings)
            if self.survey_offsets_applied and not self._survey_start_is_chosen():
                self._survey_begin_pick_start()
            entry = getattr(self, "survey_export_name_entry", None)
            if entry is not None:
                entry.setText(self._survey_name_from_filename(file_path))
            if hasattr(self, "_plot_survey_plan"):
                self._plot_survey_plan(preserve_view_limits=False)
            if hasattr(self, "_zoom_to_survey_plan"):
                self._zoom_to_survey_plan()
            self._survey_refresh_controls()
            self.survey_status_label.setText(f"Imported a survey line with {len(points)} vertices.")
            geotiff_path = settings.get("geotiff_path") if isinstance(settings, dict) else None
            if geotiff_path and hasattr(self, "_resolve_import_sidecar_path"):
                geotiff_path = self._resolve_import_sidecar_path(geotiff_path, import_dir)
            if hasattr(self, "_maybe_prompt_gmrt_download_after_import"):
                self._maybe_prompt_gmrt_download_after_import(
                    geotiff_path=geotiff_path,
                    download_callback=self._download_and_load_gmrt_after_survey_import,
                )
        except Exception as exc:
            self._show_message("error", "Import Error", f"Failed to import the survey.\n{exc}")

    def _survey_name_from_filename(self, file_path):
        name = os.path.splitext(os.path.basename(file_path))[0]
        for suffix in ("_DDD", "_DMM", "_DMS", "_DD", "_DM", "_params", "_UTM"):
            if suffix == "_UTM" and "_UTM" in name:
                name = name.split("_UTM", 1)[0]
                break
            if name.endswith(suffix):
                name = name[: -len(suffix)]
        return name or f"Survey_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"

    def _survey_points_from_file(self, file_path):
        """Return (points, settings) for a survey or line-plan file."""
        file_name = os.path.basename(file_path).lower()
        extension = os.path.splitext(file_path)[1].lower()
        directory = os.path.dirname(file_path)
        settings = self._survey_sidecar_settings(file_path)
        points = None
        if file_name.endswith("_ddd.txt") and hasattr(self, "_parse_ddd_txt_file_as_polyline"):
            points = self._parse_ddd_txt_file_as_polyline(file_path)
        elif file_name.endswith("_dms.txt") and hasattr(self, "_parse_dms_txt_file_as_polyline"):
            points = self._parse_dms_txt_file_as_polyline(file_path)
        elif file_name.endswith("_dmm.txt") and hasattr(self, "_parse_dmm_txt_file_as_polyline"):
            points = self._parse_dmm_txt_file_as_polyline(file_path)
        elif extension == ".lnw" and UTMZoneDialog is not None and hasattr(self, "_parse_lnw_file_as_polyline"):
            dialog = UTMZoneDialog.for_file(self, file_path)
            if dialog.exec() != QDialog.DialogCode.Accepted:
                return None, settings
            zone, hemisphere = dialog.get_utm_info()
            points = self._parse_lnw_file_as_polyline(file_path, zone, hemisphere)
        elif extension == ".gpx" and hasattr(self, "_parse_gpx_file_as_polyline"):
            points = self._parse_gpx_file_as_polyline(file_path)
        elif extension in (".shp", ".gpkg") and hasattr(self, "_parse_vector_file_as_polyline"):
            points = self._parse_vector_file_as_polyline(file_path)
        elif extension == ".csv":
            points = self._survey_points_from_csv(file_path)
        elif extension in (".geojson", ".json"):
            points, geo_settings = self._survey_points_from_geojson(file_path)
            settings.update(geo_settings)
        else:
            self._show_message("error", "Import Error", f"Unsupported file format: {extension}")
            return None, settings
        if isinstance(settings, dict) and settings.get("geotiff_path") and hasattr(self, "_load_geotiff_from_path"):
            stored = settings.get("geotiff_path")
            if hasattr(self, "_resolve_import_sidecar_path"):
                stored = self._resolve_import_sidecar_path(stored, directory)
            if stored and os.path.isfile(stored):
                self._load_geotiff_from_path(stored)
            if hasattr(self, "_maybe_load_topography_geotiff"):
                self._maybe_load_topography_geotiff(settings.get("topography_geotiff_path"), directory)
            if hasattr(self, "_apply_geotiff_viz_params_from_params"):
                self._apply_geotiff_viz_params_from_params(settings)
        return points or [], settings

    def _survey_sidecar_settings(self, file_path):
        directory = os.path.dirname(file_path)
        stem = os.path.splitext(os.path.basename(file_path))[0]
        candidates = [stem]
        for suffix in ("_DDD", "_DMM", "_DMS", "_DD", "_DM"):
            if stem.endswith(suffix):
                candidates.append(stem[: -len(suffix)])
        if "_UTM" in stem:
            candidates.append(stem.split("_UTM", 1)[0])
        for candidate in candidates:
            path = os.path.join(directory, f"{candidate}_params.json")
            if not os.path.isfile(path):
                continue
            try:
                with open(path, encoding="utf-8") as handle:
                    payload = json.load(handle)
            except (OSError, json.JSONDecodeError):
                continue
            if isinstance(payload, dict):
                return payload
        return {}

    def _survey_points_from_csv(self, file_path):
        groups = {}
        order = []
        with open(file_path, encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                try:
                    lat, lon = lat_lon_decimal_from_survey_csv_row(row)
                except (ValueError, TypeError):
                    continue
                name = (row.get("Line Name") or row.get("line_name") or "ML").strip() or "ML"
                if name not in groups:
                    order.append(name)
                    groups[name] = []
                groups[name].append((lat, lon))
        for preferred in ("ML", "Line"):
            if preferred in groups and len(groups[preferred]) >= 2:
                return groups[preferred]
        for name in order:
            if len(groups[name]) >= 2:
                return groups[name]
        return []

    def _survey_points_from_geojson(self, file_path):
        with open(file_path, encoding="utf-8") as handle:
            payload = json.load(handle)
        if not isinstance(payload, dict):
            return [], {}
        settings = dict(payload.get("properties") or {})
        features = payload.get("features") or []
        if payload.get("type") == "Feature":
            features = [payload]
        chosen = None
        fallback = None
        for feature in features:
            if not isinstance(feature, dict) or feature.get("type") != "Feature":
                continue
            geometry = feature.get("geometry") or {}
            if geometry.get("type") != "LineString":
                continue
            coords = []
            for coord in geometry.get("coordinates") or []:
                if len(coord) >= 2:
                    coords.append((float(coord[1]), float(coord[0])))
            if len(coords) < 2:
                continue
            props = feature.get("properties") or {}
            if fallback is None:
                fallback = coords
            role = props.get("role")
            name = props.get("name")
            if role == "line" or name in ("ML", "Line"):
                chosen = coords
                break
            if role == "mask" or name == "mask":
                continue
        for feature in features:
            if not isinstance(feature, dict):
                continue
            props = feature.get("properties") or {}
            geometry = feature.get("geometry") or {}
            if props.get("role") != "mask" and props.get("name") != "mask":
                continue
            if geometry.get("type") != "Polygon":
                continue
            rings = geometry.get("coordinates") or []
            if not rings or not isinstance(rings[0], list):
                continue
            mask = []
            for coord in rings[0]:
                if len(coord) >= 2:
                    mask.append([float(coord[1]), float(coord[0])])
            if len(mask) >= 2 and mask[0] == mask[-1]:
                mask = mask[:-1]
            if len(mask) >= 3 and "mask_latlon" not in settings:
                settings["mask_latlon"] = mask
            break
        return chosen or fallback or [], settings

    def _survey_apply_saved_settings(self, settings):
        if not isinstance(settings, dict):
            return
        if settings.get("offset_m") is not None:
            try:
                self.survey_offset_m = float(settings.get("offset_m"))
                if hasattr(self, "survey_offset_spin"):
                    self.survey_offset_spin.blockSignals(True)
                    self.survey_offset_spin.setValue(self.survey_offset_m)
                    self.survey_offset_spin.blockSignals(False)
            except (TypeError, ValueError):
                pass
        if settings.get("offset_count") is not None:
            try:
                self.survey_offset_count = max(1, int(settings.get("offset_count")))
                if hasattr(self, "survey_count_spin"):
                    self.survey_count_spin.blockSignals(True)
                    self.survey_count_spin.setValue(self.survey_offset_count)
                    self.survey_count_spin.blockSignals(False)
            except (TypeError, ValueError):
                pass
        side = settings.get("offset_side")
        if side in ("port", "starboard", "both"):
            self.survey_offset_side = side
            combo = getattr(self, "survey_side_combo", None)
            if combo is not None:
                index = combo.findData(side)
                if index >= 0:
                    combo.blockSignals(True)
                    combo.setCurrentIndex(index)
                    combo.blockSignals(False)
        if "offsets_applied" in settings:
            self.survey_offsets_applied = bool(settings.get("offsets_applied"))
        edge = settings.get("start_edge")
        vertex = settings.get("start_vertex")
        if edge in ("port", "starboard"):
            self.survey_start_edge = edge
        if vertex in ("start", "finish"):
            self.survey_start_vertex = vertex
        self.survey_mask_points = self._survey_mask_from_settings(settings)
        self.survey_mask_draft = []

    def _survey_mask_from_settings(self, settings):
        if not isinstance(settings, dict):
            return []
        raw = settings.get("mask_latlon")
        if not isinstance(raw, list):
            return []
        points = []
        for item in raw:
            if not isinstance(item, (list, tuple)) or len(item) < 2:
                continue
            try:
                points.append((float(item[0]), float(item[1])))
            except (TypeError, ValueError):
                continue
        if len(points) >= 2 and points[0] == points[-1]:
            points = points[:-1]
        return points if len(points) >= 3 else []

    def _survey_export_mask(self, export_dir, export_name, export_shapefile, export_gpkg, note):
        ring = self._survey_mask_ring()
        if len(ring) < 3 or _shapely_mapping is None or fiona is None:
            return
        try:
            from shapely.geometry import Polygon as ShapelyPolygon
        except ImportError:
            return
        coordinates = [(float(lon), float(lat)) for lat, lon in ring]
        if coordinates[0] != coordinates[-1]:
            coordinates.append(coordinates[0])
        try:
            polygon = ShapelyPolygon(coordinates)
            if not polygon.is_valid:
                polygon = polygon.buffer(0)
            if polygon.geom_type == "MultiPolygon" and polygon.geoms:
                polygon = max(polygon.geoms, key=lambda geom: geom.area)
            if polygon.is_empty or polygon.geom_type != "Polygon":
                return
        except Exception:
            return
        schema = {"geometry": "Polygon", "properties": {"name": "str"}}
        feature = {"geometry": _shapely_mapping(polygon), "properties": {"name": "mask"}}
        if export_shapefile:
            path = os.path.join(export_dir, f"{export_name}_mask.shp")
            try:
                export_utils.remove_export_file(path)
                with fiona.open(path, "w", driver="ESRI Shapefile", crs="EPSG:4326", schema=schema) as collection:
                    collection.writerecords([feature])
                note(path)
            except Exception:
                pass
        if export_gpkg:
            path = os.path.join(export_dir, f"{export_name}.gpkg")
            try:
                existed = os.path.isfile(path)
                mode = "a" if existed else "w"
                with fiona.open(
                    path, mode, driver="GPKG", crs="EPSG:4326", schema=schema, layer="mask",
                ) as collection:
                    collection.writerecords([feature])
                if not existed:
                    note(path)
            except Exception:
                pass

    def _survey_centerline_from_settings(self, settings):
        """The drawn line saved with an export. The numbered plan does not identify it."""
        if not isinstance(settings, dict):
            return None
        raw = settings.get("centerline_latlon")
        if not isinstance(raw, list):
            return None
        points = []
        for item in raw:
            if not isinstance(item, (list, tuple)) or len(item) < 2:
                continue
            try:
                points.append((float(item[0]), float(item[1])))
            except (TypeError, ValueError):
                continue
        return points if len(points) >= 2 else None

    def _download_and_load_gmrt_after_survey_import(self):
        points = getattr(self, "survey_center_points", None) or []
        if len(points) < 2 or not hasattr(self, "_gmrt_download_extent_from_points"):
            return
        west, east, south, north = self._gmrt_download_extent_from_points(points, prefix="survey")
        split_topo = True
        if hasattr(self, "_gmrt_split_topo_depths_for_prefix"):
            split_topo = self._gmrt_split_topo_depths_for_prefix("survey")
        self._download_gmrt_and_load(
            west, east, south, north,
            resolution=self._gmrt_import_resolution_meters("survey"),
            layer="topo",
            default_filename_prefix="GMRT_Bathy",
            log_func=lambda message, append=True: self.survey_status_label.setText(message),
            default_directory=getattr(self, "last_survey_import_dir", None),
            split_topo_depths=split_topo,
            gmrt_button=self._active_import_button() if hasattr(self, "_active_import_button") else None,
        )

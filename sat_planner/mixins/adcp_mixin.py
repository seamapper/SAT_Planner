"""
ADCP tab: two-circle calibration planning (36 segments per circle, opposite directions).
"""

import csv
import datetime
import json
import os

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QTextCursor
from PyQt6.QtWidgets import QFileDialog, QDialog

from sat_planner import export_utils
from sat_planner.constants import GEOSPATIAL_LIBS_AVAILABLE, LineString, fiona, pyproj
from sat_planner.utils_geo import decimal_degrees_to_ddm, lat_lon_decimal_from_survey_csv_row
from sat_planner.utils_ui import show_statistics_dialog

try:
    from shapely.geometry import mapping
except ImportError:
    mapping = None

try:
    from sat_planner.mixins.survey_parsers_mixin import UTMZoneDialog
except ImportError:
    UTMZoneDialog = None

ADCP_NUM_SEGMENTS = 36
ADCP_DEFAULT_DIAMETER_M = 300.0
ADCP_ACTIVE_BTN_STYLE = "QPushButton { color: rgb(255, 165, 0); font-weight: bold; }"


class AdcpMixin:
    """Mixin for ADCP calibration circle planning, import/export, and map interaction."""

    def _adcp_get_diameter_m(self):
        try:
            d = float(self.adcp_circle_diameter_entry.text().strip())
            if d > 0 and np.isfinite(d):
                return float(d)
        except Exception:
            pass
        return ADCP_DEFAULT_DIAMETER_M

    def _adcp_get_radius_m(self):
        return self._adcp_get_diameter_m() / 2.0

    def _adcp_get_survey_speed_kts(self):
        try:
            v = float(self.adcp_survey_speed_entry.text().strip())
            if v >= 0 and np.isfinite(v):
                return float(v)
        except Exception:
            pass
        return 6.0

    def _adcp_get_turn_time_min(self):
        try:
            raw = self.adcp_turn_time_entry.text().strip()
            if not raw:
                return 5.0
            v = float(raw)
            if v >= 0 and np.isfinite(v):
                return float(v)
        except Exception:
            pass
        return 5.0

    def _build_adcp_export_basename(self):
        d_int = int(round(self._adcp_get_diameter_m()))
        date_str = datetime.datetime.now().strftime("%Y%m%d")
        return f"ADCP_Cal_Circle{d_int}m_{date_str}"

    def _update_adcp_export_name(self):
        if not hasattr(self, "adcp_export_name_entry"):
            return
        name = self._build_adcp_export_basename()
        if hasattr(self, "_deferred_set_line_edit"):
            self._deferred_set_line_edit(self.adcp_export_name_entry, name)
        else:
            self.adcp_export_name_entry.setText(name)

    def _adcp_set_deferred_field(self, attr_name, value):
        """Set an ADCP QLineEdit from import/code without pending (amber) styling."""
        if value is None:
            return
        widget = getattr(self, attr_name, None)
        if widget is None:
            return
        text = str(value).strip()
        if hasattr(self, "_deferred_set_line_edit"):
            self._deferred_set_line_edit(widget, text)
        else:
            widget.setText(text)

    def _adcp_geod(self):
        if not GEOSPATIAL_LIBS_AVAILABLE or pyproj is None:
            return None
        return pyproj.Geod(ellps="WGS84")

    def _adcp_snap_point_to_circle(self, center, click_lat, click_lon):
        """Return (lat, lon) on the circle through center at the current plan radius."""
        geod = self._adcp_geod()
        if geod is None or center is None:
            return None
        clat, clon = center
        radius_m = self._adcp_get_radius_m()
        az, _, _ = geod.inv(clon, clat, click_lon, click_lat)
        lon, lat, _ = geod.fwd(clon, clat, az, radius_m)
        return (lat, lon)

    def _adcp_build_circle_vertices(self, center, start_point, clockwise=True, num_segments=ADCP_NUM_SEGMENTS):
        """Return list of num_segments vertices (lat, lon) starting at start_point."""
        geod = self._adcp_geod()
        if geod is None or center is None or start_point is None:
            return []
        clat, clon = center
        radius_m = self._adcp_get_radius_m()
        snapped = self._adcp_snap_point_to_circle(center, start_point[0], start_point[1])
        if snapped is None:
            return []
        slat, slon = snapped
        az_start, _, _ = geod.inv(clon, clat, slon, slat)
        step = 360.0 / float(num_segments)
        vertices = []
        for i in range(num_segments):
            az = az_start + i * step if clockwise else az_start - i * step
            lon, lat, _ = geod.fwd(clon, clat, az % 360.0, radius_m)
            vertices.append((lat, lon))
        return vertices

    def _adcp_vertices_to_segments(self, vertices):
        segments = []
        n = len(vertices)
        if n < 2:
            return segments
        for i in range(n):
            segments.append([vertices[i], vertices[(i + 1) % n]])
        return segments

    def _adcp_build_circle_segments(self, center, start_point, clockwise=True):
        vertices = self._adcp_build_circle_vertices(center, start_point, clockwise=clockwise)
        if len(vertices) != ADCP_NUM_SEGMENTS:
            return []
        snapped = vertices[0]
        return self._adcp_vertices_to_segments(vertices), snapped

    def _adcp_rebuild_circle(self, circle_num):
        center = getattr(self, f"adcp_circle{circle_num}_center", None)
        start = getattr(self, f"adcp_circle{circle_num}_start", None)
        if center is None or start is None:
            setattr(self, f"adcp_circle{circle_num}_segments", [])
            return
        clockwise = circle_num == 1
        segments, snapped_start = self._adcp_build_circle_segments(center, start, clockwise=clockwise)
        setattr(self, f"adcp_circle{circle_num}_start", snapped_start)
        setattr(self, f"adcp_circle{circle_num}_segments", segments)

    def _adcp_rebuild_all_circles(self):
        self._adcp_rebuild_circle(1)
        self._adcp_rebuild_circle(2)

    def _adcp_on_diameter_changed(self, *_args):
        """Legacy name; use ``_apply_adcp_params_commit`` (deferred Enter / blur)."""
        self._apply_adcp_params_commit()

    def _apply_adcp_params_commit(self, _widget=None):
        """Apply committed ADCP parameter edits: rebuild circles and refresh export name."""
        self._update_adcp_export_name()
        if getattr(self, "adcp_circle1_start", None) is not None:
            self._adcp_rebuild_circle(1)
        if getattr(self, "adcp_circle2_start", None) is not None:
            self._adcp_rebuild_circle(2)
        if self._adcp_plan_has_any_geometry():
            self._plot_survey_plan(preserve_view_limits=True)
            if hasattr(self, "_draw_adcp_profile"):
                self._draw_adcp_profile()

    def _adcp_circle_complete(self, circle_num):
        segs = getattr(self, f"adcp_circle{circle_num}_segments", None) or []
        return len(segs) == ADCP_NUM_SEGMENTS

    def _adcp_plan_complete(self):
        return self._adcp_circle_complete(1) and self._adcp_circle_complete(2)

    def _adcp_plan_has_any_geometry(self):
        return bool(
            getattr(self, "adcp_circle1_center", None)
            or getattr(self, "adcp_circle2_center", None)
            or getattr(self, "adcp_circle1_segments", None)
            or getattr(self, "adcp_circle2_segments", None)
        )

    def _adcp_all_segment_points(self):
        points = []
        for n in (1, 2):
            for seg in getattr(self, f"adcp_circle{n}_segments", None) or []:
                points.extend(seg)
            center = getattr(self, f"adcp_circle{n}_center", None)
            if center is not None:
                points.append(center)
        return points

    def _adcp_pick_button_for_mode(self, mode):
        return {
            "c1_center": "adcp_pick_circle1_center_btn",
            "c1_start": "adcp_pick_circle1_start_btn",
            "c2_center": "adcp_pick_circle2_center_btn",
            "c2_start": "adcp_pick_circle2_start_btn",
        }.get(mode)

    def _adcp_clear_pick_button_styles(self):
        for attr in (
            "adcp_pick_circle1_center_btn",
            "adcp_pick_circle1_start_btn",
            "adcp_pick_circle2_center_btn",
            "adcp_pick_circle2_start_btn",
        ):
            btn = getattr(self, attr, None)
            if btn is not None:
                btn.setStyleSheet("")

    def _adcp_cancel_other_map_modes(self):
        self.pick_center_mode = False
        if hasattr(self, "pick_pitch_line_mode"):
            self.pick_pitch_line_mode = False
        if hasattr(self, "pick_roll_line_mode"):
            self.pick_roll_line_mode = False
        if hasattr(self, "line_planning_mode"):
            self.line_planning_mode = False
        if hasattr(self, "measurement_tool_mode"):
            self.measurement_tool_mode = False
        if hasattr(self, "pick_center_btn"):
            self.pick_center_btn.setStyleSheet("")
            self.pick_center_btn.setText("Pick Center from GeoTIFF")
        if hasattr(self, "performance_pick_center_btn"):
            self.performance_pick_center_btn.setStyleSheet("")
            self.performance_pick_center_btn.setText("Pick Center from GeoTIFF")

    def _adcp_set_pick_mode(self, mode):
        self._adcp_cancel_other_map_modes()
        self._adcp_clear_pick_button_styles()
        self.adcp_pick_mode = mode
        if mode:
            btn_name = self._adcp_pick_button_for_mode(mode)
            btn = getattr(self, btn_name, None) if btn_name else None
            if btn is not None:
                btn.setStyleSheet(ADCP_ACTIVE_BTN_STYLE)
            if hasattr(self, "canvas_widget"):
                self.canvas_widget.setCursor(Qt.CursorShape.CrossCursor)
        else:
            if hasattr(self, "canvas_widget"):
                self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)

    def _toggle_adcp_pick_circle1_center(self):
        if self.adcp_pick_mode == "c1_center":
            self._adcp_set_pick_mode(None)
            return
        self._adcp_set_pick_mode("c1_center")

    def _toggle_adcp_pick_circle1_start(self):
        if self.adcp_pick_mode == "c1_start":
            self._adcp_set_pick_mode(None)
            return
        if getattr(self, "adcp_circle1_center", None) is None:
            self._show_message("info", "ADCP Planning", "Pick Circle 1 center first.")
            return
        self._adcp_set_pick_mode("c1_start")

    def _toggle_adcp_pick_circle2_center(self):
        if self.adcp_pick_mode == "c2_center":
            self._adcp_set_pick_mode(None)
            return
        if not self._adcp_circle_complete(1) and getattr(self, "adcp_circle2_center", None) is None:
            self._show_message("info", "ADCP Planning", "Complete Circle 1 first (center and start point).")
            return
        self._adcp_set_pick_mode("c2_center")

    def _toggle_adcp_pick_circle2_start(self):
        if self.adcp_pick_mode == "c2_start":
            self._adcp_set_pick_mode(None)
            return
        if getattr(self, "adcp_circle2_center", None) is None:
            self._show_message("info", "ADCP Planning", "Pick Circle 2 center first.")
            return
        self._adcp_set_pick_mode("c2_start")

    def _update_adcp_button_states(self):
        has_geotiff = bool(getattr(self, "geotiff_dataset_original", None))
        c1_center = getattr(self, "adcp_circle1_center", None) is not None
        c2_center = getattr(self, "adcp_circle2_center", None) is not None
        c1_done = self._adcp_circle_complete(1)
        c2_done = self._adcp_circle_complete(2)

        for btn_name in (
            "adcp_pick_circle1_center_btn",
            "adcp_pick_circle1_start_btn",
            "adcp_pick_circle2_center_btn",
            "adcp_pick_circle2_start_btn",
            "adcp_zoom_btn",
            "adcp_clear_btn",
            "adcp_show_info_btn",
            "adcp_export_btn",
        ):
            btn = getattr(self, btn_name, None)
            if btn is not None:
                btn.setEnabled(has_geotiff)

        if hasattr(self, "adcp_pick_circle1_start_btn"):
            self.adcp_pick_circle1_start_btn.setEnabled(has_geotiff and c1_center)
        if hasattr(self, "adcp_pick_circle2_center_btn"):
            self.adcp_pick_circle2_center_btn.setEnabled(has_geotiff and (c1_done or c2_center))
        if hasattr(self, "adcp_pick_circle2_start_btn"):
            self.adcp_pick_circle2_start_btn.setEnabled(has_geotiff and c2_center)
        if hasattr(self, "adcp_zoom_btn"):
            self.adcp_zoom_btn.setEnabled(has_geotiff and self._adcp_plan_has_any_geometry())
        if hasattr(self, "adcp_clear_btn"):
            self.adcp_clear_btn.setEnabled(self._adcp_plan_has_any_geometry())
        if hasattr(self, "adcp_show_info_btn"):
            self.adcp_show_info_btn.setEnabled(self._adcp_plan_complete())
        if hasattr(self, "adcp_export_btn"):
            self.adcp_export_btn.setEnabled(self._adcp_plan_complete())
        if hasattr(self, "adcp_show_direction_checkbox"):
            has_segments = bool(
                getattr(self, "adcp_circle1_segments", None)
                or getattr(self, "adcp_circle2_segments", None)
            )
            self.adcp_show_direction_checkbox.setEnabled(has_geotiff and has_segments)

    def _on_adcp_show_direction_changed(self, *_args):
        if hasattr(self, "_update_travel_direction_arrows"):
            self._update_travel_direction_arrows()

    def _handle_adcp_plot_click(self, event):
        mode = getattr(self, "adcp_pick_mode", None)
        if not mode:
            return False
        if not hasattr(self, "param_notebook") or self.param_notebook.currentIndex() != 5:
            return False
        if event.inaxes != self.ax:
            return False
        clicked_lon = event.xdata
        clicked_lat = event.ydata
        if clicked_lon is None or clicked_lat is None:
            self._show_message(
                "warning",
                "Click Error",
                "Could not get valid coordinates from click. Try clicking within the displayed GeoTIFF area.",
            )
            return True

        if mode == "c1_center":
            self.adcp_circle1_center = (clicked_lat, clicked_lon)
            self.adcp_circle1_start = None
            self.adcp_circle1_segments = []
            self._adcp_set_pick_mode(None)
            self.set_adcp_activity_text(
                "Circle 1 center set. Pick Circle 1 start point on the circle (clockwise path).",
                append=False,
            )
        elif mode == "c1_start":
            snapped = self._adcp_snap_point_to_circle(self.adcp_circle1_center, clicked_lat, clicked_lon)
            if snapped is None:
                return True
            segments, snapped_start = self._adcp_build_circle_segments(
                self.adcp_circle1_center, snapped, clockwise=True
            )
            self.adcp_circle1_start = snapped_start
            self.adcp_circle1_segments = segments
            self._adcp_set_pick_mode(None)
            self.set_adcp_activity_text(
                "Circle 1 planned (clockwise). Pick Circle 2 center on the GeoTIFF.",
                append=False,
            )
        elif mode == "c2_center":
            self.adcp_circle2_center = (clicked_lat, clicked_lon)
            self.adcp_circle2_start = None
            self.adcp_circle2_segments = []
            self._adcp_set_pick_mode(None)
            self.set_adcp_activity_text(
                "Circle 2 center set. Pick Circle 2 start point on the circle (counter-clockwise path).",
                append=False,
            )
        elif mode == "c2_start":
            snapped = self._adcp_snap_point_to_circle(self.adcp_circle2_center, clicked_lat, clicked_lon)
            if snapped is None:
                return True
            segments, snapped_start = self._adcp_build_circle_segments(
                self.adcp_circle2_center, snapped, clockwise=False
            )
            self.adcp_circle2_start = snapped_start
            self.adcp_circle2_segments = segments
            self._adcp_set_pick_mode(None)
            self.set_adcp_activity_text("ADCP calibration plan complete (both circles).", append=False)
        else:
            return False

        self._update_adcp_button_states()
        self._plot_survey_plan(preserve_view_limits=True)
        if hasattr(self, "_draw_adcp_profile"):
            self._draw_adcp_profile()
        return True

    def _zoom_to_adcp_cal(self):
        points = self._adcp_all_segment_points()
        if not points:
            self._show_message("info", "Zoom to ADCP Cal", "No ADCP calibration plan to zoom to.")
            return
        lats = [p[0] for p in points]
        lons = [p[1] for p in points]
        min_lat, max_lat = min(lats), max(lats)
        min_lon, max_lon = min(lons), max(lons)
        buffer_lat = (max_lat - min_lat) * 0.08 if (max_lat - min_lat) != 0 else 0.01
        buffer_lon = (max_lon - min_lon) * 0.08 if (max_lon - min_lon) != 0 else 0.01
        self._apply_map_zoom_limits_and_reload_geotiff(
            (min_lon - buffer_lon, max_lon + buffer_lon),
            (min_lat - buffer_lat, max_lat + buffer_lat),
        )

    def _clear_adcp_plan(self):
        self.adcp_circle1_center = None
        self.adcp_circle1_start = None
        self.adcp_circle1_segments = []
        self.adcp_circle2_center = None
        self.adcp_circle2_start = None
        self.adcp_circle2_segments = []
        self._adcp_set_pick_mode(None)
        self._update_adcp_button_states()
        self._plot_survey_plan(preserve_view_limits=True)
        if hasattr(self, "_draw_adcp_profile"):
            self._draw_adcp_profile()
        self.set_adcp_activity_text("ADCP calibration plan cleared.", append=False)

    def _adcp_segment_distance_m(self, seg):
        geod = self._adcp_geod()
        if geod is None or len(seg) != 2:
            return 0.0
        _, _, d = geod.inv(seg[0][1], seg[0][0], seg[1][1], seg[1][0])
        return float(d)

    def _build_adcp_cal_info_text(self, speed_kts, turn_min):
        geod = self._adcp_geod()
        if geod is None:
            raise ValueError("Geospatial libraries not available.")
        if not self._adcp_plan_complete():
            raise ValueError("Both circles must be planned before building ADCP info.")

        speed_mps = speed_kts * 0.514444
        diameter_m = self._adcp_get_diameter_m()
        radius_m = diameter_m / 2.0

        def fmt_time(seconds_value):
            minutes_value = float(seconds_value) / 60.0
            hours_value = minutes_value / 60.0
            return f"{minutes_value:.2f} min ({hours_value:.2f} hr)"

        def leg_time_sec(d_m):
            if speed_mps <= 0 or d_m <= 0:
                return 0.0
            return d_m / speed_mps

        t = []
        t.append("ADCP CALIBRATION SUMMARY\n")
        t.append("=" * 26 + "\n")
        t.append(f"Circle diameter: {diameter_m:.1f} m (radius {radius_m:.1f} m)\n")
        t.append(f"Segments per circle: {ADCP_NUM_SEGMENTS}\n")
        t.append(f"Survey speed: {speed_kts:.2f} kts")
        if speed_mps <= 0:
            t.append(" (underway time at speed is not computed)\n")
        else:
            t.append(f" ({speed_mps:.3f} m/s)\n")
        t.append(f"Turn time (each turn): {turn_min:.2f} min\n\n")

        circle_rows = []
        for circle_num, direction_label in ((1, "clockwise"), (2, "counter-clockwise")):
            segs = getattr(self, f"adcp_circle{circle_num}_segments", []) or []
            total_dist = sum(self._adcp_segment_distance_m(s) for s in segs)
            total_time = sum(leg_time_sec(self._adcp_segment_distance_m(s)) for s in segs)
            circle_rows.append(
                {
                    "num": circle_num,
                    "direction": direction_label,
                    "segments": segs,
                    "total_dist": total_dist,
                    "total_time": total_time,
                }
            )

        t.append("CIRCLE PATTERNS\n")
        t.append("-" * 15 + "\n")
        t.append("Circle 1: 36 straight segments, clockwise, starting at the picked start point.\n")
        t.append("Circle 2: 36 straight segments, counter-clockwise, starting at the picked start point.\n\n")

        grand_dist = 0.0
        grand_time = 0.0
        for row in circle_rows:
            grand_dist += row["total_dist"]
            grand_time += row["total_time"]

        c1_start = self.adcp_circle1_start
        c2_start = self.adcp_circle2_start
        transit_dist = 0.0
        transit_time = 0.0
        if c1_start and c2_start:
            _, _, transit_dist = geod.inv(c1_start[1], c1_start[0], c2_start[1], c2_start[0])
            transit_dist = float(transit_dist)
            transit_time = leg_time_sec(transit_dist)

        t.append("TRANSIT (Circle 1 start → Circle 2 start)\n")
        t.append("-" * 38 + "\n")
        t.append(f"  Distance: {transit_dist:.1f} m")
        if speed_mps > 0:
            t.append(f"  |  {fmt_time(transit_time)}\n")
        else:
            t.append("\n")
        t.append("\n")

        num_turns = 3
        total_turn_time_sec = num_turns * turn_min * 60.0
        grand_dist += transit_dist
        grand_time += transit_time + total_turn_time_sec

        t.append("TURNS\n")
        t.append("-" * 5 + "\n")
        t.append(
            f"  1 turn entering Circle 1, 1 turn at Circle 1→2 transit, "
            f"1 turn entering Circle 2 = {num_turns} × {turn_min:.2f} min\n\n"
        )

        t.append("TOTALS\n")
        t.append("-" * 6 + "\n")
        t.append(f"  Distance underway: {grand_dist:.1f} m ({grand_dist/1000.0:.3f} km)\n")
        if speed_mps > 0:
            t.append(f"  Time at speed:     {fmt_time(grand_time - total_turn_time_sec)}\n")
        t.append(f"  Turn time:         {fmt_time(total_turn_time_sec)}\n")
        if speed_mps > 0:
            t.append(f"  Grand total:       {fmt_time(grand_time)}\n")
        t.append("\n")

        t.append("ORDERED WAYPOINTS (segment endpoints)\n")
        t.append("-" * 35 + "\n")
        labels = []
        for i in range(1, ADCP_NUM_SEGMENTS + 1):
            labels.append(f"C1S{i:02d}")
            labels.append(f"C1E{i:02d}")
        labels.append("Transit")
        for i in range(1, ADCP_NUM_SEGMENTS + 1):
            labels.append(f"C2S{i:02d}")
            labels.append(f"C2E{i:02d}")
        t.append(" → ".join(labels[:12]) + " → … → " + " → ".join(labels[-4:]) + "\n\n")

        def append_block(title, use_ddm):
            t.append(f"{title}\n")
            t.append("-" * len(title) + "\n")
            for circle_num in (1, 2):
                segs = getattr(self, f"adcp_circle{circle_num}_segments", []) or []
                center = getattr(self, f"adcp_circle{circle_num}_center", None)
                start = getattr(self, f"adcp_circle{circle_num}_start", None)
                if center:
                    if use_ddm:
                        t.append(
                            f"C{circle_num} Center: {decimal_degrees_to_ddm(center[0], is_latitude=True)}, "
                            f"{decimal_degrees_to_ddm(center[1], is_latitude=False)}\n"
                        )
                    else:
                        t.append(f"C{circle_num} Center: {center[0]:.6f}, {center[1]:.6f}\n")
                if start:
                    if use_ddm:
                        t.append(
                            f"C{circle_num} Start: {decimal_degrees_to_ddm(start[0], is_latitude=True)}, "
                            f"{decimal_degrees_to_ddm(start[1], is_latitude=False)}\n"
                        )
                    else:
                        t.append(f"C{circle_num} Start: {start[0]:.6f}, {start[1]:.6f}\n")
                for i, seg in enumerate(segs, start=1):
                    p0, p1 = seg[0], seg[1]
                    if use_ddm:
                        t.append(
                            f"C{circle_num} Seg {i:02d} start: "
                            f"{decimal_degrees_to_ddm(p0[0], is_latitude=True)}, "
                            f"{decimal_degrees_to_ddm(p0[1], is_latitude=False)}\n"
                        )
                        t.append(
                            f"C{circle_num} Seg {i:02d} end:   "
                            f"{decimal_degrees_to_ddm(p1[0], is_latitude=True)}, "
                            f"{decimal_degrees_to_ddm(p1[1], is_latitude=False)}\n"
                        )
                    else:
                        t.append(f"C{circle_num} Seg {i:02d} start: {p0[0]:.6f}, {p0[1]:.6f}\n")
                        t.append(f"C{circle_num} Seg {i:02d} end:   {p1[0]:.6f}, {p1[1]:.6f}\n")
                t.append("\n")

        append_block("Waypoints (DMM)", True)
        append_block("Waypoints (decimal degrees)", False)
        return "".join(t)

    def _show_adcp_cal_info(self):
        if not GEOSPATIAL_LIBS_AVAILABLE or pyproj is None:
            self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded.")
            return
        if not self._adcp_plan_complete():
            self._show_message("info", "No ADCP Plan", "Plan both calibration circles first.")
            return
        try:
            stats_text = self._build_adcp_cal_info_text(
                self._adcp_get_survey_speed_kts(),
                self._adcp_get_turn_time_min(),
            )
        except ValueError as e:
            self._show_message("warning", "ADCP Info", str(e))
            return
        show_statistics_dialog(self, "ADCP Calibration Info", stats_text)

    def set_adcp_activity_text(self, message, append=False):
        if not hasattr(self, "activity_log_text"):
            return
        prefixed_message = f"[ADCP] {message}"
        self.activity_log_text.setReadOnly(False)
        if append:
            current_text = self.activity_log_text.toPlainText().rstrip("\n")
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

    def _adcp_export_rows(self):
        rows = []
        line_num = 1
        for circle_num in (1, 2):
            segs = getattr(self, f"adcp_circle{circle_num}_segments", []) or []
            for i, seg in enumerate(segs, start=1):
                lname = f"ADCP_C{circle_num}_{i:02d}"
                rows.append((line_num, lname, f"C{circle_num}S{i:02d}", seg[0][0], seg[0][1]))
                rows.append((line_num, lname, f"C{circle_num}E{i:02d}", seg[1][0], seg[1][1]))
                line_num += 1
        return rows

    def _export_adcp_cal_files(self):
        if hasattr(self, "_commit_all_deferred_line_edits"):
            self._commit_all_deferred_line_edits()
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded. Cannot export.")
            return
        if not self._adcp_plan_complete():
            self._show_message("warning", "No Data", "Plan both ADCP calibration circles before exporting.")
            return

        export_shapefile = self._export_type_enabled("esri_shapefile")
        export_gpkg = self._export_type_enabled("gpkg")
        export_sis = self._export_type_enabled("sis_asciiplan")
        export_gpx = self._export_type_enabled("gpx")
        export_text_csv = self._export_type_enabled("text_csv")
        export_text_txt = self._export_type_enabled("text_txt")
        export_hypack = self._export_type_enabled("hypack_lnw")
        export_map_png = self._export_map_png_enabled()
        export_profiles_png = self._export_profiles_png_enabled()
        if (export_shapefile or export_gpkg) and mapping is None:
            self._show_message("warning", "Export Error", "Shapely is required for shapefile/GeoPackage export.")
            return

        export_name = ""
        if hasattr(self, "adcp_export_name_entry"):
            export_name = self.adcp_export_name_entry.text().strip()
        if not export_name:
            export_name = self._build_adcp_export_basename()
        bad = '<>:"/\\|?*'
        for c in bad:
            export_name = export_name.replace(c, "_")
        export_name = export_name.strip().strip(".")

        export_dir = QFileDialog.getExistingDirectory(
            self,
            "Select Export Directory",
            getattr(self, "last_adcp_import_dir", os.path.expanduser("~")),
        )
        if not export_dir:
            return
        if export_dir and os.path.isdir(export_dir):
            self.last_adcp_import_dir = export_dir
            if hasattr(self, "_save_last_adcp_import_dir"):
                self._save_last_adcp_import_dir()

        speed_kts = self._adcp_get_survey_speed_kts()
        turn_min = self._adcp_get_turn_time_min()
        diameter_m = self._adcp_get_diameter_m()
        geotiff_path = self.current_geotiff_path if hasattr(self, "current_geotiff_path") and self.current_geotiff_path else None

        try:
            rows = self._adcp_export_rows()
            csv_file_path = os.path.join(export_dir, f"{export_name}_DDD.csv")
            ddm_file_path = os.path.join(export_dir, f"{export_name}_DMM.csv")
            dms_file_path = os.path.join(export_dir, f"{export_name}_DMS.csv")
            ddm_txt_file_path = os.path.join(export_dir, f"{export_name}_DMM.txt")
            dms_txt_file_path = os.path.join(export_dir, f"{export_name}_DMS.txt")
            txt_file_path = os.path.join(export_dir, f"{export_name}_DDD.txt")
            if export_text_csv:
                export_utils.write_ddd_csv(csv_file_path, rows, newline="")
                export_utils.write_dmm_csv(ddm_file_path, rows)
                export_utils.write_dms_csv(dms_file_path, rows)
            if export_text_txt:
                export_utils.write_dmm_txt(ddm_txt_file_path, rows)
                export_utils.write_dms_txt(dms_txt_file_path, rows)
                export_utils.write_ddd_txt(txt_file_path, rows)

            shapefile_path = os.path.join(export_dir, f"{export_name}.shp")
            features = []
            line_num = 1
            for circle_num in (1, 2):
                segs = getattr(self, f"adcp_circle{circle_num}_segments", []) or []
                direction = "clockwise" if circle_num == 1 else "counter_clockwise"
                for i, seg in enumerate(segs, start=1):
                    shapely_line = LineString([(p[1], p[0]) for p in seg])
                    features.append(
                        {
                            "geometry": mapping(shapely_line),
                            "properties": {
                                "line_num": line_num,
                                "line_name": f"ADCP_C{circle_num}_{i:02d}",
                                "circle_num": circle_num,
                                "segment_num": i,
                                "direction": direction,
                            },
                        }
                    )
                    line_num += 1
            if export_shapefile or export_gpkg:
                schema = {
                    "geometry": "LineString",
                    "properties": {
                        "line_num": "int",
                        "line_name": "str",
                        "circle_num": "int",
                        "segment_num": "int",
                        "direction": "str",
                    },
                }
                if export_shapefile:
                    export_utils.remove_export_file(shapefile_path)
                    with fiona.open(
                        shapefile_path, "w", driver="ESRI Shapefile", crs="EPSG:4326", schema=schema
                    ) as collection:
                        collection.writerecords(features)
                self._write_gpkg_if_enabled(shapefile_path, schema, features, crs="EPSG:4326")

            geojson_features = []
            line_num = 1
            for circle_num in (1, 2):
                segs = getattr(self, f"adcp_circle{circle_num}_segments", []) or []
                direction = "clockwise" if circle_num == 1 else "counter_clockwise"
                center = getattr(self, f"adcp_circle{circle_num}_center", None)
                start = getattr(self, f"adcp_circle{circle_num}_start", None)
                for i, seg in enumerate(segs, start=1):
                    geojson_features.append(
                        {
                            "type": "Feature",
                            "geometry": {
                                "type": "LineString",
                                "coordinates": [[seg[0][1], seg[0][0]], [seg[1][1], seg[1][0]]],
                            },
                            "properties": {
                                "line_num": line_num,
                                "line_name": f"ADCP_C{circle_num}_{i:02d}",
                                "circle_num": circle_num,
                                "segment_num": i,
                                "direction": direction,
                                "survey_speed": speed_kts,
                                "circle_diameter_m": diameter_m,
                                "geotiff_path": geotiff_path,
                                "points": [
                                    {"point_num": 1, "lat": seg[0][0], "lon": seg[0][1]},
                                    {"point_num": 2, "lat": seg[1][0], "lon": seg[1][1]},
                                ],
                                "circle_center": (
                                    {"lat": center[0], "lon": center[1]} if center else None
                                ),
                                "circle_start": (
                                    {"lat": start[0], "lon": start[1]} if start else None
                                ),
                            },
                        }
                    )
                    line_num += 1
            geojson_collection = {
                "type": "FeatureCollection",
                "properties": {
                    "geotiff_path": geotiff_path,
                    "geotiff_nan_value": float(getattr(self, "geotiff_nan_value", -11000.0)),
                    "plan_type": "adcp",
                    "circle_diameter_m": diameter_m,
                    "survey_speed_kts": speed_kts,
                    "turn_time_min": turn_min,
                    "circle1_center": list(self.adcp_circle1_center) if self.adcp_circle1_center else None,
                    "circle1_start": list(self.adcp_circle1_start) if self.adcp_circle1_start else None,
                    "circle2_center": list(self.adcp_circle2_center) if self.adcp_circle2_center else None,
                    "circle2_start": list(self.adcp_circle2_start) if self.adcp_circle2_start else None,
                },
                "features": geojson_features,
            }
            geojson_file_path = os.path.join(export_dir, f"{export_name}.geojson")
            export_utils.remove_export_file(geojson_file_path)
            with open(geojson_file_path, "w", encoding="utf-8") as f:
                json.dump(geojson_collection, f, indent=2)

            lnw_lines = []
            ascii_lines = []
            gpx_lines = []
            for circle_num in (1, 2):
                segs = getattr(self, f"adcp_circle{circle_num}_segments", []) or []
                for i, seg in enumerate(segs, start=1):
                    name = f"ADCP_C{circle_num}_{i:02d}"
                    lnw_lines.append((name, seg))
                    ascii_lines.append((name, seg))
                    gpx_lines.append((name, seg))
            lnw_file_path = None
            if export_hypack and lnw_lines:
                all_pts = [p for _name, pts in lnw_lines for p in pts]
                zone, hem = export_utils.compute_utm_zone_from_points(all_pts)
                utm_suffix = f"_UTM{zone}{'N' if hem == 'North' else 'S'}"
                lnw_file_path = os.path.join(export_dir, f"{export_name}{utm_suffix}.lnw")
                if not export_utils.write_lnw(lnw_file_path, lnw_lines):
                    lnw_file_path = None

            sis_file_path = os.path.join(export_dir, f"{export_name}.asciiplan")
            if export_sis:
                export_utils.write_asciiplan(sis_file_path, ascii_lines)

            gpx_file_path = os.path.join(export_dir, f"{export_name}.gpx")
            gpx_written = False
            gpx_per_test_names = []
            if export_gpx:
                gpx_written = export_utils.write_gpx(
                    gpx_file_path, gpx_lines, creator="SAT Planner ADCP"
                )
                gpx_per_test_tests = [
                    (name, name, seg)
                    for name, seg in gpx_lines
                ]
                gpx_per_test_names = export_utils.write_gpx_per_test_files(
                    export_dir, export_name, gpx_per_test_tests, creator="SAT Planner ADCP"
                )

            stats_file_path = os.path.join(export_dir, f"{export_name}_info.txt")
            info_text = self._build_adcp_cal_info_text(speed_kts, turn_min)
            export_utils.remove_export_file(stats_file_path)
            with open(stats_file_path, "w", encoding="utf-8") as f:
                f.write("ADCP CALIBRATION SURVEY INFORMATION\n")
                f.write("=" * 50 + "\n\n")
                f.write(f"ADCP Calibration: {export_name}\n")
                f.write(f"Export Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                f.write(info_text)

            json_metadata_path = os.path.join(export_dir, f"{export_name}_adcp_params.json")
            meta = {
                "export_name": export_name,
                "plan_type": "adcp",
                "geotiff_path": geotiff_path,
                "geotiff_nan_value": float(getattr(self, "geotiff_nan_value", -11000.0)),
                "circle_diameter_m": diameter_m,
                "survey_speed_kts": speed_kts,
                "turn_time_min": turn_min,
                "circle1_center": list(self.adcp_circle1_center),
                "circle1_start": list(self.adcp_circle1_start),
                "circle2_center": list(self.adcp_circle2_center),
                "circle2_start": list(self.adcp_circle2_start),
                "visualization_shapefile_paths": list(getattr(self, "visualization_shapefile_paths", []) or []),
                "show_contours_var": bool(getattr(self, "show_contours_var", False)),
                "contour_interval_m": (
                    float(self.contour_interval_entry.text())
                    if hasattr(self, "contour_interval_entry") and self.contour_interval_entry.text()
                    else 200.0
                ),
            }
            export_utils.remove_export_file(json_metadata_path)
            with open(json_metadata_path, "w", encoding="utf-8") as f:
                json.dump(meta, f, indent=2)

            map_png_path = os.path.join(export_dir, f"{export_name}_map.png")
            if export_map_png:
                if hasattr(self, "_hide_map_hover_tooltip_for_export"):
                    self._hide_map_hover_tooltip_for_export()
                self._save_export_map_png(map_png_path, dpi=300, bbox_inches="tight", facecolor="white")

            profile_png_path = os.path.join(export_dir, f"{export_name}_profile.png")
            if export_profiles_png and hasattr(self, "profile_fig") and self.profile_fig is not None:
                self._save_export_profile_png(
                    profile_png_path, dpi=300, bbox_inches="tight", facecolor="white"
                )

            status_lines = []

            def _add_status(path):
                if not path:
                    return
                status_lines.append(
                    f"{'OK' if os.path.exists(path) else 'FAILED'}: {os.path.basename(path)}"
                )

            _add_status(geojson_file_path)
            if export_text_csv:
                _add_status(csv_file_path)
                _add_status(ddm_file_path)
                _add_status(dms_file_path)
            if export_text_txt:
                _add_status(txt_file_path)
                _add_status(ddm_txt_file_path)
                _add_status(dms_txt_file_path)
            if export_shapefile:
                _add_status(shapefile_path)
            _add_status(lnw_file_path)
            if export_sis:
                _add_status(sis_file_path)
            if gpx_written:
                _add_status(gpx_file_path)
            _add_status(stats_file_path)
            if export_map_png:
                _add_status(map_png_path)
            _add_status(json_metadata_path)
            if export_profiles_png and hasattr(self, "profile_fig") and self.profile_fig is not None:
                _add_status(profile_png_path)

            if status_lines:
                self.set_adcp_activity_text("ADCP export results:\n" + "\n".join(status_lines), append=True)

            self._show_message(
                "info",
                "Export Complete",
                f"ADCP calibration exported to:\n{export_dir}\n\nPrimary file:\n{os.path.basename(geojson_file_path)}",
            )
        except Exception as e:
            self._show_message("error", "Export Error", f"Failed to export ADCP calibration: {e}")
            import traceback

            traceback.print_exc()

    def _adcp_preview_circle_vertices(self, center, clockwise=True):
        """Vertices for dashed preview circle when center is set but start is not."""
        geod = self._adcp_geod()
        if geod is None or center is None:
            return []
        clat, clon = center
        radius_m = self._adcp_get_radius_m()
        step = 360.0 / float(ADCP_NUM_SEGMENTS)
        vertices = []
        for i in range(ADCP_NUM_SEGMENTS):
            az = i * step if clockwise else -i * step
            lon, lat, _ = geod.fwd(clon, clat, az % 360.0, radius_m)
            vertices.append((lat, lon))
        return vertices

    def _adcp_parse_csv_line_list(self, file_path):
        """Parse exported DDD CSV into list of two-point line segments."""
        try:
            with open(file_path, "r", encoding="utf-8") as csvfile:
                reader = csv.DictReader(csvfile)
                lines_data = {}
                for row in reader:
                    try:
                        line_num = int(row.get("Line Number", -1))
                        point_label = row.get("Point Label", "").strip()
                        lat, lon = lat_lon_decimal_from_survey_csv_row(row)
                    except (ValueError, TypeError):
                        continue
                    if line_num <= 0:
                        continue
                    if line_num not in lines_data:
                        lines_data[line_num] = {}
                    lines_data[line_num][point_label] = (lat, lon)
                imported_lines = []
                for line_num in sorted(lines_data.keys()):
                    points = lines_data[line_num]
                    start_keys = [k for k in points if k.endswith("S") or k.endswith("S01")]
                    end_keys = [k for k in points if k.endswith("E") or k.endswith("E01")]
                    start_key = start_keys[0] if start_keys else None
                    end_key = end_keys[0] if end_keys else None
                    if not start_key:
                        for k in sorted(points.keys()):
                            if "S" in k:
                                start_key = k
                                break
                    if not end_key:
                        for k in sorted(points.keys()):
                            if "E" in k:
                                end_key = k
                                break
                    if start_key and end_key and start_key in points and end_key in points:
                        imported_lines.append([points[start_key], points[end_key]])
                return imported_lines
        except Exception as e:
            self._show_message("error", "Import Error", f"Failed to parse CSV: {e}")
            return None

    def _try_parse_adcp_geojson(self, geojson_data):
        props = geojson_data.get("properties") or {}
        plan_type = props.get("plan_type")
        features = geojson_data.get("features") or []
        if plan_type != "adcp" and not features:
            return None

        c1_segs = [None] * ADCP_NUM_SEGMENTS
        c2_segs = [None] * ADCP_NUM_SEGMENTS
        found_plan = plan_type == "adcp"
        fallback_lines = []
        for feat in features:
            geom = feat.get("geometry") or {}
            if geom.get("type") != "LineString":
                continue
            coords = geom.get("coordinates") or []
            if len(coords) < 2:
                continue
            p0 = (coords[0][1], coords[0][0])
            p1 = (coords[-1][1], coords[-1][0])
            seg = [p0, p1]
            fprops = feat.get("properties") or {}
            circle_num = fprops.get("circle_num")
            segment_num = fprops.get("segment_num")
            if circle_num in (1, 2) and segment_num is not None:
                try:
                    idx = int(segment_num) - 1
                except (TypeError, ValueError):
                    idx = -1
                if 0 <= idx < ADCP_NUM_SEGMENTS:
                    found_plan = True
                    if circle_num == 1:
                        c1_segs[idx] = seg
                    else:
                        c2_segs[idx] = seg
            else:
                fallback_lines.append(seg)

        if found_plan and all(c1_segs) and all(c2_segs):
            meta = dict(props)
            return c1_segs, c2_segs, meta

        if len(fallback_lines) == ADCP_NUM_SEGMENTS * 2:
            c1 = fallback_lines[:ADCP_NUM_SEGMENTS]
            c2 = fallback_lines[ADCP_NUM_SEGMENTS:]
            return c1, c2, props
        return None

    def _adcp_apply_imported_segments(self, c1_segs, c2_segs, meta=None):
        meta = meta or {}
        self.adcp_circle1_segments = c1_segs
        self.adcp_circle2_segments = c2_segs
        if meta.get("circle1_center"):
            self.adcp_circle1_center = tuple(meta["circle1_center"])
        elif c1_segs:
            self.adcp_circle1_center = c1_segs[0][0]
        if meta.get("circle1_start"):
            self.adcp_circle1_start = tuple(meta["circle1_start"])
        elif c1_segs:
            self.adcp_circle1_start = c1_segs[0][0]
        if meta.get("circle2_center"):
            self.adcp_circle2_center = tuple(meta["circle2_center"])
        elif c2_segs:
            self.adcp_circle2_center = c2_segs[0][0]
        if meta.get("circle2_start"):
            self.adcp_circle2_start = tuple(meta["circle2_start"])
        elif c2_segs:
            self.adcp_circle2_start = c2_segs[0][0]

        if meta.get("circle_diameter_m") is not None:
            self._adcp_set_deferred_field("adcp_circle_diameter_entry", meta.get("circle_diameter_m"))
        elif meta.get("circle1_center") and meta.get("circle1_start"):
            self._adcp_rebuild_all_circles()
        if meta.get("survey_speed_kts") is not None:
            self._adcp_set_deferred_field("adcp_survey_speed_entry", meta.get("survey_speed_kts"))
        if meta.get("turn_time_min") is not None:
            self._adcp_set_deferred_field("adcp_turn_time_entry", meta.get("turn_time_min"))
        en = meta.get("export_name")
        if en:
            self._adcp_set_deferred_field("adcp_export_name_entry", en)
        else:
            self._update_adcp_export_name()

    def _adcp_import_post_import(self, file_path):
        if not self._adcp_plan_complete():
            self._show_message("warning", "Import Warning", "Need 36 segments per circle for a complete ADCP plan.")
            return
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        dir_name = os.path.dirname(file_path)
        params = None
        candidates = [base_name]
        for suffix in ("_DDD", "_DMM", "_DMS", "_DD", "_DM"):
            if base_name.endswith(suffix):
                candidates.append(base_name[: -len(suffix)])
        for candidate in candidates:
            meta_path = os.path.join(dir_name, f"{candidate}_adcp_params.json")
            if not os.path.exists(meta_path):
                continue
            try:
                with open(meta_path, "r", encoding="utf-8") as f:
                    params = json.load(f)
                break
            except Exception as e:
                print(f"Warning: could not load ADCP metadata: {e}")
        if params:
            if params.get("circle_diameter_m") is not None:
                self._adcp_set_deferred_field("adcp_circle_diameter_entry", params.get("circle_diameter_m"))
            if params.get("circle1_center"):
                self.adcp_circle1_center = tuple(params["circle1_center"])
            if params.get("circle1_start"):
                self.adcp_circle1_start = tuple(params["circle1_start"])
            if params.get("circle2_center"):
                self.adcp_circle2_center = tuple(params["circle2_center"])
            if params.get("circle2_start"):
                self.adcp_circle2_start = tuple(params["circle2_start"])
            self._adcp_rebuild_all_circles()
            if params.get("survey_speed_kts") is not None:
                self._adcp_set_deferred_field("adcp_survey_speed_entry", params.get("survey_speed_kts"))
            if params.get("turn_time_min") is not None:
                self._adcp_set_deferred_field("adcp_turn_time_entry", params.get("turn_time_min"))
            if params.get("export_name"):
                self._adcp_set_deferred_field("adcp_export_name_entry", params.get("export_name"))
            if params.get("geotiff_nan_value") is not None and hasattr(self, "_set_geotiff_nan_cutoff"):
                self._set_geotiff_nan_cutoff(params.get("geotiff_nan_value"), update_entry=True)
            gtp = params.get("geotiff_path")
            if gtp and hasattr(self, "_load_geotiff_from_path") and os.path.exists(gtp):
                self._load_geotiff_from_path(gtp)

        self._update_adcp_button_states()
        self._plot_survey_plan(preserve_view_limits=True)
        if hasattr(self, "_draw_adcp_profile"):
            self._draw_adcp_profile()
        if hasattr(self, "_zoom_to_adcp_cal"):
            self._zoom_to_adcp_cal()
        msg = f"Imported ADCP calibration from {os.path.basename(file_path)}."
        self.set_adcp_activity_text(msg, append=False)
        if getattr(self, "adcp_download_gmrt_checkbox", None) and self.adcp_download_gmrt_checkbox.isChecked():
            self._download_and_load_gmrt_after_adcp_import()

    def _download_and_load_gmrt_after_adcp_import(self):
        points = self._adcp_all_segment_points()
        if not points:
            self._show_message("warning", "GMRT Download", "No ADCP plan points to compute extent.")
            return
        lats = [p[0] for p in points]
        lons = [p[1] for p in points]
        mid_lat = (min(lats) + max(lats)) / 2.0
        mid_lon = (min(lons) + max(lons)) / 2.0
        buffer_deg = 0.5
        if hasattr(self, "adcp_gmrt_buffer_spin"):
            try:
                buffer_deg = float(self.adcp_gmrt_buffer_spin.value())
            except (ValueError, TypeError):
                pass
        split_topo_depths = True
        if hasattr(self, "adcp_split_topo_depths_checkbox"):
            split_topo_depths = bool(self.adcp_split_topo_depths_checkbox.isChecked())
        self._download_gmrt_and_load(
            mid_lon - buffer_deg,
            mid_lon + buffer_deg,
            mid_lat - buffer_deg,
            mid_lat + buffer_deg,
            resolution=100,
            layer="topo",
            default_filename_prefix="GMRT_Bathy",
            log_func=lambda msg, append=True: self.set_adcp_activity_text(msg, append=append),
            default_directory=getattr(self, "last_adcp_import_dir", None),
            split_topo_depths=split_topo_depths,
            gmrt_button=getattr(self, "adcp_import_btn", None),
        )

    def _import_adcp_cal(self):
        if hasattr(self, "_gmrt_is_downloading") and self._gmrt_is_downloading():
            self._gmrt_cancel_active_download()
            return
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded.")
            return

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select ADCP Calibration File to Import",
            getattr(self, "last_adcp_import_dir", os.path.expanduser("~")),
            "Known Survey Files (*_DMS.txt *_DMM.txt *_DDD.txt *_DDD.csv *_DMM.csv *_DMS.csv *.csv *.geojson *.json *.gpx *.lnw *.shp *.gpkg);;"
            "Hypack LNW files (*.lnw);;Shapefile (*.shp);;GeoPackage (*.gpkg);;GeoJSON (*.geojson);;JSON (*.json);;"
            "Decimal Degree CSV (*_DDD.csv);;CSV (*.csv);;GPX (*.gpx)",
        )
        if not file_path:
            return
        import_dir = os.path.dirname(file_path)
        if import_dir and os.path.isdir(import_dir):
            self.last_adcp_import_dir = import_dir
            if hasattr(self, "_save_last_adcp_import_dir"):
                self._save_last_adcp_import_dir()

        try:
            file_ext = os.path.splitext(file_path)[1].lower()
            file_basename = os.path.basename(file_path)

            if file_ext in (".geojson", ".json"):
                with open(file_path, "r", encoding="utf-8") as f:
                    geojson_data = json.load(f)
                parsed = self._try_parse_adcp_geojson(geojson_data)
                if parsed is not None:
                    c1, c2, meta = parsed
                    self._adcp_apply_imported_segments(c1, c2, meta)
                    nan_cutoff = (geojson_data.get("properties") or {}).get("geotiff_nan_value")
                    if nan_cutoff is not None and hasattr(self, "_set_geotiff_nan_cutoff"):
                        self._set_geotiff_nan_cutoff(nan_cutoff, update_entry=True)
                    gtp = (geojson_data.get("properties") or {}).get("geotiff_path")
                    if gtp and hasattr(self, "_load_geotiff_from_path") and os.path.exists(gtp):
                        self._load_geotiff_from_path(gtp)
                    self._adcp_import_post_import(file_path)
                    return
                self._show_message(
                    "warning",
                    "Import Warning",
                    "GeoJSON is not a recognized ADCP calibration plan (need 72 segments or plan_type adcp).",
                )
                return

            if file_ext in (".shp", ".gpkg"):
                imported_lines = self._parse_vector_file_as_line_list(file_path)
                if imported_lines is None:
                    return
                if len(imported_lines) == ADCP_NUM_SEGMENTS * 2:
                    c1 = imported_lines[:ADCP_NUM_SEGMENTS]
                    c2 = imported_lines[ADCP_NUM_SEGMENTS:]
                    self._adcp_apply_imported_segments(c1, c2, {})
                    self._adcp_import_post_import(file_path)
                    return
                self._show_message(
                    "warning",
                    "Import Warning",
                    f"Expected {ADCP_NUM_SEGMENTS * 2} line segments for ADCP import; found {len(imported_lines)}.",
                )
                return

            if file_ext == ".lnw":
                if UTMZoneDialog is None:
                    self._show_message("error", "Import Error", "UTM zone dialog not available.")
                    return
                utm_dialog = UTMZoneDialog.for_file(self, file_path)
                if utm_dialog.exec() != QDialog.DialogCode.Accepted:
                    return
                utm_zone, hemisphere = utm_dialog.get_utm_info()
                imported_lines = self._parse_lnw_file(file_path, utm_zone, hemisphere)
                if imported_lines is None:
                    return
                if len(imported_lines) == ADCP_NUM_SEGMENTS * 2:
                    c1 = imported_lines[:ADCP_NUM_SEGMENTS]
                    c2 = imported_lines[ADCP_NUM_SEGMENTS:]
                    self._adcp_apply_imported_segments(c1, c2, {})
                    self._adcp_import_post_import(file_path)
                    return
                self._show_message(
                    "warning",
                    "Import Warning",
                    f"Expected {ADCP_NUM_SEGMENTS * 2} LNW segments; found {len(imported_lines)}.",
                )
                return

            if file_ext == ".csv" or file_basename.lower().endswith(("_ddd.txt", "_dmm.txt", "_dms.txt")):
                if file_ext == ".csv":
                    imported_lines = self._adcp_parse_csv_line_list(file_path)
                elif file_basename.lower().endswith("_ddd.txt"):
                    imported_lines = self._parse_ddd_txt_file(file_path)
                elif file_basename.lower().endswith("_dmm.txt"):
                    imported_lines = self._parse_dmm_txt_file(file_path)
                else:
                    imported_lines = self._parse_dms_txt_file(file_path)
                if imported_lines is None:
                    return
                if len(imported_lines) == ADCP_NUM_SEGMENTS * 2:
                    c1 = imported_lines[:ADCP_NUM_SEGMENTS]
                    c2 = imported_lines[ADCP_NUM_SEGMENTS:]
                    self._adcp_apply_imported_segments(c1, c2, {})
                    self._adcp_import_post_import(file_path)
                    return
                self._show_message(
                    "warning",
                    "Import Warning",
                    f"Expected {ADCP_NUM_SEGMENTS * 2} survey lines; found {len(imported_lines)}.",
                )
                return

            if file_ext == ".gpx":
                imported_lines = self._gpx_to_unassigned_segments(file_path)
                if imported_lines is None:
                    return
                if len(imported_lines) == ADCP_NUM_SEGMENTS * 2:
                    c1 = imported_lines[:ADCP_NUM_SEGMENTS]
                    c2 = imported_lines[ADCP_NUM_SEGMENTS:]
                    self._adcp_apply_imported_segments(c1, c2, {})
                    self._adcp_import_post_import(file_path)
                    return
                self._show_message(
                    "warning",
                    "Import Warning",
                    f"Expected {ADCP_NUM_SEGMENTS * 2} GPX segments; found {len(imported_lines)}.",
                )
                return

            meta_path = file_path
            if file_path.endswith("_adcp_params.json"):
                with open(file_path, "r", encoding="utf-8") as f:
                    params = json.load(f)
                if params.get("circle1_center") and params.get("circle1_start"):
                    self.adcp_circle1_center = tuple(params["circle1_center"])
                    self.adcp_circle1_start = tuple(params["circle1_start"])
                if params.get("circle2_center") and params.get("circle2_start"):
                    self.adcp_circle2_center = tuple(params["circle2_center"])
                    self.adcp_circle2_start = tuple(params["circle2_start"])
                if params.get("circle_diameter_m") is not None:
                    self._adcp_set_deferred_field("adcp_circle_diameter_entry", params.get("circle_diameter_m"))
                self._adcp_rebuild_all_circles()
                self._adcp_import_post_import(meta_path)
                return

            self._show_message("error", "Import Error", f"Unsupported or unrecognized ADCP file format: {file_ext}")

        except Exception as e:
            self._show_message("error", "Import Error", f"Failed to import ADCP calibration: {e}")
            import traceback

            traceback.print_exc()

    def _draw_adcp_profile(self):
        if not hasattr(self, "profile_ax") or not hasattr(self, "profile_canvas"):
            return
        for ax in self.profile_fig.get_axes():
            if ax != self.profile_ax:
                ax.remove()
        self.profile_ax.clear()
        self.profile_ax.set_xlabel("Distance (m)", fontsize=8)
        self.profile_ax.set_ylabel("Elevation (m)", fontsize=8)
        self.profile_ax.tick_params(axis="both", which="major", labelsize=7)

        if not self._adcp_plan_complete():
            self.profile_ax.set_title("ADCP — Elevation Profile", fontsize=8)
            self.profile_fig.tight_layout(pad=1.0)
            self.profile_canvas.draw_idle()
            return

        geod = self._adcp_geod()
        if geod is None:
            return

        segments = []
        for circle_num, color, circle_label in (
            (1, "dodgerblue", "Circle 1"),
            (2, "darkorange", "Circle 2"),
        ):
            first_in_circle = True
            for seg in getattr(self, f"adcp_circle{circle_num}_segments", []) or []:
                label = circle_label if first_in_circle else "_nolegend_"
                segments.append((seg[0], seg[1], color, label))
                first_in_circle = False

        self.profile_ax.set_title("ADCP Calibration — Elevation Profile", fontsize=8)
        cumulative = 0.0
        elev_arrays = []
        label_used = set()
        for seg_start, seg_end, color, label in segments:
            lat1, lon1 = seg_start
            lat2, lon2 = seg_end
            lats = np.linspace(lat1, lat2, 25)
            lons = np.linspace(lon1, lon2, 25)
            elevations, slopes, _ = self._get_profile_data_from_geotiff(lats, lons)
            if elevations is None:
                continue
            dists = [0.0]
            for i in range(1, len(lats)):
                _, _, dist = geod.inv(lons[i - 1], lats[i - 1], lons[i], lats[i])
                dists.append(dists[-1] + dist)
            dists_plot = np.array(dists, dtype=float) + cumulative
            plot_label = label if label not in label_used and label != "_nolegend_" else "_nolegend_"
            if plot_label != "_nolegend_":
                label_used.add(plot_label)
            self.profile_ax.plot(dists_plot, elevations, color=color, lw=1.2, label=plot_label)
            elev_arrays.append(elevations)
            cumulative = float(dists_plot[-1])

        if elev_arrays:
            ymin = min(float(np.nanmin(e)) for e in elev_arrays)
            ymax = max(float(np.nanmax(e)) for e in elev_arrays)
            if np.isfinite(ymin) and np.isfinite(ymax):
                pad = max(1.0, (ymax - ymin) * 0.05)
                self.profile_ax.set_ylim(ymin - pad, ymax + pad)
        handles, labels = self.profile_ax.get_legend_handles_labels()
        if handles:
            self.profile_ax.legend(handles, labels, fontsize=9, loc="upper right")
        self.profile_fig.tight_layout(pad=1.0)
        self.profile_canvas.draw_idle()

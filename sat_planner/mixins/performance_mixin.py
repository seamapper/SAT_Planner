"""
Performance tab behavior: swath ping-time calculations and pick-center field updates.
"""

import csv
import json
import os
import time

import numpy as np
from PyQt6.QtCore import QTimer
from PyQt6.QtGui import QTextCursor
from PyQt6.QtWidgets import QFileDialog, QDialog

from sat_planner.constants import GEOSPATIAL_LIBS_AVAILABLE, pyproj
from sat_planner.performance_import_dialog import PerformanceLineAssignmentDialog
from sat_planner.utils_geo import decimal_degrees_to_ddm
from sat_planner.utils_ui import show_statistics_dialog

try:
    from sat_planner.mixins.survey_parsers_mixin import UTMZoneDialog
except ImportError:
    UTMZoneDialog = None


class PerformanceMixin:
    """Mixin for Performance tab calculation and field helpers."""

    def _performance_get_text(self, attr_name):
        widget = getattr(self, attr_name, None)
        if widget is None:
            return ""
        text_getter = getattr(widget, "text", None)
        if callable(text_getter):
            try:
                return str(text_getter()).strip()
            except Exception:
                return ""
        return ""

    def _performance_set_value(self, attr_name, text):
        widget = getattr(self, attr_name, None)
        if widget is None:
            return
        if hasattr(widget, "setText"):
            widget.setText(str(text))

    def _performance_clear_value(self, attr_name):
        widget = getattr(self, attr_name, None)
        if widget is None:
            return
        if hasattr(widget, "clear"):
            widget.clear()
        elif hasattr(widget, "setText"):
            widget.setText("-")

    def _debug_performance_log(self, hypothesis_id, location, message, data):
        """Append an NDJSON debug event for this debug session."""
        payload = {
            "sessionId": "f0b845",
            "runId": "pre-fix",
            "hypothesisId": hypothesis_id,
            "location": location,
            "message": message,
            "data": data,
            "timestamp": int(time.time() * 1000),
        }
        try:
            with open("debug-f0b845.log", "a", encoding="utf-8") as f:
                f.write(json.dumps(payload) + "\n")
        except Exception:
            pass

    def _performance_ping_time_sec(self):
        """Return two-way ping time in seconds, or None if inputs are invalid."""
        if not hasattr(self, "performance_test_depth_entry"):
            return None
        try:
            depth = float(self.performance_test_depth_entry.text().strip())
            swath_angle_deg = float(self.performance_swath_angle_entry.text().strip())
            sound_velocity = float(self.performance_sound_velocity_entry.text().strip())
        except Exception:
            return None

        if depth <= 0 or sound_velocity <= 0:
            return None
        if swath_angle_deg < 0 or swath_angle_deg >= 89.9:
            return None

        cos_theta = np.cos(np.radians(swath_angle_deg))
        if abs(cos_theta) < 1e-9:
            return None

        ping_time_sec = (2.0 * depth) / (sound_velocity * cos_theta)
        if ping_time_sec > 0 and np.isfinite(ping_time_sec):
            return float(ping_time_sec)
        return None

    def _update_performance_ping_time(self):
        """Compute two-way ping time using depth, swath angle, and sound velocity."""
        if not hasattr(self, "performance_ping_time_entry"):
            return
        ping_time_sec = self._performance_ping_time_sec()
        if ping_time_sec is not None:
            self._performance_set_value("performance_ping_time_entry", f"{ping_time_sec:.3f}")
        else:
            self._performance_clear_value("performance_ping_time_entry")
        self._update_performance_total_test_time()

    def _update_performance_total_test_time(self):
        """Total time = (# pings × ping time) + BIST (minutes as seconds)."""
        if not hasattr(self, "performance_total_test_time_min_entry"):
            return
        ping_time_sec = self._performance_ping_time_sec()
        if ping_time_sec is None:
            # region agent log
            self._debug_performance_log(
                "H1",
                "performance_mixin.py:_update_performance_total_test_time",
                "Ping time unavailable, clearing total test time",
                {"ping_time_sec": None},
            )
            # endregion
            self._performance_clear_value("performance_total_test_time_min_entry")
            self._performance_clear_value("performance_total_test_time_sec_entry")
            if hasattr(self, "performance_line_length_m_entry"):
                self._performance_clear_value("performance_line_length_m_entry")
                self._performance_clear_value("performance_line_length_km_entry")
            return
        try:
            num_pings = float(self.performance_num_pings_entry.text().strip())
        except Exception:
            # region agent log
            self._debug_performance_log(
                "H2",
                "performance_mixin.py:_update_performance_total_test_time",
                "Invalid number of pings",
                {"raw_num_pings": self.performance_num_pings_entry.text().strip()},
            )
            # endregion
            self._performance_clear_value("performance_total_test_time_min_entry")
            self._performance_clear_value("performance_total_test_time_sec_entry")
            if hasattr(self, "performance_line_length_m_entry"):
                self._performance_clear_value("performance_line_length_m_entry")
                self._performance_clear_value("performance_line_length_km_entry")
            return
        if num_pings < 0 or not np.isfinite(num_pings):
            # region agent log
            self._debug_performance_log(
                "H2",
                "performance_mixin.py:_update_performance_total_test_time",
                "Out-of-range number of pings",
                {"num_pings": num_pings},
            )
            # endregion
            self._performance_clear_value("performance_total_test_time_min_entry")
            self._performance_clear_value("performance_total_test_time_sec_entry")
            if hasattr(self, "performance_line_length_m_entry"):
                self._performance_clear_value("performance_line_length_m_entry")
                self._performance_clear_value("performance_line_length_km_entry")
            return

        bist_raw = self.performance_bist_time_entry.text().strip()
        if not bist_raw:
            bist_min = 0.0
        else:
            try:
                bist_min = float(bist_raw)
            except Exception:
                # region agent log
                self._debug_performance_log(
                    "H2",
                    "performance_mixin.py:_update_performance_total_test_time",
                    "Invalid BIST time",
                    {"raw_bist_min": bist_raw},
                )
                # endregion
                self._performance_clear_value("performance_total_test_time_min_entry")
                self._performance_clear_value("performance_total_test_time_sec_entry")
                if hasattr(self, "performance_line_length_m_entry"):
                    self._performance_clear_value("performance_line_length_m_entry")
                    self._performance_clear_value("performance_line_length_km_entry")
                return
        if bist_min < 0 or not np.isfinite(bist_min):
            # region agent log
            self._debug_performance_log(
                "H2",
                "performance_mixin.py:_update_performance_total_test_time",
                "Out-of-range BIST time",
                {"bist_min": bist_min},
            )
            # endregion
            self._performance_clear_value("performance_total_test_time_min_entry")
            self._performance_clear_value("performance_total_test_time_sec_entry")
            if hasattr(self, "performance_line_length_m_entry"):
                self._performance_clear_value("performance_line_length_m_entry")
                self._performance_clear_value("performance_line_length_km_entry")
            return

        collect_sec = num_pings * ping_time_sec
        total_sec = collect_sec + bist_min * 60.0
        if not np.isfinite(total_sec):
            # region agent log
            self._debug_performance_log(
                "H2",
                "performance_mixin.py:_update_performance_total_test_time",
                "Computed non-finite total seconds",
                {"collect_sec": collect_sec, "bist_min": bist_min, "total_sec": total_sec},
            )
            # endregion
            self._performance_clear_value("performance_total_test_time_min_entry")
            self._performance_clear_value("performance_total_test_time_sec_entry")
            if hasattr(self, "performance_line_length_m_entry"):
                self._performance_clear_value("performance_line_length_m_entry")
                self._performance_clear_value("performance_line_length_km_entry")
            return

        total_min = total_sec / 60.0
        self._performance_set_value("performance_total_test_time_min_entry", f"{total_min:.2f}")
        self._performance_set_value("performance_total_test_time_sec_entry", f"{total_sec:.1f}")
        # region agent log
        self._debug_performance_log(
            "H3",
            "performance_mixin.py:_update_performance_total_test_time",
            "Computed total test time",
            {"ping_time_sec": ping_time_sec, "num_pings": num_pings, "bist_min": bist_min, "total_sec": total_sec},
        )
        # endregion
        self._update_performance_line_length()

    def _update_performance_line_length(self):
        """Compute line length from test speed and ping-collection time (excludes BIST)."""
        if not hasattr(self, "performance_line_length_m_entry"):
            return

        ping_time_sec = self._performance_ping_time_sec()
        if ping_time_sec is None:
            self._performance_clear_value("performance_line_length_m_entry")
            self._performance_clear_value("performance_line_length_km_entry")
            return
        try:
            num_pings = float(self.performance_num_pings_entry.text().strip())
        except Exception:
            self._performance_clear_value("performance_line_length_m_entry")
            self._performance_clear_value("performance_line_length_km_entry")
            return
        if num_pings < 0 or not np.isfinite(num_pings):
            self._performance_clear_value("performance_line_length_m_entry")
            self._performance_clear_value("performance_line_length_km_entry")
            return
        total_sec_raw = f"{num_pings * ping_time_sec}"
        try:
            total_sec = float(total_sec_raw)
        except Exception:
            # region agent log
            self._debug_performance_log(
                "H4",
                "performance_mixin.py:_update_performance_line_length",
                "Total seconds unavailable for line length",
                {"raw_total_sec": total_sec_raw},
            )
            # endregion
            self._performance_clear_value("performance_line_length_m_entry")
            self._performance_clear_value("performance_line_length_km_entry")
            return
        try:
            speed_kts = float(self.performance_test_speed_entry.text().strip())
        except Exception:
            # region agent log
            self._debug_performance_log(
                "H4",
                "performance_mixin.py:_update_performance_line_length",
                "Invalid test speed",
                {"raw_speed_kts": self.performance_test_speed_entry.text().strip()},
            )
            # endregion
            self._performance_clear_value("performance_line_length_m_entry")
            self._performance_clear_value("performance_line_length_km_entry")
            return
        if total_sec < 0 or speed_kts < 0 or not np.isfinite(total_sec) or not np.isfinite(speed_kts):
            # region agent log
            self._debug_performance_log(
                "H4",
                "performance_mixin.py:_update_performance_line_length",
                "Out-of-range values for line length",
                {"total_sec": total_sec, "speed_kts": speed_kts},
            )
            # endregion
            self._performance_clear_value("performance_line_length_m_entry")
            self._performance_clear_value("performance_line_length_km_entry")
            return

        speed_mps = speed_kts * 0.514444
        line_length_m = speed_mps * total_sec
        line_length_km = line_length_m / 1000.0
        self._performance_set_value("performance_line_length_m_entry", f"{line_length_m:.1f}")
        self._performance_set_value("performance_line_length_km_entry", f"{line_length_km:.3f}")
        # region agent log
        self._debug_performance_log(
            "H5",
            "performance_mixin.py:_update_performance_line_length",
            "Computed line length from speed and ping collection time",
            {
                "speed_kts": speed_kts,
                "speed_mps": speed_mps,
                "total_sec": total_sec,
                "line_length_m": line_length_m,
                "line_length_km": line_length_km,
            },
        )
        # endregion
        self._schedule_autoplot_performance_test_lines(400)

    def _update_performance_center_from_pick(self, clicked_lat, clicked_lon):
        """Update Performance tab center fields from map pick-center click."""
        if hasattr(self, "performance_central_lat_entry"):
            self.performance_central_lat_entry.setText(f"{clicked_lat:.6f}")
        if hasattr(self, "performance_central_lon_entry"):
            self.performance_central_lon_entry.setText(f"{clicked_lon:.6f}")

    def _schedule_autoplot_performance_test_lines(self, delay_ms=400):
        """Debounce auto-plot: replot after the user pauses editing or after pick (delay_ms=0 or small)."""
        if not hasattr(self, "_performance_autoplot_timer"):
            self._performance_autoplot_timer = QTimer(self)
            self._performance_autoplot_timer.setSingleShot(True)
            self._performance_autoplot_timer.timeout.connect(self._run_autoplot_performance_test_lines)
        self._performance_autoplot_timer.stop()
        self._performance_autoplot_timer.start(int(delay_ms))

    def _run_autoplot_performance_test_lines(self):
        self._plot_performance_test_lines(quiet=True)

    def _update_performance_depth_from_pick(self, z_value):
        """Update Performance tab depth field from picked GeoTIFF z-value."""
        if not hasattr(self, "performance_test_depth_entry"):
            return
        try:
            if z_value is not None and not np.isnan(z_value):
                self.performance_test_depth_entry.setText(f"{abs(float(z_value)):.1f}")
            else:
                self.performance_test_depth_entry.clear()
        except Exception:
            self.performance_test_depth_entry.clear()

    def _plot_performance_test_lines(self, quiet=False):
        """Plot 4 performance test lines from center/length/swell and draw profile for line 1.

        If quiet is True, skip message boxes (used for debounced auto-plot after parameter edits).
        """
        if not GEOSPATIAL_LIBS_AVAILABLE or pyproj is None:
            if not quiet:
                self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded. Cannot plot test lines.")
            return
        if getattr(self, "geotiff_data_array", None) is None:
            if not quiet:
                self._show_message("warning", "No GeoTIFF", "Load a GeoTIFF first to plot performance test lines.")
            return

        if hasattr(self, "performance_swell_direction_entry"):
            if not self.performance_swell_direction_entry.text().strip():
                self.performance_swell_direction_entry.blockSignals(True)
                self.performance_swell_direction_entry.setText("0")
                self.performance_swell_direction_entry.blockSignals(False)

        try:
            center_lat = float(self.performance_central_lat_entry.text().strip())
            center_lon = float(self.performance_central_lon_entry.text().strip())
            line_length_m = float(self._performance_get_text("performance_line_length_m_entry"))
            swell_raw = self.performance_swell_direction_entry.text().strip()
            swell_direction = float(swell_raw) if swell_raw else 0.0
            bist_min = float(self.performance_bist_time_entry.text().strip()) if self.performance_bist_time_entry.text().strip() else 0.0
            speed_kts = float(self.performance_test_speed_entry.text().strip())
        except Exception:
            if not quiet:
                self._show_message(
                    "warning",
                    "Invalid Inputs",
                    "Provide valid Central Latitude, Central Longitude, Swell Direction, Test Speed, and calculated Line Length (m).",
                )
            return

        if not np.isfinite(line_length_m) or line_length_m <= 0:
            if not quiet:
                self._show_message("warning", "Invalid Line Length", "Line Length (m) must be greater than zero.")
            return
        if not np.isfinite(center_lat) or not np.isfinite(center_lon):
            if not quiet:
                self._show_message("warning", "Invalid Center", "Central Latitude/Longitude values are not valid.")
            return
        if bist_min < 0 or not np.isfinite(bist_min):
            if not quiet:
                self._show_message("warning", "Invalid BIST Time", "BIST Time (min) must be zero or greater.")
            return
        if speed_kts < 0 or not np.isfinite(speed_kts):
            if not quiet:
                self._show_message("warning", "Invalid Test Speed", "Test Speed (kts) must be zero or greater.")
            return

        geod = pyproj.Geod(ellps="WGS84")
        half_length = line_length_m / 2.0
        bist_length_m = speed_kts * 0.514444 * bist_min * 60.0
        headings = [((swell_direction + offset) % 360.0) for offset in (0.0, 45.0, 90.0, 135.0)]

        lines = []
        bist_segments = []
        for heading in headings:
            lon1, lat1, _ = geod.fwd(center_lon, center_lat, heading + 180.0, half_length)
            lon2, lat2, _ = geod.fwd(center_lon, center_lat, heading, half_length)
            lines.append([(lat1, lon1), (lat2, lon2)])
            if bist_length_m > 0:
                # Extend only from the line end (P?E) in the line heading direction.
                bist_end_lon, bist_end_lat, _ = geod.fwd(lon2, lat2, heading, bist_length_m)
                bist_segments.append([(lat2, lon2), (bist_end_lat, bist_end_lon)])

        self.performance_test_lines_data = lines
        self.performance_bist_segments_data = bist_segments
        self.performance_profile_line = lines[0]
        self.performance_central_point_coords = (center_lat, center_lon)
        self.central_point_coords = (center_lat, center_lon)

        self._plot_survey_plan(preserve_view_limits=True)
        if hasattr(self, "_draw_performance_profile"):
            self._draw_performance_profile()

    def _show_performance_test_info(self):
        """Open a dialog summarizing the performance test pattern, legs, transits, and timing."""
        if not GEOSPATIAL_LIBS_AVAILABLE or pyproj is None:
            self._show_message(
                "warning",
                "Disabled Feature",
                "Geospatial libraries not loaded. Cannot summarize performance test geometry.",
            )
            return
        lines = getattr(self, "performance_test_lines_data", None) or []
        if len(lines) != 4:
            self._show_message(
                "info",
                "No Performance Plan",
                "Plot performance test lines first (Plot Performance Lines).",
            )
            return

        bist_segs = getattr(self, "performance_bist_segments_data", None) or []
        has_bist = len(bist_segs) == 4
        if len(bist_segs) not in (0, 4):
            self._show_message(
                "warning",
                "BIST Geometry",
                "BIST segments are incomplete; remove lines and plot again, or ignore BIST legs below.",
            )
            has_bist = False

        try:
            speed_kts = float(self.performance_test_speed_entry.text().strip())
        except Exception:
            self._show_message("warning", "Invalid Speed", "Test Speed (kts) is not a valid number.")
            return
        if not np.isfinite(speed_kts) or speed_kts < 0:
            self._show_message("warning", "Invalid Speed", "Test Speed (kts) must be zero or greater.")
            return

        try:
            turn_min = float(self.performance_turn_time_entry.text().strip())
        except Exception:
            turn_min = 10.0
        if not np.isfinite(turn_min) or turn_min < 0:
            turn_min = 0.0

        geod = pyproj.Geod(ellps="WGS84")
        speed_mps = speed_kts * 0.514444

        def dist_m(lat1, lon1, lat2, lon2):
            _, _, d = geod.inv(lon1, lat1, lon2, lat2)
            return float(d)

        def leg_time_sec(d_m):
            if speed_mps <= 0 or d_m <= 0:
                return 0.0
            return d_m / speed_mps

        lines = [list(line) for line in lines]
        for i, seg in enumerate(lines):
            if len(seg) != 2:
                self._show_message("warning", "Invalid Lines", "Each performance line must have a start and end.")
                return
        if has_bist:
            for i, seg in enumerate(bist_segs):
                if len(seg) != 2:
                    self._show_message("warning", "Invalid BIST", "Each BIST segment must have two endpoints.")
                    return

        per_line_rows = []
        total_line_dist_m = 0.0
        total_line_time_sec = 0.0
        for idx in range(4):
            n = idx + 1
            p_s = lines[idx][0]
            p_e = lines[idx][1]
            sw_out = dist_m(p_s[0], p_s[1], p_e[0], p_e[1])
            sw_back = dist_m(p_e[0], p_e[1], p_s[0], p_s[1])
            if has_bist:
                b_s = bist_segs[idx][0]
                b_e = bist_segs[idx][1]
                b_out = dist_m(b_s[0], b_s[1], b_e[0], b_e[1])
                b_back = dist_m(b_e[0], b_e[1], b_s[0], b_s[1])
            else:
                b_out = b_back = 0.0
            line_dist = sw_out + sw_back + b_out + b_back
            line_time = leg_time_sec(sw_out) + leg_time_sec(sw_back) + leg_time_sec(b_out) + leg_time_sec(b_back)
            total_line_dist_m += line_dist
            total_line_time_sec += line_time
            per_line_rows.append(
                {
                    "n": n,
                    "sw_out": sw_out,
                    "sw_back": sw_back,
                    "b_out": b_out,
                    "b_back": b_back,
                    "line_dist": line_dist,
                    "line_time_sec": line_time,
                }
            )

        transit_rows = []
        total_transit_dist_m = 0.0
        total_transit_time_sec = 0.0
        for idx in range(3):
            n_from = idx + 1
            n_to = idx + 2
            a = lines[idx][0]
            b = lines[idx + 1][0]
            d = dist_m(a[0], a[1], b[0], b[1])
            t = leg_time_sec(d)
            transit_rows.append({"from": n_from, "to": n_to, "dist_m": d, "time_sec": t})
            total_transit_dist_m += d
            total_transit_time_sec += t

        num_turns = 8
        total_turn_time_sec = num_turns * turn_min * 60.0
        grand_dist_m = total_line_dist_m + total_transit_dist_m
        grand_time_sec = total_line_time_sec + total_transit_time_sec + total_turn_time_sec

        t = []
        t.append("PERFORMANCE TEST SUMMARY\n")
        t.append("=" * 24 + "\n")
        t.append(f"Test speed: {speed_kts:.2f} kts")
        if speed_mps <= 0:
            t.append(" (underway time at speed is not computed)\n")
        else:
            t.append(f" ({speed_mps:.3f} m/s)\n")
        t.append(f"Turn time (each of {num_turns} turns): {turn_min:.2f} min\n")
        t.append(f"BIST legs included: {'yes' if has_bist else 'no (swath only)'}\n\n")

        t.append("PATTERN (per line i = 1…4)\n")
        t.append("-" * 22 + "\n")
        t.append("1. Swath outbound:   PiS → PiE along the plotted performance line.\n")
        if has_bist:
            t.append(
                "2. BIST outbound:  BiS → BiE (BiS coincides with PiE, same heading as the line; RXnoise BIST).\n"
            )
            t.append("3. BIST inbound:     BiE → BiS.\n")
            t.append("4. Swath inbound:    PiE → PiS.\n")
        else:
            t.append("2. Swath inbound:    PiE → PiS.\n")
        t.append(
            "After finishing line i, the vessel is at PiS; transit along a geodesic to P(i+1)S before starting line i+1.\n"
        )
        t.append("\n")
        if has_bist:
            t.append(
                "BIST geometry matches the gold map segments: each segment runs from PiE to the outer BIST endpoint.\n\n"
            )
        else:
            t.append(
                "BIST Time was zero when lines were plotted, so this summary includes only swath legs (out and back).\n\n"
            )

        t.append("LEGS — distance (m) and time at test speed\n")
        t.append("-" * 42 + "\n")
        for row in per_line_rows:
            n = row["n"]
            t.append(f"Line {n}:\n")
            t.append(f"  P{n}S → P{n}E (swath out):     {row['sw_out']:.1f} m")
            if speed_mps > 0:
                t.append(f"  |  {leg_time_sec(row['sw_out']):.1f} s\n")
            else:
                t.append("\n")
            if has_bist:
                t.append(f"  B{n}S → B{n}E (BIST out):     {row['b_out']:.1f} m")
                if speed_mps > 0:
                    t.append(f"  |  {leg_time_sec(row['b_out']):.1f} s\n")
                else:
                    t.append("\n")
                t.append(f"  B{n}E → B{n}S (BIST return):   {row['b_back']:.1f} m")
                if speed_mps > 0:
                    t.append(f"  |  {leg_time_sec(row['b_back']):.1f} s\n")
                else:
                    t.append("\n")
            t.append(f"  P{n}E → P{n}S (swath return):  {row['sw_back']:.1f} m")
            if speed_mps > 0:
                t.append(f"  |  {leg_time_sec(row['sw_back']):.1f} s\n")
            else:
                t.append("\n")
            t.append(
                f"  Line {n} subtotal (pattern only): {row['line_dist']:.1f} m"
                + (f"  |  {row['line_time_sec']:.1f} s ({row['line_time_sec']/60.0:.2f} min)\n" if speed_mps > 0 else "\n")
            )
            t.append("\n")

        t.append("TRANSITS (geodesic PiS → P(i+1)S after completing line i)\n")
        t.append("-" * 52 + "\n")
        for row in transit_rows:
            t.append(
                f"  P{row['from']}S → P{row['to']}S:  {row['dist_m']:.1f} m"
                + (f"  |  {row['time_sec']:.1f} s\n" if speed_mps > 0 else "\n")
            )
        t.append("\n")

        t.append("TURNS\n")
        t.append("-" * 5 + "\n")
        t.append(f"  2 turns per line × 4 lines = {num_turns} turns × {turn_min:.2f} min = {total_turn_time_sec/60.0:.2f} min\n\n")

        t.append("TOTALS (all four line patterns + transits + turns)\n")
        t.append("-" * 48 + "\n")
        t.append(f"  Distance underway (legs + transits): {grand_dist_m:.1f} m ({grand_dist_m/1000.0:.3f} km)\n")
        if speed_mps > 0:
            t.append(
                f"  Time at speed (legs + transits):     {total_line_time_sec + total_transit_time_sec:.1f} s "
                f"({(total_line_time_sec + total_transit_time_sec)/60.0:.2f} min)\n"
            )
        t.append(f"  Turn time (fixed):                   {total_turn_time_sec:.1f} s ({total_turn_time_sec/60.0:.2f} min)\n")
        if speed_mps > 0:
            t.append(
                f"  Grand total (speed + turns):         {grand_time_sec:.1f} s ({grand_time_sec/60.0:.2f} min)\n"
            )
        t.append("\n")

        t.append("ORDERED WAYPOINTS (full sequence)\n")
        t.append("-" * 31 + "\n")
        order_labels = []
        for idx in range(4):
            n = idx + 1
            order_labels.append(f"P{n}S")
            order_labels.append(f"P{n}E")
            if has_bist:
                order_labels.append(f"B{n}E")
                order_labels.append(f"P{n}E")
            order_labels.append(f"P{n}S")
            if idx < 3:
                order_labels.append(f"P{idx+2}S")
        t.append(" → ".join(order_labels) + "\n\n")

        def append_waypoints_block(title, use_ddm):
            t.append(f"{title}\n")
            t.append("-" * len(title) + "\n")
            for idx in range(4):
                n = idx + 1
                p_s = lines[idx][0]
                p_e = lines[idx][1]
                if use_ddm:
                    t.append(
                        f"P{n}S: {decimal_degrees_to_ddm(p_s[0], is_latitude=True)}, "
                        f"{decimal_degrees_to_ddm(p_s[1], is_latitude=False)}\n"
                    )
                    t.append(
                        f"P{n}E: {decimal_degrees_to_ddm(p_e[0], is_latitude=True)}, "
                        f"{decimal_degrees_to_ddm(p_e[1], is_latitude=False)}\n"
                    )
                    if has_bist:
                        b_e = bist_segs[idx][1]
                        t.append(
                            f"B{n}E: {decimal_degrees_to_ddm(b_e[0], is_latitude=True)}, "
                            f"{decimal_degrees_to_ddm(b_e[1], is_latitude=False)}\n"
                        )
                else:
                    t.append(f"P{n}S: {p_s[0]:.6f}, {p_s[1]:.6f}\n")
                    t.append(f"P{n}E: {p_e[0]:.6f}, {p_e[1]:.6f}\n")
                    if has_bist:
                        b_e = bist_segs[idx][1]
                        t.append(f"B{n}E: {b_e[0]:.6f}, {b_e[1]:.6f}\n")
            t.append("\n")

        append_waypoints_block("Waypoints (DMM)", True)
        append_waypoints_block("Waypoints (decimal degrees)", False)

        show_statistics_dialog(self, "Performance Test Info", "".join(t))

    def _zoom_to_performance_lines(self):
        """Zoom map view to all currently plotted performance test lines."""
        if not hasattr(self, "performance_test_lines_data") or len(self.performance_test_lines_data) != 4:
            self._show_message("info", "Zoom to Performance Lines", "No performance lines to zoom to.")
            return
        points = []
        for line in self.performance_test_lines_data:
            points.extend(line)
        if hasattr(self, "performance_bist_segments_data"):
            for segment in self.performance_bist_segments_data:
                points.extend(segment)
        if not points:
            self._show_message("info", "Zoom to Performance Lines", "No performance lines to zoom to.")
            return
        lats = [p[0] for p in points]
        lons = [p[1] for p in points]
        min_lat, max_lat = min(lats), max(lats)
        min_lon, max_lon = min(lons), max(lons)
        buffer_lat = (max_lat - min_lat) * 0.05 if (max_lat - min_lat) != 0 else 0.01
        buffer_lon = (max_lon - min_lon) * 0.05 if (max_lon - min_lon) != 0 else 0.01
        xlim = (min_lon - buffer_lon, max_lon + buffer_lon)
        ylim = (min_lat - buffer_lat, max_lat + buffer_lat)
        self._apply_map_zoom_limits_and_reload_geotiff(xlim, ylim)

    def _clear_performance_lines(self):
        """Remove performance test lines from map and clear performance profile line."""
        self.performance_test_lines_data = []
        self.performance_bist_segments_data = []
        self.performance_profile_line = []
        self.performance_central_point_coords = (None, None)
        self._plot_survey_plan(preserve_view_limits=True)
        if hasattr(self, "_draw_performance_profile"):
            self._draw_performance_profile()

    def set_performance_activity_text(self, message, append=False):
        """Append or replace a line in the Activity Log with a [Performance] prefix."""
        if not hasattr(self, "activity_log_text"):
            return
        prefixed_message = f"[Performance] {message}"
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

    def _update_performance_export_name(self):
        if not hasattr(self, "performance_export_name_entry"):
            return
        try:
            swath_a = self.performance_swath_angle_entry.text().strip() or "0"
            spd = self.performance_test_speed_entry.text().strip() or "0"
            self.performance_export_name_entry.setText(f"Performance_{swath_a}deg_{spd}_kts")
        except Exception:
            pass

    def _apply_performance_assignments(self, assignments, imported_lines):
        """Build performance_test_lines_data and performance_bist_segments_data from dialog."""
        swath = {}
        bist = {}
        for line_idx, role in assignments.items():
            if role.startswith("swath_"):
                swath[int(role.split("_")[1])] = imported_lines[line_idx]
            elif role.startswith("bist_"):
                bist[int(role.split("_")[1])] = imported_lines[line_idx]
        self.performance_test_lines_data = [swath[i] for i in range(1, 5)]
        if len(bist) == 4:
            self.performance_bist_segments_data = [bist[i] for i in range(1, 5)]
        else:
            self.performance_bist_segments_data = []
        self.performance_profile_line = self.performance_test_lines_data[0]
        self.central_point_coords = self.performance_test_lines_data[0][0]

    def _try_parse_performance_csv(self, file_path):
        """Return (lines4, bist_segs) if CSV uses performance P/B labels; else (None, None)."""
        try:
            with open(file_path, "r", encoding="utf-8") as csvfile:
                reader = csv.DictReader(csvfile)
                lines_data = {}
                for row in reader:
                    try:
                        line_num = int(row.get("Line Number", -1))
                        point_label = (row.get("Point Label") or "").strip()
                        lat = float(row.get("Latitude", 0))
                        lon = float(row.get("Longitude", 0))
                    except (ValueError, TypeError):
                        continue
                    if line_num < 0:
                        continue
                    lines_data.setdefault(line_num, {})[point_label] = (lat, lon)
            swath_lines = []
            for i in range(1, 5):
                pts = lines_data.get(i, {})
                s = pts.get(f"P{i}S")
                e = pts.get(f"P{i}E")
                if s is None or e is None:
                    return None, None
                swath_lines.append([s, e])
            bist_lines = []
            for i in range(1, 5):
                key = 10 + i
                pts = lines_data.get(key, {})
                s = pts.get(f"B{i}S")
                e = pts.get(f"B{i}E")
                if s is None or e is None:
                    bist_lines = []
                    break
                bist_lines.append([s, e])
            if len(bist_lines) == 4:
                return swath_lines, bist_lines
            return swath_lines, []
        except Exception:
            return None, None

    def _try_parse_performance_geojson(self, geojson_data):
        """Return (lines4, bist_segs) if unambiguous; else (None, None)."""
        if not isinstance(geojson_data, dict):
            return None, None
        if geojson_data.get("type") == "FeatureCollection":
            features = geojson_data.get("features", [])
        elif geojson_data.get("type") == "Feature":
            features = [geojson_data]
        else:
            features = []
        swath = {}
        bist = {}
        for feature in features:
            if not isinstance(feature, dict) or feature.get("type") != "Feature":
                continue
            geom = feature.get("geometry") or {}
            if geom.get("type") != "LineString":
                continue
            coords = geom.get("coordinates") or []
            if len(coords) < 2:
                continue
            p0 = (coords[0][1], coords[0][0])
            p1 = (coords[-1][1], coords[-1][0])
            props = feature.get("properties") or {}
            line_num = props.get("line_num")
            try:
                line_num = int(line_num)
            except (TypeError, ValueError):
                line_num = None
            if line_num is not None and 1 <= line_num <= 4:
                swath[line_num] = [p0, p1]
            elif line_num is not None and 11 <= line_num <= 14:
                bist[line_num - 10] = [p0, p1]
        if len(swath) != 4:
            return None, None
        lines = [swath[i] for i in range(1, 5)]
        if len(bist) == 0:
            return lines, []
        if len(bist) == 4:
            return lines, [bist[i] for i in range(1, 5)]
        return None, None

    def _geojson_to_unassigned_segments(self, geojson_data):
        """Each LineString feature -> one 2-point segment (first/last), for assignment dialog."""
        if not isinstance(geojson_data, dict):
            return []
        if geojson_data.get("type") == "FeatureCollection":
            features = geojson_data.get("features", [])
        elif geojson_data.get("type") == "Feature":
            features = [geojson_data]
        else:
            features = []
        out = []
        for feature in features:
            if not isinstance(feature, dict) or feature.get("type") != "Feature":
                continue
            geom = feature.get("geometry") or {}
            if geom.get("type") != "LineString":
                continue
            coords = geom.get("coordinates") or []
            if len(coords) < 2:
                continue
            p0 = (coords[0][1], coords[0][0])
            p1 = (coords[-1][1], coords[-1][0])
            out.append([p0, p1])
        return out

    def _download_and_load_gmrt_after_perf_import(self):
        """Download GMRT GeoTIFF covering performance lines + BIST with user buffer (degrees)."""
        all_points = []
        for line in getattr(self, "performance_test_lines_data", []) or []:
            all_points.extend(line)
        for seg in getattr(self, "performance_bist_segments_data", []) or []:
            all_points.extend(seg)
        if not all_points:
            self._show_message("warning", "GMRT Download", "No performance line points to compute extent.")
            return
        lats = [p[0] for p in all_points]
        lons = [p[1] for p in all_points]
        mid_lat = (min(lats) + max(lats)) / 2.0
        mid_lon = (min(lons) + max(lons)) / 2.0
        buffer_deg = 0.5
        if hasattr(self, "performance_gmrt_buffer_spin"):
            try:
                buffer_deg = float(self.performance_gmrt_buffer_spin.value())
            except (ValueError, TypeError):
                pass
        west = mid_lon - buffer_deg
        east = mid_lon + buffer_deg
        south = mid_lat - buffer_deg
        north = mid_lat + buffer_deg
        self._download_gmrt_and_load(
            west,
            east,
            south,
            north,
            resolution=100,
            layer="topo",
            default_filename_prefix="GMRT_Bathy",
            log_func=lambda msg, append=True: self.set_performance_activity_text(msg, append=append),
            default_directory=getattr(self, "last_perf_import_dir", None),
        )

    def _perf_import_post_import(self, file_path):
        """Populate performance fields from optional metadata; refresh map and profile."""
        lines = getattr(self, "performance_test_lines_data", []) or []
        if len(lines) != 4:
            self._show_message("warning", "Import Warning", "Need four performance swath lines.")
            return
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        dir_name = os.path.dirname(file_path)
        meta_path = os.path.join(dir_name, f"{base_name}_performance_params.json")
        params = None
        if os.path.exists(meta_path):
            try:
                with open(meta_path, "r", encoding="utf-8") as f:
                    params = json.load(f)
            except Exception as e:
                print(f"Warning: could not load performance metadata: {e}")
        all_pts = []
        for ln in lines:
            all_pts.extend(ln)
        for seg in getattr(self, "performance_bist_segments_data", []) or []:
            all_pts.extend(seg)
        try:
            if params:
                if params.get("central_lat") is not None and hasattr(self, "performance_central_lat_entry"):
                    self.performance_central_lat_entry.setText(str(params["central_lat"]))
                if params.get("central_lon") is not None and hasattr(self, "performance_central_lon_entry"):
                    self.performance_central_lon_entry.setText(str(params["central_lon"]))
                if params.get("swell_direction_deg") is not None and hasattr(self, "performance_swell_direction_entry"):
                    self.performance_swell_direction_entry.setText(str(params["swell_direction_deg"]))
                if params.get("swath_angle_deg") is not None and hasattr(self, "performance_swath_angle_entry"):
                    self.performance_swath_angle_entry.setText(str(params["swath_angle_deg"]))
                if params.get("bist_time_min") is not None and hasattr(self, "performance_bist_time_entry"):
                    self.performance_bist_time_entry.setText(str(params["bist_time_min"]))
                if params.get("turn_time_min") is not None and hasattr(self, "performance_turn_time_entry"):
                    self.performance_turn_time_entry.setText(str(params["turn_time_min"]))
                if params.get("test_speed_kts") is not None and hasattr(self, "performance_test_speed_entry"):
                    self.performance_test_speed_entry.setText(str(params["test_speed_kts"]))
                if params.get("num_pings") is not None and hasattr(self, "performance_num_pings_entry"):
                    self.performance_num_pings_entry.setText(str(params["num_pings"]))
                if params.get("sound_velocity") is not None and hasattr(self, "performance_sound_velocity_entry"):
                    self.performance_sound_velocity_entry.setText(str(params["sound_velocity"]))
                en = params.get("export_name")
                if en and hasattr(self, "performance_export_name_entry"):
                    self.performance_export_name_entry.setText(str(en))
                gtp = params.get("geotiff_path")
                if gtp and hasattr(self, "_load_geotiff_from_path") and os.path.exists(gtp):
                    self._load_geotiff_from_path(gtp)
            elif all_pts and pyproj is not None:
                geod = pyproj.Geod(ellps="WGS84")
                all_lats = [p[0] for p in all_pts]
                all_lons = [p[1] for p in all_pts]
                c_lat = (min(all_lats) + max(all_lats)) / 2.0
                c_lon = (min(all_lons) + max(all_lons)) / 2.0
                if hasattr(self, "performance_central_lat_entry"):
                    self.performance_central_lat_entry.setText(f"{c_lat:.6f}")
                if hasattr(self, "performance_central_lon_entry"):
                    self.performance_central_lon_entry.setText(f"{c_lon:.6f}")
                fl = lines[0]
                lat1, lon1 = fl[0]
                lat2, lon2 = fl[1]
                try:
                    fwd_az, _, dist = geod.inv(lon1, lat1, lon2, lat2)
                    if hasattr(self, "performance_swell_direction_entry") and not self.performance_swell_direction_entry.text().strip():
                        self.performance_swell_direction_entry.setText(f"{(fwd_az % 360):.1f}")
                except Exception:
                    pass
                if hasattr(self, "performance_export_name_entry") and not self.performance_export_name_entry.text().strip():
                    self.performance_export_name_entry.setText(base_name)
        except Exception as e:
            print(f"Warning: performance post-import field update: {e}")

        self.performance_profile_line = lines[0]
        self.performance_central_point_coords = lines[0][0]
        self.central_point_coords = lines[0][0]
        if hasattr(self, "_update_performance_line_length"):
            self._update_performance_line_length()
        self._plot_survey_plan(preserve_view_limits=True)
        if hasattr(self, "_draw_performance_profile"):
            self._draw_performance_profile()

        parts = ["4 performance swath lines"]
        if len(getattr(self, "performance_bist_segments_data", []) or []) == 4:
            parts.append("4 BIST segments")
        msg = f"Successfully imported: {', '.join(parts)}"
        if params:
            msg += " (parameters from performance metadata JSON)"
        else:
            msg += " (geometry only; fill swell direction if needed)"
        self.set_performance_activity_text(msg, append=False)
        if getattr(self, "performance_download_gmrt_checkbox", None) and self.performance_download_gmrt_checkbox.isChecked():
            self._download_and_load_gmrt_after_perf_import()

    def _import_performance_survey(self):
        """Import performance plan from the same file types as Accuracy; assign lines when needed."""
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded.")
            return

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Performance Survey File to Import",
            getattr(self, "last_perf_import_dir", os.path.expanduser("~")),
            "Known Survey Files (*_DMS.txt *_DMM.txt *_DDD.txt *_DDD.csv *.csv *.geojson *.json *.lnw);;"
            "Hypack LNW files (*.lnw);;Degrees Minutes Seconds (*_DMS.txt);;Degrees Decimal Minutes (*_DMM.txt);;"
            "Decimal Degrees (*_DDD.txt);;Decimal Degree CSV (*_DDD.csv);;CSV (*.csv);;GeoJSON (*.geojson);;JSON (*.json)",
        )
        if not file_path:
            return
        import_dir = os.path.dirname(file_path)
        if import_dir and os.path.isdir(import_dir):
            self.last_perf_import_dir = import_dir
            self._save_last_perf_import_dir()

        try:
            file_ext = os.path.splitext(file_path)[1].lower()
            file_basename = os.path.basename(file_path)
            file_processed = False

            if file_ext == ".lnw":
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
                dialog = PerformanceLineAssignmentDialog(self, imported_lines)
                if dialog.exec() != QDialog.DialogCode.Accepted:
                    return
                self._apply_performance_assignments(dialog.get_assignments(), imported_lines)
                file_processed = True

            if not file_processed and file_basename.lower().endswith("_dms.txt"):
                imported_lines = self._parse_dms_txt_file(file_path)
                if imported_lines is None:
                    return
                dialog = PerformanceLineAssignmentDialog(self, imported_lines)
                if dialog.exec() != QDialog.DialogCode.Accepted:
                    return
                self._apply_performance_assignments(dialog.get_assignments(), imported_lines)
                file_processed = True

            if not file_processed and file_basename.lower().endswith("_dmm.txt"):
                imported_lines = self._parse_dmm_txt_file(file_path)
                if imported_lines is None:
                    return
                dialog = PerformanceLineAssignmentDialog(self, imported_lines)
                if dialog.exec() != QDialog.DialogCode.Accepted:
                    return
                self._apply_performance_assignments(dialog.get_assignments(), imported_lines)
                file_processed = True

            if not file_processed and file_basename.lower().endswith("_ddd.txt"):
                imported_lines = self._parse_ddd_txt_file(file_path)
                if imported_lines is None:
                    return
                dialog = PerformanceLineAssignmentDialog(self, imported_lines)
                if dialog.exec() != QDialog.DialogCode.Accepted:
                    return
                self._apply_performance_assignments(dialog.get_assignments(), imported_lines)
                file_processed = True

            if file_processed:
                self._perf_import_post_import(file_path)
                return

            if file_ext == ".csv":
                swath, bist = self._try_parse_performance_csv(file_path)
                if swath is not None:
                    self.performance_test_lines_data = swath
                    self.performance_bist_segments_data = bist if bist else []
                    self._perf_import_post_import(file_path)
                    return
                self._show_message(
                    "warning",
                    "CSV Import",
                    "Could not parse as performance CSV (expected Line 1–4 with P1S/P1E…P4E; "
                    "optional lines 11–14 with B1S/B1E…B4E). Use GeoJSON/LNW/text formats or export from this app.",
                )
                return

            if file_ext in (".geojson", ".json"):
                with open(file_path, "r", encoding="utf-8") as f:
                    geojson_data = json.load(f)
                swath, bist = self._try_parse_performance_geojson(geojson_data)
                if swath is not None:
                    self.performance_test_lines_data = swath
                    self.performance_bist_segments_data = bist if bist else []
                    gtp = (geojson_data.get("properties") or {}).get("geotiff_path")
                    if not gtp:
                        for feat in geojson_data.get("features") or []:
                            p = (feat.get("properties") or {}).get("geotiff_path")
                            if p:
                                gtp = p
                                break
                    if gtp and hasattr(self, "_load_geotiff_from_path") and os.path.exists(gtp):
                        self._load_geotiff_from_path(gtp)
                    spd = None
                    for feat in geojson_data.get("features") or []:
                        pr = feat.get("properties") or {}
                        for key in ("test_speed_kts", "survey_speed"):
                            if pr.get(key) is not None:
                                try:
                                    spd = float(pr[key])
                                    break
                                except (TypeError, ValueError):
                                    pass
                        if spd is not None:
                            break
                    if spd is not None and hasattr(self, "performance_test_speed_entry"):
                        self.performance_test_speed_entry.setText(str(spd))
                    self._perf_import_post_import(file_path)
                    return
                imported_lines = self._geojson_to_unassigned_segments(geojson_data)
                if not imported_lines:
                    self._show_message("warning", "Import Warning", "No LineString features found in GeoJSON.")
                    return
                dialog = PerformanceLineAssignmentDialog(self, imported_lines)
                if dialog.exec() != QDialog.DialogCode.Accepted:
                    return
                self._apply_performance_assignments(dialog.get_assignments(), imported_lines)
                self._perf_import_post_import(file_path)
                return

            self._show_message("error", "Import Error", f"Unsupported file format: {file_ext}")

        except Exception as e:
            self._show_message("error", "Import Error", f"Failed to import performance survey: {e}")
            import traceback

            traceback.print_exc()

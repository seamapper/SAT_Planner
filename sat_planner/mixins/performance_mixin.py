"""
Performance tab behavior: swath ping-time calculations and pick-center field updates.
"""

import json
import time

import numpy as np
from sat_planner.constants import GEOSPATIAL_LIBS_AVAILABLE, pyproj


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

    def _update_performance_center_from_pick(self, clicked_lat, clicked_lon):
        """Update Performance tab center fields from map pick-center click."""
        if hasattr(self, "performance_central_lat_entry"):
            self.performance_central_lat_entry.setText(f"{clicked_lat:.6f}")
        if hasattr(self, "performance_central_lon_entry"):
            self.performance_central_lon_entry.setText(f"{clicked_lon:.6f}")

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

    def _plot_performance_test_lines(self):
        """Plot 4 performance test lines from center/length/swell and draw profile for line 1."""
        if not GEOSPATIAL_LIBS_AVAILABLE or pyproj is None:
            self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded. Cannot plot test lines.")
            return
        if getattr(self, "geotiff_data_array", None) is None:
            self._show_message("warning", "No GeoTIFF", "Load a GeoTIFF first to plot performance test lines.")
            return

        try:
            center_lat = float(self.performance_central_lat_entry.text().strip())
            center_lon = float(self.performance_central_lon_entry.text().strip())
            line_length_m = float(self._performance_get_text("performance_line_length_m_entry"))
            swell_direction = float(self.performance_swell_direction_entry.text().strip())
            bist_min = float(self.performance_bist_time_entry.text().strip()) if self.performance_bist_time_entry.text().strip() else 0.0
            speed_kts = float(self.performance_test_speed_entry.text().strip())
        except Exception:
            self._show_message(
                "warning",
                "Invalid Inputs",
                "Provide valid Central Latitude, Central Longitude, Swell Direction, Test Speed, and calculated Line Length (m).",
            )
            return

        if not np.isfinite(line_length_m) or line_length_m <= 0:
            self._show_message("warning", "Invalid Line Length", "Line Length (m) must be greater than zero.")
            return
        if not np.isfinite(center_lat) or not np.isfinite(center_lon):
            self._show_message("warning", "Invalid Center", "Central Latitude/Longitude values are not valid.")
            return
        if bist_min < 0 or not np.isfinite(bist_min):
            self._show_message("warning", "Invalid BIST Time", "BIST Time (min) must be zero or greater.")
            return
        if speed_kts < 0 or not np.isfinite(speed_kts):
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
        self.ax.set_xlim(min_lon - buffer_lon, max_lon + buffer_lon)
        self.ax.set_ylim(min_lat - buffer_lat, max_lat + buffer_lat)
        self.canvas.draw_idle()

    def _clear_performance_lines(self):
        """Remove performance test lines from map and clear performance profile line."""
        self.performance_test_lines_data = []
        self.performance_bist_segments_data = []
        self.performance_profile_line = []
        self.performance_central_point_coords = (None, None)
        self._plot_survey_plan(preserve_view_limits=True)
        if hasattr(self, "_draw_performance_profile"):
            self._draw_performance_profile()

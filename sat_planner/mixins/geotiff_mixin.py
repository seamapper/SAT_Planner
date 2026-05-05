"""
GeoTIFF loading, display mode, dynamic resolution, and contour controls.
"""
import os
import csv
import json
import datetime
import xml.etree.ElementTree as ET
import numpy as np
from PyQt6.QtWidgets import (
    QApplication,
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QProgressBar,
    QPushButton,
    QVBoxLayout,
    QWidget,
)
from PyQt6.QtCore import Qt, QTimer
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from sat_planner.constants import (
    GEOSPATIAL_LIBS_AVAILABLE,
    LineString,
    RasterioIOError,
    Resampling,
    Window,
    window_bounds,
    rowcol,
    reproject,
    transform,
    pyproj,
    rasterio,
)
from sat_planner import export_utils
from sat_planner.utils_geo import decimal_degrees_to_ddm
from sat_planner.utils_geo import lat_lon_decimal_from_survey_csv_row
from sat_planner.utils_ui import show_statistics_dialog


class GeoTIFFMixin:
    """Mixin providing GeoTIFF load/remove, display mode, dynamic resolution, and contours."""

    def _dynamic_resolution_reference_extent(self):
        """Return the full GeoTIFF extent to use for zoom/resolution calculations."""
        if hasattr(self, "geotiff_original_extent") and self.geotiff_original_extent is not None:
            return self.geotiff_original_extent
        return getattr(self, "geotiff_extent", None)

    def _get_geotiff_nan_cutoff(self):
        """Return current lower-bound cutoff for setting GeoTIFF values to NaN."""
        try:
            return float(getattr(self, "geotiff_nan_value", -11000.0))
        except (TypeError, ValueError):
            return -11000.0

    def _detect_geotiff_nan_cutoff_from_dataset(self):
        """Choose default NaN cutoff from loaded dataset nodata when available."""
        ds = getattr(self, "geotiff_dataset_original", None)
        if ds is not None:
            nodata_val = getattr(ds, "nodata", None)
            try:
                if nodata_val is not None and np.isfinite(float(nodata_val)):
                    return float(nodata_val)
            except (TypeError, ValueError):
                pass
        return -11000.0

    def _calculate_backscatter_line_statistics(self):
        """Compute survey statistics for the backscatter centerline (with lead-in/out)."""
        segment = getattr(self, "_get_backscatter_centerline_with_lead_in", lambda: None)()
        centerline = getattr(self, "backscatter_box_centerline", None)
        if (
            not segment
            or len(segment) != 2
            or segment[0] is None
            or segment[1] is None
            or not centerline
            or len(centerline) != 2
            or centerline[0] is None
            or centerline[1] is None
        ):
            return None
        try:
            (lat1, lon1), (lat2, lon2) = segment
            (a_lat, a_lon), (b_lat, b_lon) = centerline
            geod = pyproj.Geod(ellps="WGS84")
            fwd_az, back_az, total_distance = geod.inv(lon1, lat1, lon2, lat2)
            _, _, lead_in_m = geod.inv(lon1, lat1, a_lon, a_lat)
            _, _, lead_out_m = geod.inv(b_lon, b_lat, lon2, lat2)
            try:
                speed_knots = float(self.backscatter_survey_speed_entry.text()) if self.backscatter_survey_speed_entry.text() else 8.0
            except Exception:
                speed_knots = 8.0
            speed_ms = speed_knots * 0.514444
            total_time_seconds = total_distance / speed_ms if speed_ms > 0 else 0.0
            total_time_minutes = total_time_seconds / 60.0
            total_time_hours = total_time_minutes / 60.0

            dists, elevations, slopes = (None, None, None)
            if hasattr(self, "_profile_arrays_along_segment_endpoints"):
                dists, elevations, slopes = self._profile_arrays_along_segment_endpoints(lat1, lon1, lat2, lon2, n=250)
            depth_info = ""
            if elevations is not None and np.size(elevations) > 0 and np.any(np.isfinite(elevations)):
                valid_e = elevations[np.isfinite(elevations)]
                shallowest_depth_m = float(np.min(np.abs(valid_e)))
                max_depth_m = float(np.max(np.abs(valid_e)))
                depth_info += f"Shallowest Depth: {shallowest_depth_m:.1f} m\n"
                depth_info += f"Maximum Depth: {max_depth_m:.1f} m\n"
            if slopes is not None and np.size(slopes) > 0 and np.any(np.isfinite(slopes)):
                valid_s = slopes[np.isfinite(slopes)]
                depth_info += f"Minimum Slope: {float(np.min(valid_s)):.2f} deg\n"
                depth_info += f"Maximum Slope: {float(np.max(valid_s)):.2f} deg\n"

            return {
                "total_distance_m": float(total_distance),
                "total_distance_km": float(total_distance / 1000.0),
                "total_distance_nm": float(total_distance / 1852.0),
                "heading": float(fwd_az % 360.0),
                "reciprocal_heading": float(back_az % 360.0),
                "lead_in_m": float(lead_in_m),
                "lead_out_m": float(lead_out_m),
                "speed_knots": float(speed_knots),
                "total_time_minutes": float(total_time_minutes),
                "total_time_hours": float(total_time_hours),
                "depth_info": depth_info.strip(),
                "line_waypoints": [(lat1, lon1), (a_lat, a_lon), (b_lat, b_lon), (lat2, lon2)],
                "area_corners": list(getattr(self, "backscatter_box_vertices", []) or []),
            }
        except Exception as e:
            print(f"Error calculating backscatter line statistics: {e}")
            return None

    def _show_backscatter_line_information(self):
        """Show survey info dialog for backscatter line (with lead-in/out)."""
        stats = self._calculate_backscatter_line_statistics()
        if not stats:
            self._show_message("warning", "No Backscatter Line", "No backscatter line is defined. Select an area/line first.")
            return
        stats_text = "BACKSCATTER LINE STATISTICS\n" + "=" * 31 + "\n\n"
        stats_text += (
            f"Total Distance: {stats['total_distance_m']:.1f} m "
            f"({stats['total_distance_km']:.3f} km, {stats['total_distance_nm']:.3f} nm)\n"
        )
        stats_text += f"Heading: {stats['heading']:.1f}° (Reciprocal {stats['reciprocal_heading']:.1f}°)\n"
        stats_text += (
            f"Lead-in/out Length: {stats['lead_in_m']:.1f} m / {stats['lead_out_m']:.1f} m "
            f"({stats['lead_in_m']/1852.0:.3f} / {stats['lead_out_m']/1852.0:.3f} nm)\n"
        )
        stats_text += f"Survey Speed: {stats['speed_knots']:.1f} knots\n"
        stats_text += (
            f"Estimated Time: {stats['total_time_minutes']:.1f} min "
            f"({stats['total_time_hours']:.2f} hr)\n\n"
        )
        if stats["depth_info"]:
            stats_text += stats["depth_info"] + "\n\n"

        line_waypoints = stats.get("line_waypoints", [])
        if len(line_waypoints) == 4:
            dmm_heading = "Backscatter Line Waypoints (DMM)"
            stats_text += f"{dmm_heading}\n" + "-" * len(dmm_heading) + "\n"
            labels = ["Lead-in Start", "Area Start", "Area End", "Lead-out End"]
            for i, (lat, lon) in enumerate(line_waypoints):
                stats_text += f"{labels[i]}: {decimal_degrees_to_ddm(lat, True)}, {decimal_degrees_to_ddm(lon, False)}\n"
            ddd_heading = "Backscatter Line Waypoints (DDD)"
            stats_text += f"\n{ddd_heading}\n" + "-" * len(ddd_heading) + "\n"
            for i, (lat, lon) in enumerate(line_waypoints):
                stats_text += f"{labels[i]}: {lat:.6f}, {lon:.6f}\n"

        area_corners = stats.get("area_corners", [])
        if len(area_corners) == 4:
            dmm_heading = "Backscatter Area Corners (DMM)"
            stats_text += f"\n{dmm_heading}\n" + "-" * len(dmm_heading) + "\n"
            for i, (lon, lat) in enumerate(area_corners):
                stats_text += f"C{i+1}: {decimal_degrees_to_ddm(lat, True)}, {decimal_degrees_to_ddm(lon, False)}\n"
            ddd_heading = "Backscatter Area Corners (DDD)"
            stats_text += f"\n{ddd_heading}\n" + "-" * len(ddd_heading) + "\n"
            for i, (lon, lat) in enumerate(area_corners):
                stats_text += f"C{i+1}: {lat:.6f}, {lon:.6f}\n"

        show_statistics_dialog(self, "Survey Info", stats_text)

    def _update_backscatter_export_name_default(self):
        """Set default backscatter export name: BS_YYYYMMDD_<mean depth along line>."""
        if not hasattr(self, "backscatter_export_name_entry"):
            return
        date_str = datetime.datetime.now().strftime("%Y%m%d")
        mean_depth_m = 0.0
        try:
            segment = getattr(self, "_get_backscatter_centerline_with_lead_in", lambda: None)()
            if segment and len(segment) == 2 and segment[0] is not None and segment[1] is not None:
                (lat1, lon1), (lat2, lon2) = segment
                dists, elevations, _ = self._profile_arrays_along_segment_endpoints(lat1, lon1, lat2, lon2, n=250)
                if elevations is not None and np.size(elevations) > 0 and np.any(np.isfinite(elevations)):
                    mean_depth_m = float(np.nanmean(np.abs(elevations[np.isfinite(elevations)])))
        except Exception:
            mean_depth_m = 0.0
        self.backscatter_export_name_entry.setText(f"BS_{date_str}_{mean_depth_m:.0f}")

    def _download_and_load_gmrt_after_backscatter_import(self):
        """Download and load GMRT data around imported backscatter line extent."""
        segment = getattr(self, "_get_backscatter_centerline_with_lead_in", lambda: None)()
        if not segment or len(segment) != 2 or segment[0] is None or segment[1] is None:
            self._show_message("warning", "GMRT Download", "No backscatter line to compute extent.")
            return
        (lat1, lon1), (lat2, lon2) = segment
        mid_lat = (lat1 + lat2) / 2.0
        mid_lon = (lon1 + lon2) / 2.0
        buffer_deg = 0.5
        if hasattr(self, "backscatter_gmrt_buffer_spin"):
            try:
                buffer_deg = float(self.backscatter_gmrt_buffer_spin.value())
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
            log_func=lambda msg, append=True: self.set_line_info_text(msg, append=append) if hasattr(self, "set_line_info_text") else None,
            default_directory=getattr(self, "last_backscatter_import_dir", None),
        )

    def _import_backscatter_line(self):
        """Import backscatter line/area from supported survey files + optional params sidecar."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select Backscatter Line File to Import", self.last_backscatter_import_dir,
            "Known Line Files (*_DMS.txt *_DMM.txt *_DDD.txt *.lnw *_DDD.csv *_DMM.csv *_DMS.csv *.csv *.geojson *.json *.gpx);;"
            "All files (*.*)"
        )
        if not file_path:
            return
        import_dir = os.path.dirname(file_path)
        if import_dir and os.path.isdir(import_dir):
            self.last_backscatter_import_dir = import_dir
            if hasattr(self, "_save_last_backscatter_import_dir"):
                self._save_last_backscatter_import_dir()

        points = []
        ext = os.path.splitext(file_path)[1].lower()
        try:
            if file_path.lower().endswith("_ddd.txt") and hasattr(self, "_parse_ddd_txt_file_as_polyline"):
                points = self._parse_ddd_txt_file_as_polyline(file_path) or []
            elif file_path.lower().endswith("_dms.txt") and hasattr(self, "_parse_dms_txt_file_as_polyline"):
                points = self._parse_dms_txt_file_as_polyline(file_path) or []
            elif file_path.lower().endswith("_dmm.txt") and hasattr(self, "_parse_dmm_txt_file_as_polyline"):
                points = self._parse_dmm_txt_file_as_polyline(file_path) or []
            elif ext == ".lnw" and hasattr(self, "_parse_lnw_file_as_polyline"):
                from .survey_parsers_mixin import UTMZoneDialog
                utm_dialog = UTMZoneDialog(self)
                if utm_dialog.exec() != QDialog.DialogCode.Accepted:
                    return
                zone, hem = utm_dialog.get_utm_info()
                points = self._parse_lnw_file_as_polyline(file_path, zone, hem) or []
            elif ext == ".gpx" and hasattr(self, "_parse_gpx_file_as_polyline"):
                points = self._parse_gpx_file_as_polyline(file_path) or []
            elif ext == ".csv":
                with open(file_path, "r", encoding="utf-8") as csvfile:
                    reader = csv.DictReader(csvfile)
                    for row in reader:
                        try:
                            lat, lon = lat_lon_decimal_from_survey_csv_row(row)
                            points.append((lat, lon))
                        except Exception:
                            continue
            elif ext in [".geojson", ".json"]:
                with open(file_path, "r", encoding="utf-8") as f:
                    gj = json.load(f)
                feats = gj.get("features", []) if gj.get("type") == "FeatureCollection" else []
                for feat in feats:
                    geom = feat.get("geometry", {})
                    if geom.get("type") == "LineString":
                        for coord in geom.get("coordinates", []):
                            if len(coord) >= 2:
                                points.append((coord[1], coord[0]))
            if len(points) < 2:
                self._show_message("warning", "Import Warning", "Imported backscatter line must contain at least 2 points.")
                return

            # Backscatter exports write "<export_name>_params.json", but imported geometry might be
            # "<export_name>_DDD.csv" / "_DMM.csv" / etc. Strip common suffixes so we can still find metadata.
            import_base_name = os.path.splitext(os.path.basename(file_path))[0]
            base_name_candidates = [import_base_name]
            for suffix in ("_DDD", "_DMM", "_DMS", "_DD", "_DM"):
                if import_base_name.endswith(suffix):
                    base_name_candidates.append(import_base_name[: -len(suffix)])

            params = {}
            for candidate in base_name_candidates:
                params_path = os.path.join(import_dir, f"{candidate}_params.json")
                if not os.path.exists(params_path):
                    continue
                try:
                    with open(params_path, "r", encoding="utf-8") as pf:
                        loaded = json.load(pf)
                    if isinstance(loaded, dict):
                        params = loaded
                    break
                except Exception:
                    params = {}

            if "backscatter_centerline" in params and "backscatter_half_width_m" in params and "backscatter_vertices" in params:
                cl = params.get("backscatter_centerline") or []
                if len(cl) == 2:
                    self.backscatter_box_centerline = ((float(cl[0][0]), float(cl[0][1])), (float(cl[1][0]), float(cl[1][1])))
                self.backscatter_box_half_width_m = float(params.get("backscatter_half_width_m", 0.0))
                verts = params.get("backscatter_vertices") or []
                if len(verts) == 4:
                    self.backscatter_box_vertices = [(float(v[0]), float(v[1])) for v in verts]
                wp = params.get("backscatter_width_point")
                if wp and len(wp) == 2:
                    self.backscatter_box_width_point = (float(wp[0]), float(wp[1]))
                self.backscatter_lead_in_m = float(params.get("backscatter_lead_in_m", self.backscatter_lead_in_m))
                self.backscatter_survey_speed_kn = float(params.get("backscatter_survey_speed_kn", self.backscatter_survey_speed_kn))
                if hasattr(self, "backscatter_lead_in_entry"):
                    self.backscatter_lead_in_entry.setText(f"{self.backscatter_lead_in_m:g}")
                if hasattr(self, "backscatter_survey_speed_entry"):
                    self.backscatter_survey_speed_entry.setText(f"{self.backscatter_survey_speed_kn:g}")
            else:
                self.backscatter_box_centerline = (points[0], points[-1])
                self.backscatter_box_half_width_m = 0.0
                self.backscatter_box_vertices = None

            # Restore saved NaN sentinel and backscatter criteria (when present).
            if params.get("geotiff_nan_value") is not None and hasattr(self, "_set_geotiff_nan_cutoff"):
                self._set_geotiff_nan_cutoff(params.get("geotiff_nan_value"), update_entry=True)

            try:
                if "backscatter_slope_min_deg" in params:
                    self.backscatter_slope_min_deg = float(params.get("backscatter_slope_min_deg"))
                if "backscatter_slope_max_deg" in params:
                    self.backscatter_slope_max_deg = float(params.get("backscatter_slope_max_deg"))
                if hasattr(self, "backscatter_slope_min_spin"):
                    self.backscatter_slope_min_spin.blockSignals(True)
                    self.backscatter_slope_min_spin.setValue(float(self.backscatter_slope_min_deg))
                    self.backscatter_slope_min_spin.blockSignals(False)
                if hasattr(self, "backscatter_slope_max_spin"):
                    self.backscatter_slope_max_spin.blockSignals(True)
                    self.backscatter_slope_max_spin.setValue(float(self.backscatter_slope_max_deg))
                    self.backscatter_slope_max_spin.blockSignals(False)
            except (TypeError, ValueError):
                pass

            if "backscatter_slope_areas_show_var" in params:
                self.backscatter_slope_areas_show_var = bool(params.get("backscatter_slope_areas_show_var"))
                if hasattr(self, "backscatter_show_slope_areas_checkbox"):
                    self.backscatter_show_slope_areas_checkbox.blockSignals(True)
                    self.backscatter_show_slope_areas_checkbox.setChecked(self.backscatter_slope_areas_show_var)
                    self.backscatter_show_slope_areas_checkbox.blockSignals(False)

            if "backscatter_depth_filter_enabled_var" in params:
                self.backscatter_depth_filter_enabled_var = bool(params.get("backscatter_depth_filter_enabled_var"))
            if params.get("backscatter_depth_min_m") is not None:
                try:
                    self.backscatter_depth_min_m = float(params.get("backscatter_depth_min_m"))
                except (TypeError, ValueError):
                    pass
            if params.get("backscatter_depth_max_m") is not None:
                try:
                    self.backscatter_depth_max_m = float(params.get("backscatter_depth_max_m"))
                except (TypeError, ValueError):
                    pass
            if hasattr(self, "backscatter_depth_filter_checkbox"):
                self.backscatter_depth_filter_checkbox.blockSignals(True)
                self.backscatter_depth_filter_checkbox.setChecked(bool(self.backscatter_depth_filter_enabled_var))
                self.backscatter_depth_filter_checkbox.blockSignals(False)
            if hasattr(self, "backscatter_depth_min_m_entry"):
                self.backscatter_depth_min_m_entry.setText(f"{self.backscatter_depth_min_m:g}")
            if hasattr(self, "backscatter_depth_max_m_entry"):
                self.backscatter_depth_max_m_entry.setText(f"{self.backscatter_depth_max_m:g}")

            if "backscatter_min_area_enabled_var" in params:
                self.backscatter_min_area_enabled_var = bool(params.get("backscatter_min_area_enabled_var"))
            if params.get("backscatter_min_area_m2") is not None:
                try:
                    self.backscatter_min_area_m2 = float(params.get("backscatter_min_area_m2"))
                except (TypeError, ValueError):
                    pass
            if hasattr(self, "backscatter_min_area_checkbox"):
                self.backscatter_min_area_checkbox.blockSignals(True)
                self.backscatter_min_area_checkbox.setChecked(bool(self.backscatter_min_area_enabled_var))
                self.backscatter_min_area_checkbox.blockSignals(False)
            if hasattr(self, "backscatter_min_area_m2_entry"):
                self.backscatter_min_area_m2_entry.setText(f"{self.backscatter_min_area_m2:g}")

            if "backscatter_extent_filter_enabled_var" in params:
                self.backscatter_extent_filter_enabled_var = bool(params.get("backscatter_extent_filter_enabled_var"))
            if params.get("backscatter_min_width_m") is not None:
                try:
                    self.backscatter_min_width_m = float(params.get("backscatter_min_width_m"))
                except (TypeError, ValueError):
                    pass
            if params.get("backscatter_min_height_m") is not None:
                try:
                    self.backscatter_min_height_m = float(params.get("backscatter_min_height_m"))
                except (TypeError, ValueError):
                    pass
            if hasattr(self, "backscatter_extent_checkbox"):
                self.backscatter_extent_checkbox.blockSignals(True)
                self.backscatter_extent_checkbox.setChecked(bool(self.backscatter_extent_filter_enabled_var))
                self.backscatter_extent_checkbox.blockSignals(False)
            if hasattr(self, "backscatter_min_width_m_entry"):
                self.backscatter_min_width_m_entry.setText(f"{self.backscatter_min_width_m:g}")
            if hasattr(self, "backscatter_min_height_m_entry"):
                self.backscatter_min_height_m_entry.setText(f"{self.backscatter_min_height_m:g}")

            if "backscatter_percent_clip_enabled" in params:
                self.backscatter_percent_clip_enabled = bool(params.get("backscatter_percent_clip_enabled"))
            if params.get("backscatter_percent_clip_min") is not None:
                try:
                    self.backscatter_percent_clip_min = float(params.get("backscatter_percent_clip_min"))
                except (TypeError, ValueError):
                    pass
            if params.get("backscatter_percent_clip_max") is not None:
                try:
                    self.backscatter_percent_clip_max = float(params.get("backscatter_percent_clip_max"))
                except (TypeError, ValueError):
                    pass
            if hasattr(self, "backscatter_percent_clip_checkbox"):
                self.backscatter_percent_clip_checkbox.blockSignals(True)
                self.backscatter_percent_clip_checkbox.setChecked(bool(self.backscatter_percent_clip_enabled))
                self.backscatter_percent_clip_checkbox.blockSignals(False)
            if hasattr(self, "backscatter_percent_clip_min_entry"):
                self.backscatter_percent_clip_min_entry.setText(f"{self.backscatter_percent_clip_min:g}")
            if hasattr(self, "backscatter_percent_clip_max_entry"):
                self.backscatter_percent_clip_max_entry.setText(f"{self.backscatter_percent_clip_max:g}")

            if "show_backscatter_var" in params:
                self.show_backscatter_var = bool(params.get("show_backscatter_var"))
                if hasattr(self, "show_backscatter_checkbox"):
                    self.show_backscatter_checkbox.blockSignals(True)
                    self.show_backscatter_checkbox.setChecked(self.show_backscatter_var)
                    self.show_backscatter_checkbox.blockSignals(False)

            # Restore bathymetry/backscatter rasters when paths were saved.
            geotiff_path = params.get("geotiff_path")
            if geotiff_path and os.path.isfile(geotiff_path):
                try:
                    self._load_geotiff_from_path(geotiff_path, use_background_loading=False)
                except Exception:
                    pass
            bs_geotiff_path = params.get("backscatter_geotiff_path")
            if bs_geotiff_path and os.path.isfile(bs_geotiff_path):
                try:
                    self._load_backscatter_geotiff_from_path(bs_geotiff_path)
                except Exception:
                    pass

            self.backscatter_box_stats = self._compute_backscatter_box_stats()
            self._plot_survey_plan(preserve_view_limits=True)
            # Fit view to imported geometry so imported line/area is immediately visible.
            try:
                lon_vals = []
                lat_vals = []
                if self.backscatter_box_centerline and len(self.backscatter_box_centerline) == 2:
                    for lat, lon in self.backscatter_box_centerline:
                        if lat is not None and lon is not None:
                            lat_vals.append(float(lat))
                            lon_vals.append(float(lon))
                if self.backscatter_box_vertices and len(self.backscatter_box_vertices) == 4:
                    for lon, lat in self.backscatter_box_vertices:
                        lat_vals.append(float(lat))
                        lon_vals.append(float(lon))
                if lon_vals and lat_vals:
                    lon_min, lon_max = min(lon_vals), max(lon_vals)
                    lat_min, lat_max = min(lat_vals), max(lat_vals)
                    lon_pad = max((lon_max - lon_min) * 0.2, 0.001)
                    lat_pad = max((lat_max - lat_min) * 0.2, 0.001)
                    self.ax.set_xlim(lon_min - lon_pad, lon_max + lon_pad)
                    self.ax.set_ylim(lat_min - lat_pad, lat_max + lat_pad)
                    self.canvas.draw_idle()
            except Exception:
                pass
            if hasattr(self, "_draw_current_profile"):
                self._draw_current_profile()
            self._update_backscatter_export_name_default()
            if getattr(self, "backscatter_download_gmrt_checkbox", None) and self.backscatter_download_gmrt_checkbox.isChecked():
                self._download_and_load_gmrt_after_backscatter_import()
        except Exception as e:
            self._show_message("error", "Import Error", f"Failed to import backscatter line: {e}")

    def _export_backscatter_line(self):
        """Export backscatter line/area using Export Types and include backscatter statistics image."""
        centerline = getattr(self, "backscatter_box_centerline", None)
        if not centerline or centerline[0] is None or centerline[1] is None:
            self._show_message("warning", "No Backscatter Line", "Select an area/line first.")
            return
        export_name = "BS_export"
        if hasattr(self, "backscatter_export_name_entry"):
            export_name = self.backscatter_export_name_entry.text().strip() or export_name
            self.backscatter_export_name_entry.setText(export_name)
        export_dir = QFileDialog.getExistingDirectory(
            self,
            "Select Export Directory",
            getattr(self, "last_backscatter_export_dir", None) or self.last_export_dir,
        )
        if not export_dir:
            return
        self.last_backscatter_export_dir = export_dir
        if hasattr(self, "_save_last_backscatter_export_dir"):
            self._save_last_backscatter_export_dir()

        try:
            line_points = list(getattr(self, "_calculate_backscatter_line_statistics", lambda: {})().get("line_waypoints", []))
            if len(line_points) < 2:
                line_points = [centerline[0], centerline[1]]
            rows = [(1, "BackscatterLine", f"P{i+1}", lat, lon) for i, (lat, lon) in enumerate(line_points)]

            export_shapefile = self._export_type_enabled("esri_shapefile") if hasattr(self, "_export_type_enabled") else True
            export_sis = self._export_type_enabled("sis_asciiplan") if hasattr(self, "_export_type_enabled") else True
            export_gpx = self._export_type_enabled("gpx") if hasattr(self, "_export_type_enabled") else True
            export_text_csv = self._export_type_enabled("text_csv") if hasattr(self, "_export_type_enabled") else True
            export_text_txt = self._export_type_enabled("text_txt") if hasattr(self, "_export_type_enabled") else True
            export_hypack = self._export_type_enabled("hypack_lnw") if hasattr(self, "_export_type_enabled") else True
            export_map_png = self._export_type_enabled("map_png") if hasattr(self, "_export_type_enabled") else True
            export_profiles_png = self._export_type_enabled("profiles_png") if hasattr(self, "_export_type_enabled") else True

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

            geojson_file_path = os.path.join(export_dir, f"{export_name}.geojson")
            gj_features = [{
                "type": "Feature",
                "geometry": {"type": "LineString", "coordinates": [[lon, lat] for lat, lon in line_points]},
                "properties": {
                    "line_num": 1,
                    "survey_speed": float(getattr(self, "backscatter_survey_speed_kn", 8.0)),
                },
            }]
            geojson_collection = {
                "type": "FeatureCollection",
                "properties": {"geotiff_nan_value": float(getattr(self, "geotiff_nan_value", -11000.0))},
                "features": gj_features,
            }
            with open(geojson_file_path, "w", encoding="utf-8") as f:
                json.dump(geojson_collection, f, indent=2)

            lnw_file_path = None
            if export_hypack:
                lnw_lines = [("BACKSCATTER", [line_points[0], line_points[-1]])]
                zone, hem = export_utils.compute_utm_zone_from_points([line_points[0], line_points[-1]])
                lnw_file_path = os.path.join(export_dir, f"{export_name}_UTM{zone}{'N' if hem == 'North' else 'S'}.lnw")
                if not export_utils.write_lnw(lnw_file_path, lnw_lines):
                    lnw_file_path = None

            sis_file_path = os.path.join(export_dir, f"{export_name}.asciiplan")
            if export_sis:
                export_utils.write_asciiplan(sis_file_path, [("BackscatterLine", [line_points[0], line_points[-1]])])

            gpx_file_path = os.path.join(export_dir, f"{export_name}.gpx")
            if export_gpx:
                export_utils.write_gpx(gpx_file_path, [("BackscatterLine", line_points)])

            shapefile_path = os.path.join(export_dir, f"{export_name}.shp")
            if export_shapefile and LineString is not None:
                try:
                    from shapely.geometry import mapping
                    import fiona
                    schema = {"geometry": "LineString", "properties": {"line_num": "int", "line_name": "str"}}
                    with fiona.open(shapefile_path, "w", driver="ESRI Shapefile", crs="EPSG:4326", schema=schema) as collection:
                        collection.writerecord({
                            "geometry": mapping(LineString([(lon, lat) for lat, lon in line_points])),
                            "properties": {"line_num": 1, "line_name": "BackscatterLine"},
                        })
                except Exception:
                    shapefile_path = None

            map_png_path = os.path.join(export_dir, f"{export_name}_map.png")
            if export_map_png and hasattr(self, "figure"):
                self.figure.savefig(map_png_path, dpi=300)
            profile_png_path = os.path.join(export_dir, f"{export_name}_profiles.png")
            if export_profiles_png and hasattr(self, "profile_fig"):
                self.profile_fig.savefig(profile_png_path, dpi=300)

            if getattr(self, "backscatter_box_stats", None) is None:
                self.backscatter_box_stats = self._compute_backscatter_box_stats()
            self._show_backscatter_stats_dialog()
            stats_png_path = os.path.join(export_dir, f"{export_name}_backscatter_stats.png")
            if hasattr(self, "backscatter_stats_canvas") and self.backscatter_stats_canvas is not None:
                self._save_backscatter_stats_png_without_hover(stats_png_path, dpi=300)

            params_json_path = os.path.join(export_dir, f"{export_name}_params.json")
            payload = {
                "backscatter_centerline": list(getattr(self, "backscatter_box_centerline", []) or []),
                "backscatter_half_width_m": float(getattr(self, "backscatter_box_half_width_m", 0.0) or 0.0),
                "backscatter_vertices": list(getattr(self, "backscatter_box_vertices", []) or []),
                "backscatter_width_point": list(getattr(self, "backscatter_box_width_point", []) or []) if getattr(self, "backscatter_box_width_point", None) else None,
                "backscatter_lead_in_m": float(getattr(self, "backscatter_lead_in_m", 0.0) or 0.0),
                "backscatter_survey_speed_kn": float(getattr(self, "backscatter_survey_speed_kn", 8.0) or 8.0),
                "backscatter_slope_areas_show_var": bool(getattr(self, "backscatter_slope_areas_show_var", False)),
                "backscatter_slope_min_deg": float(getattr(self, "backscatter_slope_min_deg", 0.0) or 0.0),
                "backscatter_slope_max_deg": float(getattr(self, "backscatter_slope_max_deg", 2.0) or 2.0),
                "backscatter_depth_filter_enabled_var": bool(getattr(self, "backscatter_depth_filter_enabled_var", False)),
                "backscatter_depth_min_m": float(getattr(self, "backscatter_depth_min_m", 0.0) or 0.0),
                "backscatter_depth_max_m": float(getattr(self, "backscatter_depth_max_m", 11000.0) or 11000.0),
                "backscatter_min_area_enabled_var": bool(getattr(self, "backscatter_min_area_enabled_var", False)),
                "backscatter_min_area_m2": float(getattr(self, "backscatter_min_area_m2", 0.0) or 0.0),
                "backscatter_extent_filter_enabled_var": bool(getattr(self, "backscatter_extent_filter_enabled_var", False)),
                "backscatter_min_width_m": float(getattr(self, "backscatter_min_width_m", 0.0) or 0.0),
                "backscatter_min_height_m": float(getattr(self, "backscatter_min_height_m", 0.0) or 0.0),
                "show_backscatter_var": bool(getattr(self, "show_backscatter_var", False)),
                "backscatter_percent_clip_enabled": bool(getattr(self, "backscatter_percent_clip_enabled", True)),
                "backscatter_percent_clip_min": float(getattr(self, "backscatter_percent_clip_min", 0.5) or 0.5),
                "backscatter_percent_clip_max": float(getattr(self, "backscatter_percent_clip_max", 0.5) or 0.5),
                "geotiff_path": self.current_geotiff_path if hasattr(self, "current_geotiff_path") else None,
                "backscatter_geotiff_path": self.backscatter_geotiff_path if hasattr(self, "backscatter_geotiff_path") else None,
                "geotiff_nan_value": float(getattr(self, "geotiff_nan_value", -11000.0)),
            }
            with open(params_json_path, "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2)

            if hasattr(self, "set_line_info_text"):
                self.set_line_info_text(f"Backscatter line exported to {export_dir}", append=False)
        except Exception as e:
            self._show_message("error", "Export Error", f"Failed to export backscatter line: {e}")

    def _set_geotiff_nan_cutoff(self, value, update_entry=True):
        """Persist NaN cutoff and optionally synchronize the UI entry."""
        try:
            cutoff = float(value)
        except (TypeError, ValueError):
            return False
        self.geotiff_nan_value = cutoff
        if update_entry and hasattr(self, "geotiff_nan_entry"):
            self.geotiff_nan_entry.blockSignals(True)
            self.geotiff_nan_entry.setText(f"{cutoff:g}")
            self.geotiff_nan_entry.blockSignals(False)
        return True

    def _apply_geotiff_nan_filter(self, data_array):
        """Apply configured NaN sentinel filtering in place and return the array."""
        if data_array is None:
            return data_array
        cutoff = self._get_geotiff_nan_cutoff()
        # Treat only the configured sentinel value as NaN (not a threshold).
        data_array[np.isclose(data_array, cutoff, rtol=0.0, atol=1e-9)] = np.nan
        return data_array

    def _on_geotiff_nan_value_changed(self):
        """Debounce NaN threshold edits before applying to the loaded raster."""
        if not hasattr(self, "_geotiff_nan_update_timer"):
            self._geotiff_nan_update_timer = QTimer()
            self._geotiff_nan_update_timer.setSingleShot(True)
            self._geotiff_nan_update_timer.timeout.connect(self._apply_geotiff_nan_value_changed)
        self._geotiff_nan_update_timer.stop()
        self._geotiff_nan_update_timer.start(450)

    def _apply_geotiff_nan_value_changed(self):
        """Apply user-entered NaN threshold and refresh loaded GeoTIFF if needed."""
        if not hasattr(self, "geotiff_nan_entry"):
            return
        raw = self.geotiff_nan_entry.text().strip()
        if not raw:
            return
        try:
            cutoff = float(raw)
        except ValueError:
            return
        if not self._set_geotiff_nan_cutoff(cutoff, update_entry=True):
            return
        if self.geotiff_dataset_original is not None:
            self._reload_geotiff_at_current_zoom()

    def _cancel_geotiff_loading(self):
        """Cancel the current GeoTIFF loading operation."""
        self.loading_cancelled = True
        if hasattr(self, 'set_cal_info_text'):
            self.set_cal_info_text("GeoTIFF loading cancelled.")
        elif hasattr(self, 'set_ref_info_text'):
            self.set_ref_info_text("GeoTIFF loading cancelled.")

    def _load_geotiff_at_resolution(self, zoom_level=None):
        """Load GeoTIFF data at the specified resolution based on zoom level."""
        if not GEOSPATIAL_LIBS_AVAILABLE or self.geotiff_dataset_original is None:
            return False

        # If dynamic resolution is disabled, always use full resolution
        if not self.dynamic_resolution_enabled:
            downsample_factor = 1
            print("DEBUG: Dynamic resolution disabled, using full resolution")
        else:
            if zoom_level is None:
                zoom_level = self.geotiff_zoom_level

            # Clamp zoom level to prevent extreme values that could cause issues
            # Allow zoom levels up to 10.0 for zoom out operations
            zoom_level = max(0.01, min(10.0, zoom_level))

            # Determine downsampling factor based on zoom level
            # For zoom out (zoom_level > 1.0), we want higher downsampling
            if zoom_level <= 0.5:  # Zoomed in close - use full resolution
                downsample_factor = 1
            elif zoom_level <= 1.0:  # Medium zoom - use 2x downsampling
                downsample_factor = 2
            elif zoom_level <= 2.0:  # Far zoom - use 4x downsampling
                downsample_factor = 4
            elif zoom_level <= 5.0:  # Very far zoom - use 8x downsampling
                downsample_factor = 8
            else:  # Extremely far zoom - use 16x downsampling
                downsample_factor = 16

        try:
            # Use current axes limits (call canvas.draw() before this if limits just changed)
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()

            # Get the full extent of the original GeoTIFF
            if self.geotiff_dataset_original.crs != "EPSG:4326":
                # If not WGS84, we need to transform the bounds
                bounds = self.geotiff_dataset_original.bounds
                transformer = pyproj.Transformer.from_crs(
                    self.geotiff_dataset_original.crs, "EPSG:4326", always_xy=True
                )
                left, bottom = transformer.transform(bounds.left, bounds.bottom)
                right, top = transformer.transform(bounds.right, bounds.top)
                full_extent = [left, right, bottom, top]
            else:
                bounds = self.geotiff_dataset_original.bounds
                full_extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]

            # Calculate the visible region in the original dataset
            left, right, bottom, top = full_extent

            # Add a buffer to the view extent to ensure full coverage (10% of view size)
            view_width = xlim[1] - xlim[0]
            view_height = ylim[1] - ylim[0]
            buffer_x = view_width * 0.1
            buffer_y = view_height * 0.1

            # Expand view limits with buffer
            expanded_xlim = (xlim[0] - buffer_x, xlim[1] + buffer_x)
            expanded_ylim = (ylim[0] - buffer_y, ylim[1] + buffer_y)

            # Convert expanded plot limits to dataset coordinates
            if self.geotiff_dataset_original.crs != "EPSG:4326":
                transformer = pyproj.Transformer.from_crs("EPSG:4326", self.geotiff_dataset_original.crs, always_xy=True)
                plot_left, plot_bottom = transformer.transform(expanded_xlim[0], expanded_ylim[0])
                plot_right, plot_top = transformer.transform(expanded_xlim[1], expanded_ylim[1])
            else:
                plot_left, plot_bottom = expanded_xlim[0], expanded_ylim[0]
                plot_right, plot_top = expanded_xlim[1], expanded_ylim[1]

            # Calculate pixel indices for the visible region
            transform_obj = self.geotiff_dataset_original.transform
            row_min, col_min = rowcol(transform_obj, plot_left, plot_top)
            row_max, col_max = rowcol(transform_obj, plot_right, plot_bottom)

            # Ensure bounds are within dataset and handle edge cases
            row_min = max(0, int(row_min))
            row_max = min(self.geotiff_dataset_original.height, int(row_max))
            col_min = max(0, int(col_min))
            col_max = min(self.geotiff_dataset_original.width, int(col_max))

            # Ensure we have a valid region (at least 1 pixel)
            if row_max <= row_min:
                row_max = min(row_min + 1, self.geotiff_dataset_original.height)
            if col_max <= col_min:
                col_max = min(col_min + 1, self.geotiff_dataset_original.width)

            # Add additional padding around the visible region to ensure full coverage
            # Use larger padding to ensure full coverage of the view
            padding = max(100, min(200, (row_max - row_min) // 2))  # Increased adaptive padding
            row_min = max(0, row_min - padding)
            row_max = min(self.geotiff_dataset_original.height, row_max + padding)
            col_min = max(0, col_min - padding)
            col_max = min(self.geotiff_dataset_original.width, col_max + padding)

            # Final validation - ensure we have a reasonable region size
            if (row_max - row_min) < 10 or (col_max - col_min) < 10:
                print(f"Warning: Region too small ({row_max-row_min}x{col_max-col_min}), using full dataset")
                row_min, row_max = 0, self.geotiff_dataset_original.height
                col_min, col_max = 0, self.geotiff_dataset_original.width

            # Read the visible region
            window = Window.from_slices((row_min, row_max), (col_min, col_max))
            data = self.geotiff_dataset_original.read(1, window=window)

            # Apply downsampling if needed
            if downsample_factor > 1:
                data = data[::downsample_factor, ::downsample_factor]

            # Calculate the extent of the loaded region
            # Get the bounds of the window in the original CRS
            bounds = window_bounds(window, self.geotiff_dataset_original.transform)
            left, bottom, right, top = bounds

            if self.geotiff_dataset_original.crs != "EPSG:4326":
                # Transform the bounds to WGS84
                transformer = pyproj.Transformer.from_crs(
                    self.geotiff_dataset_original.crs, "EPSG:4326", always_xy=True
                )
                left, bottom = transformer.transform(left, bottom)
                right, top = transformer.transform(right, top)
                region_extent = [left, right, bottom, top]
            else:
                # Already in WGS84
                region_extent = [left, right, bottom, top]

            # Filter Z-values (user-configurable low cutoff and non-negative values)
            self._apply_geotiff_nan_filter(data)

            # Validate the extent to prevent flipping
            if region_extent[1] <= region_extent[0] or region_extent[3] <= region_extent[2]:
                print(f"Warning: Invalid extent calculated: {region_extent}. Skipping resolution update.")
                return False

            # Update the current data and extent
            self.geotiff_data_array = data.astype(float)
            self.geotiff_extent = region_extent
            self.geotiff_current_resolution = downsample_factor

            if hasattr(self, "_rewarp_backscatter_align_to_bathy"):
                self._rewarp_backscatter_align_to_bathy()

            return True

        except Exception as e:
            print(f"Error loading GeoTIFF at resolution: {e}")
            return False

    def _reload_geotiff_at_current_zoom(self):
        """Reload GeoTIFF data at the current zoom level."""
        if self.geotiff_dataset_original is not None:
            # Use the last user-set view limits (from zoom/pan operations) instead of current limits
            # This ensures we preserve the user's intended zoom level
            if hasattr(self, '_last_user_xlim') and hasattr(self, '_last_user_ylim'):
                current_xlim = self._last_user_xlim
                current_ylim = self._last_user_ylim
            else:
                # Fallback to current limits if not set
                current_xlim = self.ax.get_xlim()
                current_ylim = self.ax.get_ylim()

            # Store the current view limits for restoration
            self._current_view_limits = (current_xlim, current_ylim)

            success = self._load_geotiff_at_resolution()
            if success:
                # Re-plot with new resolution (preserve_view_limits=True will use _last_user_xlim/_last_user_ylim)
                self._plot_survey_plan(preserve_view_limits=True)
                if hasattr(self, "_on_geotiff_loaded_line_info_depth"):
                    self._on_geotiff_loaded_line_info_depth()
                if hasattr(self, "_update_roll_line_info_labels"):
                    self._update_roll_line_info_labels()

                # Ensure limits are set correctly (in case _plot_survey_plan didn't preserve them due to aspect ratio)
                # Use the stored limits to maintain user's zoom level
                self.ax.set_xlim(current_xlim)
                self.ax.set_ylim(current_ylim)

                self.canvas.draw_idle()

                # Update info text to show current resolution
                if not self.dynamic_resolution_enabled:
                    resolution_text = "GeoTIFF resolution: Full resolution (dynamic resolution disabled)"
                else:
                    resolution_text = f"GeoTIFF resolution: {self.geotiff_current_resolution}x downsampling"

                if hasattr(self, 'set_cal_info_text'):
                    self.set_cal_info_text(resolution_text)
                elif hasattr(self, 'set_ref_info_text'):
                    self.set_ref_info_text(resolution_text)
            else:
                # Fallback: keep current view and show warning
                if hasattr(self, 'set_cal_info_text'):
                    self.set_cal_info_text("Warning: Could not update GeoTIFF resolution. Using current view.")
                elif hasattr(self, 'set_ref_info_text'):
                    self.set_ref_info_text("Warning: Could not update GeoTIFF resolution. Using current view.")

    def _apply_map_zoom_limits_and_reload_geotiff(self, xlim, ylim):
        """Set map limits, sync zoom/user state, and reload GeoTIFF for the view.

        Programmatic zoom (toolbar buttons) must update ``geotiff_zoom_level`` and
        redraw before ``_load_geotiff_at_resolution`` runs; otherwise dynamic
        resolution stays on the old (e.g. zoomed-out) downsampling and the
        raster looks blocky until a manual zoom refreshes it.
        """
        if hasattr(self, "_zoom_timer"):
            self._zoom_timer.stop()
        if hasattr(self, "_limit_change_timer"):
            self._limit_change_timer.stop()

        self.current_xlim = xlim
        self.current_ylim = ylim
        self._last_user_xlim = xlim
        self._last_user_ylim = ylim
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)

        ref_extent = self._dynamic_resolution_reference_extent()
        if ref_extent is not None and self.geotiff_dataset_original is not None:
            if getattr(self, "dynamic_resolution_enabled", True):
                full_width = ref_extent[1] - ref_extent[0]
                full_height = ref_extent[3] - ref_extent[2]
                current_width = xlim[1] - xlim[0]
                current_height = ylim[1] - ylim[0]
                width_ratio = current_width / full_width if full_width > 0 else 1.0
                height_ratio = current_height / full_height if full_height > 0 else 1.0
                new_zoom_level = max(width_ratio, height_ratio)
                new_zoom_level = max(0.01, min(10.0, new_zoom_level))
                self.geotiff_zoom_level = new_zoom_level

            self.canvas.draw()
            self.current_xlim = self.ax.get_xlim()
            self.current_ylim = self.ax.get_ylim()
            self._last_user_xlim = self.current_xlim
            self._last_user_ylim = self.current_ylim

            success = self._load_geotiff_at_resolution()
            if success:
                self._plot_survey_plan(preserve_view_limits=True)
                if hasattr(self, "_on_geotiff_loaded_line_info_depth"):
                    self._on_geotiff_loaded_line_info_depth()
                if hasattr(self, "_update_roll_line_info_labels"):
                    self._update_roll_line_info_labels()
            else:
                self.canvas.draw_idle()
        else:
            self.canvas.draw_idle()

    def _clear_rapid_zoom_mode(self):
        """Clear the rapid zoom mode flag."""
        if hasattr(self, '_rapid_zoom_mode'):
            del self._rapid_zoom_mode
        # Reset cumulative zoom change after rapid zooming stops
        if hasattr(self, '_cumulative_zoom_change'):
            self._cumulative_zoom_change = 0.0
        if hasattr(self, '_zoom_operation_count'):
            self._zoom_operation_count = 0

        # Trigger a resolution update after rapid zooming stops to ensure proper resolution
        if (self.geotiff_dataset_original is not None and
            hasattr(self, 'geotiff_zoom_level')):
            # Use a short delay to ensure the view has settled
            if not hasattr(self, '_zoom_timer'):
                self._zoom_timer = QTimer()
                self._zoom_timer.setSingleShot(True)
                self._zoom_timer.timeout.connect(self._reload_geotiff_at_current_zoom)
            self._zoom_timer.stop()
            self._zoom_timer.start(500)

    def _clear_panning_mode(self):
        """Clear the panning mode flag."""
        if hasattr(self, '_panning_mode'):
            del self._panning_mode

    def _toggle_dynamic_resolution(self):
        """Toggle dynamic resolution feature on/off."""
        self.dynamic_resolution_enabled = not self.dynamic_resolution_enabled
        status = "ON" if self.dynamic_resolution_enabled else "OFF"

        # Update dynamic resolution button
        if hasattr(self, 'dynamic_resolution_btn'):
            self.dynamic_resolution_btn.setText(f"Dynamic Resolution: {status}")

        # Update info text
        if hasattr(self, 'set_cal_info_text'):
            self.set_cal_info_text(f"Dynamic resolution: {status}")
        elif hasattr(self, 'set_ref_info_text'):
            self.set_ref_info_text(f"Dynamic resolution: {status}")

        # Reload GeoTIFF with new resolution settings if a GeoTIFF is loaded
        if self.geotiff_dataset_original is not None:
            print(f"DEBUG: Dynamic resolution toggled to {status}, reloading GeoTIFF")
            # Use a short delay to ensure the toggle is complete
            if hasattr(self, '_zoom_timer'):
                self._zoom_timer.stop()
            self._zoom_timer = QTimer()
            self._zoom_timer.setSingleShot(True)
            self._zoom_timer.timeout.connect(self._reload_geotiff_at_current_zoom)
            self._zoom_timer.start(100)
        else:
            print(f"DEBUG: Dynamic resolution toggled to {status}, but no GeoTIFF loaded")

    def _reset_to_consistent_view(self):
        """Reset the plot view to the consistent limits that maintain window size."""
        if self.geotiff_extent is None:
            # No GeoTIFF loaded, use global limits
            self.ax.set_xlim(self.fixed_xlim)
            self.ax.set_ylim(self.fixed_ylim)
        else:
            # Use the calculated consistent limits
            self.ax.set_xlim(self.current_xlim)
            self.ax.set_ylim(self.current_ylim)
        self.canvas.draw_idle()

    def _load_geotiff(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded. Cannot load GeoTIFF.")
            return

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select GeoTIFF File",
            self.last_geotiff_dir,
            "GeoTIFF files (*.tif *.tiff);;All files (*.*)"
        )
        if not file_path:
            return
        self.last_geotiff_dir = os.path.dirname(file_path)
        self._save_last_geotiff_dir()

        # Check file size before loading
        file_size_mb = os.path.getsize(file_path) / (1024 * 1024)
        use_background_loading = False

        if file_size_mb > 1000:  # 1GB threshold - use background loading
            response = self._ask_yes_no(
                "Very Large File Warning",
                f"GeoTIFF file is {file_size_mb:.1f} MB. This is a very large file that may take several minutes to load.\n\n"
                "Would you like to load it in the background? (Recommended)\n"
                "Or cancel and use a smaller file for better performance."
            )
            if not response:
                return
            use_background_loading = True
        elif file_size_mb > 500:  # 500MB threshold
            response = self._ask_yes_no(
                "Large File Warning",
                f"GeoTIFF file is {file_size_mb:.1f} MB. Large files may take a long time to load and could cause the application to freeze.\n\n"
                "Would you like to continue? Consider using a smaller or downsampled file for better performance."
            )
            if not response:
                return

        self._load_geotiff_from_path(file_path, use_background_loading=use_background_loading)

    def _load_geotiff_from_path(self, file_path, use_background_loading=False):
        """Load and display a GeoTIFF from the given path. Used by Load GeoTIFF and by GMRT download."""
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded. Cannot load GeoTIFF.")
            return
        if not file_path or not os.path.isfile(file_path):
            return
        self.last_geotiff_dir = os.path.dirname(file_path)
        if hasattr(self, '_save_last_geotiff_dir'):
            self._save_last_geotiff_dir()

        file_size_mb = os.path.getsize(file_path) / (1024 * 1024)

        # Show loading progress
        progress_window = QDialog(self)
        progress_window.setWindowTitle("Loading GeoTIFF")
        progress_window.setModal(True)
        progress_layout = QVBoxLayout(progress_window)

        progress_label = QLabel("Opening GeoTIFF file...")
        progress_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        progress_layout.addWidget(progress_label)

        progress_bar = QProgressBar()
        progress_bar.setRange(0, 0)  # Indeterminate mode
        progress_layout.addWidget(progress_bar)

        # Add cancel button for background loading
        if use_background_loading:
            cancel_button = QPushButton("Cancel Loading")
            cancel_button.clicked.connect(lambda: setattr(self, 'loading_cancelled', True))
            cancel_button.clicked.connect(progress_window.reject)
            progress_layout.addWidget(cancel_button)
            self.loading_cancelled = False

        progress_window.show()
        QApplication.processEvents()

        try:
            if self.geotiff_dataset_original:  # Close previously opened dataset
                self.geotiff_dataset_original.close()

            # Update progress
            progress_label.setText("Opening raster dataset...")
            QApplication.processEvents()

            # Add timeout for opening large files
            timeout_seconds = 30 if file_size_mb > 500 else 10

            # For very large files, just proceed with normal loading but show progress
            # The timeout will be handled by the progress window and user can cancel
            self.geotiff_dataset_original = rasterio.open(file_path)

            # Check if we need to downsample for large files
            total_pixels = self.geotiff_dataset_original.width * self.geotiff_dataset_original.height
            downsample_factor = 1

            if total_pixels > 10000000:  # 10 million pixels threshold
                downsample_factor = 4
                progress_label.setText(f"Large file detected. Downsampling by factor {downsample_factor} for performance...")
                QApplication.processEvents()
            elif total_pixels > 5000000:  # 5 million pixels threshold
                downsample_factor = 2
                progress_label.setText(f"Downsampling by factor {downsample_factor} for performance...")
                QApplication.processEvents()

            # Check if loading was cancelled
            if hasattr(self, 'loading_cancelled') and self.loading_cancelled:
                progress_window.close()
                return

            # Reproject GeoTIFF data to WGS84 (EPSG:4326) for consistent plotting
            src_crs = self.geotiff_dataset_original.crs
            dst_crs = "EPSG:4326"

            if src_crs != dst_crs:
                progress_label.setText(f"Reprojecting from {src_crs} to {dst_crs}...")
                QApplication.processEvents()

                # Calculate destination transform and dimensions
                left, bottom, right, top = self.geotiff_dataset_original.bounds.left, self.geotiff_dataset_original.bounds.bottom, \
                    self.geotiff_dataset_original.bounds.right, self.geotiff_dataset_original.bounds.top

                # Transform bounds to destination CRS
                transformer_bounds = pyproj.Transformer.from_crs(src_crs, dst_crs, always_xy=True)
                dst_left, dst_bottom = transformer_bounds.transform(left, bottom)
                dst_right, dst_top = transformer_bounds.transform(right, top)

                # Apply downsampling if needed
                target_height = self.geotiff_dataset_original.height // downsample_factor
                target_width = self.geotiff_dataset_original.width // downsample_factor

                # Create a new array for the reprojected data
                reprojected_data = np.empty((target_height, target_width),
                                           dtype=self.geotiff_dataset_original.dtypes[0])

                # Perform the reprojection with downsampling
                if downsample_factor > 1:
                    # Read with downsampling
                    window = Window.from_slices((0, self.geotiff_dataset_original.height), (0, self.geotiff_dataset_original.width))
                    data = self.geotiff_dataset_original.read(1, window=window)
                    # Simple downsampling by taking every nth pixel
                    data = data[::downsample_factor, ::downsample_factor]
                    source_data = data
                else:
                    source_data = rasterio.band(self.geotiff_dataset_original, 1)

                # Check if loading was cancelled before reprojection
                if hasattr(self, 'loading_cancelled') and self.loading_cancelled:
                    progress_window.close()
                    return

                reproject(
                    source=source_data,
                    destination=reprojected_data,
                    src_transform=self.geotiff_dataset_original.transform,
                    src_crs=src_crs,
                    dst_transform=transform.from_bounds(dst_left, dst_bottom, dst_right, dst_top, target_width, target_height),
                    dst_crs=dst_crs,
                    resampling=Resampling.bilinear
                )
                self.geotiff_data_array = reprojected_data.astype(float)
                self.geotiff_extent = [dst_left, dst_right, dst_bottom, dst_top]
                # Store original full extent for profile calculations (always use full extent for profiles)
                self.geotiff_original_extent = [dst_left, dst_right, dst_bottom, dst_top]

                # Create transformers for picking if reprojection happened
                self.wgs84_to_geotiff_transformer = pyproj.Transformer.from_crs(dst_crs, src_crs, always_xy=True)
                self.geotiff_to_wgs84_transformer = pyproj.Transformer.from_crs(src_crs, dst_crs, always_xy=True)

            else:  # Already WGS84
                progress_label.setText("Reading GeoTIFF data...")
                QApplication.processEvents()

                if downsample_factor > 1:
                    # Read with downsampling
                    window = Window.from_slices((0, self.geotiff_dataset_original.height), (0, self.geotiff_dataset_original.width))
                    data = self.geotiff_dataset_original.read(1, window=window)
                    # Simple downsampling
                    data = data[::downsample_factor, ::downsample_factor]
                    self.geotiff_data_array = data.astype(float)
                else:
                    self.geotiff_data_array = self.geotiff_dataset_original.read(1).astype(float)

                bounds = self.geotiff_dataset_original.bounds
                self.geotiff_extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]
                # Store original full extent for profile calculations (always use full extent for profiles)
                self.geotiff_original_extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]
                self.wgs84_to_geotiff_transformer = None
                self.geotiff_to_wgs84_transformer = None

            # Filter Z-values: user-configurable low cutoff and non-negative values
            progress_label.setText("Processing elevation data...")
            QApplication.processEvents()
            self._set_geotiff_nan_cutoff(self._detect_geotiff_nan_cutoff_from_dataset(), update_entry=True)
            self._apply_geotiff_nan_filter(self.geotiff_data_array)

            # Store the full resolution data for dynamic loading
            self.geotiff_full_resolution = self.geotiff_data_array.copy()
            self.geotiff_current_resolution = downsample_factor
            self.geotiff_zoom_level = 1.0  # Start at full resolution

            # Close progress window
            progress_window.close()

            # Remember path for suggested export names (e.g. after calibration import)
            if hasattr(self, 'current_geotiff_path'):
                self.current_geotiff_path = file_path

            # Enable Remove GeoTIFF, Pick Center, and Zoom to GeoTIFF buttons
            if hasattr(self, 'remove_geotiff_btn'):
                self.remove_geotiff_btn.setEnabled(True)
            if hasattr(self, 'pick_center_btn'):
                self.pick_center_btn.setEnabled(True)
            if hasattr(self, 'performance_pick_center_btn'):
                self.performance_pick_center_btn.setEnabled(True)
            if hasattr(self, 'zoom_to_geotiff_btn'):
                self.zoom_to_geotiff_btn.setEnabled(True)
                # Make "Pick Center from GeoTIFF" button bold and orange after loading GeoTIFF
                self.pick_center_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")
                if hasattr(self, 'performance_pick_center_btn'):
                    self.performance_pick_center_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")
            # Reset Load GeoTIFF button to normal style after successful loading
            if hasattr(self, 'load_geotiff_btn'):
                self.load_geotiff_btn.setStyleSheet("")  # Reset to default style
            # Make "Draw a Pitch Line" button bold and orange after loading GeoTIFF
            if hasattr(self, 'pick_pitch_line_btn'):
                self.pick_pitch_line_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")
            # Make "Start Drawing Line" button bold and orange after loading GeoTIFF
            if hasattr(self, 'line_start_draw_btn'):
                self.line_start_draw_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")

            # Show success message
            filename = os.path.basename(file_path)
            if downsample_factor > 1:
                success_msg = f"GeoTIFF '{filename}' loaded successfully (downsampled {downsample_factor}x for performance). Dynamic resolution enabled - zoom in for higher detail."
            else:
                success_msg = f"GeoTIFF '{filename}' loaded successfully. Dynamic resolution enabled - zoom in for higher detail."

            # Show message in activity log based on current tab
            current_tab = self.param_notebook.currentIndex() if hasattr(self, 'param_notebook') else 0
            if current_tab == 0 and hasattr(self, 'set_cal_info_text'):
                self.set_cal_info_text(success_msg)
            elif current_tab == 1 and hasattr(self, 'set_ref_info_text'):
                self.set_ref_info_text(success_msg)
            elif current_tab == 2 and hasattr(self, 'set_line_info_text'):
                self.set_line_info_text(success_msg)
            else:
                # Fallback: try to use any available method
                if hasattr(self, 'set_cal_info_text'):
                    self.set_cal_info_text(success_msg)
                elif hasattr(self, 'set_ref_info_text'):
                    self.set_ref_info_text(success_msg)
                elif hasattr(self, 'set_line_info_text'):
                    self.set_line_info_text(success_msg)

            if hasattr(self, "_rewarp_backscatter_align_to_bathy"):
                self._rewarp_backscatter_align_to_bathy()

            # Update the plot and zoom after loading
            self._plot_survey_plan()
            # Store the initial limits after first plot
            self.initial_geotiff_xlim = self.current_xlim
            self.initial_geotiff_ylim = self.current_ylim
            self._zoom_to_geotiff()

            if hasattr(self, "_on_geotiff_loaded_performance_depth"):
                self._on_geotiff_loaded_performance_depth()
            if hasattr(self, "_on_geotiff_loaded_line_info_depth"):
                self._on_geotiff_loaded_line_info_depth()
            if hasattr(self, "_update_roll_line_info_labels"):
                self._update_roll_line_info_labels()

        except RasterioIOError as e:
            progress_window.destroy()
            self._show_message("error", "GeoTIFF Error", f"Could not open GeoTIFF file: {e}")
            self._clear_plot(full_clear=True)
        except Exception as e:
            progress_window.destroy()
            self._show_message("error", "GeoTIFF Error", f"An unexpected error occurred while loading GeoTIFF: {e}")
            self._clear_plot(full_clear=True)

        # Update profile window for active tab and refresh canvas
        if hasattr(self, '_draw_current_profile'):
            self._draw_current_profile()
        else:
            self._draw_crossline_profile()
        if hasattr(self, 'profile_canvas'):
            self.profile_canvas.draw()

    def _on_geotiff_display_mode_changed(self, mode_text):
        """Handle display mode selection from combo box."""
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded. Cannot change display mode.")
            # Reset to previous mode
            if hasattr(self, 'elevation_slope_combo'):
                self._update_display_mode_combo()
            return

        if self.geotiff_data_array is None:
            self._show_message("warning", "No GeoTIFF", "Load a GeoTIFF first to change display mode.")
            # Reset to previous mode
            if hasattr(self, 'elevation_slope_combo'):
                self._update_display_mode_combo()
            return

        # Map combo box text to internal display modes
        if mode_text == "Shaded Relief":
            self.geotiff_display_mode = "elevation"
            self.hillshade_vs_slope_viz_mode = "hillshade"
        elif mode_text == "Shaded Slope":
            self.geotiff_display_mode = "slope"
            self.hillshade_vs_slope_viz_mode = "hillshade"
        elif mode_text == "Hillshade":
            self.geotiff_display_mode = "hillshade_only"
            self.hillshade_vs_slope_viz_mode = "hillshade"
        elif mode_text == "Slope":
            self.geotiff_display_mode = "slope_viz"
            self.hillshade_vs_slope_viz_mode = "slope_viz"

        # --- Preserve current plot area ---
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        self._plot_survey_plan()  # Re-plot with new display mode
        # Restore previous plot area if valid
        if xlim and ylim:
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)
        self.canvas.draw_idle()

        # Report the current layer in the info/error box
        msg = f"Layer shown: {mode_text}"
        if hasattr(self, 'set_cal_info_text') and self.param_notebook.currentIndex() == 0:
            self.set_cal_info_text(msg)
        elif hasattr(self, 'set_ref_info_text') and self.param_notebook.currentIndex() == 1:
            self.set_ref_info_text(msg)

    def _update_display_mode_combo(self):
        """Update the combo box to reflect the current display mode."""
        if not hasattr(self, 'elevation_slope_combo'):
            return

        # Map internal display mode to combo box text
        if self.geotiff_display_mode == "elevation":
            mode_text = "Shaded Relief"
        elif self.geotiff_display_mode == "slope":
            mode_text = "Shaded Slope"
        elif self.geotiff_display_mode == "hillshade_only":
            mode_text = "Hillshade"
        elif self.geotiff_display_mode == "slope_viz":
            mode_text = "Slope"
        else:
            mode_text = "Shaded Relief"  # Default

        # Block signals to prevent triggering the change handler
        self.elevation_slope_combo.blockSignals(True)
        self.elevation_slope_combo.setCurrentText(mode_text)
        self.elevation_slope_combo.blockSignals(False)

    def _on_contour_checkbox_changed(self):
        """Handle checkbox change for showing/hiding contours."""
        # Get the checkbox state
        if hasattr(self, 'show_contours_checkbox'):
            is_checked = self.show_contours_checkbox.isChecked()
        else:
            is_checked = False

        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded. Cannot show contours.")
            self.show_contours_var = False
            if hasattr(self, 'show_contours_checkbox'):
                self.show_contours_checkbox.setChecked(False)
            return

        if self.geotiff_data_array is None:
            self._show_message("warning", "No GeoTIFF", "Load a GeoTIFF first to show contours.")
            self.show_contours_var = False
            if hasattr(self, 'show_contours_checkbox'):
                self.show_contours_checkbox.setChecked(False)
            return

        # Update the variable based on checkbox state
        self.show_contours_var = is_checked

        # Preserve current plot area
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self._plot_survey_plan(preserve_view_limits=True)
            # Restore previous plot area if valid
            if xlim and ylim:
                self.ax.set_xlim(xlim)
                self.ax.set_ylim(ylim)
            self.canvas.draw_idle()
        except Exception:
            pass  # Silently handle errors

    def _on_contour_interval_changed(self):
        """Handle contour interval entry change with a short typing debounce."""
        if not hasattr(self, '_contour_interval_update_timer'):
            self._contour_interval_update_timer = QTimer()
            self._contour_interval_update_timer.setSingleShot(True)
            self._contour_interval_update_timer.timeout.connect(self._apply_contour_interval_changed)
        self._contour_interval_update_timer.start(450)

    def _apply_contour_interval_changed(self):
        """Apply contour interval updates after debounce delay."""
        # Keep all entry fields synchronized.
        # Sync all entry fields (calibration, reference, and line planning)
        try:
            if getattr(self, '_syncing_contour_interval', False):
                return  # Prevent infinite loop
            self._syncing_contour_interval = True
            try:
                # Try to get the value from whichever entry field exists and was changed
                val = None
                if hasattr(self, 'contour_interval_entry') and self.contour_interval_entry is not None:
                    val = self.contour_interval_entry.text()
            except Exception:
                pass
            finally:
                self._syncing_contour_interval = False
        except Exception:
            pass

        if not (hasattr(self, 'show_contours_var') and self.show_contours_var):
            return  # Don't update if contours are not enabled

        if not GEOSPATIAL_LIBS_AVAILABLE or self.geotiff_data_array is None:
            return

        # Preserve current plot area and update contours
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self._plot_survey_plan(preserve_view_limits=True)
            # Restore previous plot area if valid
            if xlim and ylim:
                self.ax.set_xlim(xlim)
                self.ax.set_ylim(ylim)
            self.canvas.draw_idle()
        except Exception:
            pass  # Silently handle errors

    def _on_slope_overlay_checkbox_changed(self):
        """Handle checkbox change for showing/hiding slope overlay."""
        # Get the checkbox state
        if hasattr(self, 'show_slope_overlay_checkbox'):
            is_checked = self.show_slope_overlay_checkbox.isChecked()
        else:
            is_checked = False

        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded. Cannot show slope overlay.")
            self.show_slope_overlay_var = False
            if hasattr(self, 'show_slope_overlay_checkbox'):
                self.show_slope_overlay_checkbox.setChecked(False)
            return

        if self.geotiff_data_array is None:
            self._show_message("warning", "No GeoTIFF", "Load a GeoTIFF first to show slope overlay.")
            self.show_slope_overlay_var = False
            if hasattr(self, 'show_slope_overlay_checkbox'):
                self.show_slope_overlay_checkbox.setChecked(False)
            return

        # Update the variable based on checkbox state
        self.show_slope_overlay_var = is_checked

        # Preserve current plot area
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self._plot_survey_plan(preserve_view_limits=True)
            # Restore previous plot area if valid
            if xlim and ylim:
                self.ax.set_xlim(xlim)
                self.ax.set_ylim(ylim)
            self.canvas.draw_idle()
        except Exception:
            pass  # Silently handle errors

    def _on_slope_overlay_min_changed(self):
        """Handle slope overlay min entry change with a short typing debounce."""
        if not hasattr(self, '_slope_overlay_min_update_timer'):
            self._slope_overlay_min_update_timer = QTimer()
            self._slope_overlay_min_update_timer.setSingleShot(True)
            self._slope_overlay_min_update_timer.timeout.connect(self._apply_slope_overlay_min_changed)
        self._slope_overlay_min_update_timer.start(450)

    def _apply_slope_overlay_min_changed(self):
        """Apply slope overlay min updates after debounce delay."""
        try:
            if hasattr(self, 'slope_overlay_min_entry') and self.slope_overlay_min_entry:
                min_val = float(self.slope_overlay_min_entry.text())
                self.slope_overlay_min_var = min_val
        except (ValueError, AttributeError):
            pass  # Silently handle invalid input

        # Only trigger redraw if slope overlay is currently visible.
        if not (hasattr(self, 'show_slope_overlay_var') and self.show_slope_overlay_var):
            return
        if not GEOSPATIAL_LIBS_AVAILABLE or self.geotiff_data_array is None:
            return

        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self._plot_survey_plan(preserve_view_limits=True)
            if xlim and ylim:
                self.ax.set_xlim(xlim)
                self.ax.set_ylim(ylim)
            self.canvas.draw_idle()
        except Exception:
            pass

    def _on_slope_overlay_max_changed(self):
        """Handle slope overlay max entry change with a short typing debounce."""
        if not hasattr(self, '_slope_overlay_max_update_timer'):
            self._slope_overlay_max_update_timer = QTimer()
            self._slope_overlay_max_update_timer.setSingleShot(True)
            self._slope_overlay_max_update_timer.timeout.connect(self._apply_slope_overlay_max_changed)
        self._slope_overlay_max_update_timer.start(450)

    def _apply_slope_overlay_max_changed(self):
        """Apply slope overlay max updates after debounce delay."""
        try:
            if hasattr(self, 'slope_overlay_max_entry') and self.slope_overlay_max_entry:
                max_val = float(self.slope_overlay_max_entry.text())
                self.slope_overlay_max_var = max_val
        except (ValueError, AttributeError):
            pass  # Silently handle invalid input

        # Only trigger redraw if slope overlay is currently visible.
        if not (hasattr(self, 'show_slope_overlay_var') and self.show_slope_overlay_var):
            return
        if not GEOSPATIAL_LIBS_AVAILABLE or self.geotiff_data_array is None:
            return

        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self._plot_survey_plan(preserve_view_limits=True)
            # Restore previous plot area if valid
            if xlim and ylim:
                self.ax.set_xlim(xlim)
                self.ax.set_ylim(ylim)
            self.canvas.draw_idle()
        except Exception:
            pass  # Silently handle errors

    def _on_backscatter_slope_areas_checkbox_changed(self):
        """Show/hide Backscatter bathymetry slope band overlay (magenta)."""
        if hasattr(self, 'backscatter_show_slope_areas_checkbox'):
            is_checked = self.backscatter_show_slope_areas_checkbox.isChecked()
        else:
            is_checked = False

        if is_checked:
            if not GEOSPATIAL_LIBS_AVAILABLE:
                self._show_message(
                    "warning",
                    "Disabled Feature",
                    "Geospatial libraries are not loaded. Cannot show slope areas.",
                )
                self.backscatter_slope_areas_show_var = False
                if hasattr(self, 'backscatter_show_slope_areas_checkbox'):
                    self.backscatter_show_slope_areas_checkbox.blockSignals(True)
                    self.backscatter_show_slope_areas_checkbox.setChecked(False)
                    self.backscatter_show_slope_areas_checkbox.blockSignals(False)
                return
            if self.geotiff_data_array is None:
                self._show_message("warning", "No GeoTIFF", "Load a GeoTIFF first to show slope areas.")
                self.backscatter_slope_areas_show_var = False
                if hasattr(self, 'backscatter_show_slope_areas_checkbox'):
                    self.backscatter_show_slope_areas_checkbox.blockSignals(True)
                    self.backscatter_show_slope_areas_checkbox.setChecked(False)
                    self.backscatter_show_slope_areas_checkbox.blockSignals(False)
                return

        self.backscatter_slope_areas_show_var = is_checked
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self._plot_survey_plan(preserve_view_limits=True)
            if xlim and ylim:
                self.ax.set_xlim(xlim)
                self.ax.set_ylim(ylim)
            self.canvas.draw_idle()
        except Exception:
            pass

    def _on_backscatter_slope_min_max_changed(self, *_args):
        """Update Backscatter slope band from spin boxes and refresh overlay if visible."""
        if hasattr(self, 'backscatter_slope_min_spin'):
            self.backscatter_slope_min_deg = float(self.backscatter_slope_min_spin.value())
        if hasattr(self, 'backscatter_slope_max_spin'):
            self.backscatter_slope_max_deg = float(self.backscatter_slope_max_spin.value())
        if not getattr(self, 'backscatter_slope_areas_show_var', False):
            return
        if not GEOSPATIAL_LIBS_AVAILABLE or self.geotiff_data_array is None:
            return
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self._plot_survey_plan(preserve_view_limits=True)
            if xlim and ylim:
                self.ax.set_xlim(xlim)
                self.ax.set_ylim(ylim)
            self.canvas.draw_idle()
        except Exception:
            pass

    def _on_backscatter_min_area_checkbox_changed(self, *_args):
        """Enable/disable minimum contiguous area filter for the magenta Backscatter overlay."""
        if hasattr(self, "backscatter_min_area_checkbox"):
            self.backscatter_min_area_enabled_var = self.backscatter_min_area_checkbox.isChecked()
        if hasattr(self, "backscatter_min_area_m2_entry"):
            raw = self.backscatter_min_area_m2_entry.text().strip()
            if raw:
                try:
                    v = float(raw)
                    if v >= 0:
                        self.backscatter_min_area_m2 = v
                except ValueError:
                    pass
        if not getattr(self, "backscatter_slope_areas_show_var", False):
            return
        if not GEOSPATIAL_LIBS_AVAILABLE or self.geotiff_data_array is None:
            return
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self._plot_survey_plan(preserve_view_limits=True)
            if xlim and ylim:
                self.ax.set_xlim(xlim)
                self.ax.set_ylim(ylim)
            self.canvas.draw_idle()
        except Exception:
            pass

    def _on_backscatter_min_area_m2_text_changed(self, *_args):
        """Debounce minimum area (m²) typing before applying and redrawing."""
        if not hasattr(self, "_backscatter_min_area_m2_debounce_timer"):
            self._backscatter_min_area_m2_debounce_timer = QTimer()
            self._backscatter_min_area_m2_debounce_timer.setSingleShot(True)
            self._backscatter_min_area_m2_debounce_timer.timeout.connect(self._apply_backscatter_min_area_m2_text)
        self._backscatter_min_area_m2_debounce_timer.stop()
        self._backscatter_min_area_m2_debounce_timer.start(450)

    def _apply_backscatter_min_area_m2_text(self):
        """Parse minimum area from the line edit after debounce; replot when Show Areas + Min Area are on."""
        if not hasattr(self, "backscatter_min_area_m2_entry"):
            return
        raw = self.backscatter_min_area_m2_entry.text().strip()
        if not raw:
            return
        try:
            v = float(raw)
        except ValueError:
            return
        if v < 0:
            return
        self.backscatter_min_area_m2 = v
        if not getattr(self, "backscatter_slope_areas_show_var", False):
            return
        if not getattr(self, "backscatter_min_area_enabled_var", False):
            return
        if not GEOSPATIAL_LIBS_AVAILABLE or self.geotiff_data_array is None:
            return
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self._plot_survey_plan(preserve_view_limits=True)
            if xlim and ylim:
                self.ax.set_xlim(xlim)
                self.ax.set_ylim(ylim)
            self.canvas.draw_idle()
        except Exception:
            pass

    def _on_backscatter_extent_m_text_changed(self, *_args):
        """Debounce min width/height (m) typing before applying and redrawing."""
        if not hasattr(self, "_backscatter_extent_m_debounce_timer"):
            self._backscatter_extent_m_debounce_timer = QTimer()
            self._backscatter_extent_m_debounce_timer.setSingleShot(True)
            self._backscatter_extent_m_debounce_timer.timeout.connect(self._apply_backscatter_extent_m_text)
        self._backscatter_extent_m_debounce_timer.stop()
        self._backscatter_extent_m_debounce_timer.start(450)

    def _apply_backscatter_extent_m_text(self):
        """Parse W/H from line edits after debounce; replot when Show Areas + extent filter are on."""
        if hasattr(self, "backscatter_min_width_m_entry"):
            tw = self.backscatter_min_width_m_entry.text().strip()
            if tw:
                try:
                    wv = float(tw)
                    if wv >= 0:
                        self.backscatter_min_width_m = wv
                except ValueError:
                    pass
        if hasattr(self, "backscatter_min_height_m_entry"):
            th = self.backscatter_min_height_m_entry.text().strip()
            if th:
                try:
                    hv = float(th)
                    if hv >= 0:
                        self.backscatter_min_height_m = hv
                except ValueError:
                    pass
        if not getattr(self, "backscatter_slope_areas_show_var", False):
            return
        if not getattr(self, "backscatter_extent_filter_enabled_var", False):
            return
        if not GEOSPATIAL_LIBS_AVAILABLE or self.geotiff_data_array is None:
            return
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self._plot_survey_plan(preserve_view_limits=True)
            if xlim and ylim:
                self.ax.set_xlim(xlim)
                self.ax.set_ylim(ylim)
            self.canvas.draw_idle()
        except Exception:
            pass

    def _on_backscatter_extent_checkbox_changed(self, *_args):
        """Enable/disable minimum bounding width/height filter (meters) for the magenta overlay.
        W/H fields stay editable while the filter is off so thresholds can be set first."""
        if hasattr(self, "backscatter_extent_checkbox"):
            self.backscatter_extent_filter_enabled_var = self.backscatter_extent_checkbox.isChecked()
        if hasattr(self, "backscatter_min_width_m_entry"):
            tw = self.backscatter_min_width_m_entry.text().strip()
            if tw:
                try:
                    wv = float(tw)
                    if wv >= 0:
                        self.backscatter_min_width_m = wv
                except ValueError:
                    pass
        if hasattr(self, "backscatter_min_height_m_entry"):
            th = self.backscatter_min_height_m_entry.text().strip()
            if th:
                try:
                    hv = float(th)
                    if hv >= 0:
                        self.backscatter_min_height_m = hv
                except ValueError:
                    pass
        if not getattr(self, "backscatter_slope_areas_show_var", False):
            return
        if not GEOSPATIAL_LIBS_AVAILABLE or self.geotiff_data_array is None:
            return
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self._plot_survey_plan(preserve_view_limits=True)
            if xlim and ylim:
                self.ax.set_xlim(xlim)
                self.ax.set_ylim(ylim)
            self.canvas.draw_idle()
        except Exception:
            pass

    def _on_backscatter_depth_filter_changed(self, *_args):
        """Enable/disable depth-range filter for magenta overlay and refresh."""
        if hasattr(self, "backscatter_depth_filter_checkbox"):
            self.backscatter_depth_filter_enabled_var = self.backscatter_depth_filter_checkbox.isChecked()
        if hasattr(self, "backscatter_depth_min_m_entry"):
            raw_min = self.backscatter_depth_min_m_entry.text().strip()
            if raw_min:
                try:
                    v = float(raw_min)
                    if v >= 0.0:
                        self.backscatter_depth_min_m = v
                except ValueError:
                    pass
        if hasattr(self, "backscatter_depth_max_m_entry"):
            raw_max = self.backscatter_depth_max_m_entry.text().strip()
            if raw_max:
                try:
                    v = float(raw_max)
                    if v >= 0.0:
                        self.backscatter_depth_max_m = v
                except ValueError:
                    pass
        if not getattr(self, "backscatter_slope_areas_show_var", False):
            return
        if not GEOSPATIAL_LIBS_AVAILABLE or self.geotiff_data_array is None:
            return
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self._plot_survey_plan(preserve_view_limits=True)
            if xlim and ylim:
                self.ax.set_xlim(xlim)
                self.ax.set_ylim(ylim)
            self.canvas.draw_idle()
        except Exception:
            pass

    def _on_backscatter_depth_text_changed(self, *_args):
        """Debounce depth-range edits before refreshing magenta overlay."""
        if not hasattr(self, "_backscatter_depth_debounce_timer"):
            self._backscatter_depth_debounce_timer = QTimer()
            self._backscatter_depth_debounce_timer.setSingleShot(True)
            self._backscatter_depth_debounce_timer.timeout.connect(self._apply_backscatter_depth_text)
        self._backscatter_depth_debounce_timer.stop()
        self._backscatter_depth_debounce_timer.start(350)

    def _apply_backscatter_depth_text(self):
        """Parse depth min/max and redraw if magenta overlay is enabled."""
        if hasattr(self, "backscatter_depth_min_m_entry"):
            raw_min = self.backscatter_depth_min_m_entry.text().strip()
            if raw_min:
                try:
                    v = float(raw_min)
                    if v >= 0.0:
                        self.backscatter_depth_min_m = v
                except ValueError:
                    pass
        if hasattr(self, "backscatter_depth_max_m_entry"):
            raw_max = self.backscatter_depth_max_m_entry.text().strip()
            if raw_max:
                try:
                    v = float(raw_max)
                    if v >= 0.0:
                        self.backscatter_depth_max_m = v
                except ValueError:
                    pass
        if not getattr(self, "backscatter_slope_areas_show_var", False):
            return
        if not GEOSPATIAL_LIBS_AVAILABLE or self.geotiff_data_array is None:
            return
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self._plot_survey_plan(preserve_view_limits=True)
            if xlim and ylim:
                self.ax.set_xlim(xlim)
                self.ax.set_ylim(ylim)
            self.canvas.draw_idle()
        except Exception:
            pass

    def _on_backscatter_line_info_text_changed(self, *_args):
        """Debounce Lead-in / Survey Speed typing, then apply and redraw centerline if needed."""
        if not hasattr(self, "_backscatter_line_info_debounce_timer"):
            self._backscatter_line_info_debounce_timer = QTimer()
            self._backscatter_line_info_debounce_timer.setSingleShot(True)
            self._backscatter_line_info_debounce_timer.timeout.connect(self._apply_backscatter_line_info_text)
        self._backscatter_line_info_debounce_timer.stop()
        self._backscatter_line_info_debounce_timer.start(350)

    def _apply_backscatter_line_info_text(self):
        """Parse Lead-in / Survey Speed values; refresh map if centerline is present."""
        if hasattr(self, "backscatter_lead_in_entry"):
            raw_lead = self.backscatter_lead_in_entry.text().strip()
            if raw_lead:
                try:
                    lead_m = float(raw_lead)
                    if lead_m >= 0.0:
                        self.backscatter_lead_in_m = lead_m
                except ValueError:
                    pass
        if hasattr(self, "backscatter_survey_speed_entry"):
            raw_speed = self.backscatter_survey_speed_entry.text().strip()
            if raw_speed:
                try:
                    speed_kn = float(raw_speed)
                    if speed_kn > 0.0:
                        self.backscatter_survey_speed_kn = speed_kn
                except ValueError:
                    pass
        centerline = getattr(self, "backscatter_box_centerline", None)
        if not centerline or centerline[0] is None or centerline[1] is None:
            return
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self._plot_survey_plan(preserve_view_limits=True)
            if xlim and ylim:
                self.ax.set_xlim(xlim)
                self.ax.set_ylim(ylim)
            self.canvas.draw_idle()
        except Exception:
            pass
        if hasattr(self, "_draw_current_profile"):
            self._draw_current_profile()
        if hasattr(self, "_update_backscatter_export_name_default"):
            self._update_backscatter_export_name_default()
        if hasattr(self, "_update_backscatter_export_name_default"):
            self._update_backscatter_export_name_default()

    def _toggle_slope_visualization(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning", "Disabled Feature",
                               "Geospatial libraries not loaded. Cannot toggle slope visualization.")
            return

        if self.geotiff_data_array is None:
            self._show_message("warning", "No GeoTIFF", "Load a GeoTIFF first to toggle slope visualization.")
            return

        # When toggling slope visualization, if we were in the "elevation" overlay mode,
        # we don't need to change anything specifically because we're not toggling *overlay* anymore.
        # The relevant part is just switching the main hillshade vs. slope_viz mode.

        if self.hillshade_vs_slope_viz_mode == "hillshade":
            self.hillshade_vs_slope_viz_mode = "slope_viz"
            # Disable combo box in this mode
            if hasattr(self, 'elevation_slope_combo'):
                self.elevation_slope_combo.setEnabled(False)
        else:
            self.hillshade_vs_slope_viz_mode = "hillshade"
            # Re-enable combo box
            self.geotiff_display_mode = "elevation"  # Reset to elevation when going back to hillshade view
            if hasattr(self, 'elevation_slope_combo'):
                self.elevation_slope_combo.setEnabled(True)
                self._update_display_mode_combo()

        self._plot_survey_plan()

    def _remove_geotiff(self):
        """Remove the loaded GeoTIFF and clear it from the plot."""
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded.")
            return

        if self.geotiff_dataset_original is None:
            self._show_message("warning", "No GeoTIFF", "No GeoTIFF is currently loaded.")
            return

        # Clear the GeoTIFF using the existing clear plot method with full_clear=True
        self._clear_plot(full_clear=True)

        # Re-plot without GeoTIFF
        self._plot_survey_plan(preserve_view_limits=False)
        if hasattr(self, "_on_geotiff_loaded_line_info_depth"):
            self._on_geotiff_loaded_line_info_depth()
        if hasattr(self, "_update_roll_line_info_labels"):
            self._update_roll_line_info_labels()

        # Show message in activity log
        if hasattr(self, 'set_cal_info_text'):
            self.set_cal_info_text("GeoTIFF removed.")
        elif hasattr(self, 'set_ref_info_text'):
            self.set_ref_info_text("GeoTIFF removed.")

    def _reproject_backscatter_band_to_bathy_grid(self, src_ds, band_index=1):
        """Reproject one band of an open rasterio dataset onto the current bathymetry WGS84 grid."""
        nrows, ncols = self.geotiff_data_array.shape
        left, right, bottom, top = self.geotiff_extent
        dst_transform = transform.from_bounds(left, bottom, right, top, ncols, nrows)
        dst = np.full((nrows, ncols), np.nan, dtype=np.float32)
        reproject(
            source=rasterio.band(src_ds, int(band_index)),
            destination=dst,
            src_transform=src_ds.transform,
            src_crs=src_ds.crs,
            dst_transform=dst_transform,
            dst_crs="EPSG:4326",
            resampling=Resampling.bilinear,
        )
        dst[np.isnan(self.geotiff_data_array)] = np.nan
        return dst

    def _set_backscatter_band_controls_visibility(self, is_rgb):
        """Show/hide RGB band picker controls based on loaded backscatter type."""
        if hasattr(self, "backscatter_band_combo"):
            self.backscatter_band_combo.setVisible(bool(is_rgb))
            self.backscatter_band_combo.setEnabled(bool(is_rgb))
            if is_rgb:
                self.backscatter_band_combo.blockSignals(True)
                self.backscatter_band_combo.setCurrentIndex(max(0, min(2, int(getattr(self, "backscatter_selected_band", 1)) - 1)))
                self.backscatter_band_combo.blockSignals(False)

    def _on_backscatter_band_changed(self, *_args):
        """Switch displayed backscatter band for RGB rasters."""
        if not getattr(self, "backscatter_is_rgb", False):
            return
        if not hasattr(self, "backscatter_band_combo"):
            return
        idx = int(self.backscatter_band_combo.currentIndex())
        self.backscatter_selected_band = idx + 1
        bands = getattr(self, "backscatter_raster_bands", []) or []
        if 0 <= idx < len(bands):
            self.backscatter_raster_array = bands[idx]
        if self.geotiff_data_array is None:
            return
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self._plot_survey_plan(preserve_view_limits=True)
            if xlim and ylim:
                self.ax.set_xlim(xlim)
                self.ax.set_ylim(ylim)
            self.canvas.draw_idle()
        except Exception:
            pass

    def _rewarp_backscatter_align_to_bathy(self):
        """Recompute backscatter_raster_array from backscatter_geotiff_path to match current bathymetry grid."""
        path = getattr(self, "backscatter_geotiff_path", None)
        if not path or not os.path.isfile(path):
            self.backscatter_raster_array = None
            self.backscatter_raster_bands = []
            self.backscatter_is_rgb = False
            self._set_backscatter_band_controls_visibility(False)
            return
        if self.geotiff_data_array is None or self.geotiff_extent is None:
            self.backscatter_raster_array = None
            self.backscatter_raster_bands = []
            self.backscatter_is_rgb = False
            self._set_backscatter_band_controls_visibility(False)
            return
        try:
            with rasterio.open(path) as src:
                if int(getattr(src, "count", 1)) >= 3:
                    self.backscatter_is_rgb = True
                    self.backscatter_raster_bands = [
                        self._reproject_backscatter_band_to_bathy_grid(src, 1),
                        self._reproject_backscatter_band_to_bathy_grid(src, 2),
                        self._reproject_backscatter_band_to_bathy_grid(src, 3),
                    ]
                    self.backscatter_selected_band = max(1, min(3, int(getattr(self, "backscatter_selected_band", 1))))
                    self.backscatter_raster_array = self.backscatter_raster_bands[self.backscatter_selected_band - 1]
                    self._set_backscatter_band_controls_visibility(True)
                else:
                    self.backscatter_is_rgb = False
                    self.backscatter_raster_bands = []
                    self.backscatter_selected_band = 1
                    self.backscatter_raster_array = self._reproject_backscatter_band_to_bathy_grid(src, 1)
                    self._set_backscatter_band_controls_visibility(False)
        except Exception as e:
            self.backscatter_raster_array = None
            self.backscatter_raster_bands = []
            self.backscatter_is_rgb = False
            self._set_backscatter_band_controls_visibility(False)
            self._show_message(
                "warning",
                "Backscatter GeoTIFF",
                f"Could not align backscatter grid to bathymetry: {e}",
            )

    def _load_backscatter_geotiff(self):
        """Prompt for a backscatter GeoTIFF and load it aligned to the bathymetry grid."""
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message(
                "warning",
                "Disabled Feature",
                "Geospatial libraries are not loaded. Cannot load backscatter GeoTIFF.",
            )
            return
        if self.geotiff_data_array is None or self.geotiff_extent is None:
            self._show_message(
                "warning",
                "No bathymetry grid",
                "Load the bathymetry GeoTIFF first. Backscatter is resampled to match its current grid size and extent.",
            )
            return
        start_dir = getattr(self, "last_backscatter_dir", os.path.expanduser("~"))
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Load Backscatter GeoTIFF",
            start_dir,
            "GeoTIFF (*.tif *.tiff *.gtiff);;All files (*.*)",
        )
        if not path:
            return
        self.last_backscatter_dir = os.path.dirname(path)
        if hasattr(self, "_save_last_backscatter_dir"):
            self._save_last_backscatter_dir()
        self._load_backscatter_geotiff_from_path(path)

    def _load_backscatter_geotiff_from_path(self, path):
        """Open path, warp band 1 to the current bathymetry grid, store path and array."""
        if not GEOSPATIAL_LIBS_AVAILABLE:
            return
        if self.geotiff_data_array is None or self.geotiff_extent is None:
            self._show_message("warning", "No bathymetry grid", "Load the bathymetry GeoTIFF first.")
            return
        try:
            self.last_backscatter_dir = os.path.dirname(path)
            if hasattr(self, "_save_last_backscatter_dir"):
                self._save_last_backscatter_dir()
        except Exception:
            pass
        try:
            with rasterio.open(path) as src:
                self.backscatter_geotiff_path = path
                if int(getattr(src, "count", 1)) >= 3:
                    self.backscatter_is_rgb = True
                    self.backscatter_raster_bands = [
                        self._reproject_backscatter_band_to_bathy_grid(src, 1),
                        self._reproject_backscatter_band_to_bathy_grid(src, 2),
                        self._reproject_backscatter_band_to_bathy_grid(src, 3),
                    ]
                    self.backscatter_selected_band = 1  # default Red
                    self.backscatter_raster_array = self.backscatter_raster_bands[0]
                    self._set_backscatter_band_controls_visibility(True)
                else:
                    self.backscatter_is_rgb = False
                    self.backscatter_raster_bands = []
                    self.backscatter_selected_band = 1
                    self.backscatter_raster_array = self._reproject_backscatter_band_to_bathy_grid(src, 1)
                    self._set_backscatter_band_controls_visibility(False)
        except Exception as e:
            self.backscatter_geotiff_path = None
            self.backscatter_raster_array = None
            self.backscatter_raster_bands = []
            self.backscatter_is_rgb = False
            self._set_backscatter_band_controls_visibility(False)
            self._show_message("error", "Backscatter GeoTIFF", f"Could not load file: {e}")
            return

        nrows, ncols = self.geotiff_data_array.shape
        msg = (
            f"Loaded backscatter GeoTIFF ({os.path.basename(path)}), "
            f"resampled to bathymetry grid {nrows} × {ncols}."
        )
        if hasattr(self, "set_cal_info_text"):
            self.set_cal_info_text(msg)
        elif hasattr(self, "set_ref_info_text"):
            self.set_ref_info_text(msg)
        elif hasattr(self, "set_line_info_text"):
            self.set_line_info_text(msg)

        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self._plot_survey_plan(preserve_view_limits=True)
            if xlim and ylim:
                self.ax.set_xlim(xlim)
                self.ax.set_ylim(ylim)
            self.canvas.draw_idle()
        except Exception:
            pass

    def _on_show_backscatter_checkbox_changed(self, *_args):
        """Toggle visibility of the backscatter raster overlay."""
        if hasattr(self, "show_backscatter_checkbox"):
            self.show_backscatter_var = self.show_backscatter_checkbox.isChecked()
        if self.geotiff_data_array is None:
            return
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self._plot_survey_plan(preserve_view_limits=True)
            if xlim and ylim:
                self.ax.set_xlim(xlim)
                self.ax.set_ylim(ylim)
            self.canvas.draw_idle()
        except Exception:
            pass

    def _on_backscatter_percent_clip_changed(self, *_args):
        """Debounce backscatter percent-clip control edits and refresh overlay."""
        if not hasattr(self, "_backscatter_percent_clip_debounce_timer"):
            self._backscatter_percent_clip_debounce_timer = QTimer()
            self._backscatter_percent_clip_debounce_timer.setSingleShot(True)
            self._backscatter_percent_clip_debounce_timer.timeout.connect(self._apply_backscatter_percent_clip_changed)
        self._backscatter_percent_clip_debounce_timer.stop()
        self._backscatter_percent_clip_debounce_timer.start(350)

    def _apply_backscatter_percent_clip_changed(self):
        """Apply percent-clip settings and redraw map when possible."""
        if hasattr(self, "backscatter_percent_clip_checkbox"):
            self.backscatter_percent_clip_enabled = self.backscatter_percent_clip_checkbox.isChecked()
        if hasattr(self, "backscatter_percent_clip_min_entry"):
            raw_min = self.backscatter_percent_clip_min_entry.text().strip()
            if raw_min:
                try:
                    v = float(raw_min)
                    if v >= 0.0:
                        self.backscatter_percent_clip_min = v
                except ValueError:
                    pass
        if hasattr(self, "backscatter_percent_clip_max_entry"):
            raw_max = self.backscatter_percent_clip_max_entry.text().strip()
            if raw_max:
                try:
                    v = float(raw_max)
                    if v >= 0.0:
                        self.backscatter_percent_clip_max = v
                except ValueError:
                    pass

        if self.geotiff_data_array is None:
            return
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self._plot_survey_plan(preserve_view_limits=True)
            if xlim and ylim:
                self.ax.set_xlim(xlim)
                self.ax.set_ylim(ylim)
            self.canvas.draw_idle()
        except Exception:
            pass

    def _on_backscatter_draw_box_toggled(self, checked):
        """Enable/disable map box drawing mode for backscatter statistics."""
        self.backscatter_box_draw_mode = bool(checked)
        if self.backscatter_box_draw_mode:
            self.backscatter_box_edit_width_mode = False
            if hasattr(self, "backscatter_edit_width_btn"):
                self.backscatter_edit_width_btn.blockSignals(True)
                self.backscatter_edit_width_btn.setChecked(False)
                self.backscatter_edit_width_btn.blockSignals(False)
            self.backscatter_box_stage = 0
            self.backscatter_box_centerline = None
            self.backscatter_box_width_point = None
            if hasattr(self, "backscatter_draw_box_btn"):
                self.backscatter_draw_box_btn.setText("Draw Center Line")
                self.backscatter_draw_box_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")
            self.canvas_widget.setCursor(Qt.CursorShape.CrossCursor)
        else:
            had_partial = int(getattr(self, "backscatter_box_stage", 0)) > 0
            self.backscatter_box_stage = 0
            if had_partial:
                self.backscatter_box_centerline = None
                self.backscatter_box_width_point = None
            if hasattr(self, "_clear_backscatter_draw_tooltip"):
                self._clear_backscatter_draw_tooltip()
            if hasattr(self, "backscatter_draw_box_btn"):
                self.backscatter_draw_box_btn.setText("Select Area/Line")
                self.backscatter_draw_box_btn.setStyleSheet("")
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
            if had_partial:
                try:
                    self._plot_survey_plan(preserve_view_limits=True)
                except Exception:
                    self.canvas.draw_idle()

    def _on_backscatter_edit_width_toggled(self, checked):
        """Enable/disable width edit mode for an existing backscatter area centerline."""
        self.backscatter_box_edit_width_mode = bool(checked)
        if self.backscatter_box_edit_width_mode:
            centerline = getattr(self, "backscatter_box_centerline", None)
            if not centerline or centerline[0] is None or centerline[1] is None:
                self.backscatter_box_edit_width_mode = False
                if hasattr(self, "backscatter_edit_width_btn"):
                    self.backscatter_edit_width_btn.blockSignals(True)
                    self.backscatter_edit_width_btn.setChecked(False)
                    self.backscatter_edit_width_btn.setText("Edit Area Width")
                    self.backscatter_edit_width_btn.setStyleSheet("")
                    self.backscatter_edit_width_btn.blockSignals(False)
                self._show_message("info", "Edit Area Width", "Select an area first, then edit width.")
                return
            self.backscatter_box_draw_mode = False
            self.backscatter_box_stage = 0
            if hasattr(self, "backscatter_draw_box_btn"):
                self.backscatter_draw_box_btn.blockSignals(True)
                self.backscatter_draw_box_btn.setChecked(False)
                self.backscatter_draw_box_btn.blockSignals(False)
            if hasattr(self, "backscatter_edit_width_btn"):
                self.backscatter_edit_width_btn.setText("Set Width")
                self.backscatter_edit_width_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")
            self.canvas_widget.setCursor(Qt.CursorShape.CrossCursor)
        else:
            if hasattr(self, "backscatter_edit_width_btn"):
                self.backscatter_edit_width_btn.setText("Edit Area Width")
                self.backscatter_edit_width_btn.setStyleSheet("")
            if hasattr(self, "_clear_backscatter_draw_tooltip"):
                self._clear_backscatter_draw_tooltip()
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)

    def _on_backscatter_show_stats_clicked(self):
        """Open backscatter statistics dialog for the currently saved box."""
        stats = getattr(self, "backscatter_box_stats", None)
        has_saved_geometry = bool(
            (
                getattr(self, "backscatter_box_vertices", None)
                and len(getattr(self, "backscatter_box_vertices", []) or []) == 4
            )
            or (
                getattr(self, "backscatter_box_centerline", None)
                and len(getattr(self, "backscatter_box_centerline", []) or []) == 2
            )
        )
        if not stats:
            if not has_saved_geometry:
                self._show_message("info", "Backscatter Statistics", "Draw a box first to compute statistics.")
                return

            # Imported plans can have saved geometry without cached stats; compute on demand.
            if getattr(self, "geotiff_extent", None) is None:
                self._show_message(
                    "info",
                    "Backscatter Statistics",
                    "Load bathymetry first so area statistics can be computed.",
                )
                return
            try:
                self.backscatter_box_stats = self._compute_backscatter_box_stats()
            except Exception:
                self.backscatter_box_stats = None
            stats = getattr(self, "backscatter_box_stats", None)
        if not stats:
            self._show_message(
                "info",
                "Backscatter Statistics",
                "Could not compute statistics for the current area. Try reloading the GeoTIFF(s) and reopening stats.",
            )
            return
        self._show_backscatter_stats_dialog()

    def _on_backscatter_zoom_to_line_clicked(self):
        """Zoom map to selected backscatter area bounds (box, not centerline)."""
        vertices = getattr(self, "backscatter_box_vertices", None)
        if not vertices or len(vertices) != 4:
            self._show_message("info", "Backscatter Area", "Select an area first.")
            return
        verts = np.array(vertices, dtype=float)
        if verts.shape != (4, 2):
            self._show_message("info", "Backscatter Area", "Select an area first.")
            return
        lon_min = float(np.min(verts[:, 0]))
        lon_max = float(np.max(verts[:, 0]))
        lat_min = float(np.min(verts[:, 1]))
        lat_max = float(np.max(verts[:, 1]))
        dlon = max(1e-6, lon_max - lon_min)
        dlat = max(1e-6, lat_max - lat_min)
        xlim = (lon_min - 0.05 * dlon, lon_max + 0.05 * dlon)
        ylim = (lat_min - 0.05 * dlat, lat_max + 0.05 * dlat)
        if getattr(self, "geotiff_dataset_original", None) is not None and hasattr(self, "_apply_map_zoom_limits_and_reload_geotiff"):
            self._apply_map_zoom_limits_and_reload_geotiff(xlim, ylim)
        else:
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)
            self.canvas.draw_idle()

    def _on_backscatter_remove_line_clicked(self):
        """Remove backscatter normalization centerline and area overlay."""
        self._on_backscatter_clear_box_clicked()

    def _on_backscatter_clear_box_clicked(self):
        """Remove backscatter stats box overlay and clear cached statistics."""
        self.backscatter_box_draw_mode = False
        self.backscatter_box_edit_width_mode = False
        self.backscatter_box_stage = 0
        self.backscatter_box_centerline = None
        self.backscatter_box_width_point = None
        if hasattr(self, "_clear_backscatter_draw_tooltip"):
            self._clear_backscatter_draw_tooltip()
        self.backscatter_box_vertices = None
        self.backscatter_box_half_width_m = None
        self.backscatter_box_stats = None
        if hasattr(self, "backscatter_draw_box_btn"):
            self.backscatter_draw_box_btn.blockSignals(True)
            self.backscatter_draw_box_btn.setChecked(False)
            self.backscatter_draw_box_btn.setText("Select Area/Line")
            self.backscatter_draw_box_btn.setStyleSheet("")
            self.backscatter_draw_box_btn.blockSignals(False)
        if hasattr(self, "backscatter_edit_width_btn"):
            self.backscatter_edit_width_btn.blockSignals(True)
            self.backscatter_edit_width_btn.setChecked(False)
            self.backscatter_edit_width_btn.setText("Edit Area Width")
            self.backscatter_edit_width_btn.setStyleSheet("")
            self.backscatter_edit_width_btn.blockSignals(False)
        if hasattr(self, "backscatter_box_patch") and self.backscatter_box_patch is not None:
            try:
                self.backscatter_box_patch.remove()
            except Exception:
                pass
            self.backscatter_box_patch = None
        if hasattr(self, "backscatter_box_centerline_line") and self.backscatter_box_centerline_line is not None:
            try:
                self.backscatter_box_centerline_line.remove()
            except Exception:
                pass
            self.backscatter_box_centerline_line = None
        if hasattr(self, "backscatter_stats_dialog") and self.backscatter_stats_dialog is not None:
            try:
                self.backscatter_stats_dialog.close()
            except Exception:
                pass
            self.backscatter_stats_dialog = None
            self.backscatter_stats_canvas = None
            self.backscatter_stats_img_ax = None
            self.backscatter_stats_ax = None
            self.backscatter_stats_text_label = None
        self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
        try:
            self._plot_survey_plan(preserve_view_limits=True)
        except Exception:
            self.canvas.draw_idle()
        if hasattr(self, "_draw_current_profile"):
            self._draw_current_profile()

    def _update_backscatter_box_patch(self, vertices):
        """Create/update temporary cyan oriented rectangle while drawing."""
        if not vertices or len(vertices) != 4:
            return
        if hasattr(self, "backscatter_box_patch") and self.backscatter_box_patch is not None:
            try:
                self.backscatter_box_patch.set_xy(vertices)
            except Exception:
                self.backscatter_box_patch = None
        if self.backscatter_box_patch is None:
            try:
                from matplotlib.patches import Polygon
                self.backscatter_box_patch = Polygon(
                    vertices,
                    closed=True,
                    fill=False,
                    edgecolor="cyan",
                    linewidth=2.0,
                    linestyle="-",
                    zorder=40,
                )
                self.ax.add_patch(self.backscatter_box_patch)
            except Exception:
                self.backscatter_box_patch = None

    def _backscatter_oriented_box_vertices(self, a_lat, a_lon, b_lat, b_lon, half_width_m):
        """Return 4 oriented rectangle vertices around centerline AB with symmetric half-width."""
        center_lat = (a_lat + b_lat) / 2.0
        m_per_deg_lat = 111320.0
        m_per_deg_lon = 111320.0 * np.cos(np.radians(center_lat))
        if m_per_deg_lon == 0:
            return None
        ax, ay = a_lon * m_per_deg_lon, a_lat * m_per_deg_lat
        bx, by = b_lon * m_per_deg_lon, b_lat * m_per_deg_lat
        vx, vy = bx - ax, by - ay
        norm = np.hypot(vx, vy)
        if norm <= 0:
            return None
        px, py = -vy / norm, vx / norm
        ox, oy = px * half_width_m, py * half_width_m
        p1l = ((ax + ox) / m_per_deg_lon, (ay + oy) / m_per_deg_lat)
        p2l = ((bx + ox) / m_per_deg_lon, (by + oy) / m_per_deg_lat)
        p2r = ((bx - ox) / m_per_deg_lon, (by - oy) / m_per_deg_lat)
        p1r = ((ax - ox) / m_per_deg_lon, (ay - oy) / m_per_deg_lat)
        return [p1l, p2l, p2r, p1r]

    def _backscatter_half_width_from_point(self, a_lat, a_lon, b_lat, b_lon, p_lat, p_lon):
        """Perpendicular distance (meters) from point P to centerline AB."""
        center_lat = (a_lat + b_lat + p_lat) / 3.0
        m_per_deg_lat = 111320.0
        m_per_deg_lon = 111320.0 * np.cos(np.radians(center_lat))
        if m_per_deg_lon == 0:
            return 0.0
        ax, ay = a_lon * m_per_deg_lon, a_lat * m_per_deg_lat
        bx, by = b_lon * m_per_deg_lon, b_lat * m_per_deg_lat
        px, py = p_lon * m_per_deg_lon, p_lat * m_per_deg_lat
        vx, vy = bx - ax, by - ay
        wx, wy = px - ax, py - ay
        norm = np.hypot(vx, vy)
        if norm <= 0:
            return 0.0
        return abs(vx * wy - vy * wx) / norm

    def _get_backscatter_centerline_with_lead_in(self):
        """Return saved centerline endpoints extended by Lead-in at both ends."""
        centerline = getattr(self, "backscatter_box_centerline", None)
        if not centerline or len(centerline) != 2 or centerline[0] is None or centerline[1] is None:
            return None
        (a_lat, a_lon), (b_lat, b_lon) = centerline
        lead_in_m = float(max(0.0, getattr(self, "backscatter_lead_in_m", 0.0) or 0.0))
        if lead_in_m <= 0.0:
            return (a_lat, a_lon), (b_lat, b_lon)
        geod = pyproj.Geod(ellps="WGS84")
        az_ab, az_ba, _ = geod.inv(a_lon, a_lat, b_lon, b_lat)
        a_lon_ext, a_lat_ext, _ = geod.fwd(a_lon, a_lat, az_ba, lead_in_m)
        b_lon_ext, b_lat_ext, _ = geod.fwd(b_lon, b_lat, az_ab, lead_in_m)
        return (a_lat_ext, a_lon_ext), (b_lat_ext, b_lon_ext)

    def _finalize_backscatter_box(self, vertices, centerline, half_width_m):
        """Finalize oriented box polygon, compute stats, redraw, and open stats dialog."""
        if not vertices or len(vertices) != 4:
            self._show_message("warning", "Backscatter Box", "Could not create oriented box.")
            return
        self.backscatter_box_vertices = vertices
        self.backscatter_box_centerline = centerline
        self.backscatter_box_half_width_m = float(max(0.0, half_width_m))
        self.backscatter_box_stats = self._compute_backscatter_box_stats()
        self._plot_survey_plan(preserve_view_limits=True)
        if hasattr(self, "_draw_current_profile"):
            self._draw_current_profile()
        if hasattr(self, "_update_backscatter_export_name_default"):
            self._update_backscatter_export_name_default()
        self._show_backscatter_stats_dialog()

    def _compute_backscatter_box_stats(self):
        """Compute histogram input and dimensions for current backscatter box."""
        if self.backscatter_box_vertices is None:
            return None
        extent = getattr(self, "geotiff_extent", None)
        if extent is None:
            return None
        vertices = np.array(self.backscatter_box_vertices, dtype=float)
        if vertices.shape != (4, 2):
            return None
        lon_min = float(np.min(vertices[:, 0]))
        lon_max = float(np.max(vertices[:, 0]))
        lat_min = float(np.min(vertices[:, 1]))
        lat_max = float(np.max(vertices[:, 1]))
        left, right, bottom, top = extent
        if right == left or top == bottom:
            return None
        geod = pyproj.Geod(ellps="WGS84")
        centerline = getattr(self, "backscatter_box_centerline", None)
        if centerline and len(centerline) == 2:
            (a_lat, a_lon), (b_lat, b_lon) = centerline
            _, _, height_m = geod.inv(a_lon, a_lat, b_lon, b_lat)
        else:
            height_m = 0.0
        width_m = 2.0 * float(getattr(self, "backscatter_box_half_width_m", 0.0) or 0.0)
        rectified = self._compute_backscatter_box_rectified_image()
        rectified_img = rectified.get("image")
        if rectified_img is not None and np.size(rectified_img) > 0:
            valid = np.array(rectified_img[np.isfinite(rectified_img)], dtype=float)
        else:
            valid = np.array([], dtype=float)
        bathy_rect = self._compute_rectified_slope_aspect_from_bathy()
        bathy_rectified = self._compute_rectified_box_image_from_source("bathy").get("image")
        return {
            "values": valid,
            "width_m": float(abs(width_m)),
            "height_m": float(abs(height_m)),
            "rectified_image": rectified_img,
            "rectified_long_m": float(rectified.get("long_m", 0.0)),
            "rectified_short_m": float(rectified.get("short_m", 0.0)),
            "bathy_rectified": bathy_rectified,
            "slope_rectified": bathy_rect.get("slope"),
        }

    def _get_backscatter_display_vmin_vmax(self, values):
        """Compute display scaling for backscatter using the same percent-clip rules as main map."""
        valid = np.isfinite(values)
        if not np.any(valid):
            return 0.0, 1.0
        arr = values[valid]
        if getattr(self, "backscatter_percent_clip_enabled", True):
            min_pct = float(getattr(self, "backscatter_percent_clip_min", 0.5))
            max_pct = float(getattr(self, "backscatter_percent_clip_max", 0.5))
            min_pct = max(0.0, min(49.9, min_pct))
            max_pct = max(0.0, min(49.9, max_pct))
            low_q = min_pct
            high_q = 100.0 - max_pct
            if high_q <= low_q:
                low_q = 0.5
                high_q = 99.5
            vmin = float(np.nanpercentile(arr, low_q))
            vmax = float(np.nanpercentile(arr, high_q))
        else:
            vmin = float(np.nanmin(arr))
            vmax = float(np.nanmax(arr))
        if vmin >= vmax:
            vmin = float(np.nanmin(arr))
            vmax = float(np.nanmax(arr))
        if vmin >= vmax:
            vmax = vmin + 1e-6
        return vmin, vmax

    def _sample_dataset_band_on_wgs84_grid(self, src_ds, band_index, sample_lon, sample_lat, apply_bathy_nan=False):
        """Sample a rasterio dataset band at WGS84 lon/lat points and return image-shaped array."""
        out = np.full(sample_lon.shape, np.nan, dtype=float)
        lon_flat = np.asarray(sample_lon, dtype=float).ravel()
        lat_flat = np.asarray(sample_lat, dtype=float).ravel()
        if lon_flat.size == 0:
            return out

        if src_ds.crs is not None and str(src_ds.crs) != "EPSG:4326":
            transformer = pyproj.Transformer.from_crs("EPSG:4326", src_ds.crs, always_xy=True)
            xs, ys = transformer.transform(lon_flat, lat_flat)
        else:
            xs, ys = lon_flat, lat_flat

        rows, cols = rowcol(src_ds.transform, xs, ys)
        rows = np.asarray(rows, dtype=int)
        cols = np.asarray(cols, dtype=int)
        valid = (
            np.isfinite(xs)
            & np.isfinite(ys)
            & (rows >= 0)
            & (rows < int(src_ds.height))
            & (cols >= 0)
            & (cols < int(src_ds.width))
        )
        if not np.any(valid):
            return out

        valid_rows = rows[valid]
        valid_cols = cols[valid]
        rmin = int(np.min(valid_rows))
        rmax = int(np.max(valid_rows))
        cmin = int(np.min(valid_cols))
        cmax = int(np.max(valid_cols))
        window = Window.from_slices((rmin, rmax + 1), (cmin, cmax + 1))
        data = src_ds.read(int(band_index), window=window).astype(float)

        local_rows = valid_rows - rmin
        local_cols = valid_cols - cmin
        sampled = data[local_rows, local_cols]
        out_flat = out.ravel()
        out_flat[np.where(valid)[0]] = sampled

        if src_ds.nodata is not None:
            out[np.isclose(out, float(src_ds.nodata), rtol=0.0, atol=1e-9)] = np.nan
        if apply_bathy_nan:
            self._apply_geotiff_nan_filter(out)
        return out

    def _get_source_pixel_size_m(self, source_key, center_lat):
        """Estimate source pixel size in meters for full-resolution stats sampling."""
        if source_key == "bathy":
            ds = getattr(self, "geotiff_dataset_original", None)
            if ds is None:
                return 5.0
            tr = ds.transform
            res_x = abs(float(tr.a))
            res_y = abs(float(tr.e))
            if ds.crs is not None and str(ds.crs) != "EPSG:4326":
                return max(0.25, min(res_x, res_y))
            meters_per_deg_lon = max(1.0, 111320.0 * np.cos(np.radians(center_lat)))
            meters_per_deg_lat = 111320.0
            return max(0.25, min(res_x * meters_per_deg_lon, res_y * meters_per_deg_lat))

        path = getattr(self, "backscatter_geotiff_path", None)
        if not path or not os.path.isfile(path):
            return 5.0
        try:
            with rasterio.open(path) as ds:
                tr = ds.transform
                res_x = abs(float(tr.a))
                res_y = abs(float(tr.e))
                if ds.crs is not None and str(ds.crs) != "EPSG:4326":
                    return max(0.25, min(res_x, res_y))
                meters_per_deg_lon = max(1.0, 111320.0 * np.cos(np.radians(center_lat)))
                meters_per_deg_lat = 111320.0
                return max(0.25, min(res_x * meters_per_deg_lon, res_y * meters_per_deg_lat))
        except Exception:
            return 5.0

    def _compute_rectified_box_image_from_source(self, source_key):
        """Return de-rotated image using full-resolution source raster(s), independent of display grid."""
        centerline = getattr(self, "backscatter_box_centerline", None)
        half_width = float(getattr(self, "backscatter_box_half_width_m", 0.0) or 0.0)
        if (
            centerline is None
            or len(centerline) != 2
            or centerline[0] is None
            or centerline[1] is None
            or half_width <= 0.0
        ):
            return {"image": None, "long_m": 0.0, "short_m": 0.0}

        (a_lat, a_lon), (b_lat, b_lon) = centerline
        geod = pyproj.Geod(ellps="WGS84")
        _, _, centerline_len_m = geod.inv(a_lon, a_lat, b_lon, b_lat)
        width_m = 2.0 * half_width
        long_m = float(max(centerline_len_m, width_m))
        short_m = float(min(centerline_len_m, width_m))
        if long_m <= 0.0 or short_m <= 0.0:
            return {"image": None, "long_m": 0.0, "short_m": 0.0}

        center_lat = 0.5 * (a_lat + b_lat)
        center_lon = 0.5 * (a_lon + b_lon)
        m_per_deg_lat = 111320.0
        m_per_deg_lon = 111320.0 * np.cos(np.radians(center_lat))
        if abs(m_per_deg_lon) < 1e-9:
            return {"image": None, "long_m": long_m, "short_m": short_m}

        ax_m = np.array([a_lon * m_per_deg_lon, a_lat * m_per_deg_lat], dtype=float)
        bx_m = np.array([b_lon * m_per_deg_lon, b_lat * m_per_deg_lat], dtype=float)
        c_m = 0.5 * (ax_m + bx_m)
        ab_vec = bx_m - ax_m
        ab_norm = float(np.hypot(ab_vec[0], ab_vec[1]))
        if ab_norm <= 0.0:
            return {"image": None, "long_m": long_m, "short_m": short_m}
        u_ab = ab_vec / ab_norm
        v_ab = np.array([-u_ab[1], u_ab[0]], dtype=float)
        width_point = getattr(self, "backscatter_box_width_point", None)
        if width_point is not None and len(width_point) == 2:
            wp_lat, wp_lon = width_point
            wp_m = np.array([wp_lon * m_per_deg_lon, wp_lat * m_per_deg_lat], dtype=float)
            side_sign = float(np.dot(wp_m - c_m, v_ab))
            if side_sign < 0.0:
                v_ab = -v_ab

        # Deterministic convention:
        # - If centerline is longer, +X is first centerline click -> second centerline click.
        # - If width is longer, +X is toward the side used for width-click definition.
        if centerline_len_m >= width_m:
            u_long = u_ab
            v_short = v_ab
        else:
            u_long = v_ab
            v_short = u_ab

        base_res_m = self._get_source_pixel_size_m(source_key, center_lat)
        nx = int(np.clip(np.ceil(long_m / base_res_m) + 1, 32, 2400))
        ny = int(np.clip(np.ceil(short_m / base_res_m) + 1, 16, 1200))

        x_vals = np.linspace(-0.5 * long_m, 0.5 * long_m, nx)
        y_vals = np.linspace(-0.5 * short_m, 0.5 * short_m, ny)
        xx, yy = np.meshgrid(x_vals, y_vals)
        xy_m = c_m.reshape(1, 1, 2) + xx[..., None] * u_long.reshape(1, 1, 2) + yy[..., None] * v_short.reshape(1, 1, 2)
        sample_lon = xy_m[..., 0] / m_per_deg_lon
        sample_lat = xy_m[..., 1] / m_per_deg_lat

        if source_key == "bathy":
            ds = getattr(self, "geotiff_dataset_original", None)
            if ds is None:
                return {"image": None, "long_m": long_m, "short_m": short_m}
            image = self._sample_dataset_band_on_wgs84_grid(
                ds,
                1,
                sample_lon,
                sample_lat,
                apply_bathy_nan=True,
            )
            return {"image": image, "long_m": long_m, "short_m": short_m}

        path = getattr(self, "backscatter_geotiff_path", None)
        if not path or not os.path.isfile(path):
            return {"image": None, "long_m": long_m, "short_m": short_m}
        band_index = int(getattr(self, "backscatter_selected_band", 1) or 1)
        try:
            with rasterio.open(path) as src:
                band_index = max(1, min(int(getattr(src, "count", 1)), band_index))
                image = self._sample_dataset_band_on_wgs84_grid(
                    src,
                    band_index,
                    sample_lon,
                    sample_lat,
                    apply_bathy_nan=False,
                )
        except Exception:
            return {"image": None, "long_m": long_m, "short_m": short_m}
        return {"image": image, "long_m": long_m, "short_m": short_m}

    def _compute_backscatter_box_rectified_image(self):
        """Return a de-rotated backscatter image from full-resolution source raster."""
        return self._compute_rectified_box_image_from_source("backscatter")

    def _compute_rectified_slope_aspect_from_bathy(self):
        """Compute slope and aspect from bathymetry in the selected rectified box frame."""
        rect = self._compute_rectified_box_image_from_source("bathy")
        bathy = rect.get("image")
        long_m = float(rect.get("long_m", 0.0))
        short_m = float(rect.get("short_m", 0.0))
        if bathy is None or np.size(bathy) == 0:
            return {"slope": None, "aspect": None}
        ny, nx = bathy.shape
        if nx < 2 or ny < 2 or long_m <= 0.0 or short_m <= 0.0:
            return {"slope": None, "aspect": None}

        dx_m = long_m / max(1, nx - 1)
        dy_m = short_m / max(1, ny - 1)
        z = np.array(bathy, dtype=float)
        z_filled = np.where(np.isfinite(z), z, 0.0)
        dz_dy, dz_dx = np.gradient(z_filled, dy_m, dx_m)
        slope_deg = np.degrees(np.arctan(np.sqrt(dz_dx ** 2 + dz_dy ** 2)))
        aspect_deg = (np.degrees(np.arctan2(dz_dy, -dz_dx)) + 360.0) % 360.0
        invalid = ~np.isfinite(z)
        slope_deg[invalid] = np.nan
        aspect_deg[invalid] = np.nan
        return {"slope": slope_deg, "aspect": aspect_deg}

    def _show_backscatter_stats_dialog(self):
        """Create/show Backscatter Statistics dialog with histogram and box dimensions."""
        stats = getattr(self, "backscatter_box_stats", None)
        if not stats:
            return
        if self.backscatter_stats_dialog is None:
            self.backscatter_stats_dialog = QDialog(self)
            self.backscatter_stats_dialog.setWindowTitle("Backscatter Statistics")
            self.backscatter_stats_dialog.resize(1100, 760)
            layout = QVBoxLayout(self.backscatter_stats_dialog)
            top_row = QWidget()
            top_layout = QVBoxLayout(top_row)
            top_layout.setContentsMargins(0, 0, 0, 0)
            fig = Figure(figsize=(6.8, 6.4))
            self.backscatter_stats_img_ax = fig.add_subplot(221)
            self.backscatter_stats_ax = fig.add_subplot(222)
            self.backscatter_stats_aspect_ax = fig.add_subplot(223)
            self.backscatter_stats_slope_ax = fig.add_subplot(224)
            fig.subplots_adjust(wspace=0.26, hspace=0.32, bottom=0.09, top=0.94, left=0.07, right=0.98)
            self.backscatter_stats_canvas = FigureCanvas(fig)
            self.backscatter_stats_motion_cid = self.backscatter_stats_canvas.mpl_connect(
                "motion_notify_event", self._on_backscatter_stats_motion
            )
            top_layout.addWidget(self.backscatter_stats_canvas)
            layout.addWidget(top_row)
            self.backscatter_stats_text_label = QLabel("")
            self.backscatter_stats_text_label.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft)
            self.backscatter_stats_text_label.setWordWrap(True)
            layout.addWidget(self.backscatter_stats_text_label)

        vals = stats.get("values", np.array([]))
        rectified_img = stats.get("rectified_image", None)
        long_m = float(stats.get("rectified_long_m", 0.0))
        short_m = float(stats.get("rectified_short_m", 0.0))
        slope_rectified = stats.get("slope_rectified", None)
        bathy_rectified = stats.get("bathy_rectified", None)
        self._backscatter_stats_rectified_img = rectified_img
        self._backscatter_stats_slope_rectified = slope_rectified
        self._backscatter_stats_bathy_rectified = bathy_rectified
        self._backscatter_stats_long_m = long_m
        self._backscatter_stats_short_m = short_m
        self.backscatter_stats_img_ax.clear()
        if rectified_img is not None and np.size(rectified_img) > 0 and np.any(np.isfinite(rectified_img)):
            vmin, vmax = self._get_backscatter_display_vmin_vmax(rectified_img)
            self.backscatter_stats_img_ax.imshow(
                rectified_img,
                cmap="gray",
                origin="lower",
                aspect="auto",
                extent=(0.0, long_m, 0.0, short_m),
                vmin=vmin,
                vmax=vmax,
                interpolation="nearest",
            )
            self.backscatter_stats_img_ax.set_title("Backscatter Area (Long Axis × Short Axis)")
            self.backscatter_stats_img_ax.set_xlabel("Long Direction (m)")
            self.backscatter_stats_img_ax.set_ylabel("Short Direction (m)")
            self.backscatter_stats_img_hover_text = self.backscatter_stats_img_ax.text(
                0.02, 0.98, "Value: -",
                transform=self.backscatter_stats_img_ax.transAxes,
                ha="left", va="top", fontsize=8,
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
                zorder=20,
            )
        else:
            self.backscatter_stats_img_ax.text(
                0.5,
                0.5,
                "No valid backscatter image for selected area.",
                transform=self.backscatter_stats_img_ax.transAxes,
                ha="center",
                va="center",
            )
            self.backscatter_stats_img_ax.set_title("Backscatter Area (Long Axis × Short Axis)")
            self.backscatter_stats_img_ax.set_xlabel("Long Direction (m)")
            self.backscatter_stats_img_ax.set_ylabel("Short Direction (m)")
            self.backscatter_stats_img_hover_text = None

        self.backscatter_stats_slope_ax.clear()
        if slope_rectified is not None and np.size(slope_rectified) > 0 and np.any(np.isfinite(slope_rectified)):
            self.backscatter_stats_slope_ax.imshow(
                slope_rectified,
                cmap="inferno",
                origin="lower",
                aspect="auto",
                extent=(0.0, long_m, 0.0, short_m),
                interpolation="nearest",
            )
            self.backscatter_stats_slope_ax.set_title("Bathymetry Slope (deg)")
            self.backscatter_stats_slope_ax.set_xlabel("Long Direction (m)")
            self.backscatter_stats_slope_ax.set_ylabel("Short Direction (m)")
            self.backscatter_stats_slope_hover_text = self.backscatter_stats_slope_ax.text(
                0.02, 0.98, "Value: -",
                transform=self.backscatter_stats_slope_ax.transAxes,
                ha="left", va="top", fontsize=8,
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
                zorder=20,
            )
        else:
            self.backscatter_stats_slope_ax.text(
                0.5, 0.5, "No valid slope map.",
                transform=self.backscatter_stats_slope_ax.transAxes,
                ha="center", va="center",
            )
            self.backscatter_stats_slope_ax.set_title("Bathymetry Slope (deg)")
            self.backscatter_stats_slope_ax.set_xlabel("Long Direction (m)")
            self.backscatter_stats_slope_ax.set_ylabel("Short Direction (m)")
            self.backscatter_stats_slope_hover_text = None

        self.backscatter_stats_aspect_ax.clear()
        if bathy_rectified is not None and np.size(bathy_rectified) > 0 and np.any(np.isfinite(bathy_rectified)):
            bmin = float(np.nanmin(bathy_rectified))
            bmax = float(np.nanmax(bathy_rectified))
            if bmin >= bmax:
                bmax = bmin + 1e-6
            self.backscatter_stats_aspect_ax.imshow(
                bathy_rectified,
                cmap="terrain",
                origin="lower",
                aspect="auto",
                extent=(0.0, long_m, 0.0, short_m),
                vmin=bmin,
                vmax=bmax,
                interpolation="nearest",
            )
            self.backscatter_stats_aspect_ax.set_title("Bathymetry (m)")
            self.backscatter_stats_aspect_ax.set_xlabel("Long Direction (m)")
            self.backscatter_stats_aspect_ax.set_ylabel("Short Direction (m)")
            self.backscatter_stats_bathy_hover_text = self.backscatter_stats_aspect_ax.text(
                0.02, 0.98, "Value: -",
                transform=self.backscatter_stats_aspect_ax.transAxes,
                ha="left", va="top", fontsize=8,
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
                zorder=20,
            )
        else:
            self.backscatter_stats_aspect_ax.text(
                0.5, 0.5, "No valid bathymetry map.",
                transform=self.backscatter_stats_aspect_ax.transAxes,
                ha="center", va="center",
            )
            self.backscatter_stats_aspect_ax.set_title("Bathymetry (m)")
            self.backscatter_stats_aspect_ax.set_xlabel("Long Direction (m)")
            self.backscatter_stats_aspect_ax.set_ylabel("Short Direction (m)")
            self.backscatter_stats_bathy_hover_text = None

        self.backscatter_stats_ax.clear()
        if vals.size > 0:
            counts, bin_edges, _ = self.backscatter_stats_ax.hist(
                vals, bins=40, color="steelblue", edgecolor="black", alpha=0.85
            )
            self._backscatter_hist_counts = np.array(counts, dtype=float)
            self._backscatter_hist_edges = np.array(bin_edges, dtype=float)
            self.backscatter_stats_ax.set_title("Backscatter Value Distribution")
            self.backscatter_stats_ax.set_xlabel("Backscatter Value")
            self.backscatter_stats_ax.set_ylabel("Count")
            self.backscatter_stats_hist_hover_text = self.backscatter_stats_ax.text(
                0.02, 0.98, "Bin: -\nCount: -",
                transform=self.backscatter_stats_ax.transAxes,
                ha="left", va="top", fontsize=8,
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
                zorder=20,
            )
        else:
            self._backscatter_hist_counts = None
            self._backscatter_hist_edges = None
            self.backscatter_stats_ax.text(
                0.5,
                0.5,
                "No valid backscatter values in selected box.",
                transform=self.backscatter_stats_ax.transAxes,
                ha="center",
                va="center",
            )
            self.backscatter_stats_ax.set_title("Backscatter Value Distribution")
            self.backscatter_stats_hist_hover_text = None
        self.backscatter_stats_ax.figure.subplots_adjust(wspace=0.26, hspace=0.32, bottom=0.09, top=0.94, left=0.07, right=0.98)
        self.backscatter_stats_canvas.draw_idle()

        width_m = float(stats.get("width_m", 0.0))
        height_m = float(stats.get("height_m", 0.0))
        width_km = width_m / 1000.0
        height_km = height_m / 1000.0
        width_nm = width_m / 1852.0
        height_nm = height_m / 1852.0
        text = (
            f"Width: {width_m:.2f} m | {width_km:.4f} km | {width_nm:.4f} nm\n"
            f"Height: {height_m:.2f} m | {height_km:.4f} km | {height_nm:.4f} nm"
        )
        self.backscatter_stats_text_label.setText(text)
        self.backscatter_stats_dialog.show()
        self.backscatter_stats_dialog.raise_()
        self.backscatter_stats_dialog.activateWindow()

    def _save_backscatter_stats_png_without_hover(self, path, dpi=300):
        """Write stats figure to PNG without hover readout boxes (matplotlib text artists)."""
        canvas = getattr(self, "backscatter_stats_canvas", None)
        if canvas is None:
            return
        hover_artists = [
            getattr(self, "backscatter_stats_img_hover_text", None),
            getattr(self, "backscatter_stats_slope_hover_text", None),
            getattr(self, "backscatter_stats_bathy_hover_text", None),
            getattr(self, "backscatter_stats_hist_hover_text", None),
        ]
        prev_vis = []
        for t in hover_artists:
            if t is not None:
                prev_vis.append((t, t.get_visible()))
                t.set_visible(False)
            else:
                prev_vis.append(None)
        try:
            fig = canvas.figure
            canvas.draw()
            fig.savefig(path, dpi=dpi)
        finally:
            for entry in prev_vis:
                if entry is not None:
                    t, vis = entry
                    t.set_visible(vis)
            canvas.draw_idle()

    def _on_backscatter_stats_motion(self, event):
        """Show upper-left value readout for backscatter/slope/bathymetry stats maps."""
        long_m = float(getattr(self, "_backscatter_stats_long_m", 0.0) or 0.0)
        short_m = float(getattr(self, "_backscatter_stats_short_m", 0.0) or 0.0)
        if long_m <= 0.0 or short_m <= 0.0:
            return

        def _update_for_axis(target_ax, arr, text_artist, fmt):
            if target_ax is None or text_artist is None or arr is None:
                return False
            if event.inaxes != target_ax or event.xdata is None or event.ydata is None:
                return False
            if np.size(arr) == 0:
                text_artist.set_text("Value: -")
                return True
            ny, nx = arr.shape
            x = float(np.clip(event.xdata, 0.0, long_m))
            y = float(np.clip(event.ydata, 0.0, short_m))
            ix = int(np.clip(round((x / max(long_m, 1e-12)) * (nx - 1)), 0, nx - 1))
            iy = int(np.clip(round((y / max(short_m, 1e-12)) * (ny - 1)), 0, ny - 1))
            val = arr[iy, ix]
            text_artist.set_text("Value: -" if not np.isfinite(val) else fmt.format(float(val)))
            return True

        handled = False
        handled |= _update_for_axis(
            getattr(self, "backscatter_stats_img_ax", None),
            getattr(self, "_backscatter_stats_rectified_img", None),
            getattr(self, "backscatter_stats_img_hover_text", None),
            "Value: {:.3f}",
        )
        handled |= _update_for_axis(
            getattr(self, "backscatter_stats_slope_ax", None),
            getattr(self, "_backscatter_stats_slope_rectified", None),
            getattr(self, "backscatter_stats_slope_hover_text", None),
            "Value: {:.2f} deg",
        )
        handled |= _update_for_axis(
            getattr(self, "backscatter_stats_aspect_ax", None),
            getattr(self, "_backscatter_stats_bathy_rectified", None),
            getattr(self, "backscatter_stats_bathy_hover_text", None),
            "Value: {:.2f} m",
        )
        if (
            event.inaxes == getattr(self, "backscatter_stats_ax", None)
            and hasattr(self, "backscatter_stats_hist_hover_text")
            and self.backscatter_stats_hist_hover_text is not None
            and getattr(self, "_backscatter_hist_counts", None) is not None
            and getattr(self, "_backscatter_hist_edges", None) is not None
            and event.xdata is not None
            and event.ydata is not None
        ):
            counts = self._backscatter_hist_counts
            edges = self._backscatter_hist_edges
            idx = int(np.searchsorted(edges, event.xdata, side="right") - 1)
            if 0 <= idx < len(counts):
                lo = float(edges[idx])
                hi = float(edges[idx + 1])
                count = int(round(float(counts[idx])))
                if 0.0 <= event.ydata <= float(max(counts[idx], 0.0)):
                    self.backscatter_stats_hist_hover_text.set_text(
                        f"Bin: {lo:.3f} to {hi:.3f}\nCount: {count}"
                    )
                else:
                    self.backscatter_stats_hist_hover_text.set_text("Bin: -\nCount: -")
            else:
                self.backscatter_stats_hist_hover_text.set_text("Bin: -\nCount: -")
            handled = True
        if handled and hasattr(self, "backscatter_stats_canvas") and self.backscatter_stats_canvas is not None:
            self.backscatter_stats_canvas.draw_idle()

    def _zoom_to_geotiff(self):
        """Zooms the plot to the initial GeoTIFF bounds (same as when first loaded)."""
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded. Cannot zoom to GeoTIFF.")
            return

        if self.geotiff_extent is None:
            self._show_message("warning", "No GeoTIFF", "No GeoTIFF underlay is loaded to zoom to.")
            return

        # Use stored initial limits if available, otherwise calculate them
        if hasattr(self, 'initial_geotiff_xlim') and self.initial_geotiff_xlim is not None and \
           hasattr(self, 'initial_geotiff_ylim') and self.initial_geotiff_ylim is not None:
            # Restore to the exact initial limits when GeoTIFF was first loaded
            self.current_xlim = self.initial_geotiff_xlim
            self.current_ylim = self.initial_geotiff_ylim
        else:
            # Fallback: calculate consistent plot limits (same method used when GeoTIFF is first loaded)
            self.current_xlim, self.current_ylim = self._calculate_consistent_plot_limits()
            # Store them for future use
            self.initial_geotiff_xlim = self.current_xlim
            self.initial_geotiff_ylim = self.current_ylim

        # Store these as the user's intended view limits BEFORE setting axes
        self._last_user_xlim = self.current_xlim
        self._last_user_ylim = self.current_ylim

        # Set axis limits FIRST so _load_geotiff_at_resolution can see them
        self.ax.set_xlim(self.current_xlim)
        self.ax.set_ylim(self.current_ylim)

        # Calculate zoom level based on the view extent vs original full GeoTIFF extent
        if self.dynamic_resolution_enabled:
            # Use original extent for comparison, not current extent (which may be a subset)
            if hasattr(self, 'geotiff_original_extent') and self.geotiff_original_extent is not None:
                full_extent = self.geotiff_original_extent
            elif hasattr(self, 'geotiff_extent') and self.geotiff_extent is not None:
                full_extent = self.geotiff_extent
            else:
                full_extent = None

            if full_extent is not None:
                full_width = full_extent[1] - full_extent[0]
                full_height = full_extent[3] - full_extent[2]
                current_width = self.current_xlim[1] - self.current_xlim[0]
                current_height = self.current_ylim[1] - self.current_ylim[0]

                width_ratio = current_width / full_width if full_width > 0 else 1.0
                height_ratio = current_height / full_height if full_height > 0 else 1.0
                # Use the larger ratio (more zoomed in) to determine resolution
                new_zoom_level = max(width_ratio, height_ratio)
                new_zoom_level = max(0.01, min(10.0, new_zoom_level))
                self.geotiff_zoom_level = new_zoom_level
            else:
                # If no extent available, use full resolution
                self.geotiff_zoom_level = 1.0
        else:
            # If dynamic resolution is disabled, use full resolution
            self.geotiff_zoom_level = 1.0

        # Reload GeoTIFF data at appropriate resolution for the view
        # This now uses the axis limits we just set
        if self.geotiff_dataset_original is not None:
            success = self._load_geotiff_at_resolution()
            if not success:
                self._show_message("warning", "GeoTIFF Reload", "Could not reload GeoTIFF data for the new view.")

        # Replot everything with preserve_view_limits=True to keep our limits
        # This will properly display the GeoTIFF and preserve all existing lines
        self._plot_survey_plan(preserve_view_limits=True)
        if hasattr(self, "_on_geotiff_loaded_line_info_depth"):
            self._on_geotiff_loaded_line_info_depth()
        if hasattr(self, "_update_roll_line_info_labels"):
            self._update_roll_line_info_labels()

    def _calculate_slope_at_point(self, lat, lon):
        """Calculate slope at a given lat/lon point."""
        if (self.geotiff_data_array is None or self.geotiff_extent is None):
            return None, None

        try:
            left, right, bottom, top = tuple(self.geotiff_extent)
            nrows, ncols = self.geotiff_data_array.shape

            col = int((lon - left) / (right - left) * (ncols - 1))
            row = int((top - lat) / (top - bottom) * (nrows - 1))

            if 0 <= row < nrows and 0 <= col < ncols:
                elevation = self.geotiff_data_array[row, col]

                center_lat_geotiff = (self.geotiff_extent[2] + self.geotiff_extent[3]) / 2
                m_per_deg_lat = 111320.0
                m_per_deg_lon = 111320.0 * np.cos(np.radians(center_lat_geotiff))

                res_lat_deg = (self.geotiff_extent[3] - self.geotiff_extent[2]) / self.geotiff_data_array.shape[0]
                res_lon_deg = (self.geotiff_extent[1] - self.geotiff_extent[0]) / self.geotiff_data_array.shape[1]

                dx_m = res_lon_deg * m_per_deg_lon
                dy_m = res_lat_deg * m_per_deg_lat

                window_size = 3
                r_start = max(0, row - window_size // 2)
                r_end = min(nrows, row + window_size // 2 + 1)
                c_start = max(0, col - window_size // 2)
                c_end = min(ncols, col + window_size // 2 + 1)

                window_data = self.geotiff_data_array[r_start:r_end, c_start:c_end]

                if window_data.size > 1 and not np.all(np.isnan(window_data)):
                    dz_dy, dz_dx = np.gradient(window_data, dy_m, dx_m)
                    slope_rad = np.arctan(np.sqrt(dz_dx**2 + dz_dy**2))
                    slope_degrees = np.degrees(slope_rad)
                    center_idx = window_size // 2
                    if center_idx < slope_degrees.shape[0] and center_idx < slope_degrees.shape[1]:
                        return elevation, slope_degrees[center_idx, center_idx]
                    return elevation, slope_degrees[0, 0]
                return elevation, None
            return None, None
        except Exception as e:
            print(f"Error calculating slope: {e}")
            return None, None

    def _get_depth_at_point(self, lat, lon):
        """Get depth/elevation at a given lat/lon point."""
        if (self.geotiff_data_array is None or self.geotiff_extent is None):
            return 0.0

        try:
            left, right, bottom, top = tuple(self.geotiff_extent)
            nrows, ncols = self.geotiff_data_array.shape

            col = int((lon - left) / (right - left) * (ncols - 1))
            row = int((top - lat) / (top - bottom) * (nrows - 1))

            if 0 <= row < nrows and 0 <= col < ncols:
                elevation = self.geotiff_data_array[row, col]
                return elevation if not np.isnan(elevation) else 0.0
            return 0.0
        except Exception as e:
            print(f"Error getting depth: {e}")
            return 0.0

    def _on_draw_event_update_colormap(self, event=None):
        """Update elevation colormap to visible region when view changes."""
        if not hasattr(self, 'geotiff_data_array') or self.geotiff_data_array is None:
            return
        if not hasattr(self, 'geotiff_extent') or self.geotiff_extent is None:
            return
        if not hasattr(self, 'geotiff_display_mode') or self.geotiff_display_mode != 'elevation':
            return
        if not hasattr(self, 'geotiff_image_plot') or self.geotiff_image_plot is None:
            return
        display_data = self.geotiff_data_array
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        left, right, bottom, top = tuple(self.geotiff_extent)
        nrows, ncols = display_data.shape
        col_min = int(np.clip((min(xlim) - left) / (right - left) * (ncols - 1), 0, ncols - 1))
        col_max = int(np.clip((max(xlim) - left) / (right - left) * (ncols - 1), 0, ncols - 1))
        row_min = int(np.clip((top - max(ylim)) / (top - bottom) * (nrows - 1), 0, nrows - 1))
        row_max = int(np.clip((top - min(ylim)) / (top - bottom) * (nrows - 1), 0, nrows - 1))
        r0, r1 = sorted([row_min, row_max])
        c0, c1 = sorted([col_min, col_max])
        visible_region = display_data[r0:r1+1, c0:c1+1]
        if visible_region.size > 0 and not np.all(np.isnan(visible_region)):
            vmin_elev = np.nanmin(visible_region)
            vmax_elev = np.nanmax(visible_region)
        else:
            vmin_elev = np.nanmin(display_data) if not np.all(np.isnan(display_data)) else None
            vmax_elev = np.nanmax(display_data) if not np.all(np.nan(display_data)) else None
        if vmin_elev is not None and vmax_elev is not None and vmin_elev != vmax_elev:
            self.geotiff_image_plot.set_clim(vmin=vmin_elev, vmax=vmax_elev)

    def _on_temp_line_motion(self, event):
        """Draw temporary line from start point to current mouse; show length, azimuth, time."""
        if not hasattr(self, '_temp_line') or not hasattr(self, '_temp_line_start'):
            return
        if event.inaxes != self.ax:
            return
        x0, y0 = self._temp_line_start[1], self._temp_line_start[0]
        x1, y1 = event.xdata, event.ydata
        if x1 is None or y1 is None:
            return
        self._temp_line.set_data([x0, x1], [y0, y1])
        try:
            geod = pyproj.Geod(ellps="WGS84")
            az12, az21, dist = geod.inv(x0, y0, x1, y1)
        except Exception:
            dist = np.sqrt((x1 - x0)**2 + (y1 - y0)**2) * 111320
            az12 = np.degrees(np.arctan2(x1 - x0, y1 - y0)) % 360
        speed_knots = None
        try:
            if hasattr(self, 'cal_survey_speed_entry') and self.cal_survey_speed_entry.isVisible():
                speed_knots = float(self.cal_survey_speed_entry.text())
            elif hasattr(self, 'survey_speed_entry'):
                speed_knots = float(self.survey_speed_entry.text())
        except Exception:
            speed_knots = None
        msg = f"Length: {dist:.1f} m, Azimuth: {az12:.1f}°"
        if speed_knots and speed_knots > 0:
            speed_m_per_s = speed_knots * 1852 / 3600.0
            time_sec = dist / speed_m_per_s
            if time_sec < 60:
                time_str = f"{time_sec:.1f} s"
            elif time_sec < 3600:
                time_str = f"{time_sec/60:.1f} min"
            else:
                time_str = f"{time_sec/3600:.2f} hr"
            msg += f", Time: {time_str}"
        self.canvas.draw_idle()

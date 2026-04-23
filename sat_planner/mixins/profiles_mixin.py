"""
Profile plots: elevation/slope along crossline, pitch line, roll line, line planning, performance.
_get_profile_data_from_geotiff, _draw_crossline_profile, _draw_pitch_line_profile,
_draw_roll_line_profile, _draw_current_profile.
"""
import numpy as np

from sat_planner.constants import pyproj


class ProfilesMixin:
    """Mixin for profile plots: _get_profile_data_from_geotiff, _draw_crossline_profile,
    _draw_pitch_line_profile, _draw_roll_line_profile, _draw_current_profile."""

    def _show_slope_on_profile(self):
        """Whether elevation profile plots should include the slope (twin axis) curve."""
        cb = getattr(self, "slope_profile_checkbox", None)
        if cb is not None:
            return cb.isChecked()
        return bool(getattr(self, "show_slope_profile_var", True))

    def _calibration_profile_select_is_roll(self):
        """True if Calibration tab profile dropdown is set to Roll (index 1)."""
        combo = getattr(self, "cal_profile_select_combo", None)
        if combo is None:
            return False
        return combo.currentIndex() == 1

    def _on_calibration_profile_select_changed(self, *_args):
        if not hasattr(self, "param_notebook"):
            return
        if self.param_notebook.currentIndex() != 0:
            return
        self._draw_current_profile()

    def _on_accuracy_profile_select_changed(self, *_args):
        if not hasattr(self, "param_notebook"):
            return
        if self.param_notebook.currentIndex() != 1:
            return
        self._draw_current_profile()

    def _sync_accuracy_profile_select_options(self):
        """Refresh Accuracy profile dropdown options from current generated lines."""
        combo = getattr(self, "ref_profile_select_combo", None)
        if combo is None:
            return
        previous_text = combo.currentText().strip() if combo.count() > 0 else "Crossline"
        options = []
        if hasattr(self, "cross_line_data") and self.cross_line_data and len(self.cross_line_data) == 2:
            options.append("Crossline")
        elif not options:
            options.append("Crossline")
        main_lines = getattr(self, "survey_lines_data", []) or []
        for i in range(len(main_lines)):
            options.append(f"Main Line {i + 1}")

        combo.blockSignals(True)
        combo.clear()
        combo.addItems(options)
        if previous_text in options:
            combo.setCurrentText(previous_text)
        else:
            combo.setCurrentIndex(0)
        combo.blockSignals(False)

    def _draw_accuracy_selected_profile(self):
        """Draw accuracy profile for selected Crossline or Main Line option."""
        self._sync_accuracy_profile_select_options()
        combo = getattr(self, "ref_profile_select_combo", None)
        selection = combo.currentText().strip() if combo is not None else "Crossline"
        segment = None
        title = "Crossline Elevation Profile"
        color = "darkorchid"

        if selection.startswith("Main Line"):
            try:
                idx = int(selection.replace("Main Line", "").strip()) - 1
            except Exception:
                idx = -1
            lines = getattr(self, "survey_lines_data", []) or []
            if 0 <= idx < len(lines) and len(lines[idx]) == 2:
                segment = lines[idx]
                title = f"Main Line {idx + 1} Elevation Profile"
                color = "blue"
        else:
            cross = getattr(self, "cross_line_data", []) or []
            if len(cross) == 2:
                segment = cross

        self._draw_segment_profile(segment, title, color)

    def _draw_segment_profile(self, segment, title, color):
        """Draw profile (elevation+slope) for a 2-point segment."""
        self.profile_ax.clear()
        for ax in self.profile_fig.get_axes():
            if ax != self.profile_ax:
                ax.remove()
        self.profile_ax.set_title(title, fontsize=8)
        self.profile_ax.set_xlabel("Distance (m)", fontsize=8)
        self.profile_ax.set_ylabel("Elevation (m)", fontsize=8)
        self.profile_ax.tick_params(axis='both', which='major', labelsize=7)
        self._profile_dists = None
        self._profile_elevations = None
        self._profile_slopes = None

        if not (
            self.geotiff_data_array is not None and
            self.geotiff_extent is not None and
            segment and
            len(segment) == 2
        ):
            self.profile_canvas.draw_idle()
            return

        (lat1, lon1), (lat2, lon2) = segment
        dists, elevations, slopes = self._profile_arrays_along_segment_endpoints(lat1, lon1, lat2, lon2, n=100)
        if dists is None or elevations is None:
            self.profile_canvas.draw_idle()
            return

        self._profile_dists = dists
        self._profile_elevations = elevations
        self._profile_slopes = slopes
        self.profile_ax.plot(dists, elevations, color=color, lw=1, label='Elevation')
        slope_ax = None
        if self._show_slope_on_profile():
            slope_ax = self.profile_ax.twinx()
            slope_ax.plot(dists, slopes, color='teal', lw=1, linestyle='--', label='Slope (deg)')
            slope_ax.set_ylabel('Slope (deg)', fontsize=8)
            slope_ax.tick_params(axis='y', labelsize=7)
            slope_ax.grid(False)
        self.profile_ax.set_xlim(dists[0], dists[-1] if len(dists) > 0 else 1)
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

    def _get_profile_data_from_geotiff(self, lats, lons):
        """Get elevation and slope data from GeoTIFF for profile calculations.
        Always uses the original full dataset to ensure consistent results regardless of zoom level.
        Returns (elevations, slopes, extent) or (None, None, None) if data unavailable."""
        if self.geotiff_dataset_original is None:
            return None, None, None

        if hasattr(self, 'geotiff_original_extent') and self.geotiff_original_extent is not None:
            extent = self.geotiff_original_extent
        elif self.geotiff_extent is not None:
            extent = self.geotiff_extent
        else:
            return None, None, None

        try:
            from rasterio.warp import reproject, Resampling, calculate_default_transform
            from rasterio.transform import array_bounds

            if self.geotiff_dataset_original.crs != "EPSG:4326":
                transformer = pyproj.Transformer.from_crs("EPSG:4326", self.geotiff_dataset_original.crs, always_xy=True)
                left, bottom = transformer.transform(extent[0], extent[2])
                right, top = transformer.transform(extent[1], extent[3])
            else:
                left, right, bottom, top = extent

            data = self.geotiff_dataset_original.read(1).astype(float)

            if self.geotiff_dataset_original.crs != "EPSG:4326":
                src_crs = self.geotiff_dataset_original.crs
                dst_crs = "EPSG:4326"
                src_bounds = self.geotiff_dataset_original.bounds
                # bounds may be BoundingBox (named tuple) or plain tuple (left, bottom, right, top)
                try:
                    src_left, src_bottom, src_right, src_top = src_bounds.left, src_bounds.bottom, src_bounds.right, src_bounds.top
                except AttributeError:
                    src_left, src_bottom, src_right, src_top = src_bounds[0], src_bounds[1], src_bounds[2], src_bounds[3]
                src_width = self.geotiff_dataset_original.width
                src_height = self.geotiff_dataset_original.height
                dst_transform, dst_width, dst_height = calculate_default_transform(
                    src_crs, dst_crs, src_width, src_height,
                    left=src_left, bottom=src_bottom,
                    right=src_right, top=src_top
                )
                reprojected_data = np.zeros((dst_height, dst_width), dtype=np.float32)
                reproject(
                    source=data,
                    destination=reprojected_data,
                    src_transform=self.geotiff_dataset_original.transform,
                    src_crs=src_crs,
                    dst_transform=dst_transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.bilinear
                )
                data = reprojected_data
                # array_bounds returns (left, bottom, right, top)
                dst_left, dst_bottom, dst_right, dst_top = array_bounds(dst_height, dst_width, dst_transform)
                extent = [dst_left, dst_right, dst_bottom, dst_top]

            data[data < -11000] = np.nan
            data[data >= 0] = np.nan

            left, right, bottom, top = extent
            nrows, ncols = data.shape
            rows = ((top - lats) / (top - bottom) * (nrows - 1)).clip(0, nrows - 1)
            cols = ((lons - left) / (right - left) * (ncols - 1)).clip(0, ncols - 1)

            elevations = []
            slopes = []
            for r, c in zip(rows, cols):
                ir, ic = int(round(r)), int(round(c))
                elevations.append(data[ir, ic])
                slope = None
                if 0 < ir < nrows-1 and 0 < ic < ncols-1:
                    center_lat = (extent[2] + extent[3]) / 2
                    m_per_deg_lat = 111320.0
                    m_per_deg_lon = 111320.0 * np.cos(np.radians(center_lat))
                    res_lat_deg = (extent[3] - extent[2]) / nrows
                    res_lon_deg = (extent[1] - extent[0]) / ncols
                    dx_m = res_lon_deg * m_per_deg_lon
                    dy_m = res_lat_deg * m_per_deg_lat
                    window = data[ir-1:ir+2, ic-1:ic+2]
                    if window.shape == (3, 3) and not np.all(np.isnan(window)):
                        dz_dy, dz_dx = np.gradient(window, dy_m, dx_m)
                        slope_rad = np.arctan(np.sqrt(dz_dx[1, 1]**2 + dz_dy[1, 1]**2))
                        slope = np.degrees(slope_rad)
                slopes.append(slope if slope is not None else np.nan)

            return np.array(elevations), np.array(slopes), extent

        except Exception as e:
            print(f"Error reading profile data from GeoTIFF: {e}")
            return None, None, None

    def _profile_arrays_along_segment_endpoints(self, lat1, lon1, lat2, lon2, n=100):
        """Sample elevation/slope along a straight line in lat/lon; cumulative distance via geodesic (m).

        Same sampling as crossline/pitch profile plots. Returns (dists, elevations, slopes) or
        (None, None, None) if no elevation source is available.
        """
        lats = np.linspace(lat1, lat2, n)
        lons = np.linspace(lon1, lon2, n)
        elevations, slopes, _ = self._get_profile_data_from_geotiff(lats, lons)
        if elevations is None:
            if self.geotiff_data_array is None or self.geotiff_extent is None:
                return None, None, None
            left, right, bottom, top = tuple(self.geotiff_extent)
            nrows, ncols = self.geotiff_data_array.shape
            rows = ((top - lats) / (top - bottom) * (nrows - 1)).clip(0, nrows - 1)
            cols = ((lons - left) / (right - left) * (ncols - 1)).clip(0, ncols - 1)
            elev_list = []
            slope_list = []
            for r, c in zip(rows, cols):
                ir, ic = int(round(r)), int(round(c))
                elev_list.append(self.geotiff_data_array[ir, ic])
                slope = None
                if 0 < ir < nrows - 1 and 0 < ic < ncols - 1:
                    center_lat_geotiff = (self.geotiff_extent[2] + self.geotiff_extent[3]) / 2
                    m_per_deg_lat = 111320.0
                    m_per_deg_lon = 111320.0 * np.cos(np.radians(center_lat_geotiff))
                    res_lat_deg = (self.geotiff_extent[3] - self.geotiff_extent[2]) / nrows
                    res_lon_deg = (self.geotiff_extent[1] - self.geotiff_extent[0]) / ncols
                    dx_m = res_lon_deg * m_per_deg_lon
                    dy_m = res_lat_deg * m_per_deg_lat
                    window = self.geotiff_data_array[ir - 1 : ir + 2, ic - 1 : ic + 2]
                    if window.shape == (3, 3) and not np.all(np.isnan(window)):
                        dz_dy, dz_dx = np.gradient(window, dy_m, dx_m)
                        slope_rad = np.arctan(np.sqrt(dz_dx[1, 1] ** 2 + dz_dy[1, 1] ** 2))
                        slope = np.degrees(slope_rad)
                slope_list.append(slope if slope is not None else np.nan)
            elevations = np.array(elev_list)
            slopes = np.array(slope_list)
        geod = pyproj.Geod(ellps="WGS84")
        dists = [0.0]
        for i in range(1, n):
            try:
                _, _, seg = geod.inv(lons[i - 1], lats[i - 1], lons[i], lats[i])
                if seg is not None and not (np.isnan(seg) or np.isinf(seg)):
                    dists.append(dists[-1] + float(seg))
                else:
                    dists.append(dists[-1])
            except Exception:
                dists.append(dists[-1])
        return np.array(dists, dtype=float), elevations, slopes

    def _draw_crossline_profile(self):
        if not (
            self.geotiff_data_array is not None and
            self.geotiff_extent is not None and
            self.cross_line_data and
            len(self.cross_line_data) == 2
        ):
            self.profile_ax.clear()
            self.profile_canvas.draw_idle()
            self._profile_dists = None
            self._profile_elevations = None
            self._profile_slopes = None
            if hasattr(self, 'profile_info_text') and self.profile_info_text:
                try:
                    self.profile_info_text.set_visible(False)
                except Exception:
                    pass
                self.profile_info_text = None
            return
        self.profile_ax.clear()
        for ax in self.profile_fig.get_axes():
            if ax != self.profile_ax:
                ax.remove()
        self.profile_ax.set_title("Crossline Elevation Profile", fontsize=8)
        self.profile_ax.set_xlabel("Distance (m)", fontsize=8)
        self.profile_ax.set_ylabel("Elevation (m)", fontsize=8)
        self.profile_ax.tick_params(axis='both', which='major', labelsize=7)
        slope_ax = None
        self._profile_dists = None
        self._profile_elevations = None
        self._profile_slopes = None
        (lat1, lon1), (lat2, lon2) = self.cross_line_data
        lats = np.linspace(lat1, lat2, 100)
        lons = np.linspace(lon1, lon2, 100)

        elevations, slopes, profile_extent = self._get_profile_data_from_geotiff(lats, lons)
        if elevations is None:
            if self.geotiff_data_array is None or self.geotiff_extent is None:
                return
            left, right, bottom, top = tuple(self.geotiff_extent)
            nrows, ncols = self.geotiff_data_array.shape
            rows = ((top - lats) / (top - bottom) * (nrows - 1)).clip(0, nrows - 1)
            cols = ((lons - left) / (right - left) * (ncols - 1)).clip(0, ncols - 1)
            elevations = []
            slopes = []
            for r, c in zip(rows, cols):
                ir, ic = int(round(r)), int(round(c))
                elevations.append(self.geotiff_data_array[ir, ic])
                slope = None
                if 0 < ir < nrows-1 and 0 < ic < ncols-1:
                    center_lat_geotiff = (self.geotiff_extent[2] + self.geotiff_extent[3]) / 2
                    m_per_deg_lat = 111320.0
                    m_per_deg_lon = 111320.0 * np.cos(np.radians(center_lat_geotiff))
                    res_lat_deg = (self.geotiff_extent[3] - self.geotiff_extent[2]) / nrows
                    res_lon_deg = (self.geotiff_extent[1] - self.geotiff_extent[0]) / ncols
                    dx_m = res_lon_deg * m_per_deg_lon
                    dy_m = res_lat_deg * m_per_deg_lat
                    window = self.geotiff_data_array[ir-1:ir+2, ic-1:ic+2]
                    if window.shape == (3, 3) and not np.all(np.isnan(window)):
                        dz_dy, dz_dx = np.gradient(window, dy_m, dx_m)
                        slope_rad = np.arctan(np.sqrt(dz_dx[1, 1]**2 + dz_dy[1, 1]**2))
                        slope = np.degrees(slope_rad)
                slopes.append(slope if slope is not None else np.nan)
            elevations = np.array(elevations)
            slopes = np.array(slopes)
        if hasattr(self, 'local_proj_transformer') and self.local_proj_transformer is not None:
            eastings, northings = self.local_proj_transformer.transform(lons, lats)
            dists = np.sqrt((eastings - eastings[0])**2 + (northings - northings[0])**2)
        else:
            dists = np.linspace(0, 1, 100)
        self._profile_dists = dists
        self._profile_elevations = elevations
        self._profile_slopes = slopes
        self.profile_ax.plot(dists, elevations, color='darkorchid', lw=1, label='Elevation')
        if self._show_slope_on_profile():
            slope_ax = self.profile_ax.twinx()
            slope_ax.plot(dists, slopes, color='teal', lw=1, linestyle='--', label='Slope (deg)')
            slope_ax.set_ylabel('Slope (deg)', fontsize=8)
            slope_ax.tick_params(axis='y', labelsize=7)
            slope_ax.grid(False)
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

        if (
            hasattr(self, "set_cal_info_text")
            and self.geotiff_data_array is not None
            and self.cross_line_data
            and len(self.cross_line_data) == 2
        ):
            try:
                valid_elevs = elevations[~np.isnan(elevations)]
                if valid_elevs.size > 0:
                    mean_depth = np.mean(valid_elevs)
                    median_depth = np.median(valid_elevs)
                    min_depth = np.min(valid_elevs)
                    max_depth = np.max(valid_elevs)
                    msg = (
                        f"Crossline Depth Stats (meters):\n"
                        f"Mean: {mean_depth:.2f}\n"
                        f"Median: {median_depth:.2f}\n"
                        f"Minimum Elevation: {min_depth:.2f}\n"
                        f"Maximum Elevation: {max_depth:.2f}"
                    )
                else:
                    msg = "No valid elevation data under crossline."
            except Exception as e:
                msg = f"Error calculating crossline stats: {e}"
            self.set_cal_info_text(msg)

    def _draw_pitch_line_profile(self):
        for ax in self.profile_fig.get_axes():
            if ax != self.profile_ax:
                ax.remove()
        self.profile_ax.clear()
        self.profile_ax.set_title("Pitch Line Elevation Profile", fontsize=8)
        self.profile_ax.set_xlabel("Distance (m)", fontsize=8)
        self.profile_ax.set_ylabel("Elevation (m)", fontsize=8)
        self.profile_ax.tick_params(axis='both', which='major', labelsize=7)
        slope_ax = None

        if not (hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2):
            self.profile_fig.tight_layout(pad=1.0)
            if hasattr(self, 'profile_canvas'):
                self.profile_canvas.draw_idle()
            return

        (lat1, lon1), (lat2, lon2) = self.pitch_line_points
        lats = np.linspace(lat1, lat2, 100)
        lons = np.linspace(lon1, lon2, 100)

        elevations, slopes, profile_extent = self._get_profile_data_from_geotiff(lats, lons)
        if elevations is None:
            if self.geotiff_data_array is None or self.geotiff_extent is None:
                return
            left, right, bottom, top = tuple(self.geotiff_extent)
            nrows, ncols = self.geotiff_data_array.shape
            rows = ((top - lats) / (top - bottom) * (nrows - 1)).clip(0, nrows - 1)
            cols = ((lons - left) / (right - left) * (ncols - 1)).clip(0, ncols - 1)
            elevations = []
            slopes = []
            for r, c in zip(rows, cols):
                ir, ic = int(round(r)), int(round(c))
                elevations.append(self.geotiff_data_array[ir, ic])
                slope = None
                if 0 < ir < nrows-1 and 0 < ic < ncols-1:
                    center_lat_geotiff = (self.geotiff_extent[2] + self.geotiff_extent[3]) / 2
                    m_per_deg_lat = 111320.0
                    m_per_deg_lon = 111320.0 * np.cos(np.radians(center_lat_geotiff))
                    res_lat_deg = (self.geotiff_extent[3] - self.geotiff_extent[2]) / nrows
                    res_lon_deg = (self.geotiff_extent[1] - self.geotiff_extent[0]) / ncols
                    dx_m = res_lon_deg * m_per_deg_lon
                    dy_m = res_lat_deg * m_per_deg_lat
                    window = self.geotiff_data_array[ir-1:ir+2, ic-1:ic+2]
                    if window.shape == (3, 3) and not np.all(np.isnan(window)):
                        dz_dy, dz_dx = np.gradient(window, dy_m, dx_m)
                        slope_rad = np.arctan(np.sqrt(dz_dx[1, 1]**2 + dz_dy[1, 1]**2))
                        slope = np.degrees(slope_rad)
                slopes.append(slope if slope is not None else np.nan)
            elevations = np.array(elevations)
            slopes = np.array(slopes)
        try:
            geod = pyproj.Geod(ellps="WGS84")
            dists = [0.0]
            for i in range(1, len(lats)):
                _, _, dist = geod.inv(lons[i-1], lats[i-1], lons[i], lats[i])
                dists.append(dists[-1] + dist)
            dists = np.array(dists)
        except Exception:
            dists = np.linspace(0, 1, 100)
        self.profile_ax.plot(dists, elevations, color='red', lw=1, label='Elevation')
        if self._show_slope_on_profile():
            slope_ax = self.profile_ax.twinx()
            slope_ax.plot(dists, slopes, color='teal', lw=1, linestyle='--', label='Slope (deg)')
            slope_ax.set_ylabel('Slope (deg)', fontsize=8)
            slope_ax.tick_params(axis='y', labelsize=7)
            slope_ax.grid(False)
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

        if (
            hasattr(self, "set_cal_info_text")
            and self.geotiff_data_array is not None
            and hasattr(self, "pitch_line_points")
            and len(self.pitch_line_points) == 2
        ):
            try:
                valid_elevs = elevations[~np.isnan(elevations)]
                if valid_elevs.size > 0:
                    mean_depth = np.mean(valid_elevs)
                    median_depth = np.median(valid_elevs)
                    min_depth = np.min(valid_elevs)
                    max_depth = np.max(valid_elevs)
                    msg = (
                        f"Pitch Line Depth Stats (meters):\n"
                        f"Mean: {mean_depth:.2f}\n"
                        f"Median: {median_depth:.2f}\n"
                        f"Minimum Elevation: {min_depth:.2f}\n"
                        f"Maximum Elevation: {max_depth:.2f}"
                    )
                else:
                    msg = "No valid elevation data under pitch line."
            except Exception as e:
                msg = f"Error calculating pitch line stats: {e}"
            self.set_cal_info_text(msg)
            self.set_cal_info_text("Using Median Depth for Line Offset", append=False)

    def _draw_roll_line_profile(self):
        for ax in self.profile_fig.get_axes():
            if ax != self.profile_ax:
                ax.remove()
        self.profile_ax.clear()
        self.profile_ax.set_title("Roll Line Elevation Profile", fontsize=8)
        self.profile_ax.set_xlabel("Distance (m)", fontsize=8)
        self.profile_ax.set_ylabel("Elevation (m)", fontsize=8)
        self.profile_ax.tick_params(axis="both", which="major", labelsize=7)
        slope_ax = None

        if not (hasattr(self, "roll_line_points") and len(self.roll_line_points) == 2):
            self.profile_fig.tight_layout(pad=1.0)
            if hasattr(self, "profile_canvas"):
                self.profile_canvas.draw_idle()
            return

        (lat1, lon1), (lat2, lon2) = self.roll_line_points
        lats = np.linspace(lat1, lat2, 100)
        lons = np.linspace(lon1, lon2, 100)

        elevations, slopes, profile_extent = self._get_profile_data_from_geotiff(lats, lons)
        if elevations is None:
            if self.geotiff_data_array is None or self.geotiff_extent is None:
                return
            left, right, bottom, top = tuple(self.geotiff_extent)
            nrows, ncols = self.geotiff_data_array.shape
            rows = ((top - lats) / (top - bottom) * (nrows - 1)).clip(0, nrows - 1)
            cols = ((lons - left) / (right - left) * (ncols - 1)).clip(0, ncols - 1)
            elevations = []
            slopes = []
            for r, c in zip(rows, cols):
                ir, ic = int(round(r)), int(round(c))
                elevations.append(self.geotiff_data_array[ir, ic])
                slope = None
                if 0 < ir < nrows - 1 and 0 < ic < ncols - 1:
                    center_lat_geotiff = (self.geotiff_extent[2] + self.geotiff_extent[3]) / 2
                    m_per_deg_lat = 111320.0
                    m_per_deg_lon = 111320.0 * np.cos(np.radians(center_lat_geotiff))
                    res_lat_deg = (self.geotiff_extent[3] - self.geotiff_extent[2]) / nrows
                    res_lon_deg = (self.geotiff_extent[1] - self.geotiff_extent[0]) / ncols
                    dx_m = res_lon_deg * m_per_deg_lon
                    dy_m = res_lat_deg * m_per_deg_lat
                    window = self.geotiff_data_array[ir - 1 : ir + 2, ic - 1 : ic + 2]
                    if window.shape == (3, 3) and not np.all(np.isnan(window)):
                        dz_dy, dz_dx = np.gradient(window, dy_m, dx_m)
                        slope_rad = np.arctan(np.sqrt(dz_dx[1, 1] ** 2 + dz_dy[1, 1] ** 2))
                        slope = np.degrees(slope_rad)
                slopes.append(slope if slope is not None else np.nan)
            elevations = np.array(elevations)
            slopes = np.array(slopes)
        try:
            geod = pyproj.Geod(ellps="WGS84")
            dists = [0.0]
            for i in range(1, len(lats)):
                _, _, dist = geod.inv(lons[i - 1], lats[i - 1], lons[i], lats[i])
                dists.append(dists[-1] + dist)
            dists = np.array(dists)
        except Exception:
            dists = np.linspace(0, 1, 100)
        self.profile_ax.plot(dists, elevations, color="purple", lw=1, label="Elevation")
        if self._show_slope_on_profile():
            slope_ax = self.profile_ax.twinx()
            slope_ax.plot(dists, slopes, color="teal", lw=1, linestyle="--", label="Slope (deg)")
            slope_ax.set_ylabel("Slope (deg)", fontsize=8)
            slope_ax.tick_params(axis="y", labelsize=7)
            slope_ax.grid(False)
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
            self.profile_ax.legend(handles, labels, fontsize=10, loc="upper right")
        self.profile_fig.tight_layout(pad=1.0)
        self.profile_canvas.draw_idle()

        if (
            hasattr(self, "set_cal_info_text")
            and self.geotiff_data_array is not None
            and hasattr(self, "roll_line_points")
            and len(self.roll_line_points) == 2
        ):
            try:
                valid_elevs = elevations[~np.isnan(elevations)]
                if valid_elevs.size > 0:
                    mean_depth = np.mean(valid_elevs)
                    median_depth = np.median(valid_elevs)
                    min_depth = np.min(valid_elevs)
                    max_depth = np.max(valid_elevs)
                    msg = (
                        f"Roll Line Depth Stats (meters):\n"
                        f"Mean: {mean_depth:.2f}\n"
                        f"Median: {median_depth:.2f}\n"
                        f"Minimum Elevation: {min_depth:.2f}\n"
                        f"Maximum Elevation: {max_depth:.2f}"
                    )
                else:
                    msg = "No valid elevation data under roll line."
            except Exception as e:
                msg = f"Error calculating roll line stats: {e}"
            self.set_cal_info_text(msg)

    def _draw_performance_profile(self):
        """Draw elevation/slope profile for performance line 1 and BIST segment 1 (map colors: sienna, gold)."""
        for ax in self.profile_fig.get_axes():
            if ax != self.profile_ax:
                ax.remove()
        self.profile_ax.clear()
        self.profile_ax.set_xlabel("Distance (m)", fontsize=8)
        self.profile_ax.set_ylabel("Elevation (m)", fontsize=8)
        self.profile_ax.tick_params(axis='both', which='major', labelsize=7)
        slope_ax = None

        # Same colors as plotting_mixin performance map: line 1 swath = sienna, BIST = gold
        perf_swath_color = "sienna"
        perf_bist_color = "gold"

        segments = []
        if hasattr(self, "performance_profile_line") and self.performance_profile_line and len(self.performance_profile_line) == 2:
            segments.append(
                (
                    self.performance_profile_line[0],
                    self.performance_profile_line[1],
                    perf_swath_color,
                    "Swath (line 1)",
                )
            )
        bist_data = getattr(self, "performance_bist_segments_data", None) or []
        if len(bist_data) > 0 and len(bist_data[0]) == 2:
            segments.append((bist_data[0][0], bist_data[0][1], perf_bist_color, "BIST"))

        if not segments:
            self.profile_ax.set_title("Performance — Elevation Profile", fontsize=8)
            self.profile_fig.tight_layout(pad=1.0)
            if hasattr(self, "profile_canvas"):
                self.profile_canvas.draw_idle()
            return

        title = "Performance Line 1 + BIST — Elevation Profile" if len(segments) > 1 else "Performance Line 1 — Elevation Profile"
        self.profile_ax.set_title(title, fontsize=8)

        geod = pyproj.Geod(ellps="WGS84")
        cumulative = 0.0
        elev_arrays = []
        slope_dists_concat = []
        slopes_concat = []

        for seg_start, seg_end, color, label in segments:
            lat1, lon1 = seg_start
            lat2, lon2 = seg_end
            lats = np.linspace(lat1, lat2, 100)
            lons = np.linspace(lon1, lon2, 100)
            elevations, slopes, _ = self._get_profile_data_from_geotiff(lats, lons)
            if elevations is None:
                continue
            try:
                dists = [0.0]
                for i in range(1, len(lats)):
                    _, _, dist = geod.inv(lons[i - 1], lats[i - 1], lons[i], lats[i])
                    dists.append(dists[-1] + dist)
                dists = np.array(dists, dtype=float)
            except Exception:
                dists = np.linspace(0.0, 1.0, len(lats))
            dists_plot = dists + cumulative
            self.profile_ax.plot(dists_plot, elevations, color=color, lw=1.5, label=label)
            elev_arrays.append(elevations)
            if self._show_slope_on_profile() and slopes is not None:
                # Skip first sample when appending after first segment to avoid duplicate junction point
                start_i = 1 if slope_dists_concat else 0
                slope_dists_concat.extend(dists_plot[start_i:].tolist())
                slopes_concat.extend(slopes[start_i:].tolist())
            cumulative = float(dists_plot[-1])

        if not elev_arrays:
            self.profile_fig.tight_layout(pad=1.0)
            if hasattr(self, "profile_canvas"):
                self.profile_canvas.draw_idle()
            return

        if self._show_slope_on_profile() and slope_dists_concat:
            slope_ax = self.profile_ax.twinx()
            sd = np.array(slope_dists_concat, dtype=float)
            ss = np.array(slopes_concat, dtype=float)
            slope_ax.plot(sd, ss, color="teal", lw=1, linestyle="--", label="Slope (deg)")
            slope_ax.set_ylabel("Slope (deg)", fontsize=8)
            slope_ax.tick_params(axis="y", labelsize=7)
            slope_ax.grid(False)
            if np.any(~np.isnan(ss)):
                slope_ax.set_ylim(0, np.nanmax(ss[~np.isnan(ss)]) * 1.1)

        self.profile_ax.set_xlim(0.0, max(cumulative, 1e-6))
        stacked_e = np.concatenate(elev_arrays)
        if np.any(~np.isnan(stacked_e)):
            self.profile_ax.set_ylim(np.nanmin(stacked_e), np.nanmax(stacked_e))

        handles, labels = self.profile_ax.get_legend_handles_labels()
        if slope_ax:
            slope_handles, slope_labels = slope_ax.get_legend_handles_labels()
            handles.extend(slope_handles)
            labels.extend(slope_labels)
        if handles:
            self.profile_ax.legend(handles, labels, fontsize=10, loc="upper right")
        self.profile_fig.tight_layout(pad=1.0)
        self.profile_canvas.draw_idle()

    def _on_slope_profile_checkbox_changed(self, *_args):
        if hasattr(self, "slope_profile_checkbox"):
            self.show_slope_profile_var = self.slope_profile_checkbox.isChecked()
        self._draw_current_profile()

    def _draw_current_profile(self):
        if not hasattr(self, 'param_notebook'):
            return
        if not hasattr(self, 'profile_ax') or not hasattr(self, 'profile_canvas'):
            return
        current_tab = self.param_notebook.currentIndex()
        if current_tab == 0:
            if self._calibration_profile_select_is_roll():
                self._draw_roll_line_profile()
            else:
                self._draw_pitch_line_profile()
        elif current_tab == 1:
            self._draw_accuracy_selected_profile()
        elif current_tab == 2:
            self._draw_line_planning_profile()
        elif current_tab == 3:
            self._draw_performance_profile()
        else:
            self.profile_ax.clear()
            self.profile_ax.set_title("Elevation Profile", fontsize=8)
            self.profile_ax.set_xlabel("Distance (m)", fontsize=8)
            self.profile_ax.set_ylabel("Elevation (m)", fontsize=8)
            self.profile_fig.tight_layout(pad=1.0)
            self.profile_canvas.draw()

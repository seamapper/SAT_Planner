"""
GeoTIFF loading, display mode, dynamic resolution, and contour controls.
"""
import os
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


class GeoTIFFMixin:
    """Mixin providing GeoTIFF load/remove, display mode, dynamic resolution, and contours."""

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

        if self.geotiff_extent is not None and self.geotiff_dataset_original is not None:
            if getattr(self, "dynamic_resolution_enabled", True):
                full_width = self.geotiff_extent[1] - self.geotiff_extent[0]
                full_height = self.geotiff_extent[3] - self.geotiff_extent[2]
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
            self.backscatter_box_stage = 0
            self.backscatter_box_centerline = None
            self.canvas_widget.setCursor(Qt.CursorShape.CrossCursor)
        else:
            had_partial = int(getattr(self, "backscatter_box_stage", 0)) > 0
            self.backscatter_box_stage = 0
            self.backscatter_box_centerline = None
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
            if had_partial:
                try:
                    self._plot_survey_plan(preserve_view_limits=True)
                except Exception:
                    self.canvas.draw_idle()

    def _on_backscatter_show_stats_clicked(self):
        """Open backscatter statistics dialog for the currently saved box."""
        if not getattr(self, "backscatter_box_stats", None):
            self._show_message("info", "Backscatter Statistics", "Draw a box first to compute statistics.")
            return
        self._show_backscatter_stats_dialog()

    def _on_backscatter_clear_box_clicked(self):
        """Remove backscatter stats box overlay and clear cached statistics."""
        self.backscatter_box_draw_mode = False
        self.backscatter_box_stage = 0
        self.backscatter_box_centerline = None
        self.backscatter_box_vertices = None
        self.backscatter_box_half_width_m = None
        self.backscatter_box_stats = None
        if hasattr(self, "backscatter_draw_box_btn"):
            self.backscatter_draw_box_btn.blockSignals(True)
            self.backscatter_draw_box_btn.setChecked(False)
            self.backscatter_draw_box_btn.blockSignals(False)
        if hasattr(self, "backscatter_box_patch") and self.backscatter_box_patch is not None:
            try:
                self.backscatter_box_patch.remove()
            except Exception:
                pass
            self.backscatter_box_patch = None
        if hasattr(self, "backscatter_stats_dialog") and self.backscatter_stats_dialog is not None:
            try:
                self.backscatter_stats_dialog.close()
            except Exception:
                pass
            self.backscatter_stats_dialog = None
            self.backscatter_stats_canvas = None
            self.backscatter_stats_ax = None
            self.backscatter_stats_text_label = None
        self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
        try:
            self._plot_survey_plan(preserve_view_limits=True)
        except Exception:
            self.canvas.draw_idle()

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
        self._show_backscatter_stats_dialog()

    def _compute_backscatter_box_stats(self):
        """Compute histogram input and dimensions for current backscatter box."""
        if self.backscatter_box_vertices is None:
            return None
        bs = getattr(self, "backscatter_raster_array", None)
        extent = getattr(self, "geotiff_extent", None)
        if bs is None or extent is None:
            return None
        vertices = np.array(self.backscatter_box_vertices, dtype=float)
        if vertices.shape != (4, 2):
            return None
        lon_min = float(np.min(vertices[:, 0]))
        lon_max = float(np.max(vertices[:, 0]))
        lat_min = float(np.min(vertices[:, 1]))
        lat_max = float(np.max(vertices[:, 1]))
        left, right, bottom, top = extent
        nrows, ncols = bs.shape
        if right == left or top == bottom:
            return None
        c0 = int(np.clip((lon_min - left) / (right - left) * (ncols - 1), 0, ncols - 1))
        c1 = int(np.clip((lon_max - left) / (right - left) * (ncols - 1), 0, ncols - 1))
        r0 = int(np.clip((top - lat_max) / (top - bottom) * (nrows - 1), 0, nrows - 1))
        r1 = int(np.clip((top - lat_min) / (top - bottom) * (nrows - 1), 0, nrows - 1))
        rr0, rr1 = sorted([r0, r1])
        cc0, cc1 = sorted([c0, c1])
        subset = bs[rr0:rr1 + 1, cc0:cc1 + 1]
        lon_vals = np.linspace(left, right, ncols)
        lat_vals = np.linspace(top, bottom, nrows)
        sub_lons = lon_vals[cc0:cc1 + 1]
        sub_lats = lat_vals[rr0:rr1 + 1]
        lon_mesh, lat_mesh = np.meshgrid(sub_lons, sub_lats)
        points = np.column_stack((lon_mesh.ravel(), lat_mesh.ravel()))
        from matplotlib.path import Path
        poly_path = Path(vertices)
        inside = poly_path.contains_points(points).reshape(subset.shape)
        valid = subset[np.isfinite(subset) & inside]
        geod = pyproj.Geod(ellps="WGS84")
        centerline = getattr(self, "backscatter_box_centerline", None)
        if centerline and len(centerline) == 2:
            (a_lat, a_lon), (b_lat, b_lon) = centerline
            _, _, height_m = geod.inv(a_lon, a_lat, b_lon, b_lat)
        else:
            height_m = 0.0
        width_m = 2.0 * float(getattr(self, "backscatter_box_half_width_m", 0.0) or 0.0)
        return {
            "values": valid,
            "width_m": float(abs(width_m)),
            "height_m": float(abs(height_m)),
        }

    def _show_backscatter_stats_dialog(self):
        """Create/show Backscatter Statistics dialog with histogram and box dimensions."""
        stats = getattr(self, "backscatter_box_stats", None)
        if not stats:
            return
        if self.backscatter_stats_dialog is None:
            self.backscatter_stats_dialog = QDialog(self)
            self.backscatter_stats_dialog.setWindowTitle("Backscatter Statistics")
            self.backscatter_stats_dialog.resize(700, 620)
            layout = QVBoxLayout(self.backscatter_stats_dialog)
            top_row = QWidget()
            top_layout = QVBoxLayout(top_row)
            top_layout.setContentsMargins(0, 0, 0, 0)
            fig = Figure(figsize=(6.4, 4.2))
            self.backscatter_stats_ax = fig.add_subplot(111)
            fig.subplots_adjust(bottom=0.18)
            self.backscatter_stats_canvas = FigureCanvas(fig)
            top_layout.addWidget(self.backscatter_stats_canvas)
            layout.addWidget(top_row)
            self.backscatter_stats_text_label = QLabel("")
            self.backscatter_stats_text_label.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft)
            self.backscatter_stats_text_label.setWordWrap(True)
            layout.addWidget(self.backscatter_stats_text_label)

        vals = stats.get("values", np.array([]))
        self.backscatter_stats_ax.clear()
        if vals.size > 0:
            self.backscatter_stats_ax.hist(vals, bins=40, color="steelblue", edgecolor="black", alpha=0.85)
            self.backscatter_stats_ax.set_title("Backscatter Value Distribution")
            self.backscatter_stats_ax.set_xlabel("Backscatter Value")
            self.backscatter_stats_ax.set_ylabel("Count")
        else:
            self.backscatter_stats_ax.text(
                0.5,
                0.5,
                "No valid backscatter values in selected box.",
                transform=self.backscatter_stats_ax.transAxes,
                ha="center",
                va="center",
            )
            self.backscatter_stats_ax.set_title("Backscatter Value Distribution")
        self.backscatter_stats_ax.figure.subplots_adjust(bottom=0.18)
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

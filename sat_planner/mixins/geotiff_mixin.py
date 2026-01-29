"""
GeoTIFF loading, display mode, dynamic resolution, and contour controls.
"""
import os
import numpy as np
from PyQt6.QtWidgets import (
    QApplication,
    QDialog,
    QFileDialog,
    QLabel,
    QProgressBar,
    QPushButton,
    QVBoxLayout,
)
from PyQt6.QtCore import Qt, QTimer

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
            # Calculate the visible region based on current plot limits
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

            # Filter Z-values (match SAT_Update: < -11000 or >= 0 → NaN)
            data[data < -11000] = np.nan
            data[data >= 0] = np.nan

            # Validate the extent to prevent flipping
            if region_extent[1] <= region_extent[0] or region_extent[3] <= region_extent[2]:
                print(f"Warning: Invalid extent calculated: {region_extent}. Skipping resolution update.")
                return False

            # Update the current data and extent
            self.geotiff_data_array = data.astype(float)
            self.geotiff_extent = region_extent
            self.geotiff_current_resolution = downsample_factor

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

            # Filter Z-values: values < -11000 or >= 0 are set to NaN (match SAT_Update)
            progress_label.setText("Processing elevation data...")
            QApplication.processEvents()

            self.geotiff_data_array[self.geotiff_data_array < -11000] = np.nan
            self.geotiff_data_array[self.geotiff_data_array >= 0] = np.nan

            # Store the full resolution data for dynamic loading
            self.geotiff_full_resolution = self.geotiff_data_array.copy()
            self.geotiff_current_resolution = downsample_factor
            self.geotiff_zoom_level = 1.0  # Start at full resolution

            # Close progress window
            progress_window.close()

            # Enable Remove GeoTIFF button and Pick Center button
            if hasattr(self, 'remove_geotiff_btn'):
                self.remove_geotiff_btn.setEnabled(True)
            if hasattr(self, 'pick_center_btn'):
                self.pick_center_btn.setEnabled(True)
                # Make "Pick Center from GeoTIFF" button bold and orange after loading GeoTIFF
                self.pick_center_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")
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

            # Update the plot and zoom after loading
            self._plot_survey_plan()
            # Store the initial limits after first plot
            self.initial_geotiff_xlim = self.current_xlim
            self.initial_geotiff_ylim = self.current_ylim
            self._zoom_to_geotiff()

        except RasterioIOError as e:
            progress_window.destroy()
            self._show_message("error", "GeoTIFF Error", f"Could not open GeoTIFF file: {e}")
            self._clear_plot(full_clear=True)
        except Exception as e:
            progress_window.destroy()
            self._show_message("error", "GeoTIFF Error", f"An unexpected error occurred while loading GeoTIFF: {e}")
            self._clear_plot(full_clear=True)

        # Draw crossline profile
        self._draw_crossline_profile()

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
        """Handle contour interval entry change. Keeps all entry fields synchronized."""
        # Sync all entry fields (calibration, reference, and line planning)
        try:
            if hasattr(self, '_syncing_contour_interval'):
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

        # Show message in activity log
        if hasattr(self, 'set_cal_info_text'):
            self.set_cal_info_text("GeoTIFF removed.")
        elif hasattr(self, 'set_ref_info_text'):
            self.set_ref_info_text("GeoTIFF removed.")

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

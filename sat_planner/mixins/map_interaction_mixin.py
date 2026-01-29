"""
Map interaction: plot click, scroll zoom, middle-button pan, zoom-to helpers, pick center.
_on_plot_click, _on_scroll, _on_middle_*, _zoom_to_*, _toggle_pick_center_mode.
"""
import traceback

import numpy as np
from PyQt6.QtCore import Qt, QTimer

from sat_planner.constants import GEOSPATIAL_LIBS_AVAILABLE, rowcol


class MapInteractionMixin:
    """Mixin for map interaction: click, scroll, pan, zoom-to, pick center."""

    def _toggle_pick_center_mode(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot pick center.")
            return

        if self.geotiff_dataset_original is None:
            self._show_message("warning","No GeoTIFF", "Load a GeoTIFF first to pick center.")
            return

        self.pick_center_mode = not self.pick_center_mode
        if self.pick_center_mode:
            self.pick_center_btn.setText("Picking Enabled (Click Plot)")
            # Change cursor to cross
            self.canvas_widget.setCursor(Qt.CursorShape.CrossCursor)
            if hasattr(self, 'set_ref_info_text') and self.param_notebook.currentIndex() == 1:
                self.set_ref_info_text("Click on the Survey Plan Plot to set Central Lat/Lon, Line Length (Depth * Line Length Multiplier) and Line Spacing (Depth * Separation Multiplier).")
            # Remove info text if present
            if hasattr(self, 'pick_center_info_text') and self.pick_center_info_text is not None:
                self.pick_center_info_text.set_visible(False)
                self.pick_center_info_text = None
                self.canvas.draw_idle()
        else:
            self.pick_center_btn.setText("Pick Center from GeoTIFF")
            # Restore cursor to default
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
            # Remove info text if present
            if hasattr(self, 'pick_center_info_text') and self.pick_center_info_text is not None:
                self.pick_center_info_text.set_visible(False)
                self.pick_center_info_text = None
                self.canvas.draw_idle()

    def _on_plot_click(self, event):
        # Ignore middle button (button 2) - let it be handled by pan handlers
        if event.button == 2:
            return
        if getattr(self, '_handle_line_planning_plot_click', lambda e: False)(event):
            return
        if self.pick_pitch_line_mode:
            if event.inaxes != self.ax:
                print("DEBUG: Clicked outside plot axes. Click ignored.")
                return
            clicked_lon = event.xdata
            clicked_lat = event.ydata
            if clicked_lon is None or clicked_lat is None:
                self._show_message("warning","Click Error", "Could not get valid coordinates from click. Try clicking within the displayed GeoTIFF area.")
                return
            self.pitch_line_points.append((clicked_lat, clicked_lon))
            if len(self.pitch_line_points) == 1:
                self.pick_pitch_line_btn.setText("Drawing Pitch Line: Click End Point")
                # Start drawing temporary line
                self._temp_line = self.ax.plot([clicked_lon, clicked_lon], [clicked_lat, clicked_lat], color='orange', linestyle='--', linewidth=2)[0]
                self._temp_line_start = (clicked_lat, clicked_lon)
                self._temp_line_motion_cid = self.canvas.mpl_connect('motion_notify_event', self._on_temp_line_motion)
                self.canvas.draw_idle()
                # Remove info text if present
                if hasattr(self, 'pitch_line_info_text') and self.pitch_line_info_text is not None:
                    self.pitch_line_info_text.set_visible(False)
                    self.pitch_line_info_text = None
                    self.canvas.draw_idle()
            elif len(self.pitch_line_points) == 2:
                self.pick_pitch_line_btn.setText("Draw a Pitch Line")
                self.pick_pitch_line_mode = False
                self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
                # Remove temporary line
                if hasattr(self, '_temp_line') and self._temp_line in self.ax.lines:
                    self._temp_line.remove()
                    del self._temp_line
                self.canvas.mpl_disconnect(self._temp_line_motion_cid)
                self._plot_survey_plan(preserve_view_limits=True)  # Redraw to show pitch line
                # Make "Draw a Pitch Line" normal, "Add Heading Lines" and "Draw a Roll Line" bold and orange
                if hasattr(self, 'pick_pitch_line_btn'):
                    self.pick_pitch_line_btn.setStyleSheet("")  # Reset to normal
                if hasattr(self, 'add_heading_lines_btn'):
                    self.add_heading_lines_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")
                if hasattr(self, 'pick_roll_line_btn'):
                    self.pick_roll_line_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")
                # Remove info text if present
                if hasattr(self, 'pitch_line_info_text') and self.pitch_line_info_text is not None:
                    self.pitch_line_info_text.set_visible(False)
                    self.pitch_line_info_text = None
                    self.canvas.draw_idle()
                # Report pitch line summary in info/error box
                try:
                    import pyproj
                    geod = pyproj.Geod(ellps="WGS84")
                    (lat1, lon1), (lat2, lon2) = self.pitch_line_points
                    az12, az21, dist_m = geod.inv(lon1, lat1, lon2, lat2)
                    dist_nm = dist_m / 1852.0
                    speed_knots = float(self.cal_survey_speed_entry.text()) if self.cal_survey_speed_entry.text() else 8.0
                    speed_m_per_h = speed_knots * 1852
                    time_hours = dist_m / speed_m_per_h if speed_m_per_h > 0 else 0
                    time_minutes = time_hours * 60
                    summary = (
                        f"Pitch Line Length: {dist_m:.1f} m\n"
                        f"Pitch Line Length: {dist_nm:.3f} nautical miles\n"
                        f"Pitch Line Azimuth: {az12:.1f}°\n"
                        f"Pitch Line Survey Time: {time_minutes:.1f} min\n"
                        f"Pitch Line Survey Time: {time_hours:.2f} hr"
                    )
                    self.set_cal_info_text(summary)
                except Exception as e:
                    self.set_cal_info_text(f"Error calculating pitch line summary: {e}")
                # self._show_message("info","Pitch Line Picked", f"Pitch line defined from\nStart: {self.pitch_line_points[0]}\nEnd: {self.pitch_line_points[1]}")
                self._update_cal_line_offset_from_pitch_line()
                self._draw_pitch_line_profile()
                self._update_cal_export_name_from_pitch_line()
                self._update_cal_line_times()
                self._update_pitch_line_button_states()
            return  # Do not process as center pick if in pitch line mode
        if self.pick_roll_line_mode:
            if event.inaxes != self.ax:
                print("DEBUG: Clicked outside plot axes. Click ignored.")
                return
            clicked_lon = event.xdata
            clicked_lat = event.ydata
            if clicked_lon is None or clicked_lat is None:
                self._show_message("warning","Click Error", "Could not get valid coordinates from click. Try clicking within the displayed GeoTIFF area.")
                return
            self.roll_line_points.append((clicked_lat, clicked_lon))
            if len(self.roll_line_points) == 1:
                self.pick_roll_line_btn.setText("Drawing Roll Line: Click End Point")
                # Start drawing temporary line
                self._temp_line = self.ax.plot([clicked_lon, clicked_lon], [clicked_lat, clicked_lat], color='purple', linestyle='--', linewidth=2)[0]
                self._temp_line_start = (clicked_lat, clicked_lon)
                self._temp_line_motion_cid = self.canvas.mpl_connect('motion_notify_event', self._on_temp_line_motion)
                self.canvas.draw_idle()
            elif len(self.roll_line_points) == 2:
                self.pick_roll_line_btn.setText("Draw a Roll Line")
                self.pick_roll_line_mode = False
                self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
                # Remove temporary line
                if hasattr(self, '_temp_line') and self._temp_line in self.ax.lines:
                    self._temp_line.remove()
                    del self._temp_line
                self.canvas.mpl_disconnect(self._temp_line_motion_cid)
                self._plot_survey_plan(preserve_view_limits=True)  # Redraw to show roll line
                # Make "Draw a Roll Line" button normal after drawing
                if hasattr(self, 'pick_roll_line_btn'):
                    self.pick_roll_line_btn.setStyleSheet("")  # Reset to normal
                # Report roll line summary in info/error box
                try:
                    import pyproj
                    geod = pyproj.Geod(ellps="WGS84")
                    (lat1, lon1), (lat2, lon2) = self.roll_line_points
                    az12, az21, dist_m = geod.inv(lon1, lat1, lon2, lat2)
                    dist_nm = dist_m / 1852.0
                    speed_knots = float(self.cal_survey_speed_entry.text()) if self.cal_survey_speed_entry.text() else 8.0
                    speed_m_per_h = speed_knots * 1852
                    time_hours = dist_m / speed_m_per_h if speed_m_per_h > 0 else 0
                    time_minutes = time_hours * 60
                    summary = (
                        f"Roll Line Length: {dist_m:.1f} m\n"
                        f"Roll Line Length: {dist_nm:.3f} nautical miles\n"
                        f"Roll Line Azimuth: {az12:.1f}°\n"
                        f"Roll Line Survey Time: {time_minutes:.1f} min\n"
                        f"Roll Line Survey Time: {time_hours:.2f} hr"
                    )
                    self.set_cal_info_text(summary)
                except Exception as e:
                    self.set_cal_info_text(f"Error calculating roll line summary: {e}")
                self._update_cal_line_times()
                self._update_roll_line_button_states()
            return  # Do not process as center pick if in roll line mode
        if not self.pick_center_mode:
            print("DEBUG: Pick center mode is OFF. Click ignored.")
            return

        if event.inaxes != self.ax:
            print("DEBUG: Clicked outside plot axes. Click ignored.")
            return

        clicked_lon = event.xdata
        clicked_lat = event.ydata

        print(f"DEBUG: Clicked raw coordinates (Lon: {clicked_lon}, Lat: {clicked_lat})")

        if clicked_lon is None or clicked_lat is None:
            self._show_message("warning","Click Error",
                                   "Could not get valid coordinates from click. Try clicking within the displayed GeoTIFF area.")
            print("DEBUG: Clicked coordinates are None. Returning.")
            return

        # Show slope information if in shaded relief or slope visualization mode
        if hasattr(self, 'hillshade_vs_slope_viz_mode') and self.hillshade_vs_slope_viz_mode in ["hillshade", "slope_viz"]:
            elevation, slope = self._calculate_slope_at_point(clicked_lat, clicked_lon)
            if elevation is not None and slope is not None:
                info_msg = f"Lat: {clicked_lat:.6f}, Lon: {clicked_lon:.6f}, Slope: {slope:.1f}°, Depth: {abs(elevation):.1f}m"
                print(f"INFO: {info_msg}")

        self.central_lat_entry.setText(f"{clicked_lat:.6f}")
        self.central_lon_entry.setText(f"{clicked_lon:.6f}")
        print(f"DEBUG: Central Lat/Lon entries updated to: {clicked_lat:.6f}, {clicked_lon:.6f}")

        export_name_to_set = None

        # Read Z-value from original GeoTIFF using the appropriate transformer
        if self.geotiff_dataset_original:
            try:
                # Transform clicked WGS84 point to original GeoTIFF CRS
                if self.wgs84_to_geotiff_transformer:
                    src_x, src_y = self.wgs84_to_geotiff_transformer.transform(clicked_lon, clicked_lat)
                    print(
                        f"DEBUG: Transformed WGS84 ({clicked_lon}, {clicked_lat}) to GeoTIFF CRS: (X: {src_x}, Y: {src_y})")
                else:  # No reprojection was done, assume original is already WGS84
                    src_x, src_y = clicked_lon, clicked_lat
                    print("DEBUG: GeoTIFF is already WGS84, no transformation needed for picking original pixel.")

                # Get row, col in original dataset
                row, col = rowcol(self.geotiff_dataset_original.transform, src_x, src_y)
                print(f"DEBUG: Calculated row: {row}, col: {col} in original GeoTIFF.")
                print(
                    f"DEBUG: Original GeoTIFF dimensions: height={self.geotiff_dataset_original.height}, width={self.geotiff_dataset_original.width}")

                # Check if row, col are within bounds of the ORIGINAL GeoTIFF
                if 0 <= row < self.geotiff_dataset_original.height and 0 <= col < self.geotiff_dataset_original.width:
                    # Make sure to read from the correct band, usually band 1
                    z_value = self.geotiff_dataset_original.read(1)[row, col]
                    z_value = float(z_value)  # Ensure float type for calculations
                    print(f"DEBUG: Raw Z-value from original GeoTIFF at pixel ({row}, {col}): {z_value}")

                    # Check if Z-value is valid and represents depth (non-positive elevation)
                    # Using abs(z_value) for multipliers to ensure positive length/distance
                    if not np.isnan(z_value) and z_value <= 0:

                        # Apply Distance Between Lines Multiplier
                        dist_multiplier = self.dist_between_lines_multiplier
                        # Use absolute value of Z for calculation, as depth is positive for multipliers
                        calculated_dist_between_lines = abs(z_value) * dist_multiplier
                        print(
                            f"DEBUG: Dist Multiplier: {dist_multiplier}, Calculated Dist Between Lines: {calculated_dist_between_lines:.2f}")

                        if calculated_dist_between_lines > 0:  # Ensure distance is positive
                            self.dist_between_lines_entry.setText(f"{calculated_dist_between_lines:.2f}")
                            # Set Crossline Lead-in/out to 0.2 * Distance Between Lines
                            bisect_lead = 0.2 * calculated_dist_between_lines
                            self.bisect_lead_entry.setText(f"{bisect_lead:.2f}")
                            # Set Export Name to 'Reference_' + int(distance between lines) + 'm_' + int(heading) + 'deg'
                            export_name_to_set = f"Reference_{int(calculated_dist_between_lines)}m_{int(float(self.heading_entry.text()))}deg"
                        else:
                            self._show_message("warning","Input Warning",
                                                   "Calculated Distance Between Lines is not positive. Not setting Distance automatically.")
                            print(
                                "DEBUG: Warning: Calculated Distance Between Lines is not positive after multiplication.")

                        # Apply Line Length Multiplier
                        len_multiplier = self.line_length_multiplier
                        # Use absolute value of Z for calculation
                        calculated_line_length = abs(z_value) * len_multiplier
                        print(
                            f"DEBUG: Length Multiplier: {len_multiplier}, Calculated Line Length: {calculated_line_length:.2f}")
                        if calculated_line_length > 0:  # Ensure line length is positive
                            self.line_length_entry.setText(f"{calculated_line_length:.2f}")
                        else:
                            self._show_message("warning","Input Warning",
                                                   "Calculated Line Length is not positive. Not setting Line Length automatically.")
                            print("DEBUG: Warning: Calculated Line Length is not positive after multiplication.")
                    else:
                        self._show_message("warning","GeoTIFF Warning",
                                               f"Invalid Z-value ({z_value}) at clicked location (NaN or positive elevation). Line spacing and length not set.")
                        print(
                            f"DEBUG: Warning: Z-value is NaN or positive ({z_value}). Conditions for setting parameters not met.")
                else:
                    self._show_message("warning","GeoTIFF Warning",
                                           f"Clicked location (row={row}, col={col}) is outside GeoTIFF bounds (height={self.geotiff_dataset_original.height}, width={self.geotiff_dataset_original.width}). Line spacing and length not set.")
                    print(f"DEBUG: Warning: Clicked location ({row}, {col}) is outside GeoTIFF bounds.")

            except Exception as e:
                self._show_message("error","GeoTIFF Read Error",
                                     f"Failed to read Z-value from GeoTIFF: {e}. See console for more details.")
                print(f"DEBUG ERROR: Failed to read Z-value from GeoTIFF: {e}")
                traceback.print_exc()  # Print full traceback for deeper debugging
        else:
            print("DEBUG: self.geotiff_dataset_original is None. GeoTIFF not loaded.")

        # Set export name if calculated
        if export_name_to_set is not None:
            self.export_name_entry.clear()
            self.export_name_entry.setText(export_name_to_set)

        # Reset "Pick Center from GeoTIFF" button to normal style after center is selected
        if hasattr(self, 'pick_center_btn'):
            self.pick_center_btn.setStyleSheet("")
        
        self._toggle_pick_center_mode()  # Turn off pick mode after click
        self._generate_and_plot(show_success_dialog=False)  # Re-generate plot with new center, suppress dialog

    def _on_scroll(self, event):
        """Callback for mouse scroll event to zoom in/out from the center of the window."""
        if event.inaxes != self.ax:
            return

        # Set rapid zoom mode to prevent resolution updates during continuous zooming
        self._rapid_zoom_mode = True
        if not hasattr(self, '_rapid_zoom_timer'):
            self._rapid_zoom_timer = QTimer()
            self._rapid_zoom_timer.setSingleShot(True)
            self._rapid_zoom_timer.timeout.connect(self._clear_rapid_zoom_mode)
        self._rapid_zoom_timer.stop()
        self._rapid_zoom_timer.start(2000)
        
        # Get current view limits (only once)
        cur_xlim = self.ax.get_xlim()
        cur_ylim = self.ax.get_ylim()
        
        # Store the current zoom level before changing it (only if GeoTIFF is loaded)
        old_zoom_level = None
        if (self.geotiff_dataset_original is not None and 
            self.dynamic_resolution_enabled and
            hasattr(self, 'geotiff_extent') and self.geotiff_extent is not None):
            full_width = self.geotiff_extent[1] - self.geotiff_extent[0]
            full_height = self.geotiff_extent[3] - self.geotiff_extent[2]
            current_width = cur_xlim[1] - cur_xlim[0]
            current_height = cur_ylim[1] - cur_ylim[0]
            width_ratio = current_width / full_width if full_width > 0 else 1.0
            height_ratio = current_height / full_height if full_height > 0 else 1.0
            old_zoom_level = min(width_ratio, height_ratio)
            
            # Initialize cumulative zoom change tracking
            if not hasattr(self, '_cumulative_zoom_change'):
                self._cumulative_zoom_change = 0.0
                self._last_zoom_level = old_zoom_level
                self._zoom_operation_count = 0
        
        # Calculate center of current view
        center_x = (cur_xlim[0] + cur_xlim[1]) / 2
        center_y = (cur_ylim[0] + cur_ylim[1]) / 2
        
        # Calculate current ranges
        cur_xrange = (cur_xlim[1] - cur_xlim[0]) * 0.5
        cur_yrange = (cur_ylim[1] - cur_ylim[0]) * 0.5

        # Handle scroll events - check both event.button and event.step for compatibility
        if hasattr(event, 'step') and event.step != 0:
            # Use event.step if available (positive = up/zoom in, negative = down/zoom out)
            if event.step > 0:
                scale_factor = 0.9  # Zoom in
            else:
                scale_factor = 1.1  # Zoom out
        elif event.button == 'up':
            # Zoom in
            scale_factor = 0.9  # Smaller factor zooms in more
        elif event.button == 'down':
            # Zoom out
            scale_factor = 1.1  # Larger factor zooms out more
        else:
            return  # Not a scroll wheel event

        # Set new limits, centered at the center of the current view
        new_xlim = (center_x - cur_xrange * scale_factor, center_x + cur_xrange * scale_factor)
        new_ylim = (center_y - cur_yrange * scale_factor, center_y + cur_yrange * scale_factor)
        
        # Set new limits
        self.ax.set_xlim(new_xlim)
        self.ax.set_ylim(new_ylim)
        
        # Set aspect ratio with 'datalim' to keep axes position fixed (like original Tkinter version)
        # Don't manually set position - let matplotlib handle it automatically
        try:
            center_lat = (new_ylim[0] + new_ylim[1]) / 2.0
            center_lat_rad = np.radians(center_lat)
            aspect_ratio = 1.0 / np.cos(center_lat_rad)
            aspect_ratio = np.clip(aspect_ratio, 0.1, 10.0)
            # Use 'datalim' to keep axes position fixed - it will adjust limits if needed
            self.ax.set_aspect(aspect_ratio, adjustable='datalim')
        except:
            pass
        
        # Update zoom level for dynamic resolution loading (only if GeoTIFF is loaded and dynamic resolution enabled)
        if (old_zoom_level is not None and 
            self.geotiff_dataset_original is not None and 
            self.dynamic_resolution_enabled and
            self.geotiff_extent is not None):
            
            # Calculate new zoom level based on the updated view
            full_width = self.geotiff_extent[1] - self.geotiff_extent[0]
            full_height = self.geotiff_extent[3] - self.geotiff_extent[2]
            new_width = new_xlim[1] - new_xlim[0]
            new_height = new_ylim[1] - new_ylim[0]
            
            # New zoom level is the ratio of current view to full extent
            width_ratio = new_width / full_width if full_width > 0 else 1.0
            height_ratio = new_height / full_height if full_height > 0 else 1.0
            new_zoom_level = min(width_ratio, height_ratio)
            new_zoom_level = max(0.01, min(10.0, new_zoom_level))
            
            # Track cumulative zoom changes
            zoom_change = abs(new_zoom_level - old_zoom_level)
            self._cumulative_zoom_change += zoom_change
            self._zoom_operation_count += 1
            
            # Check if cumulative zoom change is significant enough to trigger resolution update
            if self._cumulative_zoom_change > 0.15 or self._zoom_operation_count >= 10:
                self.geotiff_zoom_level = new_zoom_level
                self._cumulative_zoom_change = 0.0
                self._zoom_operation_count = 0
                self._last_zoom_level = new_zoom_level
                # Use a timer to avoid interfering with continuous zooming
                if not hasattr(self, '_zoom_timer'):
                    self._zoom_timer = QTimer()
                    self._zoom_timer.setSingleShot(True)
                    self._zoom_timer.timeout.connect(self._reload_geotiff_at_current_zoom)
                self._zoom_timer.stop()
                self._zoom_timer.start(1000)

        self._last_user_xlim = new_xlim
        self._last_user_ylim = new_ylim

        # Use draw_idle() for non-blocking updates during zoom
        self.canvas.draw_idle()

    def _zoom_to_plan(self):
        all_lats = []
        all_lons = []

        # Add survey line points
        for line in self.survey_lines_data:
            for p in line:
                all_lats.append(p[0])
                all_lons.append(p[1])

        # Add crossline points
        if self.cross_line_data:
            for p in self.cross_line_data:
                all_lats.append(p[0])
                all_lons.append(p[1])

        if all_lats and all_lons:
            min_lat, max_lat = min(all_lats), max(all_lats)
            min_lon, max_lon = min(all_lons), max(all_lons)

            # Add a small buffer (5% of range)
            buffer_lat = (max_lat - min_lat) * 0.05 if (max_lat - min_lat) != 0 else 0.01
            buffer_lon = (max_lon - min_lon) * 0.05 if (max_lon - min_lon) != 0 else 0.01

            xlim = (min_lon - buffer_lon, max_lon + buffer_lon)
            ylim = (min_lat - buffer_lat, max_lat + buffer_lat)
            
            # Store the new limits
            self.current_xlim = xlim
            self.current_ylim = ylim
            self._last_user_xlim = xlim
            self._last_user_ylim = ylim
            
            # Set axis limits
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)
            
            # Update zoom level and reload GeoTIFF data if available
            if self.geotiff_extent is not None and self.geotiff_dataset_original is not None:
                # Calculate zoom level based on the view extent vs full GeoTIFF extent
                if self.dynamic_resolution_enabled:
                    full_width = self.geotiff_extent[1] - self.geotiff_extent[0]
                    full_height = self.geotiff_extent[3] - self.geotiff_extent[2]
                    current_width = xlim[1] - xlim[0]
                    current_height = ylim[1] - ylim[0]
                    
                    width_ratio = current_width / full_width if full_width > 0 else 1.0
                    height_ratio = current_height / full_height if full_height > 0 else 1.0
                    # Use the larger ratio (more zoomed in) to determine resolution
                    new_zoom_level = max(width_ratio, height_ratio)
                    new_zoom_level = max(0.01, min(10.0, new_zoom_level))
                    self.geotiff_zoom_level = new_zoom_level
                
                # Force a small delay to ensure axis limits are applied before loading
                self.canvas.draw()
                
                # Reload GeoTIFF data at appropriate resolution
                success = self._load_geotiff_at_resolution()
                if success:
                    # Ensure the loaded extent covers the full view - if not, we may need to reload
                    if hasattr(self, 'geotiff_extent') and self.geotiff_extent is not None:
                        loaded_left, loaded_right, loaded_bottom, loaded_top = self.geotiff_extent
                        view_left, view_right = xlim
                        view_bottom, view_top = ylim
                        
                        # Check if loaded extent covers view extent (with small tolerance)
                        tolerance = 0.0001  # degrees
                        if (loaded_left > view_left + tolerance or loaded_right < view_right - tolerance or
                            loaded_bottom > view_bottom + tolerance or loaded_top < view_top - tolerance):
                            # Extent doesn't fully cover view - this shouldn't happen, but if it does, reload
                            # The issue might be in _load_geotiff_at_resolution, but for now we'll proceed
                            pass
                    
                    # Replot with new resolution
                    self._plot_survey_plan(preserve_view_limits=True)
                else:
                    self.canvas.draw_idle()
            else:
                self.canvas.draw_idle()
        else:
            # If nothing to zoom to, just reset view
            self.ax.autoscale_view()
            self.canvas.draw_idle()

    def _on_middle_press(self, event):
        # Only handle middle button (button 2)
        if event.button != 2:
            return
        if event.inaxes != self.ax:
            return
        self._pan_start = (event.x, event.y)
        # Store data coordinates as fallback
        if event.xdata is not None and event.ydata is not None:
            self._pan_start_data = (event.xdata, event.ydata)
        self._orig_xlim = self.ax.get_xlim()
        self._orig_ylim = self.ax.get_ylim()
        # Store the correct axes position - check if colorbar exists NOW
        has_colorbar = (hasattr(self, 'elevation_colorbar') and self.elevation_colorbar is not None) or \
                       (hasattr(self, 'slope_colorbar') and self.slope_colorbar is not None)
        if has_colorbar and hasattr(self, '_axes_pos_with_colorbar'):
            # Use the position that was set when colorbar was created
            self._pan_ax_pos = list(self._axes_pos_with_colorbar)  # (left, bottom, width, height)
        else:
            # No colorbar - use full width to fill available space
            self._pan_ax_pos = [0.08, 0.08, 0.87, 0.87]  # [left, bottom, width, height]
        # Set panning mode to prevent resolution updates during pan
        self._panning_mode = True
        # Set a safety timeout to clear panning mode if it gets stuck
        if hasattr(self, '_panning_timeout'):
            self._panning_timeout.stop()
        self._panning_timeout = QTimer()
        self._panning_timeout.setSingleShot(True)
        self._panning_timeout.timeout.connect(self._clear_panning_mode)
        self._panning_timeout.start(5000)

    def _on_middle_motion(self, event):
        if not hasattr(self, '_pan_start'):
            return
        if event.inaxes != self.ax:
            return
        
        # Convert pixel movement to data coordinates
        # Match the original Tkinter version exactly
        inv = self.ax.transData.inverted()
        try:
            x0, y0 = inv.transform((self._pan_start[0], self._pan_start[1]))
            x1, y1 = inv.transform((event.x, event.y))
            dx_data = x0 - x1  # start - current (matches original)
            dy_data = y0 - y1  # start - current (matches original)
        except:
            # Fallback: use data coordinates if available
            if event.xdata is not None and event.ydata is not None and hasattr(self, '_pan_start_data'):
                dx_data = self._pan_start_data[0] - event.xdata
                dy_data = self._pan_start_data[1] - event.ydata
            else:
                return
        
        # Update limits
        # For x: add dx_data (when dragging right, dx_data is negative, so adding it increases limits, shows right)
        # For y: add dy_data (when dragging up, dy_data is negative, so adding it increases limits, shows up)
        new_xlim = (self._orig_xlim[0] + dx_data, self._orig_xlim[1] + dx_data)
        new_ylim = (self._orig_ylim[0] + dy_data, self._orig_ylim[1] + dy_data)
        self.ax.set_xlim(new_xlim)
        self.ax.set_ylim(new_ylim)
        
        # Update stored user limits so reload uses correct position
        self._last_user_xlim = new_xlim
        self._last_user_ylim = new_ylim
        
        # Set aspect ratio with 'datalim' to keep axes position fixed (like original Tkinter version)
        # Don't manually set position - let matplotlib handle it automatically
        try:
            center_lat = (new_ylim[0] + new_ylim[1]) / 2.0
            center_lat_rad = np.radians(center_lat)
            aspect_ratio = 1.0 / np.cos(center_lat_rad)
            aspect_ratio = np.clip(aspect_ratio, 0.1, 10.0)
            # Use 'datalim' to keep axes position fixed - it will adjust limits if needed
            self.ax.set_aspect(aspect_ratio, adjustable='datalim')
        except:
            pass
        
        # Use draw_idle() for non-blocking updates during panning
        self.canvas.draw_idle()

    def _on_middle_release(self, event):
        # Only handle middle button (button 2)
        if event.button != 2:
            return
        if hasattr(self, '_pan_start'):
            del self._pan_start
            del self._orig_xlim
            del self._orig_ylim
            
            # Store current view limits immediately after panning
            current_xlim = self.ax.get_xlim()
            current_ylim = self.ax.get_ylim()
            self._current_view_limits = (current_xlim, current_ylim)
            
            # Clear panning mode and timeout
            if hasattr(self, '_panning_mode'):
                del self._panning_mode
            if hasattr(self, '_panning_timeout'):
                self._panning_timeout.stop()
                del self._panning_timeout
            
            # Always reload GeoTIFF after panning to ensure data is loaded for the new view area
            # This is especially important when dynamic resolution is enabled, as the loaded data
            # might only cover the previous view area
            if self.geotiff_dataset_original is not None:
                # Get current view limits
                xlim = self.ax.get_xlim()
                ylim = self.ax.get_ylim()
                
                # Check if current view is outside the loaded GeoTIFF extent
                # If so, we definitely need to reload
                needs_reload = False
                if hasattr(self, 'geotiff_extent') and self.geotiff_extent is not None:
                    # Check if view extends beyond loaded extent (with small tolerance)
                    tolerance = 0.01  # degrees
                    if (xlim[0] < self.geotiff_extent[0] - tolerance or
                        xlim[1] > self.geotiff_extent[1] + tolerance or
                        ylim[0] < self.geotiff_extent[2] - tolerance or
                        ylim[1] > self.geotiff_extent[3] + tolerance):
                        needs_reload = True
                
                # Also reload if dynamic resolution is enabled and zoom level changed
                if self.dynamic_resolution_enabled and hasattr(self, 'geotiff_extent') and self.geotiff_extent is not None:
                    full_width = self.geotiff_extent[1] - self.geotiff_extent[0]
                    full_height = self.geotiff_extent[3] - self.geotiff_extent[2]
                    current_width = xlim[1] - xlim[0]
                    current_height = ylim[1] - ylim[0]
                    
                    width_ratio = current_width / full_width if full_width > 0 else 1.0
                    height_ratio = current_height / full_height if full_height > 0 else 1.0
                    new_zoom_level = min(width_ratio, height_ratio)
                    
                    # Update zoom level if it changed significantly
                    if abs(new_zoom_level - self.geotiff_zoom_level) > 0.2:
                        self.geotiff_zoom_level = new_zoom_level
                        needs_reload = True
                
                # Reload if needed
                if needs_reload:
                    # Use a timer to avoid reloading too frequently
                    if hasattr(self, '_pan_timer'):
                        self._pan_timer.stop()
                    self._pan_timer = QTimer()
                    self._pan_timer.setSingleShot(True)
                    self._pan_timer.timeout.connect(self._reload_geotiff_at_current_zoom)
                    self._pan_timer.start(300)  # Short delay to batch rapid panning
                else:
                    # Even if zoom level didn't change, we might need to reload if view moved
                    # outside the currently loaded data extent. Check if we have loaded data array.
                    if hasattr(self, 'geotiff_data_array') and self.geotiff_data_array is not None:
                        # Check if the current view center is significantly outside the loaded extent
                        view_center_x = (xlim[0] + xlim[1]) / 2
                        view_center_y = (ylim[0] + ylim[1]) / 2
                        if hasattr(self, 'geotiff_extent') and self.geotiff_extent is not None:
                            loaded_center_x = (self.geotiff_extent[0] + self.geotiff_extent[1]) / 2
                            loaded_center_y = (self.geotiff_extent[2] + self.geotiff_extent[3]) / 2
                            # If view center moved significantly (more than 50% of view width/height), reload
                            view_width = xlim[1] - xlim[0]
                            view_height = ylim[1] - ylim[0]
                            if (abs(view_center_x - loaded_center_x) > view_width * 0.5 or
                                abs(view_center_y - loaded_center_y) > view_height * 0.5):
                                if hasattr(self, '_pan_timer'):
                                    self._pan_timer.stop()
                                self._pan_timer = QTimer()
                                self._pan_timer.setSingleShot(True)
                                self._pan_timer.timeout.connect(self._reload_geotiff_at_current_zoom)
                                self._pan_timer.start(300)

    def _zoom_to_pitch_line(self):
        if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2:
            (lat1, lon1), (lat2, lon2) = self.pitch_line_points
            min_lat, max_lat = min(lat1, lat2), max(lat1, lat2)
            min_lon, max_lon = min(lon1, lon2), max(lon1, lon2)
            # Add a small buffer (5% of range or 0.01 if range is 0)
            buffer_lat = (max_lat - min_lat) * 0.05 if (max_lat - min_lat) != 0 else 0.01
            buffer_lon = (max_lon - min_lon) * 0.05 if (max_lon - min_lon) != 0 else 0.01
            self.ax.set_xlim(min_lon - buffer_lon, max_lon + buffer_lon)
            self.ax.set_ylim(min_lat - buffer_lat, max_lat + buffer_lat)
            self.canvas.draw_idle()

    def _zoom_to_line(self):
        """Zoom to the line planning points."""
        if not hasattr(self, 'line_planning_points') or len(self.line_planning_points) < 2:
            self._show_message("info","Zoom to Line", "No line planning points to zoom to.")
            return
        
        lats = [p[0] for p in self.line_planning_points]
        lons = [p[1] for p in self.line_planning_points]
        min_lat, max_lat = min(lats), max(lats)
        min_lon, max_lon = min(lons), max(lons)
        
        # Add a small buffer (5% of range or 0.01 if range is 0)
        buffer_lat = (max_lat - min_lat) * 0.05 if (max_lat - min_lat) != 0 else 0.01
        buffer_lon = (max_lon - min_lon) * 0.05 if (max_lon - min_lon) != 0 else 0.01
        
        xlim = (min_lon - buffer_lon, max_lon + buffer_lon)
        ylim = (min_lat - buffer_lat, max_lat + buffer_lat)
        
        # Store the new limits
        self.current_xlim = xlim
        self.current_ylim = ylim
        self._last_user_xlim = xlim
        self._last_user_ylim = ylim
        
        # Set axis limits
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        
        # Update zoom level and reload GeoTIFF data if available
        if self.geotiff_extent is not None and self.geotiff_dataset_original is not None:
            # Calculate zoom level based on the view extent vs full GeoTIFF extent
            if self.dynamic_resolution_enabled:
                full_width = self.geotiff_extent[1] - self.geotiff_extent[0]
                full_height = self.geotiff_extent[3] - self.geotiff_extent[2]
                current_width = xlim[1] - xlim[0]
                current_height = ylim[1] - ylim[0]
                
                width_ratio = current_width / full_width if full_width > 0 else 1.0
                height_ratio = current_height / full_height if full_height > 0 else 1.0
                # Use the larger ratio (more zoomed in) to determine resolution
                new_zoom_level = max(width_ratio, height_ratio)
                new_zoom_level = max(0.01, min(10.0, new_zoom_level))
                self.geotiff_zoom_level = new_zoom_level
            
            # Reload GeoTIFF data at appropriate resolution
            success = self._load_geotiff_at_resolution()
            if success:
                # Replot with new resolution
                self._plot_survey_plan(preserve_view_limits=True)
            else:
                self.canvas.draw_idle()
        else:
            self.canvas.draw_idle()

    def _zoom_to_any_lines(self):
        # Collect all points from pitch line and heading lines
        points = []
        if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2:
            points.extend(self.pitch_line_points)
        if hasattr(self, 'heading_lines') and len(self.heading_lines) == 2:
            for line in self.heading_lines:
                points.extend(line)
        if hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2:
            points.extend(self.roll_line_points)
        if not points:
            self._show_message("info","Zoom to Lines", "No lines to zoom to.")
            return
        lats = [p[0] for p in points]
        lons = [p[1] for p in points]
        min_lat, max_lat = min(lats), max(lats)
        min_lon, max_lon = min(lons), max(lons)
        # Add a small buffer (5% of range or 0.01 if range is 0)
        buffer_lat = (max_lat - min_lat) * 0.05 if (max_lat - min_lat) != 0 else 0.01
        buffer_lon = (max_lon - min_lon) * 0.05 if (max_lon - min_lon) != 0 else 0.01
        
        xlim = (min_lon - buffer_lon, max_lon + buffer_lon)
        ylim = (min_lat - buffer_lat, max_lat + buffer_lat)
        
        # Store the new limits
        self.current_xlim = xlim
        self.current_ylim = ylim
        self._last_user_xlim = xlim
        self._last_user_ylim = ylim
        
        # Set axis limits
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        
        # Update zoom level and reload GeoTIFF data if available
        if self.geotiff_extent is not None and self.geotiff_dataset_original is not None:
            # Calculate zoom level based on the view extent vs full GeoTIFF extent
            if self.dynamic_resolution_enabled:
                full_width = self.geotiff_extent[1] - self.geotiff_extent[0]
                full_height = self.geotiff_extent[3] - self.geotiff_extent[2]
                current_width = xlim[1] - xlim[0]
                current_height = ylim[1] - ylim[0]
                
                width_ratio = current_width / full_width if full_width > 0 else 1.0
                height_ratio = current_height / full_height if full_height > 0 else 1.0
                # Use the larger ratio (more zoomed in) to determine resolution
                new_zoom_level = max(width_ratio, height_ratio)
                new_zoom_level = max(0.01, min(10.0, new_zoom_level))
                self.geotiff_zoom_level = new_zoom_level
            
            # Reload GeoTIFF data at appropriate resolution
            success = self._load_geotiff_at_resolution()
            if success:
                # Replot with new resolution
                self._plot_survey_plan(preserve_view_limits=True)
            else:
                self.canvas.draw_idle()
        else:
            self.canvas.draw_idle()

    def _zoom_to_pitch_heading_lines(self):
        # Collect all points from pitch line and heading lines
        points = []
        if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2:
            points.extend(self.pitch_line_points)
        if hasattr(self, 'heading_lines') and len(self.heading_lines) == 2:
            for line in self.heading_lines:
                points.extend(line)
        if not points:
            self._show_message("info","Zoom To Pitch/Heading", "No pitch or heading lines to zoom to.")
            return
        lats = [p[0] for p in points]
        lons = [p[1] for p in points]
        min_lat, max_lat = min(lats), max(lats)
        min_lon, max_lon = min(lons), max(lons)
        # Add a small buffer (5% of range or 0.01 if range is 0)
        buffer_lat = (max_lat - min_lat) * 0.05 if (max_lat - min_lat) != 0 else 0.01
        buffer_lon = (max_lon - min_lon) * 0.05 if (max_lon - min_lon) != 0 else 0.01
        self.ax.set_xlim(min_lon - buffer_lon, max_lon + buffer_lon)
        self.ax.set_ylim(min_lat - buffer_lat, max_lat + buffer_lat)
        self.canvas.draw_idle()

    def _zoom_to_roll_line(self):
        # Zoom to the area around the roll line only
        if not hasattr(self, 'roll_line_points') or len(self.roll_line_points) != 2:
            self._show_message("info","Zoom to Roll", "No roll line to zoom to.")
            return
        lats = [p[0] for p in self.roll_line_points]
        lons = [p[1] for p in self.roll_line_points]
        min_lat, max_lat = min(lats), max(lats)
        min_lon, max_lon = min(lons), max(lons)
        # Add a small buffer (5% of range or 0.01 if range is 0)
        buffer_lat = (max_lat - min_lat) * 0.05 if (max_lat - min_lat) != 0 else 0.01
        buffer_lon = (max_lon - min_lon) * 0.05 if (max_lon - min_lon) != 0 else 0.01
        self.ax.set_xlim(min_lon - buffer_lon, max_lon + buffer_lon)
        self.ax.set_ylim(min_lat - buffer_lat, max_lat + buffer_lat)
        self.canvas.draw_idle()

"""
Survey plan plotting: generate plan, plot survey lines/GeoTIFF/contours, clear plot, colorbars.
"""
import traceback
import numpy as np
from matplotlib.colors import LightSource
from matplotlib.ticker import FuncFormatter

from sat_planner.constants import GEOSPATIAL_LIBS_AVAILABLE, pyproj, CRSError
from sat_planner.utils_geo import decimal_degrees_to_ddm


class PlottingMixin:
    """Mixin providing _generate_and_plot, _plot_survey_plan, _clear_plot, _remove_colorbar."""

    def _generate_and_plot(self, show_success_dialog=True, auto_zoom=True):
        # Clear the elevation profile immediately when generating a new plan
        self._draw_crossline_profile()
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded. Cannot generate plot.")
            return

        is_valid, values = self._validate_inputs()
        if not is_valid:
            return

        # self._update_survey_time_display(values['line_length'], values['survey_speed'])  # Moved below
        self._clear_plot(full_clear=False)  # Clear only plot data, keep geotiff if loaded

        try:
            central_lat = values['central_lat']
            central_lon = values['central_lon']
            line_length = values['line_length']
            heading = values['heading']
            dist_between_lines = values['dist_between_lines']
            num_lines = values['num_lines']
            bisect_lead = values['bisect_lead']
            offset_direction = values['offset_direction']

            self.central_point_coords = (central_lat, central_lon)

            # 1. Define a local projected CRS (e.g., UTM zone) for calculations
            # Determine UTM zone from central_lon
            utm_zone = int((central_lon + 180) / 6) + 1
            if central_lat >= 0:
                utm_crs = f"EPSG:326{utm_zone}"  # Northern Hemisphere
            else:
                utm_crs = f"EPSG:327{utm_zone}"  # Southern Hemisphere

            print(f"DEBUG: Using UTM CRS: {utm_crs}")

            # Create a transformer from WGS84 (EPSG:4326) to the local UTM CRS
            self.local_proj_transformer = pyproj.Transformer.from_crs(
                "EPSG:4326", utm_crs, always_xy=True
            )
            # Create a transformer back from local UTM to WGS84
            self.wgs84_transformer = pyproj.Transformer.from_crs(
                utm_crs, "EPSG:4326", always_xy=True
            )

            # Convert central point to local projected coordinates (easting, northing)
            central_easting, central_northing = self.local_proj_transformer.transform(central_lon, central_lat)
            print(f"DEBUG: Central point in UTM coordinates: ({central_easting}, {central_northing})")

            self.survey_lines_data = []
            self.cross_line_data = []

            geod = pyproj.Geod(ellps="WGS84")
            # Calculate the index of the center line
            center_index = (num_lines - 1) / 2
            for i in range(num_lines):
                # Calculate offset distance from the central point
                offset_distance = (i - center_index) * dist_between_lines
                # Calculate perpendicular azimuth for offset direction
                if offset_direction == "North":
                    perp_azimuth = heading + 90
                else:  # "South"
                    perp_azimuth = heading - 90
                # Offset the center point geodetically
                lon_center, lat_center, _ = geod.fwd(central_lon, central_lat, perp_azimuth, offset_distance)

                # Use geod.fwd to get endpoints at correct heading
                lon1, lat1, _ = geod.fwd(lon_center, lat_center, heading + 180, line_length / 2)
                lon2, lat2, _ = geod.fwd(lon_center, lat_center, heading, line_length / 2)

                self.survey_lines_data.append([(lat1, lon1), (lat2, lon2)])

                # Debug: Print actual azimuth between endpoints (from start to end)
                try:
                    fwd_azimuth, back_azimuth, distance = geod.inv(lon1, lat1, lon2, lat2)
                    print(f"DEBUG: Survey line {i+1}: Start=({lat1:.6f}, {lon1:.6f}), End=({lat2:.6f}, {lon2:.6f}), Calculated Azimuth: {fwd_azimuth:.2f} deg, Input Heading: {heading}")
                except Exception as e:
                    print(f"DEBUG: Could not calculate azimuth for line {i+1}: {e}")

            # Calculate crossline (perpendicular to main lines, geodetic)
            # Connect midpoints of first and last main survey lines, then extend by lead-in/out
            if len(self.survey_lines_data) >= 2:
                # Get first and last survey lines
                first_line = self.survey_lines_data[0]
                last_line = self.survey_lines_data[-1]

                # Calculate midpoints of first and last lines
                first_mid_lat = (first_line[0][0] + first_line[1][0]) / 2
                first_mid_lon = (first_line[0][1] + first_line[1][1]) / 2
                last_mid_lat = (last_line[0][0] + last_line[1][0]) / 2
                last_mid_lon = (last_line[0][1] + last_line[1][1]) / 2

                # Calculate the azimuth from first midpoint to last midpoint
                _, _, crossline_distance = geod.inv(first_mid_lon, first_mid_lat, last_mid_lon, last_mid_lat)
                crossline_azimuth, _, _ = geod.inv(first_mid_lon, first_mid_lat, last_mid_lon, last_mid_lat)

                # Extend the crossline by lead-in/out distance on each end
                # Start point: extend from last midpoint in forward direction (closer to main survey area)
                lon1, lat1, _ = geod.fwd(last_mid_lon, last_mid_lat, crossline_azimuth, bisect_lead)
                # End point: extend from first midpoint in opposite direction (farther from main survey area)
                lon2, lat2, _ = geod.fwd(first_mid_lon, first_mid_lat, crossline_azimuth + 180, bisect_lead)

                self.cross_line_data = [(lat1, lon1), (lat2, lon2)]
            else:
                # Fallback to original method if only one line
                crossline_azimuth1 = (heading + 90) % 360
                crossline_azimuth2 = (heading + 270) % 360  # Equivalent to heading - 90
                cross_line_length_total = line_length + (2 * bisect_lead)

                # Endpoint 1
                lon1, lat1, _ = geod.fwd(central_lon, central_lat, crossline_azimuth1, cross_line_length_total / 2)
                # Endpoint 2
                lon2, lat2, _ = geod.fwd(central_lon, central_lat, crossline_azimuth2, cross_line_length_total / 2)

                self.cross_line_data = [(lat1, lon1), (lat2, lon2)]

            self._plot_survey_plan()
            if auto_zoom:
                self._zoom_to_plan()  # Auto zoom after plotting

            # --- Report Main Line and Crossline summary in info/error area ---
            try:
                geod = pyproj.Geod(ellps="WGS84")
                # Main Line: use the center line (middle index)
                center_index = (num_lines - 1) // 2
                main_line = self.survey_lines_data[center_index]
                (lat1, lon1), (lat2, lon2) = main_line
                _, _, main_length_m = geod.inv(lon1, lat1, lon2, lat2)
                main_length_km = main_length_m / 1000.0
                main_length_nm = main_length_m / 1852.0
                speed_knots = float(self.survey_speed_entry.text()) if self.survey_speed_entry.text() else 8.0
                speed_m_per_h = speed_knots * 1852
                main_time_h = main_length_m / speed_m_per_h if speed_m_per_h > 0 else 0
                main_time_min = main_time_h * 60
                # Crossline
                (clat1, clon1), (clat2, clon2) = self.cross_line_data
                _, _, cross_length_m = geod.inv(clon1, clat1, clon2, clat2)
                cross_length_km = cross_length_m / 1000.0
                cross_length_nm = cross_length_m / 1852.0
                cross_time_h = cross_length_m / speed_m_per_h if speed_m_per_h > 0 else 0
                cross_time_min = cross_time_h * 60

                # Get number of crossline passes
                try:
                    num_passes = int(self.crossline_passes_entry.text()) if self.crossline_passes_entry.text() else 1
                except (ValueError, AttributeError):
                    num_passes = 1

                # Calculate total survey time including travel between lines
                total_survey_time = self._calculate_total_survey_time()

                # Calculate single line length and heading
                num_main_lines = len(self.survey_lines_data) if self.survey_lines_data else 0
                length_line = ""
                heading_lines = ""

                if num_main_lines > 0:
                    single_line_length_m = total_survey_time['main_lines_distance_m'] / num_main_lines
                    single_line_length_km = single_line_length_m / 1000
                    single_line_length_nm = single_line_length_m / 1852
                    length_line = f"Length of Main Line: {single_line_length_m:.1f} m ({single_line_length_km:.3f} km, {single_line_length_nm:.3f} nm)\n"

                    # Calculate heading for the first main line
                    try:
                        if pyproj is not None:
                            geod = pyproj.Geod(ellps="WGS84")
                            first_line = self.survey_lines_data[0]
                            lat1, lon1 = first_line[0]
                            lat2, lon2 = first_line[1]
                            fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)

                            heading = fwd_az % 360
                            reciprocal_heading = back_az % 360

                            heading_lines = f"Heading: {heading:.1f}°\nReciprocal Heading: {reciprocal_heading:.1f}°\n"
                        else:
                            heading_lines = "Heading: pyproj not available\n"
                    except Exception:
                        heading_lines = "Heading: Unable to calculate\n"

                summary = (
                    f"=== SURVEY DISTANCE BREAKDOWN ===\n"
                    f"{heading_lines}"
                    f"{length_line}"
                    f"Main Lines Total Distance: {total_survey_time['main_lines_distance_m']:.1f} m\n"
                    f"Main Lines Total Distance: {total_survey_time['main_lines_distance_km']:.3f} km\n"
                    f"Main Lines Total Distance: {total_survey_time['main_lines_distance_nm']:.3f} nm\n"
                    f"Travel Between Lines: {total_survey_time['travel_between_lines_distance_m']:.1f} m\n"
                    f"Travel Between Lines: {total_survey_time['travel_between_lines_distance_km']:.3f} km\n"
                    f"Travel Between Lines: {total_survey_time['travel_between_lines_distance_nm']:.3f} nm\n"
                    f"Travel to Crossline: {total_survey_time['travel_to_crossline_distance_m']:.1f} m\n"
                    f"Travel to Crossline: {total_survey_time['travel_to_crossline_distance_km']:.3f} km\n"
                    f"Travel to Crossline: {total_survey_time['travel_to_crossline_distance_nm']:.3f} nm\n"
                    f"Crossline (per pass): {total_survey_time['crossline_single_pass_distance_m']:.1f} m\n"
                    f"Crossline (per pass): {total_survey_time['crossline_single_pass_distance_km']:.3f} km\n"
                    f"Crossline (per pass): {total_survey_time['crossline_single_pass_distance_nm']:.3f} nm\n"
                    f"Crossline Passes: {total_survey_time['num_crossline_passes']}\n"
                    f"Crossline Total Distance: {total_survey_time['crossline_total_distance_m']:.1f} m\n"
                    f"Crossline Total Distance: {total_survey_time['crossline_total_distance_km']:.3f} km\n"
                    f"Crossline Total Distance: {total_survey_time['crossline_total_distance_nm']:.3f} nm\n"
                    f"\n=== TOTAL SURVEY DISTANCE ===\n"
                    f"Total Survey Distance: {total_survey_time['total_distance_m']:.1f} m\n"
                    f"Total Survey Distance: {total_survey_time['total_distance_km']:.3f} km\n"
                    f"Total Survey Distance: {total_survey_time['total_distance_nm']:.3f} nm\n"
                    f"\n=== SURVEY TIME BREAKDOWN ===\n"
                    f"Main Lines Survey: {total_survey_time['main_lines_minutes']:.1f} min\n"
                    f"Crossline Survey: {total_survey_time['crossline_minutes']:.1f} min\n"
                    f"Travel Between Lines: {total_survey_time['travel_minutes']:.1f} min\n"
                    f"Travel to Crossline: {total_survey_time['travel_to_crossline_minutes']:.1f} min\n"
                    f"\n=== TOTAL SURVEY TIME ===\n"
                    f"Total Survey Time: {total_survey_time['total_minutes']:.1f} min\n"
                    f"Total Survey Time: {total_survey_time['total_hours']:.2f} hr"
                )

                # Add crossline heading information to summary
                if self.cross_line_data and len(self.cross_line_data) >= 2:
                    try:
                        if pyproj is not None:
                            geod = pyproj.Geod(ellps="WGS84")
                            lat1, lon1 = self.cross_line_data[0]
                            lat2, lon2 = self.cross_line_data[1]
                            fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)

                            crossline_heading = fwd_az % 360
                            crossline_reciprocal_heading = back_az % 360

                            summary += f"Crossline Heading: {crossline_heading:.1f}°\n"
                            summary += f"Crossline Reciprocal Heading: {crossline_reciprocal_heading:.1f}°\n"
                        else:
                            summary += f"Crossline Heading: pyproj not available\n"
                    except Exception:
                        summary += f"Crossline Heading: Unable to calculate\n"

                self.set_ref_info_text(summary)
            except Exception as e:
                self.set_ref_info_text(f"Error calculating summary: {e}")

            # Success dialog removed as requested
            # if show_success_dialog:
            #     self._show_message("info","Success", "Survey plan generated and plotted.")

        except CRSError as e:
            error_msg = f"Failed to set up projection: {str(e)}"
            print(f"ERROR: {error_msg}")
            self._show_message("error", "Projection Error", error_msg)
        except Exception as e:
            error_msg = f"Failed to generate survey lines: {str(e)}"
            print(f"ERROR: {error_msg}")
            print("Full traceback:")
            traceback.print_exc()
            self._show_message("error", "Calculation Error", error_msg)

        # In _generate_and_plot and _load_geotiff, after plotting/plan generation, call:
        # Use _draw_current_profile to draw the appropriate profile based on active tab
        if hasattr(self, '_draw_current_profile'):
            self._draw_current_profile()
        else:
            self._draw_crossline_profile()

    def _remove_colorbar(self, colorbar_attr):
        colorbar = getattr(self, colorbar_attr, None)
        if colorbar is not None:
            try:
                if hasattr(colorbar, 'ax') and colorbar.ax is not None:
                    # Only remove if the axes is still in the figure
                    if colorbar.ax in self.figure.axes:
                        self.figure.delaxes(colorbar.ax)
            except Exception as e:
                print(f"Warning: Error removing {colorbar_attr}: {e}")
            setattr(self, colorbar_attr, None)
            # Clear stored axes position if no colorbars remain
            has_any_colorbar = (hasattr(self, 'elevation_colorbar') and self.elevation_colorbar is not None) or \
                               (hasattr(self, 'slope_colorbar') and self.slope_colorbar is not None)
            if not has_any_colorbar and hasattr(self, '_axes_pos_with_colorbar'):
                delattr(self, '_axes_pos_with_colorbar')

    def _calculate_adaptive_vert_exag(self, data_array):
        """
        Calculate adaptive vertical exaggeration based on elevation range.
        Returns a vert_exag value that adapts to the relief in the data.
        """
        valid_data = data_array[~np.isnan(data_array)]
        if len(valid_data) == 0:
            return 0.1  # Default if no valid data

        elevation_range = np.max(valid_data) - np.min(valid_data)

        if elevation_range < 20:
            vert_exag = 1.5 + (elevation_range / 20) * 1.5  # 1.5 at 0m, 3.0 at 20m
        elif elevation_range < 50:
            vert_exag = 0.8 + ((elevation_range - 20) / 30) * 0.7  # 0.8 at 20m, 1.5 at 50m
        elif elevation_range < 200:
            vert_exag = 0.4 + ((elevation_range - 50) / 150) * 0.4  # 0.4 at 50m, 0.8 at 200m
        else:
            vert_exag = max(0.05, 0.1 - ((elevation_range - 200) / 800) * 0.05)  # 0.1 down to 0.05

        return vert_exag

    def _calculate_consistent_plot_limits(self):
        """Calculate plot limits that maintain consistent window size regardless of GeoTIFF dimensions."""
        if self.geotiff_extent is None:
            # No GeoTIFF loaded, use global limits
            return self.fixed_xlim, self.fixed_ylim

        # Get GeoTIFF extent
        min_lon, max_lon, min_lat, max_lat = self.geotiff_extent

        # Calculate the center of the GeoTIFF
        center_lon = (min_lon + max_lon) / 2
        center_lat = (min_lat + max_lat) / 2

        # Calculate the dimensions of the GeoTIFF
        geotiff_width = max_lon - min_lon
        geotiff_height = max_lat - min_lat

        # Calculate geographic aspect ratio at center latitude
        # Longitude degrees get shorter as you move away from equator
        center_lat_rad = np.radians(center_lat)
        geographic_aspect = 1.0 / np.cos(center_lat_rad)
        geographic_aspect = np.clip(geographic_aspect, 0.1, 10.0)

        # Get figure dimensions to calculate figure aspect ratio
        fig_width, fig_height = self.figure.get_size_inches()
        plot_width = fig_width * (0.95 - 0.08)  # right - left margins
        plot_height = fig_height * (0.95 - 0.08)  # top - bottom margins
        figure_aspect = plot_width / plot_height

        # Calculate the GeoTIFF's display aspect ratio (accounting for geographic aspect)
        # When displayed, longitude is stretched by geographic_aspect
        geotiff_display_width = geotiff_width * geographic_aspect
        geotiff_display_height = geotiff_height
        geotiff_display_aspect = geotiff_display_width / geotiff_display_height if geotiff_display_height > 0 else 1.0

        # Use GeoTIFF dimensions directly with minimal buffer (1% to avoid edge clipping)
        buffer_factor = 1.01

        # Calculate limits that will fill the figure while maintaining GeoTIFF aspect ratio
        # We want the display aspect of our limits to match the figure aspect
        # Display aspect = (width * geographic_aspect) / height
        # We want: (width * geographic_aspect) / height = figure_aspect
        # So: width / height = figure_aspect / geographic_aspect

        # Start with GeoTIFF dimensions with minimal buffer
        base_width = geotiff_width * buffer_factor
        base_height = geotiff_height * buffer_factor

        # Calculate what the display aspect would be with these dimensions
        base_display_aspect = (base_width * geographic_aspect) / base_height if base_height > 0 else 1.0

        # Adjust to match figure aspect while maintaining GeoTIFF's aspect ratio
        # Only expand in one dimension to fill the figure, keeping the other at the GeoTIFF size
        if base_display_aspect > figure_aspect:
            # Display would be wider than figure - expand height to fill, keep width at GeoTIFF size
            new_height = (base_width * geographic_aspect) / figure_aspect
            new_width = base_width
        else:
            # Display would be taller than figure - expand width to fill, keep height at GeoTIFF size
            new_width = (base_height * figure_aspect) / geographic_aspect
            new_height = base_height

        # Set minimum size
        min_plot_size = 0.1  # degrees
        new_width = max(new_width, min_plot_size)
        new_height = max(new_height, min_plot_size)

        # Calculate new limits centered on the GeoTIFF
        half_width = new_width / 2
        half_height = new_height / 2
        new_xlim = (center_lon - half_width, center_lon + half_width)
        new_ylim = (center_lat - half_height, center_lat + half_height)

        # Ensure limits don't exceed global bounds
        new_xlim = (max(-180, new_xlim[0]), min(180, new_xlim[1]))
        new_ylim = (max(-90, new_ylim[0]), min(90, new_ylim[1]))

        return new_xlim, new_ylim

    def _plot_survey_plan(self, preserve_view_limits=True):
        try:
            # Remove all axes except the main one to prevent accumulation of colorbar axes
            for ax in self.figure.axes[:]:
                if ax is not self.ax:
                    self.figure.delaxes(ax)

            # Extra: fully clear the figure and recreate the main axes to guarantee no extra axes remain
            self.figure.clear()
            self.ax = self.figure.add_subplot(111)
            # Always ensure axes fill the figure area
            self.figure.subplots_adjust(left=0.08, right=0.95, top=0.95, bottom=0.08)
            self.slope_colorbar = None
            self.elevation_colorbar = None

            # After clearing, reset image plot references (do not try to remove them)
            self.geotiff_image_plot = None
            self.geotiff_hillshade_plot = None
            self.contour_plot = None
            self.basemap_image_plot = None
            self.slope_overlay_image_plot = None

            # Calculate consistent plot limits before plotting
            self.current_xlim, self.current_ylim = self._calculate_consistent_plot_limits()

            # Patch: preserve user zoom if requested
            prev_xlim = getattr(self, '_last_user_xlim', None)
            prev_ylim = getattr(self, '_last_user_ylim', None)
            if preserve_view_limits and prev_xlim is not None and prev_ylim is not None:
                # Store the preserved limits - they will be set after aspect ratio calculation
                # to ensure they're not overridden by the aspect ratio adjustment
                preserved_xlim = prev_xlim
                preserved_ylim = prev_ylim
                # Set initial limits for aspect ratio calculation
                self.ax.set_xlim(prev_xlim)
                self.ax.set_ylim(prev_ylim)
            else:
                self.ax.set_xlim(self.current_xlim)
                self.ax.set_ylim(self.current_ylim)
                preserved_xlim = None
                preserved_ylim = None

            # Plot imagery basemap if enabled (plot first so it's in the background)
            if hasattr(self, 'show_imagery_basemap_var') and self.show_imagery_basemap_var:
                self._load_and_plot_basemap(force_reload=True)

            # Plot GeoTIFF if loaded
            if self.geotiff_data_array is not None and self.geotiff_extent is not None:
                # Reset flag for logging vertical exaggeration (only log once per plot cycle)
                self._vert_exag_logged_this_cycle = False
                if self.hillshade_vs_slope_viz_mode == "hillshade":
                    # Use multidirectional hillshade: average of 4 azimuths
                    masked_data = np.ma.array(self.geotiff_data_array, mask=np.isnan(self.geotiff_data_array))
                    data_for_hillshade = np.nan_to_num(masked_data.filled(np.nanmin(self.geotiff_data_array)))
                    # Calculate adaptive vertical exaggeration based on elevation range
                    vert_exag = self._calculate_adaptive_vert_exag(self.geotiff_data_array)
                    # Log vertical exaggeration to activity log (once per plot cycle)
                    if hasattr(self, 'param_notebook'):
                        current_tab = self.param_notebook.currentIndex()
                        if current_tab == 0 and hasattr(self, 'set_cal_info_text'):
                            self.set_cal_info_text(f"Hillshade vertical exaggeration: {vert_exag:.3f}")
                        elif current_tab == 1 and hasattr(self, 'set_ref_info_text'):
                            self.set_ref_info_text(f"Hillshade vertical exaggeration: {vert_exag:.3f}")
                        elif current_tab == 2 and hasattr(self, 'set_line_info_text'):
                            self.set_line_info_text(f"Hillshade vertical exaggeration: {vert_exag:.3f}")
                    self._vert_exag_logged_this_cycle = True
                    azimuths = [45, 135, 225, 315]
                    altitude = 45
                    hillshades = []
                    for az in azimuths:
                        ls = LightSource(azdeg=az, altdeg=altitude)
                        hillshade = ls.hillshade(data_for_hillshade, vert_exag=vert_exag)
                        hillshades.append(hillshade)
                    # Average the hillshades
                    hillshade_array = np.mean(hillshades, axis=0)
                    # Mask NaN values to make them transparent
                    masked_hillshade = np.ma.masked_where(np.isnan(self.geotiff_data_array), hillshade_array)
                    self.geotiff_hillshade_plot = self.ax.imshow(masked_hillshade, extent=tuple(self.geotiff_extent),
                                                                cmap='gray', origin='upper',
                                                                zorder=-2)  # Plot beneath everything

                    # Restored: Conditional plotting for elevation or slope overlay
                    if self.geotiff_display_mode == "elevation":
                        display_data = self.geotiff_data_array
                        cmap = 'rainbow'  # Change to rainbow

                        # Calculate min/max for the displayed elevation data (ignoring NaNs) in the current plot window
                        xlim = self.ax.get_xlim()
                        ylim = self.ax.get_ylim()
                        left, right, bottom, top = tuple(self.geotiff_extent)
                        nrows, ncols = display_data.shape
                        # Find the pixel indices that correspond to the current axes limits
                        col_min = int(np.clip((min(xlim) - left) / (right - left) * (ncols - 1), 0, ncols - 1))
                        col_max = int(np.clip((max(xlim) - left) / (right - left) * (ncols - 1), 0, ncols - 1))
                        row_min = int(np.clip((top - max(ylim)) / (top - bottom) * (nrows - 1), 0, nrows - 1))
                        row_max = int(np.clip((top - min(ylim)) / (top - bottom) * (nrows - 1), 0, nrows - 1))
                        # Ensure min <= max
                        r0, r1 = sorted([row_min, row_max])
                        c0, c1 = sorted([col_min, col_max])
                        visible_region = display_data[r0:r1+1, c0:c1+1]
                        if visible_region.size > 0 and not np.all(np.isnan(visible_region)):
                            vmin_elev = np.nanmin(visible_region)
                            vmax_elev = np.nanmax(visible_region)
                        else:
                            vmin_elev = np.nanmin(display_data) if not np.all(np.isnan(display_data)) else None
                            vmax_elev = np.nanmax(display_data) if not np.all(np.isnan(display_data)) else None

                        # Only set vmin/vmax if valid range exists
                        if vmin_elev is not None and vmax_elev is not None and vmin_elev != vmax_elev:
                            self.geotiff_image_plot = self.ax.imshow(display_data, extent=tuple(self.geotiff_extent),
                                                                     cmap=cmap, origin='upper', alpha=0.5, zorder=-1,
                                                                     vmin=vmin_elev,
                                                                     vmax=vmax_elev)
                        else:
                            self.geotiff_image_plot = self.ax.imshow(display_data, extent=tuple(self.geotiff_extent),
                                                                     cmap=cmap, origin='upper', alpha=0.5, zorder=-1)

                        # Add colorbar for elevation
                        # Let matplotlib automatically adjust the axes position (like original Tkinter version)
                        self.elevation_colorbar = self.figure.colorbar(self.geotiff_image_plot, ax=self.ax,
                                                                        orientation='vertical', label='Elevation (m)',
                                                                        pad=0.02, fraction=0.1, shrink=1.0, aspect=20)

                    elif self.geotiff_display_mode == "hillshade_only":
                        # Hillshade only mode - no overlay, just the hillshade
                        # The hillshade is already plotted above, so we don't add any overlay
                        # Remove any existing colorbars
                        if hasattr(self, 'elevation_colorbar') and self.elevation_colorbar is not None:
                            self.elevation_colorbar.remove()
                            self.elevation_colorbar = None
                        if hasattr(self, 'slope_colorbar') and self.slope_colorbar is not None:
                            self.slope_colorbar.remove()
                            self.slope_colorbar = None

                    elif self.geotiff_display_mode == "slope":
                        # Calculate slope in degrees for overlay
                        center_lat_geotiff = (self.geotiff_extent[2] + self.geotiff_extent[3]) / 2
                        m_per_deg_lat = 111320.0
                        m_per_deg_lon = 111320.0 * np.cos(np.radians(center_lat_geotiff))

                        res_lat_deg = (self.geotiff_extent[3] - self.geotiff_extent[2]) / self.geotiff_data_array.shape[0]
                        res_lon_deg = (self.geotiff_extent[1] - self.geotiff_extent[0]) / self.geotiff_data_array.shape[1]

                        dx_m = res_lon_deg * m_per_deg_lon
                        dy_m = res_lat_deg * m_per_deg_lat

                        temp_data = np.nan_to_num(self.geotiff_data_array, nan=0.0)
                        dz_dy_grid, dz_dx_grid = np.gradient(temp_data, dy_m, dx_m)

                        slope_rad = np.arctan(np.sqrt(dz_dx_grid ** 2 + dz_dy_grid ** 2))
                        slope_degrees = np.degrees(slope_rad)
                        slope_degrees[np.isnan(self.geotiff_data_array)] = np.nan

                        max_slope_for_cmap = 30.0
                        display_data = np.clip(slope_degrees, 0, max_slope_for_cmap)

                        # Add hillshade underlay from a single direction (315, 45)
                        vert_exag = self._calculate_adaptive_vert_exag(self.geotiff_data_array)
                        if not getattr(self, '_vert_exag_logged_this_cycle', True):
                            if hasattr(self, 'param_notebook'):
                                current_tab = self.param_notebook.currentIndex()
                                if current_tab == 0 and hasattr(self, 'set_cal_info_text'):
                                    self.set_cal_info_text(f"Hillshade vertical exaggeration: {vert_exag:.3f}")
                                elif current_tab == 1 and hasattr(self, 'set_ref_info_text'):
                                    self.set_ref_info_text(f"Hillshade vertical exaggeration: {vert_exag:.3f}")
                                elif current_tab == 2 and hasattr(self, 'set_line_info_text'):
                                    self.set_line_info_text(f"Hillshade vertical exaggeration: {vert_exag:.3f}")
                            self._vert_exag_logged_this_cycle = True
                        ls = LightSource(azdeg=315, altdeg=45)
                        hillshade = ls.hillshade(temp_data, vert_exag=vert_exag)
                        masked_hillshade = np.ma.masked_where(np.isnan(self.geotiff_data_array), hillshade)
                        self.ax.imshow(masked_hillshade, extent=tuple(self.geotiff_extent), cmap='gray', origin='upper', alpha=0.5, zorder=-2)

                        cmap = 'plasma'
                        self.geotiff_image_plot = self.ax.imshow(display_data, extent=tuple(self.geotiff_extent),
                                                                 cmap=cmap, origin='upper', alpha=0.5,
                                                                 vmin=0, vmax=max_slope_for_cmap,
                                                                 zorder=-1)

                        self.slope_colorbar = self.figure.colorbar(self.geotiff_image_plot, ax=self.ax,
                                                                   orientation='vertical', label='Slope (degrees)',
                                                                   pad=0.02, fraction=0.1, shrink=1.0, aspect=20)

                        num_ticks = 7
                        tick_values = np.linspace(0, max_slope_for_cmap, num_ticks)
                        self.slope_colorbar.set_ticks(tick_values.tolist())

                        tick_labels = [f"{t:.0f}" for t in tick_values]
                        if np.nanmax(slope_degrees) > max_slope_for_cmap:
                            tick_labels[-1] = f"> {max_slope_for_cmap:.0f}"
                        self.slope_colorbar.set_ticklabels(tick_labels)

                else:  # slope_viz mode
                    # Calculate slope visualization (similar to above but without overlay)
                    center_lat_geotiff = (self.geotiff_extent[2] + self.geotiff_extent[3]) / 2
                    m_per_deg_lat = 111320.0
                    m_per_deg_lon = 111320.0 * np.cos(np.radians(center_lat_geotiff))

                    res_lat_deg = (self.geotiff_extent[3] - self.geotiff_extent[2]) / self.geotiff_data_array.shape[0]
                    res_lon_deg = (self.geotiff_extent[1] - self.geotiff_extent[0]) / self.geotiff_data_array.shape[1]

                    dx_m = res_lon_deg * m_per_deg_lon
                    dy_m = res_lat_deg * m_per_deg_lat

                    temp_data = np.nan_to_num(self.geotiff_data_array, nan=0.0)
                    dz_dy_grid, dz_dx_grid = np.gradient(temp_data, dy_m, dx_m)

                    slope_rad = np.arctan(np.sqrt(dz_dx_grid ** 2 + dz_dy_grid ** 2))
                    slope_degrees = np.degrees(slope_rad)
                    slope_degrees[np.isnan(self.geotiff_data_array)] = np.nan

                    max_slope_for_cmap = 30.0
                    display_slope_data = np.clip(slope_degrees, 0, max_slope_for_cmap)

                    self.geotiff_image_plot = self.ax.imshow(display_slope_data, extent=tuple(self.geotiff_extent),
                                                             cmap='inferno', origin='upper', vmin=0,
                                                             vmax=max_slope_for_cmap, zorder=-1)

                    self.slope_colorbar = self.figure.colorbar(self.geotiff_image_plot, ax=self.ax,
                                                              orientation='vertical', label='Slope (degrees)',
                                                              pad=0.02, fraction=0.1, shrink=1.0, aspect=20)

                    num_ticks = 7
                    tick_values = np.linspace(0, max_slope_for_cmap, num_ticks)
                    self.slope_colorbar.set_ticks(tick_values.tolist())

                    tick_labels = [f"{t:.0f}" for t in tick_values]
                    if np.nanmax(slope_degrees) > max_slope_for_cmap:
                        tick_labels[-1] = f"> {max_slope_for_cmap:.0f}"
                    self.slope_colorbar.set_ticklabels(tick_labels)

            # Plot contours if enabled
            if (self.geotiff_data_array is not None and self.geotiff_extent is not None and
                hasattr(self, 'show_contours_var') and self.show_contours_var):
                try:
                    # Get contour interval from entry field (check both tabs)
                    contour_interval = 200.0  # Default
                    try:
                        if hasattr(self, 'contour_interval_entry') and self.contour_interval_entry:
                            contour_interval = float(self.contour_interval_entry.text())
                        if contour_interval <= 0:
                            contour_interval = 200.0  # Default to 200 meters if invalid
                    except (ValueError, AttributeError):
                        contour_interval = 200.0  # Default to 200 meters if invalid

                    # Remove previous contour plot if it exists
                    if hasattr(self, 'contour_plot') and self.contour_plot is not None:
                        for collection in self.contour_plot.collections:
                            collection.remove()
                        self.contour_plot = None

                    # Get bathymetry data (always use elevation data, not slope)
                    bathymetry_data = self.geotiff_data_array

                    # Create coordinate grids for contour plotting
                    left, right, bottom, top = tuple(self.geotiff_extent)
                    nrows, ncols = bathymetry_data.shape
                    lon_grid = np.linspace(left, right, ncols)
                    lat_grid = np.linspace(bottom, top, nrows)
                    lon_mesh, lat_mesh = np.meshgrid(lon_grid, lat_grid)

                    # Flip data vertically because imshow uses origin='upper' (row 0 = top)
                    # but contour expects origin='lower' (row 0 = bottom)
                    bathymetry_data_flipped = np.flipud(bathymetry_data)

                    # Calculate contour levels based on data range and interval
                    valid_data = bathymetry_data[~np.isnan(bathymetry_data)]
                    if len(valid_data) > 0:
                        min_elev = np.nanmin(valid_data)
                        max_elev = np.nanmax(valid_data)

                        # Generate contour levels
                        # Start from the first level above min_elev that's a multiple of interval
                        start_level = np.ceil(min_elev / contour_interval) * contour_interval
                        end_level = np.floor(max_elev / contour_interval) * contour_interval
                        contour_levels = np.arange(start_level, end_level + contour_interval, contour_interval)

                        if len(contour_levels) > 0:
                            # Plot contours using matplotlib contour
                            # Use cyan color when in shaded slope, hillshade, or slope mode, black for shaded relief
                            if hasattr(self, 'geotiff_display_mode'):
                                contour_color = 'cyan' if self.geotiff_display_mode in ["slope", "hillshade_only", "slope_viz"] else 'black'
                            else:
                                contour_color = 'black'
                            self.contour_plot = self.ax.contour(lon_mesh, lat_mesh, bathymetry_data_flipped,
                                                                levels=contour_levels,
                                                                colors=contour_color,
                                                                linewidths=0.5,
                                                                alpha=0.7,
                                                                zorder=12)  # Above slope overlay and other features
                except Exception as e:
                    # Silently fail if contour plotting fails (e.g., invalid data, etc.)
                    if hasattr(self, 'contour_plot') and self.contour_plot is not None:
                        try:
                            for collection in self.contour_plot.collections:
                                collection.remove()
                        except Exception:
                            pass
                        self.contour_plot = None
            else:
                # Remove contours if checkbox is disabled or no GeoTIFF loaded
                if hasattr(self, 'contour_plot') and self.contour_plot is not None:
                    try:
                        for collection in self.contour_plot.collections:
                            collection.remove()
                    except Exception:
                        pass
                    self.contour_plot = None

            # Plot slope overlay if enabled
            if (self.geotiff_data_array is not None and self.geotiff_extent is not None and
                hasattr(self, 'show_slope_overlay_var') and self.show_slope_overlay_var):
                try:
                    # Get min/max values from entry fields
                    min_slope = 10.0  # Default
                    max_slope = 20.0  # Default
                    try:
                        if hasattr(self, 'slope_overlay_min_var'):
                            min_slope = float(self.slope_overlay_min_var)
                        if hasattr(self, 'slope_overlay_max_var'):
                            max_slope = float(self.slope_overlay_max_var)
                        if min_slope >= max_slope:
                            min_slope = 10.0
                            max_slope = 20.0
                    except (ValueError, AttributeError):
                        min_slope = 10.0
                        max_slope = 20.0

                    # Remove previous slope overlay plot if it exists
                    if hasattr(self, 'slope_overlay_image_plot') and self.slope_overlay_image_plot is not None:
                        try:
                            self.slope_overlay_image_plot.remove()
                        except Exception:
                            pass
                        self.slope_overlay_image_plot = None

                    # Calculate slope in degrees (same method as used for slope display mode)
                    center_lat_geotiff = (self.geotiff_extent[2] + self.geotiff_extent[3]) / 2
                    m_per_deg_lat = 111320.0
                    m_per_deg_lon = 111320.0 * np.cos(np.radians(center_lat_geotiff))

                    res_lat_deg = (self.geotiff_extent[3] - self.geotiff_extent[2]) / self.geotiff_data_array.shape[0]
                    res_lon_deg = (self.geotiff_extent[1] - self.geotiff_extent[0]) / self.geotiff_data_array.shape[1]

                    dx_m = res_lon_deg * m_per_deg_lon
                    dy_m = res_lat_deg * m_per_deg_lat

                    temp_data = np.nan_to_num(self.geotiff_data_array, nan=0.0)
                    dz_dy_grid, dz_dx_grid = np.gradient(temp_data, dy_m, dx_m)

                    slope_rad = np.arctan(np.sqrt(dz_dx_grid ** 2 + dz_dy_grid ** 2))
                    slope_degrees = np.degrees(slope_rad)
                    slope_degrees[np.isnan(self.geotiff_data_array)] = np.nan

                    # Create mask for slopes between min and max
                    mask = (slope_degrees >= min_slope) & (slope_degrees <= max_slope)
                    # Create green overlay with transparency
                    overlay = np.zeros((slope_degrees.shape[0], slope_degrees.shape[1], 4), dtype=np.float32)
                    overlay[mask, 0] = 0.0  # Red channel
                    overlay[mask, 1] = 1.0  # Green channel
                    overlay[mask, 2] = 0.0  # Blue channel
                    overlay[mask, 3] = 0.4  # Alpha channel (40% transparency)
                    overlay[~mask, 3] = 0.0  # Transparent outside range

                    # Plot the overlay
                    self.slope_overlay_image_plot = self.ax.imshow(
                        overlay,
                        extent=tuple(self.geotiff_extent),
                        origin='upper',
                        zorder=11,  # Above other layers but below contours
                        interpolation='bilinear'
                    )
                except Exception as e:
                    # Silently fail if slope overlay plotting fails
                    if hasattr(self, 'slope_overlay_image_plot') and self.slope_overlay_image_plot is not None:
                        try:
                            self.slope_overlay_image_plot.remove()
                        except Exception:
                            pass
                        self.slope_overlay_image_plot = None
            else:
                # Remove slope overlay if checkbox is disabled or no GeoTIFF loaded
                if hasattr(self, 'slope_overlay_image_plot') and self.slope_overlay_image_plot is not None:
                    try:
                        self.slope_overlay_image_plot.remove()
                    except Exception:
                        pass
                    self.slope_overlay_image_plot = None

            # Plot main survey lines
            for i, line in enumerate(self.survey_lines_data):
                latitudes = [p[0] for p in line]
                longitudes = [p[1] for p in line]
                label = "Reference Line" if i == 0 else "_nolegend_"
                self.ax.plot(longitudes, latitudes, color='blue', linewidth=1.5, linestyle='--', label=label)

                # For alternating lines, flip the start/end points to show survey order
                # Even lines (0, 2, 4...) use normal order, odd lines (1, 3, 5...) use flipped order
                if i % 2 == 0:  # Even lines - normal order
                    start_lon, start_lat = longitudes[0], latitudes[0]
                    end_lon, end_lat = longitudes[1], latitudes[1]
                    start_label = f'L{i+1}S'
                    end_label = f'L{i+1}E'
                else:  # Odd lines - flipped order to show zigzag survey pattern
                    start_lon, start_lat = longitudes[1], latitudes[1]
                    end_lon, end_lat = longitudes[0], latitudes[0]
                    start_label = f'L{i+1}S'
                    end_label = f'L{i+1}E'

                # Add labels for start and end points
                self.ax.annotate(start_label, (start_lon, start_lat),
                                xytext=(5, 5), textcoords='offset points',
                                fontsize=8, color='blue', weight='bold',
                                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
                self.ax.annotate(end_label, (end_lon, end_lat),
                                xytext=(5, 5), textcoords='offset points',
                                fontsize=8, color='blue', weight='bold',
                                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))

            # Plot crossline
            if self.cross_line_data:
                latitudes = [p[0] for p in self.cross_line_data]
                longitudes = [p[1] for p in self.cross_line_data]
                self.ax.plot(longitudes, latitudes, color='darkorchid', linestyle='-', linewidth=1.5,
                             label='Crossline')

                # Add labels for crossline points
                self.ax.annotate('CLS', (longitudes[0], latitudes[0]),
                                xytext=(5, 5), textcoords='offset points',
                                fontsize=8, color='darkorchid', weight='bold',
                                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
                self.ax.annotate('CLE', (longitudes[1], latitudes[1]),
                                xytext=(5, 5), textcoords='offset points',
                                fontsize=8, color='darkorchid', weight='bold',
                                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))

            # Plot central point
            if self.central_point_coords[0] is not None:
                self.ax.plot(self.central_point_coords[1], self.central_point_coords[0],
                             'o', color='green', markersize=5, label='Central Point')
                self.ax.annotate('CP', (self.central_point_coords[1], self.central_point_coords[0]),
                                xytext=(5, 5), textcoords='offset points',
                                fontsize=8, color='green', weight='bold',
                                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))

            # Plot pitch line if defined
            if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2:
                latitudes = [p[0] for p in self.pitch_line_points]
                longitudes = [p[1] for p in self.pitch_line_points]
                self.ax.plot(longitudes, latitudes, color='red', linewidth=2.5, linestyle='-', marker='o', label='Pitch Line')

                # Add labels for pitch line points
                self.ax.annotate('PLS', (longitudes[0], latitudes[0]),
                                xytext=(5, 5), textcoords='offset points',
                                fontsize=8, color='red', weight='bold',
                                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
                self.ax.annotate('PLE', (longitudes[1], latitudes[1]),
                                xytext=(5, 5), textcoords='offset points',
                                fontsize=8, color='red', weight='bold',
                                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))

            # Plot roll line if defined
            if hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2:
                latitudes = [p[0] for p in self.roll_line_points]
                longitudes = [p[1] for p in self.roll_line_points]
                self.ax.plot(longitudes, latitudes, color='purple', linewidth=2.5, linestyle='-', marker='o', label='Roll')

                # Add labels for roll line points
                self.ax.annotate('RLS', (longitudes[0], latitudes[0]),
                                xytext=(5, 5), textcoords='offset points',
                                fontsize=8, color='purple', weight='bold',
                                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
                self.ax.annotate('RLE', (longitudes[1], latitudes[1]),
                                xytext=(5, 5), textcoords='offset points',
                                fontsize=8, color='purple', weight='bold',
                                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))

            # Plot heading lines if present
            if hasattr(self, 'heading_lines') and len(self.heading_lines) == 2:
                for i, line in enumerate(self.heading_lines):
                    latitudes = [p[0] for p in line]
                    longitudes = [p[1] for p in line]
                    label = 'Heading1' if i == 0 else 'Heading2'
                    self.ax.plot(longitudes, latitudes, color='deeppink', linewidth=2, linestyle='--', marker='x', label=label)

                    # Add labels for heading line points
                    prefix = 'H1' if i == 0 else 'H2'
                    self.ax.annotate(f'{prefix}S', (longitudes[0], latitudes[0]),
                                    xytext=(5, 5), textcoords='offset points',
                                    fontsize=8, color='deeppink', weight='bold',
                                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
                    self.ax.annotate(f'{prefix}E', (longitudes[1], latitudes[1]),
                                    xytext=(5, 5), textcoords='offset points',
                                    fontsize=8, color='deeppink', weight='bold',
                                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))

            # --- Plot the drawn line in Line Planning tab if it exists ---
            if hasattr(self, 'line_planning_points') and len(self.line_planning_points) >= 2:
                lats = [p[0] for p in self.line_planning_points]
                lons = [p[1] for p in self.line_planning_points]
                self.ax.plot(lons, lats, color='orange', linewidth=2, marker='o', label='Line Planning')

                # Add labels for line planning waypoints
                for i, (lon, lat) in enumerate(zip(lons, lats)):
                    self.ax.annotate(f'WP{i+1}', (lon, lat),
                                    xytext=(5, 5), textcoords='offset points',
                                    fontsize=8, color='orange', weight='bold',
                                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))

            # Standard plot elements
            self.ax.set_xlabel("Longitude")
            self.ax.set_ylabel("Latitude")
            self.ax.set_title("Survey Plan")
            self.ax.grid(True)

            # Format axis tick labels as DDM (degrees-decimal minutes)
            self.ax.xaxis.set_major_formatter(
                FuncFormatter(lambda x, _: decimal_degrees_to_ddm(x, is_latitude=False))
            )
            self.ax.yaxis.set_major_formatter(
                FuncFormatter(lambda y, _: decimal_degrees_to_ddm(y, is_latitude=True))
            )

            # Calculate aspect ratio based on center latitude for equal scaling
            # Longitude degrees get shorter as you move away from equator
            # Aspect ratio = 1 / cos(latitude) to maintain equal distance scaling
            try:
                # Get initial limits for center calculation
                xlim = self.ax.get_xlim()
                ylim = self.ax.get_ylim()
                center_lat = (ylim[0] + ylim[1]) / 2.0
                center_x = (xlim[0] + xlim[1]) / 2.0
                center_y = (ylim[0] + ylim[1]) / 2.0

                # Convert to radians and calculate aspect ratio
                center_lat_rad = np.radians(center_lat)
                aspect_ratio = 1.0 / np.cos(center_lat_rad)
                # Clamp aspect ratio to reasonable bounds to avoid extreme values at poles
                aspect_ratio = np.clip(aspect_ratio, 0.1, 10.0)

                # Get figure size
                fig_width, fig_height = self.figure.get_size_inches()
                # Use subplot_adjust margins: left=0.08, right=0.95, top=0.95, bottom=0.08
                # This gives us the actual plotable area
                plot_width = fig_width * (0.95 - 0.08)  # right - left
                plot_height = fig_height * (0.95 - 0.08)  # top - bottom
                figure_aspect = plot_width / plot_height

                # Calculate current data extent
                data_width = xlim[1] - xlim[0]
                data_height = ylim[1] - ylim[0]

                # Calculate how the data would display with the geographic aspect ratio
                # With aspect_ratio applied, longitude is stretched, so display aspect is:
                data_display_aspect = (data_width * aspect_ratio) / data_height

                # Only adjust limits to fill the figure if we're not preserving user's zoom level
                # When preserve_view_limits=True, skip this adjustment to maintain user's zoom
                if not preserve_view_limits:
                    # Adjust limits to fill the figure while maintaining geographic aspect
                    if data_display_aspect > figure_aspect:
                        # Data would appear wider than figure - expand height (latitude) to fill
                        new_height = (data_width * aspect_ratio) / figure_aspect
                        ylim_new = (center_y - new_height / 2.0, center_y + new_height / 2.0)
                        self.ax.set_ylim(ylim_new)
                        # Update xlim to match the new center in case it shifted
                        self.ax.set_xlim(xlim)
                    elif data_display_aspect < figure_aspect:
                        # Data would appear taller than figure - expand width (longitude) to fill
                        new_width = (data_height * figure_aspect) / aspect_ratio
                        xlim_new = (center_x - new_width / 2.0, center_x + new_width / 2.0)
                        self.ax.set_xlim(xlim_new)
                        # Update ylim to match the new center in case it shifted
                        self.ax.set_ylim(ylim)
                    # If equal, no adjustment needed

                # Now set the aspect ratio with 'datalim' so the axes box stays fixed
                # This matches the original Tkinter version - let matplotlib handle positioning automatically
                self.ax.set_aspect(aspect_ratio, adjustable='datalim')

                # If preserving view limits, restore them after setting aspect
                if preserve_view_limits and 'preserved_xlim' in locals() and preserved_xlim is not None:
                    # Restore the exact limits
                    self.ax.set_xlim(preserved_xlim)
                    self.ax.set_ylim(preserved_ylim)

            except Exception:
                # Fallback to auto if calculation fails
                self.ax.set_aspect('auto')

            # Remove any existing legend before adding a new one to prevent shrinking
            handles, labels = self.ax.get_legend_handles_labels()
            if handles and any(label and not label.startswith('_') for label in labels):
                self.ax.legend()
            else:
                if self.ax.get_legend() is not None:
                    self.ax.get_legend().remove()

            # Plot NOAA ENC Charts overlay if enabled (plot last so it overlays everything)
            if hasattr(self, 'show_noaa_charts_var') and self.show_noaa_charts_var:
                self._load_and_plot_noaa_charts(force_reload=True)

            self.canvas.draw_idle()

        except Exception as e:
            print(f"Error in _plot_survey_plan: {str(e)}")
            traceback.print_exc()
            self._show_message("error", "Plot Error", f"Failed to generate plot: {str(e)}")
            # Attempt to recover by reinitializing the plot
            self._clear_plot(full_clear=True)
            self.ax = self.figure.add_subplot(111)
            self.figure.subplots_adjust(left=0.08, right=0.95, top=0.95, bottom=0.08)
            self.canvas.draw_idle()

    def _clear_plot(self, full_clear=True):
        """Clear the plot and optionally reset all data."""
        # Clear the axes
        if self.ax:
            self.ax.clear()

        # Reset to fixed plot limits
        self.current_xlim = self.fixed_xlim
        self.current_ylim = self.fixed_ylim
        self.ax.set_xlim(self.current_xlim)
        self.ax.set_ylim(self.current_ylim)

        # Reset data
        self.survey_lines_data = []
        self.cross_line_data = []
        self.central_point_coords = (None, None)

        # Clear GeoTIFF related data only if full_clear is True
        if full_clear:
            if self.geotiff_dataset_original:
                self.geotiff_dataset_original.close()
                self.geotiff_dataset_original = None
            self.geotiff_data_array = None
            self.geotiff_extent = None
            self.geotiff_original_extent = None
            self.initial_geotiff_xlim = None
            self.initial_geotiff_ylim = None
            self.geotiff_image_plot = None
            self.geotiff_hillshade_plot = None
            self.slope_overlay_image_plot = None
            self.contour_plot = None
            # Robustly remove colorbars
            self._remove_colorbar('slope_colorbar')
            self._remove_colorbar('elevation_colorbar')

            # Restored: geotiff_display_mode reset
            self.geotiff_display_mode = "elevation"
            # Update combo box if it exists
            if hasattr(self, 'elevation_slope_combo'):
                self._update_display_mode_combo()

            self.hillshade_vs_slope_viz_mode = "hillshade"

            # Reset transformers
            self.geotiff_to_wgs84_transformer = None
            self.wgs84_to_geotiff_transformer = None

            # Disable Remove GeoTIFF button and Pick Center button
            if hasattr(self, 'remove_geotiff_btn'):
                self.remove_geotiff_btn.setEnabled(False)
            if hasattr(self, 'pick_center_btn'):
                self.pick_center_btn.setEnabled(False)
                self.pick_center_btn.setStyleSheet("")  # Reset to normal style when GeoTIFF is removed
            # Reset "Start Drawing Line" button to normal style when GeoTIFF is removed
            if hasattr(self, 'line_start_draw_btn'):
                self.line_start_draw_btn.setStyleSheet("")
            # Reset Load GeoTIFF button to orange and bold when GeoTIFF is removed
            if hasattr(self, 'load_geotiff_btn'):
                self.load_geotiff_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")

        # Reset plot elements
        if self.ax:
            self.ax.set_xlabel("")
            self.ax.set_ylabel("")
            self.ax.set_title("")
            if self.ax.legend_ is not None:  # Check if legend exists before trying to remove
                self.ax.legend().remove()

        # Redraw the canvas
        if hasattr(self, 'canvas'):
            self.canvas.draw_idle()

        # Also clear the crossline elevation profile
        self._draw_crossline_profile()

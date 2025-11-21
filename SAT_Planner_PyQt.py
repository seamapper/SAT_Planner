from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
                             QGridLayout, QLabel, QLineEdit, QPushButton, QCheckBox, QTextEdit,
                             QScrollArea, QTabWidget, QFileDialog, QMessageBox, QDialog,
                             QDialogButtonBox, QSlider, QComboBox, QFrame, QSizePolicy, QProgressBar, QGroupBox)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt6.QtGui import QTextCursor, QColor, QTextCharFormat
# Ensure matplotlib is imported (for PyInstaller)
import matplotlib
matplotlib.use('QtAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.colors import LightSource
from matplotlib.figure import Figure
import numpy as np
import os
import csv
import traceback  # Added for more detailed error logging
import json
import platform
import datetime
import re
import threading
import time

# Import pyproj for coordinate transformations and heading calculations
try:
    import pyproj
except ImportError:
    pyproj = None

# __version__ = "2025.01"  # 1st version of the app
#__version__ = "2025.02"  # Added metadata file saving/loading, contour interval synchronization
# __version__ = "2025.03"  # Added line planning, import/export of lines, and other improvements
__version__ = "2025.04"  # Converted to PyQt6

# --- Conditional Imports for Geospatial Libraries ---
GEOSPATIAL_LIBS_AVAILABLE = True  # Assume true until an import fails
try:
    import rasterio
    from rasterio import transform
    from rasterio.errors import RasterioIOError
    from rasterio.windows import Window
    from rasterio.windows import transform as window_transform
    from rasterio.windows import bounds as window_bounds
    from rasterio.transform import rowcol
    from rasterio.warp import reproject, Resampling
    import pyproj
    from pyproj.exceptions import CRSError
    from shapely.geometry import LineString
    import fiona
except ImportError as e:
    GEOSPATIAL_LIBS_AVAILABLE = False
    print(f"Warning: A geospatial library not found: {e}. GeoTIFF/Shapefile features will be disabled.")
    # Show warning only once on app startup if libs are missing
    # This messagebox will appear after the Tkinter root is initialized
    # so we'll save the message for later if this is a direct script execution.
    # If run via __name__ == "__main__", it will be shown when app starts.


class SurveyPlanApp(QMainWindow):
    CONFIG_FILENAME = os.path.join(os.path.expanduser("~"), ".cal_ref_planner_config.json")

    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"UNH/CCOM-JHC - SAT/QAT Planner - v{__version__} - pjohnson@ccom.unh.edu")
        
        # Set minimum window size
        self.setMinimumSize(1600, 1150)
        
        # Auto-regenerate timer: wait 800ms after last parameter change before regenerating
        self.auto_regenerate_timer = QTimer()
        self.auto_regenerate_timer.setSingleShot(True)
        self.auto_regenerate_timer.timeout.connect(self._auto_regenerate_survey_plan)
        
        # Central widget for main layout
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.main_layout = QVBoxLayout(self.central_widget)  # Changed to vertical for profile at bottom
        self.main_layout.setContentsMargins(5, 5, 5, 5)  # Add some margins
        
        # Create horizontal layout for params and plot
        self.content_layout = QHBoxLayout()
        self.content_layout.setContentsMargins(0, 0, 0, 0)
        
        # Initialize figure and axes with fixed size
        self.figure = Figure(figsize=(12, 10))  # Increased size for better visibility
        self.ax = self.figure.add_subplot(111)
        self.figure.subplots_adjust(left=0.08, right=0.95, top=0.95, bottom=0.08)
        
        # Set fixed plot window limits (global extent that will be maintained)
        # These will be the default limits when no GeoTIFF is loaded
        self.fixed_xlim = (-180, 180)  # Global longitude range
        self.fixed_ylim = (-90, 90)    # Global latitude range
        
        # Store the current plot limits (will be updated when GeoTIFF is loaded)
        self.current_xlim = self.fixed_xlim
        self.current_ylim = self.fixed_ylim

        self.survey_lines_data = []  # Stores [(lat1, lon1), (lat2, lon2)] for each line
        self.cross_line_data = []  # Stores [(lat1, lon1), (lat2, lon2)] for the crossline
        self.central_point_coords = (None, None)  # Stores (lat, lon) of central point

        self.geotiff_dataset_original = None  # Original rasterio dataset (keeps its CRS)
        self.geotiff_data_array = None  # NumPy array of *reprojected* GeoTIFF elevation data (WGS84)
        self.geotiff_extent = None  # [left, right, bottom, top] of *reprojected* data
        self.geotiff_image_plot = None  # Matplotlib imshow object for GeoTIFF
        self.geotiff_hillshade_plot = None  # Matplotlib imshow object for hillshade
        self.slope_colorbar = None  # To store the colorbar object for slope visualization
        self.elevation_colorbar = None  # To store the colorbar object for elevation visualization
        
        # Dynamic resolution variables
        self.geotiff_full_resolution = None  # Store full resolution data for dynamic loading
        self.geotiff_current_resolution = None  # Current resolution being displayed
        self.geotiff_zoom_level = 1.0  # Current zoom level (1.0 = full resolution)
        self.geotiff_zoom_threshold = 0.1  # Zoom threshold for switching to higher resolution
        self.dynamic_resolution_enabled = True  # Toggle for dynamic resolution feature

        # Restored: Variable for controlling elevation/slope overlay
        self.geotiff_display_mode = "elevation"  # "elevation" or "slope"

        self.hillshade_vs_slope_viz_mode = "hillshade"  # "hillshade" or "slope_viz" (for main underlay)

        # Contour display variables
        self.show_contours_var = False  # Default to off
        self.contour_plot = None  # Store contour plot object for removal/update

        # Variable for the new line length multiplier slider
        self.line_length_multiplier = 8.0  # Default multiplier
        # Variable for the new distance between lines multiplier slider
        self.dist_between_lines_multiplier = 1.0  # Default multiplier

        # Initialize the picking mode state
        self.pick_center_mode = False

        # Pyproj transformers
        self.geod = None  # Geodetic object for precise distance/bearing on ellipsoid
        self.local_proj_transformer = None  # For transforming lat/lon to local projected CRS (easting, northing)
        self.wgs84_transformer = None  # For transforming local projected CRS back to WGS84
        self.geotiff_to_wgs84_transformer = None  # For transforming GeoTIFF CRS to WGS84 for display
        self.wgs84_to_geotiff_transformer = None  # For transforming WGS84 to GeoTIFF's original CRS for picking

        self.pitch_line_points = []  # Stores [(lat1, lon1), (lat2, lon2)] for pitch line
        self.pick_pitch_line_mode = False  # Flag for pitch line picking mode
        self.edit_pitch_line_mode = False  # Flag for pitch line editing mode
        self.pitch_line_start_handle = None  # Draggable handle for start point
        self.pitch_line_end_handle = None  # Draggable handle for end point
        self.dragging_pitch_line_handle = None  # Which handle is being dragged
        self.pitch_line_temp_line = None  # Temporary line during drawing
        self.pitch_line_annotation = None  # Annotation for pitch line
        self.pitch_line_tooltip = None  # Tooltip widget
        self.pitch_line_tooltip_text = None  # Tooltip text
        self.pitch_line_edit_line = None  # Temporary line during editing
        self.heading_lines = []  # Stores two heading lines as [(lat1, lon1), (lat2, lon2)]
        self.roll_line_points = []  # Stores [(lat1, lon1), (lat2, lon2)] for roll line
        self.pick_roll_line_mode = False  # Flag for roll line picking mode
        self.edit_roll_line_mode = False  # Flag for roll line editing mode
        self.roll_line_start_handle = None  # Draggable handle for start point
        self.roll_line_end_handle = None  # Draggable handle for end point
        self.dragging_roll_line_handle = None  # Which handle is being dragged
        self.roll_line_edit_line = None  # Temporary line during editing

        # --- Add state for line planning ---
        self.line_planning_mode = False
        self.line_planning_points = []
        self.line_planning_line = None
        self.edit_line_planning_mode = False
        self.line_planning_handles = []
        self.dragging_line_planning_handle = None
        self.line_planning_edit_line = None

        # Track previous tab for clearing on tab switch
        self.previous_tab_index = 0
        self.tab_switch_initialized = False  # Flag to prevent clearing on first tab change

        # --- Patch in line planning methods and event handler BEFORE widgets are created ---
        def _toggle_line_planning_mode(self):
            if not GEOSPATIAL_LIBS_AVAILABLE or self.geotiff_data_array is None:
                self._show_message("warning", "No GeoTIFF", "Load a GeoTIFF first to draw a line.")
                return
            self.line_planning_mode = not self.line_planning_mode
            if self.line_planning_mode:
                self.line_planning_points = []
                if self.line_planning_line is not None:
                    self.line_planning_line.remove()
                    self.line_planning_line = None
                self.line_start_draw_btn.setText("Drawing: Left-click to add, Right-click to finish")
                self.canvas_widget.setCursor(Qt.CursorShape.CrossCursor)
                self.set_line_info_text("Left click to start line, left click to add waypoints, and Right click to end line")
            else:
                self.line_start_draw_btn.setText("Start Drawing Line")
                self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)

        def _clear_line_planning(self):
            self.line_planning_points = []
            if self.line_planning_line is not None:
                self.line_planning_line.remove()
                self.line_planning_line = None
            # Remove temp dashed line if present
            if hasattr(self, 'line_planning_temp_line') and self.line_planning_temp_line is not None:
                self.line_planning_temp_line.remove()
                self.line_planning_temp_line = None
            # Remove info text if present
            if hasattr(self, 'line_planning_info_text') and self.line_planning_info_text is not None:
                self.line_planning_info_text.set_visible(False)
                self.line_planning_info_text = None
            # Clear edit mode and handles
            if self.edit_line_planning_mode:
                self.edit_line_planning_mode = False
                self.line_planning_handles = []
                self.dragging_line_planning_handle = None
                if hasattr(self, 'line_planning_edit_line') and self.line_planning_edit_line is not None:
                    self.line_planning_edit_line.remove()
                    self.line_planning_edit_line = None
                # Disconnect events
                if hasattr(self, 'line_planning_pick_cid'):
                    self.canvas.mpl_disconnect(self.line_planning_pick_cid)
                if hasattr(self, 'line_planning_motion_cid'):
                    self.canvas.mpl_disconnect(self.line_planning_motion_cid)
                if hasattr(self, 'line_planning_release_cid'):
                    self.canvas.mpl_disconnect(self.line_planning_release_cid)
                self.line_edit_btn.setText("Edit Line Planning")
                self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
            self.canvas.draw_idle()
            self.line_start_draw_btn.setText("Start Drawing Line")
            # Reset "Start Drawing Line" button style - if GeoTIFF is loaded, make it bold and orange
            if hasattr(self, 'geotiff_dataset_original') and self.geotiff_dataset_original is not None:
                if hasattr(self, 'line_start_draw_btn'):
                    self.line_start_draw_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")
            else:
                if hasattr(self, 'line_start_draw_btn'):
                    self.line_start_draw_btn.setStyleSheet("")
            self.line_planning_mode = False
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
            # Clear the profile plot
            self._clear_line_planning_profile()
            # Update button states
            self._update_line_planning_button_states()

        def _clear_line_planning_profile(self):
            """Clear the line planning profile plot."""
            if hasattr(self, 'profile_ax'):
                self.profile_ax.clear()
                self.profile_ax.set_title("Line Planning Elevation Profile", fontsize=8)
                self.profile_ax.set_xlabel("Distance (m)", fontsize=8)
                self.profile_ax.set_ylabel("Elevation (m)", fontsize=8)
                self.profile_ax.tick_params(axis='both', which='major', labelsize=7)
                if hasattr(self, 'profile_canvas'):
                    self.profile_canvas.draw_idle()

        orig_on_plot_click = getattr(self, '_on_plot_click', None)
        def new_on_plot_click(self, event):
            # Only handle if Line Planning tab is selected and mode is on
            if hasattr(self, 'param_notebook') and self.param_notebook.currentIndex() == 2 and self.line_planning_mode:
                if event.inaxes != self.ax:
                    return
                if event.button == 1:  # Left click: add point
                    lat, lon = event.ydata, event.xdata
                    if lat is None or lon is None:
                        return
                    self.line_planning_points.append((lat, lon))
                    # Draw or update the line
                    lats = [p[0] for p in self.line_planning_points]
                    lons = [p[1] for p in self.line_planning_points]
                    if self.line_planning_line is None:
                        self.line_planning_line, = self.ax.plot(lons, lats, color='orange', linewidth=2, marker='o', label='Line Planning')
                    else:
                        self.line_planning_line.set_data(lons, lats)
                    # Remove any existing temp dashed line
                    if hasattr(self, 'line_planning_temp_line') and self.line_planning_temp_line is not None:
                        self.line_planning_temp_line.remove()
                        self.line_planning_temp_line = None
                    self.canvas.draw_idle()
                elif event.button == 3:  # Right click: finish
                    # Disable line planning mode first to prevent further clicks
                    self.line_planning_mode = False
                    self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
                    
                    if len(self.line_planning_points) >= 2:
                        self.line_start_draw_btn.setText("Line Finished. Start Drawing Line")
                        # Reset "Start Drawing Line" button to normal style after line is drawn
                        if hasattr(self, 'line_start_draw_btn'):
                            self.line_start_draw_btn.setStyleSheet("")
                        # Draw the elevation profile for the drawn line
                        try:
                            self._draw_line_planning_profile()
                        except Exception as e:
                            print(f"Error drawing line planning profile: {e}")
                        # Report line summary to info/error box
                        try:
                            import pyproj
                            geod = pyproj.Geod(ellps="WGS84")
                            num_segments = len(self.line_planning_points) - 1
                            total_length_m = 0.0
                            for i in range(1, len(self.line_planning_points)):
                                lat1, lon1 = self.line_planning_points[i-1]
                                lat2, lon2 = self.line_planning_points[i]
                                _, _, d = geod.inv(lon1, lat1, lon2, lat2)
                                if not np.isnan(d) and not np.isinf(d):
                                    total_length_m += d
                            total_length_km = total_length_m / 1000.0
                            total_length_nm = total_length_m / 1852.0
                            try:
                                speed_knots = float(self.line_survey_speed_entry.text())
                            except Exception:
                                speed_knots = 8.0
                            speed_m_per_h = speed_knots * 1852
                            time_hours = total_length_m / speed_m_per_h if speed_m_per_h > 0 else 0
                            time_minutes = time_hours * 60
                            summary = (
                                f"Line Summary:\n"
                                f"Segments: {num_segments}\n"
                                f"Total Length: {total_length_m:.1f} m\n"
                                f"Total Length: {total_length_km:.3f} km\n"
                                f"Total Length: {total_length_nm:.3f} nautical miles\n"
                                f"Time to Run: {time_minutes:.1f} min\n"
                                f"Time to Run: {time_hours:.2f} hr"
                            )
                            self.set_line_info_text(summary)
                        except Exception as e:
                            self.set_line_info_text(f"Error calculating line summary: {e}")
                    else:
                        self.line_start_draw_btn.setText("Start Drawing Line")
                        if len(self.line_planning_points) > 0:
                            # Clear points if less than 2
                            self.line_planning_points = []
                            if hasattr(self, 'line_planning_line') and self.line_planning_line is not None:
                                self.line_planning_line.remove()
                                self.line_planning_line = None
                    
                    # Remove temp dashed line if present
                    if hasattr(self, 'line_planning_temp_line') and self.line_planning_temp_line is not None:
                        self.line_planning_temp_line.remove()
                        self.line_planning_temp_line = None
                    # Remove info text if present
                    if hasattr(self, 'line_planning_info_text') and self.line_planning_info_text is not None:
                        self.line_planning_info_text.set_visible(False)
                        self.line_planning_info_text = None
                    self.canvas.draw_idle()
                    # Update button states
                    self._update_line_planning_button_states()
                    return  # Do not process as other modes
                return  # Do not process as other modes
            # Otherwise, call the original handler if it exists
            if orig_on_plot_click:
                orig_on_plot_click(event)



        def _export_drawn_line(self):
            if not GEOSPATIAL_LIBS_AVAILABLE:
                self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot export line.")
                return
            if not self.line_planning_points or len(self.line_planning_points) < 2:
                self._show_message("warning","No Line", "Draw a line with at least two points before exporting.")
                return
            # Ask for export directory
            export_dir = QFileDialog.getExistingDirectory(self, "Select Export Directory", self.last_used_dir)
            if not export_dir:
                return
            self.last_used_dir = export_dir
            self._save_last_used_dir()
            # Get export name from the field
            if hasattr(self, 'line_export_name_entry'):
                export_name = self.line_export_name_entry.text().strip()
                if not export_name:
                    export_name = f"Line_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
                    self.line_export_name_entry.setText(export_name)
            else:
                # Fallback if field doesn't exist
                export_name = f"LinePlanning_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
            try:
                import json
                import fiona
                from shapely.geometry import LineString, mapping
                # --- Export to CSV ---
                # CSV files for decimal degrees always use _DD suffix
                csv_file_path = os.path.join(export_dir, f"{export_name}_DD.csv")
                with open(csv_file_path, 'w', newline='') as csvfile:
                    csv_writer = csv.writer(csvfile)
                    csv_writer.writerow(['Point Label', 'Latitude', 'Longitude'])
                    for i, (lat, lon) in enumerate(self.line_planning_points):
                        csv_writer.writerow([i + 1, lat, lon])
                
                # --- Export to DDM format (Decimal Minutes) ---
                ddm_file_path = os.path.join(export_dir, f"{export_name}_DM.csv")
                with open(ddm_file_path, 'w', newline='', encoding='utf-8') as ddmfile:
                    ddm_writer = csv.writer(ddmfile)
                    ddm_writer.writerow(['Line Number', 'Point Label', 'Latitude with Decimal Minutes', 'Longitude with Decimal Minutes'])
                    for i, (lat, lon) in enumerate(self.line_planning_points):
                        lat_ddm = self._decimal_degrees_to_ddm(lat, is_latitude=True)
                        lon_ddm = self._decimal_degrees_to_ddm(lon, is_latitude=False)
                        ddm_writer.writerow([1, i + 1, lat_ddm, lon_ddm])
                
                # --- Export to DDM text format (Decimal Minutes) ---
                ddm_txt_file_path = os.path.join(export_dir, f"{export_name}_DM.txt")
                with open(ddm_txt_file_path, 'w', encoding='utf-8') as ddm_txt_file:
                    for i, (lat, lon) in enumerate(self.line_planning_points):
                        lat_ddm = self._decimal_degrees_to_ddm(lat, is_latitude=True)
                        lon_ddm = self._decimal_degrees_to_ddm(lon, is_latitude=False)
                        ddm_txt_file.write(f"{i + 1}, {lat_ddm}, {lon_ddm}\n")
                
                # --- Export to ESRI Shapefile (.shp) ---
                schema = {
                    'geometry': 'LineString',
                    'properties': {'name': 'str'},
                }
                crs_epsg = 'EPSG:4326'  # WGS 84
                shapely_line = LineString([(lon, lat) for lat, lon in self.line_planning_points])
                features = [{
                    'geometry': mapping(shapely_line),
                    'properties': {'name': export_name},
                }]
                shapefile_path = os.path.join(export_dir, f"{export_name}.shp")
                with fiona.open(shapefile_path, 'w', driver='ESRI Shapefile', crs=crs_epsg, schema=schema) as collection:
                    collection.writerecords(features)
                # --- Export to GeoJSON ---
                geojson_file_path = os.path.join(export_dir, f"{export_name}.geojson")
                geojson_feature = {
                    "type": "Feature",
                    "geometry": mapping(shapely_line),
                    "properties": {
                        "name": export_name,
                        "points": [
                            {"point_num": i + 1, "lat": lat, "lon": lon}
                            for i, (lat, lon) in enumerate(self.line_planning_points)
                        ]
                    }
                }
                geojson_collection = {
                    "type": "FeatureCollection",
                    "features": [geojson_feature]
                }
                with open(geojson_file_path, 'w') as f:
                    json.dump(geojson_collection, f, indent=2)
                self.set_line_info_text(
                    f"Line exported successfully to:\n"
                    f"- {os.path.basename(csv_file_path)}\n"
                    f"- {os.path.basename(shapefile_path)} (and associated files)\n"
                    f"- {os.path.basename(geojson_file_path)}\n"
                    f"in directory: {export_dir}", append=False)

                # --- Export to Hypack LNW format ---
                lnw_file_path = os.path.join(export_dir, f"{export_name}.lnw")
                with open(lnw_file_path, 'w') as f:
                    f.write("LNW 1.0\n")
                    # Export each waypoint as a separate line segment
                    for i, (lat, lon) in enumerate(self.line_planning_points):
                        # Get depth from GeoTIFF if available
                        depth = self._get_depth_at_point(lat, lon) if self.geotiff_data_array is not None else 0.0
                        # Get speed from line survey speed entry
                        try:
                            speed = float(self.line_survey_speed_entry.text()) if self.line_survey_speed_entry.text() else 8.0
                        except:
                            speed = 8.0
                        f.write(f"WAYPOINT_{i+1:03d},{lat:.6f},{lon:.6f},{abs(depth):.1f},{speed:.1f},50.0,1,1\n")

                self.set_line_info_text(
                    f"Line exported successfully to:\n"
                    f"- {os.path.basename(csv_file_path)}\n"
                    f"- {os.path.basename(shapefile_path)} (and associated files)\n"
                    f"- {os.path.basename(geojson_file_path)}\n"
                    f"- {os.path.basename(lnw_file_path)}\n"
                    f"in directory: {export_dir}", append=False)

                # --- Export to SIS ASCII Plan format ---
                sis_file_path = os.path.join(export_dir, f"{export_name}.asciiplan")
                with open(sis_file_path, 'w') as f:
                    f.write("SIS ASCII Plan\n")
                    # Export each waypoint as a separate line segment
                    for i, (lat, lon) in enumerate(self.line_planning_points):
                        # Get speed from line_survey_speed_entry if available, else default to 8 knots
                        try:
                            speed_knots = float(self.line_survey_speed_entry.text()) if self.line_survey_speed_entry.text() else 8.0
                        except:
                            speed_knots = 8.0
                        
                        # Get depth at this point if GeoTIFF is loaded
                        depth = self._get_depth_at_point(lat, lon)
                        
                        # Export as waypoint
                        f.write(f"WAYPOINT_{i+1:03d}, {lat:.6f}, {lon:.6f}, {depth:.1f}, {speed_knots:.1f}, 1, 1\n")

                # --- Export PNG images ---
                # Export main map
                map_png_path = os.path.join(export_dir, f"{export_name}_map.png")
                self.figure.savefig(map_png_path, dpi=300, bbox_inches='tight', facecolor='white')
                
                # Export profile if it exists
                profile_png_path = None
                if hasattr(self, 'profile_fig') and self.profile_fig is not None:
                    profile_png_path = os.path.join(export_dir, f"{export_name}_profile.png")
                    self.profile_fig.savefig(profile_png_path, dpi=300, bbox_inches='tight', facecolor='white')
                
                # --- Export statistics text file ---
                stats = self._calculate_line_planning_statistics()
                if stats:
                    stats_file_path = os.path.join(export_dir, f"{export_name}_statistics.txt")
                    with open(stats_file_path, 'w') as f:
                        f.write("LINE PLANNING STATISTICS\n")
                        f.write("=" * 30 + "\n\n")
                        f.write(f"Number of Points: {stats['num_points']}\n")
                        f.write(f"Total Distance: {stats['total_distance_m']:.1f} m ({stats['total_distance_km']:.3f} km, {stats['total_distance_nm']:.3f} nm)\n")
                        f.write(f"Survey Speed: {stats['speed_knots']:.1f} knots\n")
                        f.write(f"Estimated Time: {stats['total_time_minutes']:.1f} min ({stats['total_time_hours']:.2f} hr)\n\n")
                        
                        if stats['depth_info']:
                            f.write(stats['depth_info'])
                        
                        # Segment information
                        if len(stats['segment_distances']) > 1:
                            f.write(f"\nSEGMENT DISTANCE AND HEADINGS\n")
                            f.write("-" * 30 + "\n")
                            for i in range(len(self.line_planning_points) - 1):
                                lat1, lon1 = self.line_planning_points[i]
                                lat2, lon2 = self.line_planning_points[i + 1]
                                if pyproj is not None:
                                    geod = pyproj.Geod(ellps="WGS84")
                                    fwd_az, back_az, dist = geod.inv(lon1, lat1, lon2, lat2)
                                    heading = fwd_az % 360
                                    f.write(f"Segment {i+1}: {dist:.1f} m ({dist/1000:.3f} km, {dist/1852:.3f} nm) - Heading: {heading:.1f}°\n")
                                else:
                                    dist = stats['segment_distances'][i]
                                    f.write(f"Segment {i+1}: {dist:.1f} m ({dist/1000:.3f} km, {dist/1852:.3f} nm) - Heading: Unable to calculate\n")
                        
                        # Waypoints with labels
                        if hasattr(self, 'line_planning_points') and self.line_planning_points:
                            f.write(f"\nWAYPOINTS\n")
                            f.write("-" * 10 + "\n")
                            for i, point in enumerate(self.line_planning_points):
                                lat, lon = point
                                lat_ddm = self._decimal_degrees_to_ddm(lat, is_latitude=True)
                                lon_ddm = self._decimal_degrees_to_ddm(lon, is_latitude=False)
                                f.write(f"WP{i+1}: {lat_ddm}, {lon_ddm} ({lat:.6f}, {lon:.6f})\n")

                # Update success message to include all exports
                success_msg = f"Line exported successfully to:\n"
                success_msg += f"- {os.path.basename(csv_file_path)}\n"
                success_msg += f"- {os.path.basename(shapefile_path)} (and associated files)\n"
                success_msg += f"- {os.path.basename(geojson_file_path)}\n"
                success_msg += f"- {os.path.basename(lnw_file_path)}\n"
                success_msg += f"- {os.path.basename(sis_file_path)}\n"
                success_msg += f"- {os.path.basename(map_png_path)}\n"
                if profile_png_path:
                    success_msg += f"- {os.path.basename(profile_png_path)}\n"
                if stats:
                    success_msg += f"- {os.path.basename(stats_file_path)}\n"
                success_msg += f"in directory: {export_dir}"
                
                self.set_line_info_text(success_msg, append=False)
            except Exception as e:
                self.set_line_info_text(f"Failed to export drawn line: {e}")

        def _import_drawn_line(self):
            """Import line planning route from CSV or GeoJSON file."""
            import csv
            import json
            
            # Open file dialog to select import file
            file_path, _ = QFileDialog.getOpenFileName(
                self,
                "Select Line Plan File to Import",
                self.last_line_import_dir,
                "Decimal Degree CSV files (*_DD.csv);;CSV files (*.csv);;GeoJSON files (*.geojson);;JSON files (*.json);;All files (*.*)"
            )
            
            if not file_path:
                return
            
            # Save the directory for next time
            import_dir = os.path.dirname(file_path)
            if import_dir and os.path.isdir(import_dir):
                self.last_line_import_dir = import_dir
                self._save_last_line_import_dir()
            
            try:
                # Clear existing line planning points
                self.line_planning_points = []
                
                # Determine file type and import accordingly
                file_ext = os.path.splitext(file_path)[1].lower()
                
                if file_ext == '.csv':
                    # Import from CSV
                    with open(file_path, 'r', encoding='utf-8') as csvfile:
                        csv_reader = csv.DictReader(csvfile)
                        
                        for row in csv_reader:
                            try:
                                # CSV format: Point Label, Latitude, Longitude
                                lat = float(row.get('Latitude', 0))
                                lon = float(row.get('Longitude', 0))
                                self.line_planning_points.append((lat, lon))
                            except (ValueError, TypeError):
                                continue
                
                elif file_ext in ['.geojson', '.json']:
                    # Import from GeoJSON
                    with open(file_path, 'r', encoding='utf-8') as f:
                        geojson_data = json.load(f)
                    
                    # Handle FeatureCollection
                    if geojson_data.get('type') == 'FeatureCollection':
                        features = geojson_data.get('features', [])
                    elif geojson_data.get('type') == 'Feature':
                        features = [geojson_data]
                    else:
                        features = []
                    
                    for feature in features:
                        if feature.get('type') != 'Feature':
                            continue
                        
                        geometry = feature.get('geometry', {})
                        if geometry.get('type') != 'LineString':
                            continue
                        
                        coordinates = geometry.get('coordinates', [])
                        if len(coordinates) < 2:
                            continue
                        
                        # GeoJSON uses [lon, lat] format, we need [lat, lon]
                        for coord in coordinates:
                            if len(coord) >= 2:
                                lon, lat = coord[0], coord[1]
                                self.line_planning_points.append((lat, lon))
                
                else:
                    self._show_message("error","Import Error", f"Unsupported file format: {file_ext}")
                    return
                
                if len(self.line_planning_points) < 2:
                    self._show_message("warning","Import Warning", "Line plan must have at least 2 points. Found fewer points in the file.")
                    self.line_planning_points = []
                    return
                
                # Extract basename from the imported file and set it in the export name field
                if hasattr(self, 'line_export_name_entry'):
                    file_basename = os.path.splitext(os.path.basename(file_path))[0]
                    # Remove common suffixes that might be added during export
                    for suffix in ['_DD', '_DM', '.geojson', '.shp', '.lnw']:
                        if file_basename.endswith(suffix):
                            file_basename = file_basename[:-len(suffix)]
                    self.line_export_name_entry.setText(file_basename)
                
                # Update UI and plot
                self._update_line_planning_button_states()
                self._plot_survey_plan(preserve_view_limits=True)
                
                # Show success message
                self.set_line_info_text(f"Successfully imported line plan with {len(self.line_planning_points)} points.")
                # Reset "Start Drawing Line" button to normal style after importing line
                if hasattr(self, 'line_start_draw_btn') and len(self.line_planning_points) >= 2:
                    self.line_start_draw_btn.setStyleSheet("")
            
            except Exception as e:
                self._show_message("error","Import Error", f"Failed to import line plan: {e}")
                import traceback
                traceback.print_exc()

        # Add mouse motion event for temp dashed line in line planning mode
        def _on_line_planning_motion(self, event):
            if not (hasattr(self, 'param_notebook') and self.param_notebook.currentIndex() == 2 and self.line_planning_mode):
                # Remove info text if present
                if hasattr(self, 'line_planning_info_text') and self.line_planning_info_text is not None:
                    self.line_planning_info_text.set_visible(False)
                    self.line_planning_info_text = None
                    self.canvas.draw_idle()
                return
            if event.inaxes != self.ax:
                return
            if len(self.line_planning_points) == 0:
                # Remove info text if present
                if hasattr(self, 'line_planning_info_text') and self.line_planning_info_text is not None:
                    self.line_planning_info_text.set_visible(False)
                    self.line_planning_info_text = None
                    self.canvas.draw_idle()
                return
            last_lat, last_lon = self.line_planning_points[-1]
            cur_lat, cur_lon = event.ydata, event.xdata
            if cur_lat is None or cur_lon is None:
                # Remove temp line if present
                if hasattr(self, 'line_planning_temp_line') and self.line_planning_temp_line is not None:
                    self.line_planning_temp_line.remove()
                    self.line_planning_temp_line = None
                    self.canvas.draw_idle()
                # Remove info text if present
                if hasattr(self, 'line_planning_info_text') and self.line_planning_info_text is not None:
                    self.line_planning_info_text.set_visible(False)
                    self.line_planning_info_text = None
                    self.canvas.draw_idle()
                return
            # Draw or update the temp dashed line
            if hasattr(self, 'line_planning_temp_line') and self.line_planning_temp_line is not None:
                self.line_planning_temp_line.set_data([last_lon, cur_lon], [last_lat, cur_lat])
            else:
                self.line_planning_temp_line, = self.ax.plot([last_lon, cur_lon], [last_lat, cur_lat], color='orange', linewidth=2, linestyle='--', alpha=0.7)
            # Calculate segment length and total length
            try:
                import pyproj
                geod = pyproj.Geod(ellps="WGS84")
                # Segment length (last point to current mouse)
                _, _, seg_len = geod.inv(last_lon, last_lat, cur_lon, cur_lat)
                # Total length (sum of all segments + current segment)
                total_len = 0.0
                if len(self.line_planning_points) > 1:
                    for i in range(1, len(self.line_planning_points)):
                        _, _, d = geod.inv(self.line_planning_points[i-1][1], self.line_planning_points[i-1][0], self.line_planning_points[i][1], self.line_planning_points[i][0])
                        total_len += d
                total_len += seg_len
            except Exception:
                seg_len = float('nan')
                total_len = float('nan')
            # Show info in the top left of the plot
            info_str = f"Lat: {cur_lat:.6f}\nLon: {cur_lon:.6f}\nSegment: {seg_len:.1f} m\nTotal: {total_len:.1f} m"
            if hasattr(self, 'line_planning_info_text') and self.line_planning_info_text is not None:
                self.line_planning_info_text.set_text(info_str)
                self.line_planning_info_text.set_visible(True)
            else:
                self.line_planning_info_text = self.ax.text(0.02, 0.98, info_str, transform=self.ax.transAxes, fontsize=9, va='top', ha='left', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7), zorder=10)
            self.canvas.draw_idle()

        # Add mouse motion event for pitch line drawing
        def _on_pitch_line_motion(self, event):
            if not (hasattr(self, 'pick_pitch_line_mode') and self.pick_pitch_line_mode and hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 1):
                # Remove info text if present
                if hasattr(self, 'pitch_line_info_text') and self.pitch_line_info_text is not None:
                    self.pitch_line_info_text.set_visible(False)
                    self.pitch_line_info_text = None
                    self.canvas.draw_idle()
                return
            if event.inaxes != self.ax:
                return
            start_lat, start_lon = self.pitch_line_points[0]
            cur_lat, cur_lon = event.ydata, event.xdata
            if cur_lat is None or cur_lon is None:
                # Remove info text if present
                if hasattr(self, 'pitch_line_info_text') and self.pitch_line_info_text is not None:
                    self.pitch_line_info_text.set_visible(False)
                    self.pitch_line_info_text = None
                    self.canvas.draw_idle()
                return
            # Calculate length, azimuth, and time
            try:
                import pyproj
                geod = pyproj.Geod(ellps="WGS84")
                az12, az21, dist = geod.inv(start_lon, start_lat, cur_lon, cur_lat)
                # Get speed from cal_survey_speed_entry if available, else default to 8 knots
                try:
                    speed_knots = float(self.cal_survey_speed_entry.text())
                except Exception:
                    speed_knots = 8.0
                speed_m_per_h = speed_knots * 1852
                time_hours = dist / speed_m_per_h if speed_m_per_h > 0 else 0
                time_minutes = time_hours * 60
            except Exception:
                az12 = float('nan')
                dist = float('nan')
                time_minutes = float('nan')
            info_str = f"Lat: {cur_lat:.6f}\nLon: {cur_lon:.6f}\nLength: {dist:.1f} m\nAzimuth: {az12:.1f}°\nTime: {time_minutes:.1f} min"
            if hasattr(self, 'pitch_line_info_text') and self.pitch_line_info_text is not None:
                self.pitch_line_info_text.set_text(info_str)
                self.pitch_line_info_text.set_visible(True)
            else:
                self.pitch_line_info_text = self.ax.text(0.02, 0.98, info_str, transform=self.ax.transAxes, fontsize=9, va='top', ha='left', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7), zorder=10)
            self.canvas.draw_idle()

        # Add mouse motion event for pick center mode
        def _on_pick_center_motion(self, event):
            if not (hasattr(self, 'pick_center_mode') and self.pick_center_mode):
                # Remove info text if present
                if hasattr(self, 'pick_center_info_text') and self.pick_center_info_text is not None:
                    self.pick_center_info_text.set_visible(False)
                    self.pick_center_info_text = None
                    self.canvas.draw_idle()
                return
            if event.inaxes != self.ax:
                return
            cur_lat, cur_lon = event.ydata, event.xdata
            if cur_lat is None or cur_lon is None:
                # Remove info text if present
                if hasattr(self, 'pick_center_info_text') and self.pick_center_info_text is not None:
                    self.pick_center_info_text.set_visible(False)
                    self.pick_center_info_text = None
                    self.canvas.draw_idle()
                return
            # Get depth and slope at this point
            elevation, slope = self._calculate_slope_at_point(cur_lat, cur_lon)
            try:
                depth_str = f"{abs(elevation):.1f} m" if elevation is not None and not np.isnan(elevation) else "-"
                slope_str = f"{slope:.1f}°" if slope is not None and not np.isnan(slope) else "-"
            except Exception:
                depth_str = "-"
                slope_str = "-"
            info_str = f"Lat: {cur_lat:.6f}\nLon: {cur_lon:.6f}\nDepth: {depth_str}\nSlope: {slope_str}"
            if hasattr(self, 'pick_center_info_text') and self.pick_center_info_text is not None:
                self.pick_center_info_text.set_text(info_str)
                self.pick_center_info_text.set_visible(True)
            else:
                self.pick_center_info_text = self.ax.text(0.02, 0.98, info_str, transform=self.ax.transAxes, fontsize=9, va='top', ha='left', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7), zorder=10)
            self.canvas.draw_idle()

        # Add mouse motion event for general info over GeoTIFF
        def _on_geotiff_hover_motion(self, event):
            # This method is now replaced by _on_mouse_motion
            pass

        # Patch methods into self
        self._toggle_line_planning_mode = _toggle_line_planning_mode.__get__(self)
        self._clear_line_planning = _clear_line_planning.__get__(self)
        self._on_plot_click = new_on_plot_click.__get__(self)
        self._export_drawn_line = _export_drawn_line.__get__(self)
        self._import_drawn_line = _import_drawn_line.__get__(self)
        self._on_line_planning_motion = _on_line_planning_motion.__get__(self)
        self._on_pitch_line_motion = _on_pitch_line_motion.__get__(self)
        self._on_pick_center_motion = _on_pick_center_motion.__get__(self)
        self._on_geotiff_hover_motion = _on_geotiff_hover_motion.__get__(self)

        self._create_widgets()

        # Add info panel for profile plot (initialize before drawing profile)
        self.profile_info_text = None

        # Create the profile plot widgets BEFORE layout
        self.profile_fig = Figure(figsize=(3, 2.2))
        self.profile_ax = self.profile_fig.add_subplot(111)
        self.profile_canvas = FigureCanvas(self.profile_fig)
        self.profile_widget = self.profile_canvas
        self._draw_crossline_profile()

        # Slope profile checkbox for profile plot
        self.show_slope_profile_var = True
        self.slope_profile_checkbox = QCheckBox("Show Slope Profile")
        self.slope_profile_checkbox.setChecked(self.show_slope_profile_var)
        self.slope_profile_checkbox.stateChanged.connect(self._draw_current_profile)

        self._setup_layout()

        # Connect Matplotlib click event for 'Pick Center from GeoTIFF'
        self.cid_click = self.canvas.mpl_connect('button_press_event', self._on_plot_click)
        # Connect Matplotlib scroll event for zoom
        self.cid_scroll = self.canvas.mpl_connect('scroll_event', self._on_scroll)
        # Connect Matplotlib motion event for real-time info display
        self.cid_motion = self.canvas.mpl_connect('motion_notify_event', self._on_mouse_motion)
        # Connect middle mouse button for panning
        self.cid_middle_press = self.canvas.mpl_connect('button_press_event', self._on_middle_press)
        self.cid_middle_motion = self.canvas.mpl_connect('motion_notify_event', self._on_middle_motion)
        self.cid_middle_release = self.canvas.mpl_connect('button_release_event', self._on_middle_release)
        # Connect draw event for dynamic colormap scaling
        self.cid_draw = self.canvas.mpl_connect('draw_event', self._on_draw_event_update_colormap)

        # Initialize buttons states based on library availability
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Missing Libraries",
                                   f"Some geospatial libraries are not found. "
                                   "GeoTIFF and Shapefile export features will be disabled. "
                                   "Install with: pip install rasterio pyproj fiona shapely")
            self.load_geotiff_btn.setEnabled(False)
            # Restored: Enable/disable for elevation_slope_combo
            if hasattr(self, 'elevation_slope_combo'):
                self.elevation_slope_combo.setEnabled(False)
            if hasattr(self, 'pick_center_btn'):
                self.pick_center_btn.setEnabled(False)
            if hasattr(self, 'remove_geotiff_btn'):
                self.remove_geotiff_btn.setEnabled(False)  # Disable Remove GeoTIFF button

        # Bind parameter changes to update the plot
        self._updating_from_code = False  # Flag to prevent recursion
        param_entries = [
            self.central_lat_entry,
            self.central_lon_entry,
            self.line_length_entry,
            self.heading_entry,
            self.dist_between_lines_entry,
            self.num_lines_entry,
            self.bisect_lead_entry,
            self.survey_speed_entry
        ]
        for entry in param_entries:
            entry.textChanged.connect(self._on_parameter_change)

        self.last_used_dir = os.path.expanduser("~")  # Default to user's home directory
        self.last_geotiff_dir = os.path.expanduser("~")  # Default to user's home directory
        self.last_survey_params_dir = os.path.expanduser("~")  # Default to user's home directory
        self.last_export_dir = os.path.expanduser("~")  # Default to user's home directory
        # Separate import directories for each survey type
        self.last_cal_import_dir = os.path.expanduser("~")
        self.last_ref_import_dir = os.path.expanduser("~")
        self.last_line_import_dir = os.path.expanduser("~")
        self._load_last_used_dir()
        self._load_last_geotiff_dir()
        self._load_last_survey_params_dir()
        self._load_last_export_dir()
        self._load_last_cal_import_dir()
        self._load_last_ref_import_dir()
        self._load_last_line_import_dir()

        # After self._setup_layout(), connect the motion event for line planning, pitch line, pick center, and geotiff hover
        # Note: These are now handled by the main _on_mouse_motion method

        # Bind parameter changes to update the Export Name
        self.central_lat_entry.textChanged.connect(self._update_export_name)
        # Reset "Pick Center from GeoTIFF" button to normal style when Central Latitude is entered
        self.central_lat_entry.textChanged.connect(lambda: self.pick_center_btn.setStyleSheet("") if hasattr(self, 'pick_center_btn') and self.central_lat_entry.text().strip() else None)
        self.central_lon_entry.textChanged.connect(self._update_export_name)
        self.line_length_entry.textChanged.connect(self._update_export_name)
        self.heading_entry.textChanged.connect(self._update_export_name)
        self.dist_between_lines_entry.textChanged.connect(self._update_export_name)
        self.num_lines_entry.textChanged.connect(self._update_export_name)
        self.bisect_lead_entry.textChanged.connect(self._update_export_name)
        self.survey_speed_entry.textChanged.connect(self._update_export_name)
        self.offset_direction_combo.currentTextChanged.connect(self._update_export_name)
        # Slider connections already set up in widget creation
        
        # Show window maximized after everything is set up
        # Set window to minimum size on startup
        self.resize(1600, 1150)
        print("Initialization complete, window should be visible")
        
    # Helper methods for dialog conversions
    def _show_message(self, msg_type, title, message):
        """Helper method to show message boxes."""
        msg = QMessageBox(self)
        msg.setWindowTitle(title)
        msg.setText(message)
        if msg_type == "error":
            msg.setIcon(QMessageBox.Icon.Critical)
        elif msg_type == "warning":
            msg.setIcon(QMessageBox.Icon.Warning)
        elif msg_type == "info":
            msg.setIcon(QMessageBox.Icon.Information)
        elif msg_type == "question":
            msg.setIcon(QMessageBox.Icon.Question)
        msg.exec()
        return msg
        
    def _ask_yes_no(self, title, message):
        """Helper method for yes/no dialogs."""
        reply = QMessageBox.question(self, title, message, 
                                     QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
        return reply == QMessageBox.StandardButton.Yes
        
    def _ask_ok_cancel(self, title, message):
        """Helper method for ok/cancel dialogs."""
        reply = QMessageBox.question(self, title, message,
                                     QMessageBox.StandardButton.Ok | QMessageBox.StandardButton.Cancel)
        return reply == QMessageBox.StandardButton.Ok

    def set_cal_info_text(self, message, append=False):
        """Add a message to the Activity Log with [Calibration] prefix. Maintains up to 200 lines of history."""
        if not hasattr(self, 'activity_log_text'):
            return
        prefixed_message = f"[Calibration] {message}"
        self.activity_log_text.setReadOnly(False)
        if append:
            # Prepend new message (newest at top) - only when explicitly requested
            current_text = self.activity_log_text.toPlainText().rstrip('\n')
            new_text = prefixed_message + "\n" + current_text if current_text else prefixed_message + "\n"
            # Split into lines and limit to 200
            lines = new_text.splitlines()
            if len(lines) > 200:
                lines = lines[:200]
            self.activity_log_text.setPlainText("\n".join(lines) + "\n")
            # Move cursor to top
            cursor = self.activity_log_text.textCursor()
            cursor.movePosition(QTextCursor.MoveOperation.Start)
            self.activity_log_text.setTextCursor(cursor)
        else:
            # Append to end (default behavior - newest at bottom)
            self.activity_log_text.append(prefixed_message)
            # Maintain up to 200 lines by removing oldest from top
            num_lines = self.activity_log_text.document().blockCount()
            if num_lines > 200:
                cursor = self.activity_log_text.textCursor()
                cursor.movePosition(QTextCursor.MoveOperation.Start)
                cursor.movePosition(QTextCursor.MoveOperation.Down, QTextCursor.MoveMode.MoveAnchor, num_lines - 200)
                cursor.movePosition(QTextCursor.MoveOperation.Start, QTextCursor.MoveMode.KeepAnchor)
                cursor.removeSelectedText()
            # Move cursor to end and scroll to show newest message
            cursor = self.activity_log_text.textCursor()
            cursor.movePosition(QTextCursor.MoveOperation.End)
            self.activity_log_text.setTextCursor(cursor)
        self.activity_log_text.ensureCursorVisible()
        self.activity_log_text.setReadOnly(True)

        self.survey_lines_data = []  # Stores [(lat1, lon1), (lat2, lon2)] for each line
        self.cross_line_data = []  # Stores [(lat1, lon1), (lat2, lon2)] for the crossline
        self.central_point_coords = (None, None)  # Stores (lat, lon) of central point

        self.geotiff_dataset_original = None  # Original rasterio dataset (keeps its CRS)
        self.geotiff_data_array = None  # NumPy array of *reprojected* GeoTIFF elevation data (WGS84)
        self.geotiff_extent = None  # [left, right, bottom, top] of *reprojected* data
        self.geotiff_image_plot = None  # Matplotlib imshow object for GeoTIFF
        self.geotiff_hillshade_plot = None  # Matplotlib imshow object for hillshade
        self.slope_colorbar = None  # To store the colorbar object for slope visualization
        self.elevation_colorbar = None  # To store the colorbar object for elevation visualization
        
        # Dynamic resolution variables
        self.geotiff_full_resolution = None  # Store full resolution data for dynamic loading
        self.geotiff_current_resolution = None  # Current resolution being displayed
        self.geotiff_zoom_level = 1.0  # Current zoom level (1.0 = full resolution)
        self.geotiff_zoom_threshold = 0.1  # Zoom threshold for switching to higher resolution
        self.dynamic_resolution_enabled = True  # Toggle for dynamic resolution feature

        # Restored: Variable for controlling elevation/slope overlay
        self.geotiff_display_mode = "elevation"  # "elevation" or "slope"

        self.hillshade_vs_slope_viz_mode = "hillshade"  # "hillshade" or "slope_viz" (for main underlay)

        # Contour display variables
        self.show_contours_var = False  # Default to off
        self.contour_plot = None  # Store contour plot object for removal/update

        # Variable for the new line length multiplier slider
        self.line_length_multiplier = 8.0  # Default multiplier
        # Variable for the new distance between lines multiplier slider
        self.dist_between_lines_multiplier = 1.0  # Default multiplier

        # Initialize the picking mode state
        self.pick_center_mode = False

        # Pyproj transformers
        self.geod = None  # Geodetic object for precise distance/bearing on ellipsoid
        self.local_proj_transformer = None  # For transforming lat/lon to local projected CRS (easting, northing)
        self.wgs84_transformer = None  # For transforming local projected CRS back to WGS84
        self.geotiff_to_wgs84_transformer = None  # For transforming GeoTIFF CRS to WGS84 for display
        self.wgs84_to_geotiff_transformer = None  # For transforming WGS84 to GeoTIFF's original CRS for picking

        self.pitch_line_points = []  # Stores [(lat1, lon1), (lat2, lon2)] for pitch line
        self.pick_pitch_line_mode = False  # Flag for pitch line picking mode
        self.edit_pitch_line_mode = False  # Flag for pitch line editing mode
        self.pitch_line_start_handle = None  # Draggable handle for start point
        self.pitch_line_end_handle = None  # Draggable handle for end point
        self.dragging_pitch_line_handle = None  # Which handle is being dragged
        self.pitch_line_temp_line = None  # Temporary line during drawing
        self.pitch_line_annotation = None  # Annotation for pitch line
        self.pitch_line_tooltip = None  # Tooltip widget
        self.pitch_line_tooltip_text = None  # Tooltip text
        self.pitch_line_edit_line = None  # Temporary line during editing
        self.heading_lines = []  # Stores two heading lines as [(lat1, lon1), (lat2, lon2)]
        self.roll_line_points = []  # Stores [(lat1, lon1), (lat2, lon2)] for roll line
        self.pick_roll_line_mode = False  # Flag for roll line picking mode
        self.edit_roll_line_mode = False  # Flag for roll line editing mode
        self.roll_line_start_handle = None  # Draggable handle for start point
        self.roll_line_end_handle = None  # Draggable handle for end point
        self.dragging_roll_line_handle = None  # Which handle is being dragged
        self.roll_line_edit_line = None  # Temporary line during editing

        # --- Add state for line planning ---
        self.line_planning_mode = False
        self.line_planning_points = []
        self.line_planning_line = None
        self.edit_line_planning_mode = False
        self.line_planning_handles = []
        self.dragging_line_planning_handle = None
        self.line_planning_edit_line = None

        # Track previous tab for clearing on tab switch
        self.previous_tab_index = 0
        self.tab_switch_initialized = False  # Flag to prevent clearing on first tab change

        # --- Patch in line planning methods and event handler BEFORE widgets are created ---
        def _toggle_line_planning_mode(self):
            if not GEOSPATIAL_LIBS_AVAILABLE or self.geotiff_data_array is None:
                self._show_message("warning", "No GeoTIFF", "Load a GeoTIFF first to draw a line.")
                return
            self.line_planning_mode = not self.line_planning_mode
            if self.line_planning_mode:
                self.line_planning_points = []
                if self.line_planning_line is not None:
                    self.line_planning_line.remove()
                    self.line_planning_line = None
                self.line_start_draw_btn.setText("Drawing: Left-click to add, Right-click to finish")
                self.canvas_widget.setCursor(Qt.CursorShape.CrossCursor)
                self.set_line_info_text("Left click to start line, left click to add waypoints, and Right click to end line")
            else:
                self.line_start_draw_btn.setText("Start Drawing Line")
                self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)

        def _clear_line_planning(self):
            self.line_planning_points = []
            if self.line_planning_line is not None:
                self.line_planning_line.remove()
                self.line_planning_line = None
            # Remove temp dashed line if present
            if hasattr(self, 'line_planning_temp_line') and self.line_planning_temp_line is not None:
                self.line_planning_temp_line.remove()
                self.line_planning_temp_line = None
            # Remove info text if present
            if hasattr(self, 'line_planning_info_text') and self.line_planning_info_text is not None:
                self.line_planning_info_text.set_visible(False)
                self.line_planning_info_text = None
            # Clear edit mode and handles
            if self.edit_line_planning_mode:
                self.edit_line_planning_mode = False
                self.line_planning_handles = []
                self.dragging_line_planning_handle = None
                if hasattr(self, 'line_planning_edit_line') and self.line_planning_edit_line is not None:
                    self.line_planning_edit_line.remove()
                    self.line_planning_edit_line = None
                # Disconnect events
                if hasattr(self, 'line_planning_pick_cid'):
                    self.canvas.mpl_disconnect(self.line_planning_pick_cid)
                if hasattr(self, 'line_planning_motion_cid'):
                    self.canvas.mpl_disconnect(self.line_planning_motion_cid)
                if hasattr(self, 'line_planning_release_cid'):
                    self.canvas.mpl_disconnect(self.line_planning_release_cid)
                self.line_edit_btn.setText("Edit Line Planning")
                self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
            self.canvas.draw_idle()
            self.line_start_draw_btn.setText("Start Drawing Line")
            # Reset "Start Drawing Line" button style - if GeoTIFF is loaded, make it bold and orange
            if hasattr(self, 'geotiff_dataset_original') and self.geotiff_dataset_original is not None:
                if hasattr(self, 'line_start_draw_btn'):
                    self.line_start_draw_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")
            else:
                if hasattr(self, 'line_start_draw_btn'):
                    self.line_start_draw_btn.setStyleSheet("")
            self.line_planning_mode = False
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
            # Clear the profile plot
            self._clear_line_planning_profile()
            # Update button states
            self._update_line_planning_button_states()

        def _clear_line_planning_profile(self):
            """Clear the line planning profile plot."""
            if hasattr(self, 'profile_ax'):
                self.profile_ax.clear()
                self.profile_ax.set_title("Line Planning Elevation Profile", fontsize=8)
                self.profile_ax.set_xlabel("Distance (m)", fontsize=8)
                self.profile_ax.set_ylabel("Elevation (m)", fontsize=8)
                self.profile_ax.tick_params(axis='both', which='major', labelsize=7)
                if hasattr(self, 'profile_canvas'):
                    self.profile_canvas.draw_idle()

        orig_on_plot_click = getattr(self, '_on_plot_click', None)
        def new_on_plot_click(self, event):
            # Only handle if Line Planning tab is selected and mode is on
            if hasattr(self, 'param_notebook') and self.param_notebook.currentIndex() == 2 and self.line_planning_mode:
                if event.inaxes != self.ax:
                    return
                if event.button == 1:  # Left click: add point
                    lat, lon = event.ydata, event.xdata
                    if lat is None or lon is None:
                        return
                    self.line_planning_points.append((lat, lon))
                    # Draw or update the line
                    lats = [p[0] for p in self.line_planning_points]
                    lons = [p[1] for p in self.line_planning_points]
                    if self.line_planning_line is None:
                        self.line_planning_line, = self.ax.plot(lons, lats, color='orange', linewidth=2, marker='o', label='Line Planning')
                    else:
                        self.line_planning_line.set_data(lons, lats)
                    # Remove any existing temp dashed line
                    if hasattr(self, 'line_planning_temp_line') and self.line_planning_temp_line is not None:
                        self.line_planning_temp_line.remove()
                        self.line_planning_temp_line = None
                    self.canvas.draw_idle()
                elif event.button == 3:  # Right click: finish
                    # Disable line planning mode first to prevent further clicks
                    self.line_planning_mode = False
                    self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
                    
                    if len(self.line_planning_points) >= 2:
                        self.line_start_draw_btn.setText("Line Finished. Start Drawing Line")
                        # Reset "Start Drawing Line" button to normal style after line is drawn
                        if hasattr(self, 'line_start_draw_btn'):
                            self.line_start_draw_btn.setStyleSheet("")
                        # Draw the elevation profile for the drawn line
                        try:
                            self._draw_line_planning_profile()
                        except Exception as e:
                            print(f"Error drawing line planning profile: {e}")
                        # Report line summary to info/error box
                        try:
                            import pyproj
                            geod = pyproj.Geod(ellps="WGS84")
                            num_segments = len(self.line_planning_points) - 1
                            total_length_m = 0.0
                            for i in range(1, len(self.line_planning_points)):
                                lat1, lon1 = self.line_planning_points[i-1]
                                lat2, lon2 = self.line_planning_points[i]
                                _, _, d = geod.inv(lon1, lat1, lon2, lat2)
                                if not np.isnan(d) and not np.isinf(d):
                                    total_length_m += d
                            total_length_km = total_length_m / 1000.0
                            total_length_nm = total_length_m / 1852.0
                            try:
                                speed_knots = float(self.line_survey_speed_entry.text())
                            except Exception:
                                speed_knots = 8.0
                            speed_m_per_h = speed_knots * 1852
                            time_hours = total_length_m / speed_m_per_h if speed_m_per_h > 0 else 0
                            time_minutes = time_hours * 60
                            summary = (
                                f"Line Summary:\n"
                                f"Segments: {num_segments}\n"
                                f"Total Length: {total_length_m:.1f} m\n"
                                f"Total Length: {total_length_km:.3f} km\n"
                                f"Total Length: {total_length_nm:.3f} nautical miles\n"
                                f"Time to Run: {time_minutes:.1f} min\n"
                                f"Time to Run: {time_hours:.2f} hr"
                            )
                            self.set_line_info_text(summary)
                        except Exception as e:
                            self.set_line_info_text(f"Error calculating line summary: {e}")
                    else:
                        self.line_start_draw_btn.setText("Start Drawing Line")
                        if len(self.line_planning_points) > 0:
                            # Clear points if less than 2
                            self.line_planning_points = []
                            if hasattr(self, 'line_planning_line') and self.line_planning_line is not None:
                                self.line_planning_line.remove()
                                self.line_planning_line = None
                    
                    # Remove temp dashed line if present
                    if hasattr(self, 'line_planning_temp_line') and self.line_planning_temp_line is not None:
                        self.line_planning_temp_line.remove()
                        self.line_planning_temp_line = None
                    # Remove info text if present
                    if hasattr(self, 'line_planning_info_text') and self.line_planning_info_text is not None:
                        self.line_planning_info_text.set_visible(False)
                        self.line_planning_info_text = None
                    self.canvas.draw_idle()
                    # Update button states
                    self._update_line_planning_button_states()
                    return  # Do not process as other modes
                return  # Do not process as other modes
            # Otherwise, call the original handler if it exists
            if orig_on_plot_click:
                orig_on_plot_click(event)



        def _export_drawn_line(self):
            if not GEOSPATIAL_LIBS_AVAILABLE:
                self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot export line.")
                return
            if not self.line_planning_points or len(self.line_planning_points) < 2:
                self._show_message("warning","No Line", "Draw a line with at least two points before exporting.")
                return
            # Ask for export directory
            export_dir = QFileDialog.getExistingDirectory(self, "Select Export Directory", self.last_used_dir)
            if not export_dir:
                return
            self.last_used_dir = export_dir
            self._save_last_used_dir()
            # Get export name from the field
            if hasattr(self, 'line_export_name_entry'):
                export_name = self.line_export_name_entry.text().strip()
                if not export_name:
                    export_name = f"Line_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
                    self.line_export_name_entry.setText(export_name)
            else:
                # Fallback if field doesn't exist
                export_name = f"LinePlanning_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
            try:
                import json
                import fiona
                from shapely.geometry import LineString, mapping
                # --- Export to CSV ---
                # CSV files for decimal degrees always use _DD suffix
                csv_file_path = os.path.join(export_dir, f"{export_name}_DD.csv")
                with open(csv_file_path, 'w', newline='') as csvfile:
                    csv_writer = csv.writer(csvfile)
                    csv_writer.writerow(['Point Label', 'Latitude', 'Longitude'])
                    for i, (lat, lon) in enumerate(self.line_planning_points):
                        csv_writer.writerow([i + 1, lat, lon])
                
                # --- Export to DDM format (Decimal Minutes) ---
                ddm_file_path = os.path.join(export_dir, f"{export_name}_DM.csv")
                with open(ddm_file_path, 'w', newline='', encoding='utf-8') as ddmfile:
                    ddm_writer = csv.writer(ddmfile)
                    ddm_writer.writerow(['Line Number', 'Point Label', 'Latitude with Decimal Minutes', 'Longitude with Decimal Minutes'])
                    for i, (lat, lon) in enumerate(self.line_planning_points):
                        lat_ddm = self._decimal_degrees_to_ddm(lat, is_latitude=True)
                        lon_ddm = self._decimal_degrees_to_ddm(lon, is_latitude=False)
                        ddm_writer.writerow([1, i + 1, lat_ddm, lon_ddm])
                
                # --- Export to DDM text format (Decimal Minutes) ---
                ddm_txt_file_path = os.path.join(export_dir, f"{export_name}_DM.txt")
                with open(ddm_txt_file_path, 'w', encoding='utf-8') as ddm_txt_file:
                    for i, (lat, lon) in enumerate(self.line_planning_points):
                        lat_ddm = self._decimal_degrees_to_ddm(lat, is_latitude=True)
                        lon_ddm = self._decimal_degrees_to_ddm(lon, is_latitude=False)
                        ddm_txt_file.write(f"{i + 1}, {lat_ddm}, {lon_ddm}\n")
                
                # --- Export to ESRI Shapefile (.shp) ---
                schema = {
                    'geometry': 'LineString',
                    'properties': {'name': 'str'},
                }
                crs_epsg = 'EPSG:4326'  # WGS 84
                shapely_line = LineString([(lon, lat) for lat, lon in self.line_planning_points])
                features = [{
                    'geometry': mapping(shapely_line),
                    'properties': {'name': export_name},
                }]
                shapefile_path = os.path.join(export_dir, f"{export_name}.shp")
                with fiona.open(shapefile_path, 'w', driver='ESRI Shapefile', crs=crs_epsg, schema=schema) as collection:
                    collection.writerecords(features)
                # --- Export to GeoJSON ---
                geojson_file_path = os.path.join(export_dir, f"{export_name}.geojson")
                geojson_feature = {
                    "type": "Feature",
                    "geometry": mapping(shapely_line),
                    "properties": {
                        "name": export_name,
                        "points": [
                            {"point_num": i + 1, "lat": lat, "lon": lon}
                            for i, (lat, lon) in enumerate(self.line_planning_points)
                        ]
                    }
                }
                geojson_collection = {
                    "type": "FeatureCollection",
                    "features": [geojson_feature]
                }
                with open(geojson_file_path, 'w') as f:
                    json.dump(geojson_collection, f, indent=2)
                self.set_line_info_text(
                    f"Line exported successfully to:\n"
                    f"- {os.path.basename(csv_file_path)}\n"
                    f"- {os.path.basename(shapefile_path)} (and associated files)\n"
                    f"- {os.path.basename(geojson_file_path)}\n"
                    f"in directory: {export_dir}", append=False)

                # --- Export to Hypack LNW format ---
                lnw_file_path = os.path.join(export_dir, f"{export_name}.lnw")
                with open(lnw_file_path, 'w') as f:
                    f.write("LNW 1.0\n")
                    # Export each waypoint as a separate line segment
                    for i, (lat, lon) in enumerate(self.line_planning_points):
                        # Get depth from GeoTIFF if available
                        depth = self._get_depth_at_point(lat, lon) if self.geotiff_data_array is not None else 0.0
                        # Get speed from line survey speed entry
                        try:
                            speed = float(self.line_survey_speed_entry.text()) if self.line_survey_speed_entry.text() else 8.0
                        except:
                            speed = 8.0
                        f.write(f"WAYPOINT_{i+1:03d},{lat:.6f},{lon:.6f},{abs(depth):.1f},{speed:.1f},50.0,1,1\n")

                self.set_line_info_text(
                    f"Line exported successfully to:\n"
                    f"- {os.path.basename(csv_file_path)}\n"
                    f"- {os.path.basename(shapefile_path)} (and associated files)\n"
                    f"- {os.path.basename(geojson_file_path)}\n"
                    f"- {os.path.basename(lnw_file_path)}\n"
                    f"in directory: {export_dir}", append=False)

                # --- Export to SIS ASCII Plan format ---
                sis_file_path = os.path.join(export_dir, f"{export_name}.asciiplan")
                with open(sis_file_path, 'w') as f:
                    f.write("SIS ASCII Plan\n")
                    # Export each waypoint as a separate line segment
                    for i, (lat, lon) in enumerate(self.line_planning_points):
                        # Get speed from line_survey_speed_entry if available, else default to 8 knots
                        try:
                            speed_knots = float(self.line_survey_speed_entry.text()) if self.line_survey_speed_entry.text() else 8.0
                        except:
                            speed_knots = 8.0
                        
                        # Get depth at this point if GeoTIFF is loaded
                        depth = self._get_depth_at_point(lat, lon)
                        
                        # Export as waypoint
                        f.write(f"WAYPOINT_{i+1:03d}, {lat:.6f}, {lon:.6f}, {depth:.1f}, {speed_knots:.1f}, 1, 1\n")

                # --- Export PNG images ---
                # Export main map
                map_png_path = os.path.join(export_dir, f"{export_name}_map.png")
                self.figure.savefig(map_png_path, dpi=300, bbox_inches='tight', facecolor='white')
                
                # Export profile if it exists
                profile_png_path = None
                if hasattr(self, 'profile_fig') and self.profile_fig is not None:
                    profile_png_path = os.path.join(export_dir, f"{export_name}_profile.png")
                    self.profile_fig.savefig(profile_png_path, dpi=300, bbox_inches='tight', facecolor='white')
                
                # --- Export statistics text file ---
                stats = self._calculate_line_planning_statistics()
                if stats:
                    stats_file_path = os.path.join(export_dir, f"{export_name}_statistics.txt")
                    with open(stats_file_path, 'w') as f:
                        f.write("LINE PLANNING STATISTICS\n")
                        f.write("=" * 30 + "\n\n")
                        f.write(f"Number of Points: {stats['num_points']}\n")
                        f.write(f"Total Distance: {stats['total_distance_m']:.1f} m ({stats['total_distance_km']:.3f} km, {stats['total_distance_nm']:.3f} nm)\n")
                        f.write(f"Survey Speed: {stats['speed_knots']:.1f} knots\n")
                        f.write(f"Estimated Time: {stats['total_time_minutes']:.1f} min ({stats['total_time_hours']:.2f} hr)\n\n")
                        
                        if stats['depth_info']:
                            f.write(stats['depth_info'])
                        
                        # Segment information
                        if len(stats['segment_distances']) > 1:
                            f.write(f"\nSEGMENT DISTANCE AND HEADINGS\n")
                            f.write("-" * 30 + "\n")
                            for i in range(len(self.line_planning_points) - 1):
                                lat1, lon1 = self.line_planning_points[i]
                                lat2, lon2 = self.line_planning_points[i + 1]
                                if pyproj is not None:
                                    geod = pyproj.Geod(ellps="WGS84")
                                    fwd_az, back_az, dist = geod.inv(lon1, lat1, lon2, lat2)
                                    heading = fwd_az % 360
                                    f.write(f"Segment {i+1}: {dist:.1f} m ({dist/1000:.3f} km, {dist/1852:.3f} nm) - Heading: {heading:.1f}°\n")
                                else:
                                    dist = stats['segment_distances'][i]
                                    f.write(f"Segment {i+1}: {dist:.1f} m ({dist/1000:.3f} km, {dist/1852:.3f} nm) - Heading: Unable to calculate\n")
                        
                        # Waypoints with labels
                        if hasattr(self, 'line_planning_points') and self.line_planning_points:
                            f.write(f"\nWAYPOINTS\n")
                            f.write("-" * 10 + "\n")
                            for i, point in enumerate(self.line_planning_points):
                                lat, lon = point
                                lat_ddm = self._decimal_degrees_to_ddm(lat, is_latitude=True)
                                lon_ddm = self._decimal_degrees_to_ddm(lon, is_latitude=False)
                                f.write(f"WP{i+1}: {lat_ddm}, {lon_ddm} ({lat:.6f}, {lon:.6f})\n")

                # Update success message to include all exports
                success_msg = f"Line exported successfully to:\n"
                success_msg += f"- {os.path.basename(csv_file_path)}\n"
                success_msg += f"- {os.path.basename(shapefile_path)} (and associated files)\n"
                success_msg += f"- {os.path.basename(geojson_file_path)}\n"
                success_msg += f"- {os.path.basename(lnw_file_path)}\n"
                success_msg += f"- {os.path.basename(sis_file_path)}\n"
                success_msg += f"- {os.path.basename(map_png_path)}\n"
                if profile_png_path:
                    success_msg += f"- {os.path.basename(profile_png_path)}\n"
                if stats:
                    success_msg += f"- {os.path.basename(stats_file_path)}\n"
                success_msg += f"in directory: {export_dir}"
                
                self.set_line_info_text(success_msg, append=False)
            except Exception as e:
                self.set_line_info_text(f"Failed to export drawn line: {e}")

        def _import_drawn_line(self):
            """Import line planning route from CSV or GeoJSON file."""
            import csv
            import json
            
            # Open file dialog to select import file
            file_path, _ = QFileDialog.getOpenFileName(
                self,
                "Select Line Plan File to Import",
                self.last_line_import_dir,
                "Decimal Degree CSV files (*_DD.csv);;CSV files (*.csv);;GeoJSON files (*.geojson);;JSON files (*.json);;All files (*.*)"
            )
            
            if not file_path:
                return
            
            # Save the directory for next time
            import_dir = os.path.dirname(file_path)
            if import_dir and os.path.isdir(import_dir):
                self.last_line_import_dir = import_dir
                self._save_last_line_import_dir()
            
            try:
                # Clear existing line planning points
                self.line_planning_points = []
                
                # Determine file type and import accordingly
                file_ext = os.path.splitext(file_path)[1].lower()
                
                if file_ext == '.csv':
                    # Import from CSV
                    with open(file_path, 'r', encoding='utf-8') as csvfile:
                        csv_reader = csv.DictReader(csvfile)
                        
                        for row in csv_reader:
                            try:
                                # CSV format: Point Label, Latitude, Longitude
                                lat = float(row.get('Latitude', 0))
                                lon = float(row.get('Longitude', 0))
                                self.line_planning_points.append((lat, lon))
                            except (ValueError, TypeError):
                                continue
                
                elif file_ext in ['.geojson', '.json']:
                    # Import from GeoJSON
                    with open(file_path, 'r', encoding='utf-8') as f:
                        geojson_data = json.load(f)
                    
                    # Handle FeatureCollection
                    if geojson_data.get('type') == 'FeatureCollection':
                        features = geojson_data.get('features', [])
                    elif geojson_data.get('type') == 'Feature':
                        features = [geojson_data]
                    else:
                        features = []
                    
                    for feature in features:
                        if feature.get('type') != 'Feature':
                            continue
                        
                        geometry = feature.get('geometry', {})
                        if geometry.get('type') != 'LineString':
                            continue
                        
                        coordinates = geometry.get('coordinates', [])
                        if len(coordinates) < 2:
                            continue
                        
                        # GeoJSON uses [lon, lat] format, we need [lat, lon]
                        for coord in coordinates:
                            if len(coord) >= 2:
                                lon, lat = coord[0], coord[1]
                                self.line_planning_points.append((lat, lon))
                
                else:
                    self._show_message("error","Import Error", f"Unsupported file format: {file_ext}")
                    return
                
                if len(self.line_planning_points) < 2:
                    self._show_message("warning","Import Warning", "Line plan must have at least 2 points. Found fewer points in the file.")
                    self.line_planning_points = []
                    return
                
                # Extract basename from the imported file and set it in the export name field
                if hasattr(self, 'line_export_name_entry'):
                    file_basename = os.path.splitext(os.path.basename(file_path))[0]
                    # Remove common suffixes that might be added during export
                    for suffix in ['_DD', '_DM', '.geojson', '.shp', '.lnw']:
                        if file_basename.endswith(suffix):
                            file_basename = file_basename[:-len(suffix)]
                    self.line_export_name_entry.setText(file_basename)
                
                # Update UI and plot
                self._update_line_planning_button_states()
                self._plot_survey_plan(preserve_view_limits=True)
                
                # Show success message
                self.set_line_info_text(f"Successfully imported line plan with {len(self.line_planning_points)} points.")
                # Reset "Start Drawing Line" button to normal style after importing line
                if hasattr(self, 'line_start_draw_btn') and len(self.line_planning_points) >= 2:
                    self.line_start_draw_btn.setStyleSheet("")
            
            except Exception as e:
                self._show_message("error","Import Error", f"Failed to import line plan: {e}")
                import traceback
                traceback.print_exc()

        # Add mouse motion event for temp dashed line in line planning mode
        def _on_line_planning_motion(self, event):
            if not (hasattr(self, 'param_notebook') and self.param_notebook.currentIndex() == 2 and self.line_planning_mode):
                # Remove info text if present
                if hasattr(self, 'line_planning_info_text') and self.line_planning_info_text is not None:
                    self.line_planning_info_text.set_visible(False)
                    self.line_planning_info_text = None
                    self.canvas.draw_idle()
                return
            if event.inaxes != self.ax:
                return
            if len(self.line_planning_points) == 0:
                # Remove info text if present
                if hasattr(self, 'line_planning_info_text') and self.line_planning_info_text is not None:
                    self.line_planning_info_text.set_visible(False)
                    self.line_planning_info_text = None
                    self.canvas.draw_idle()
                return
            last_lat, last_lon = self.line_planning_points[-1]
            cur_lat, cur_lon = event.ydata, event.xdata
            if cur_lat is None or cur_lon is None:
                # Remove temp line if present
                if hasattr(self, 'line_planning_temp_line') and self.line_planning_temp_line is not None:
                    self.line_planning_temp_line.remove()
                    self.line_planning_temp_line = None
                    self.canvas.draw_idle()
                # Remove info text if present
                if hasattr(self, 'line_planning_info_text') and self.line_planning_info_text is not None:
                    self.line_planning_info_text.set_visible(False)
                    self.line_planning_info_text = None
                    self.canvas.draw_idle()
                return
            # Draw or update the temp dashed line
            if hasattr(self, 'line_planning_temp_line') and self.line_planning_temp_line is not None:
                self.line_planning_temp_line.set_data([last_lon, cur_lon], [last_lat, cur_lat])
            else:
                self.line_planning_temp_line, = self.ax.plot([last_lon, cur_lon], [last_lat, cur_lat], color='orange', linewidth=2, linestyle='--', alpha=0.7)
            # Calculate segment length and total length
            try:
                import pyproj
                geod = pyproj.Geod(ellps="WGS84")
                # Segment length (last point to current mouse)
                _, _, seg_len = geod.inv(last_lon, last_lat, cur_lon, cur_lat)
                # Total length (sum of all segments + current segment)
                total_len = 0.0
                if len(self.line_planning_points) > 1:
                    for i in range(1, len(self.line_planning_points)):
                        _, _, d = geod.inv(self.line_planning_points[i-1][1], self.line_planning_points[i-1][0], self.line_planning_points[i][1], self.line_planning_points[i][0])
                        total_len += d
                total_len += seg_len
            except Exception:
                seg_len = float('nan')
                total_len = float('nan')
            # Show info in the top left of the plot
            info_str = f"Lat: {cur_lat:.6f}\nLon: {cur_lon:.6f}\nSegment: {seg_len:.1f} m\nTotal: {total_len:.1f} m"
            if hasattr(self, 'line_planning_info_text') and self.line_planning_info_text is not None:
                self.line_planning_info_text.set_text(info_str)
                self.line_planning_info_text.set_visible(True)
            else:
                self.line_planning_info_text = self.ax.text(0.02, 0.98, info_str, transform=self.ax.transAxes, fontsize=9, va='top', ha='left', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7), zorder=10)
            self.canvas.draw_idle()

        # Add mouse motion event for pitch line drawing
        def _on_pitch_line_motion(self, event):
            if not (hasattr(self, 'pick_pitch_line_mode') and self.pick_pitch_line_mode and hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 1):
                # Remove info text if present
                if hasattr(self, 'pitch_line_info_text') and self.pitch_line_info_text is not None:
                    self.pitch_line_info_text.set_visible(False)
                    self.pitch_line_info_text = None
                    self.canvas.draw_idle()
                return
            if event.inaxes != self.ax:
                return
            start_lat, start_lon = self.pitch_line_points[0]
            cur_lat, cur_lon = event.ydata, event.xdata
            if cur_lat is None or cur_lon is None:
                # Remove info text if present
                if hasattr(self, 'pitch_line_info_text') and self.pitch_line_info_text is not None:
                    self.pitch_line_info_text.set_visible(False)
                    self.pitch_line_info_text = None
                    self.canvas.draw_idle()
                return
            # Calculate length, azimuth, and time
            try:
                import pyproj
                geod = pyproj.Geod(ellps="WGS84")
                az12, az21, dist = geod.inv(start_lon, start_lat, cur_lon, cur_lat)
                # Get speed from cal_survey_speed_entry if available, else default to 8 knots
                try:
                    speed_knots = float(self.cal_survey_speed_entry.text())
                except Exception:
                    speed_knots = 8.0
                speed_m_per_h = speed_knots * 1852
                time_hours = dist / speed_m_per_h if speed_m_per_h > 0 else 0
                time_minutes = time_hours * 60
            except Exception:
                az12 = float('nan')
                dist = float('nan')
                time_minutes = float('nan')
            info_str = f"Lat: {cur_lat:.6f}\nLon: {cur_lon:.6f}\nLength: {dist:.1f} m\nAzimuth: {az12:.1f}°\nTime: {time_minutes:.1f} min"
            if hasattr(self, 'pitch_line_info_text') and self.pitch_line_info_text is not None:
                self.pitch_line_info_text.set_text(info_str)
                self.pitch_line_info_text.set_visible(True)
            else:
                self.pitch_line_info_text = self.ax.text(0.02, 0.98, info_str, transform=self.ax.transAxes, fontsize=9, va='top', ha='left', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7), zorder=10)
            self.canvas.draw_idle()

        # Add mouse motion event for pick center mode
        def _on_pick_center_motion(self, event):
            if not (hasattr(self, 'pick_center_mode') and self.pick_center_mode):
                # Remove info text if present
                if hasattr(self, 'pick_center_info_text') and self.pick_center_info_text is not None:
                    self.pick_center_info_text.set_visible(False)
                    self.pick_center_info_text = None
                    self.canvas.draw_idle()
                return
            if event.inaxes != self.ax:
                return
            cur_lat, cur_lon = event.ydata, event.xdata
            if cur_lat is None or cur_lon is None:
                # Remove info text if present
                if hasattr(self, 'pick_center_info_text') and self.pick_center_info_text is not None:
                    self.pick_center_info_text.set_visible(False)
                    self.pick_center_info_text = None
                    self.canvas.draw_idle()
                return
            # Get depth and slope at this point
            elevation, slope = self._calculate_slope_at_point(cur_lat, cur_lon)
            try:
                depth_str = f"{abs(elevation):.1f} m" if elevation is not None and not np.isnan(elevation) else "-"
                slope_str = f"{slope:.1f}°" if slope is not None and not np.isnan(slope) else "-"
            except Exception:
                depth_str = "-"
                slope_str = "-"
            info_str = f"Lat: {cur_lat:.6f}\nLon: {cur_lon:.6f}\nDepth: {depth_str}\nSlope: {slope_str}"
            if hasattr(self, 'pick_center_info_text') and self.pick_center_info_text is not None:
                self.pick_center_info_text.set_text(info_str)
                self.pick_center_info_text.set_visible(True)
            else:
                self.pick_center_info_text = self.ax.text(0.02, 0.98, info_str, transform=self.ax.transAxes, fontsize=9, va='top', ha='left', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7), zorder=10)
            self.canvas.draw_idle()

        # Add mouse motion event for general info over GeoTIFF
        def _on_geotiff_hover_motion(self, event):
            # This method is now replaced by _on_mouse_motion
            pass

        # Patch methods into self
        self._toggle_line_planning_mode = _toggle_line_planning_mode.__get__(self)
        self._clear_line_planning = _clear_line_planning.__get__(self)
        self._on_plot_click = new_on_plot_click.__get__(self)
        self._export_drawn_line = _export_drawn_line.__get__(self)
        self._import_drawn_line = _import_drawn_line.__get__(self)
        self._on_line_planning_motion = _on_line_planning_motion.__get__(self)
        self._on_pitch_line_motion = _on_pitch_line_motion.__get__(self)
        self._on_pick_center_motion = _on_pick_center_motion.__get__(self)
        self._on_geotiff_hover_motion = _on_geotiff_hover_motion.__get__(self)

        self._create_widgets()

        # Add info panel for profile plot (initialize before drawing profile)
        self.profile_info_text = None

        # Create the profile plot widgets BEFORE layout
        self.profile_fig = Figure(figsize=(3, 2.2))
        self.profile_ax = self.profile_fig.add_subplot(111)
        self.profile_canvas = FigureCanvas(self.profile_fig)
        self.profile_widget = self.profile_canvas
        self._draw_crossline_profile()

        # Slope profile checkbox for profile plot
        self.show_slope_profile_var = True
        self.slope_profile_checkbox = QCheckBox("Show Slope Profile")
        self.slope_profile_checkbox.setChecked(self.show_slope_profile_var)
        self.slope_profile_checkbox.stateChanged.connect(self._draw_current_profile)

        self._setup_layout()

        # Connect Matplotlib click event for 'Pick Center from GeoTIFF'
        self.cid_click = self.canvas.mpl_connect('button_press_event', self._on_plot_click)
        # Connect Matplotlib scroll event for zoom
        self.cid_scroll = self.canvas.mpl_connect('scroll_event', self._on_scroll)
        # Connect Matplotlib motion event for real-time info display
        self.cid_motion = self.canvas.mpl_connect('motion_notify_event', self._on_mouse_motion)
        # Connect middle mouse button for panning
        self.cid_middle_press = self.canvas.mpl_connect('button_press_event', self._on_middle_press)
        self.cid_middle_motion = self.canvas.mpl_connect('motion_notify_event', self._on_middle_motion)
        self.cid_middle_release = self.canvas.mpl_connect('button_release_event', self._on_middle_release)
        # Connect draw event for dynamic colormap scaling
        self.cid_draw = self.canvas.mpl_connect('draw_event', self._on_draw_event_update_colormap)

        # Initialize buttons states based on library availability
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Missing Libraries",
                                   f"Some geospatial libraries are not found. "
                                   "GeoTIFF and Shapefile export features will be disabled. "
                                   "Install with: pip install rasterio pyproj fiona shapely")
            self.load_geotiff_btn.setEnabled(False)
            # Restored: Enable/disable for elevation_slope_combo
            if hasattr(self, 'elevation_slope_combo'):
                self.elevation_slope_combo.setEnabled(False)
            if hasattr(self, 'pick_center_btn'):
                self.pick_center_btn.setEnabled(False)
            if hasattr(self, 'remove_geotiff_btn'):
                self.remove_geotiff_btn.setEnabled(False)  # Disable Remove GeoTIFF button

        # Bind parameter changes to update the plot
        self._updating_from_code = False  # Flag to prevent recursion
        param_entries = [
            self.central_lat_entry,
            self.central_lon_entry,
            self.line_length_entry,
            self.heading_entry,
            self.dist_between_lines_entry,
            self.num_lines_entry,
            self.bisect_lead_entry,
            self.survey_speed_entry
        ]
        for entry in param_entries:
            entry.textChanged.connect(self._on_parameter_change)

        self.last_used_dir = os.path.expanduser("~")  # Default to user's home directory
        self.last_geotiff_dir = os.path.expanduser("~")  # Default to user's home directory
        self.last_survey_params_dir = os.path.expanduser("~")  # Default to user's home directory
        self.last_export_dir = os.path.expanduser("~")  # Default to user's home directory
        # Separate import directories for each survey type
        self.last_cal_import_dir = os.path.expanduser("~")
        self.last_ref_import_dir = os.path.expanduser("~")
        self.last_line_import_dir = os.path.expanduser("~")
        self._load_last_used_dir()
        self._load_last_geotiff_dir()
        self._load_last_survey_params_dir()
        self._load_last_export_dir()
        self._load_last_cal_import_dir()
        self._load_last_ref_import_dir()
        self._load_last_line_import_dir()

        # After self._setup_layout(), connect the motion event for line planning, pitch line, pick center, and geotiff hover
        # Note: These are now handled by the main _on_mouse_motion method

        # Bind parameter changes to update the Export Name
        self.central_lat_entry.textChanged.connect(self._update_export_name)
        # Reset "Pick Center from GeoTIFF" button to normal style when Central Latitude is entered
        self.central_lat_entry.textChanged.connect(lambda: self.pick_center_btn.setStyleSheet("") if hasattr(self, 'pick_center_btn') and self.central_lat_entry.text().strip() else None)
        self.central_lon_entry.textChanged.connect(self._update_export_name)
        self.line_length_entry.textChanged.connect(self._update_export_name)
        self.heading_entry.textChanged.connect(self._update_export_name)
        self.dist_between_lines_entry.textChanged.connect(self._update_export_name)
        self.num_lines_entry.textChanged.connect(self._update_export_name)
        self.bisect_lead_entry.textChanged.connect(self._update_export_name)
        self.survey_speed_entry.textChanged.connect(self._update_export_name)
        self.offset_direction_combo.currentTextChanged.connect(self._update_export_name)
        # Slider connections already set up in widget creation
        
        # Show window maximized after everything is set up
        # Set window to minimum size on startup
        self.resize(1600, 1150)
        print("Initialization complete, window should be visible")

    def set_cal_info_text(self, message, append=False):
        """Add a message to the Activity Log with [Calibration] prefix. Maintains up to 200 lines of history."""
        if not hasattr(self, 'activity_log_text'):
            return
        prefixed_message = f"[Calibration] {message}"
        self.activity_log_text.setReadOnly(False)
        if append:
            # Prepend new message (newest at top) - only when explicitly requested
            current_text = self.activity_log_text.toPlainText().rstrip('\n')
            new_text = prefixed_message + "\n" + current_text if current_text else prefixed_message + "\n"
            # Split into lines and limit to 200
            lines = new_text.splitlines()
            if len(lines) > 200:
                lines = lines[:200]
            self.activity_log_text.setPlainText("\n".join(lines) + "\n")
            # Move cursor to top
            cursor = self.activity_log_text.textCursor()
            cursor.movePosition(QTextCursor.MoveOperation.Start)
            self.activity_log_text.setTextCursor(cursor)
        else:
            # Append to end (default behavior - newest at bottom)
            self.activity_log_text.append(prefixed_message)
            # Maintain up to 200 lines by removing oldest from top
            num_lines = self.activity_log_text.document().blockCount()
            if num_lines > 200:
                cursor = self.activity_log_text.textCursor()
                cursor.movePosition(QTextCursor.MoveOperation.Start)
                cursor.movePosition(QTextCursor.MoveOperation.Down, QTextCursor.MoveMode.MoveAnchor, num_lines - 200)
                cursor.movePosition(QTextCursor.MoveOperation.Start, QTextCursor.MoveMode.KeepAnchor)
                cursor.removeSelectedText()
            # Move cursor to end and scroll to show newest message
            cursor = self.activity_log_text.textCursor()
            cursor.movePosition(QTextCursor.MoveOperation.End)
            self.activity_log_text.setTextCursor(cursor)
        self.activity_log_text.ensureCursorVisible()
        self.activity_log_text.setReadOnly(True)

    def _update_multiplier_label_len(self, val):
        """Updates the label next to the line length multiplier slider."""
        self.multiplier_label_len.setText(f"{float(val):.1f}")

    def _update_multiplier_label_dist(self, val):
        """Updates the label next to the distance between lines multiplier slider."""
        self.multiplier_label_dist.setText(f"{float(val):.1f}")

    def _setup_layout(self):
        try:
            # Create vertical layout for right side (map + profile)
            right_side_layout = QVBoxLayout()
            right_side_layout.setContentsMargins(0, 0, 0, 0)
            
            # Add map plot to right side layout
            if hasattr(self, 'plot_frame'):
                right_side_layout.addWidget(self.plot_frame, 1)  # Stretch factor 1
            else:
                print(f"ERROR: Missing widget - plot_frame: {hasattr(self, 'plot_frame')}")
                return
            
            # Add profile plot at bottom of right side
            if hasattr(self, 'profile_widget') and hasattr(self, 'slope_profile_checkbox'):
                profile_layout = QVBoxLayout()
                profile_layout.setContentsMargins(0, 0, 0, 0)
                profile_layout.addWidget(self.profile_widget)
                profile_layout.addWidget(self.slope_profile_checkbox)
                profile_widget = QWidget()
                profile_widget.setLayout(profile_layout)
                profile_widget.setMaximumHeight(250)  # Limit profile height
                right_side_layout.addWidget(profile_widget)
            else:
                print(f"WARNING: Missing profile widgets - profile_widget: {hasattr(self, 'profile_widget')}, slope_profile_checkbox: {hasattr(self, 'slope_profile_checkbox')}")
            
            # Create widget to hold right side layout
            right_side_widget = QWidget()
            right_side_widget.setLayout(right_side_layout)
            
            # Add horizontal layout for params (left, full height) and right side (map + profile)
            if hasattr(self, 'param_scroll'):
                self.content_layout.addWidget(self.param_scroll)
                self.content_layout.addWidget(right_side_widget, 1)  # Stretch factor 1
            else:
                print(f"ERROR: Missing widget - param_scroll: {hasattr(self, 'param_scroll')}")
                return
            
            # Add content layout to main layout
            self.main_layout.addLayout(self.content_layout, 1)  # stretch factor 1
            
            # Connect signals for line_length_entry and survey_speed_entry
            if hasattr(self, 'line_length_entry') and hasattr(self, 'survey_speed_entry'):
                self.line_length_entry.textChanged.connect(self._on_line_length_or_speed_change)
                self.survey_speed_entry.textChanged.connect(self._on_line_length_or_speed_change)
            
            # Force update and show
            self.central_widget.update()
            self.central_widget.show()
            print(f"Layout setup completed - param_scroll: {self.param_scroll is not None}, plot_frame: {self.plot_frame is not None}")
            print(f"Main layout children: {self.main_layout.count()}, Content layout children: {self.content_layout.count()}")
        except Exception as e:
            print(f"ERROR in _setup_layout: {e}")
            import traceback
            traceback.print_exc()

    def _validate_inputs(self):
        """Validates all input fields and returns a tuple (is_valid, values_dict)."""
        values = {}
        try:
            # Central Latitude
            central_lat_str = self.central_lat_entry.text().strip()
            if not central_lat_str:
                self._show_message("error","Input Error", "Central Latitude cannot be empty.")
                return False, {}
            values['central_lat'] = float(central_lat_str)
            if not (-90 <= values['central_lat'] <= 90):
                self._show_message("error","Input Error", "Central Latitude must be between -90 and 90.")
                return False, {}

            # Central Longitude
            central_lon_str = self.central_lon_entry.text().strip()
            if not central_lon_str:
                self._show_message("error","Input Error", "Central Longitude cannot be empty.")
                return False, {}
            values['central_lon'] = float(central_lon_str)
            if not (-180 <= values['central_lon'] <= 180):
                self._show_message("error","Input Error", "Central Longitude must be between -180 and 180.")
                return False, {}

            # Line Length
            line_length_str = self.line_length_entry.text().strip()
            if not line_length_str:
                self._show_message("error","Input Error", "Line Length cannot be empty.")
                return False, {}
            values['line_length'] = float(line_length_str)
            if not (values['line_length'] > 0):
                self._show_message("error","Input Error", "Line Length must be positive.")
                return False, {}

            # Heading
            heading_str = self.heading_entry.text().strip()
            if not heading_str:
                self._show_message("error","Input Error", "Heading cannot be empty.")
                return False, {}
            values['heading'] = float(heading_str)
            if not (0 <= values['heading'] <= 360):
                self._show_message("error","Input Error", "Heading must be between 0 and 360 degrees.")
                return False, {}

            # Distance Between Lines
            dist_between_lines_str = self.dist_between_lines_entry.text().strip()
            if not dist_between_lines_str:
                self._show_message("error","Input Error", "Distance Between Lines cannot be empty.")
                return False, {}
            values['dist_between_lines'] = float(dist_between_lines_str)
            if not (values['dist_between_lines'] > 0):
                self._show_message("error","Input Error", "Distance Between Lines must be positive.")
                return False, {}

            # Number of Lines
            num_lines_str = self.num_lines_entry.text().strip()
            if not num_lines_str:
                self._show_message("error","Input Error", "Number of Lines cannot be empty.")
                return False, {}
            values['num_lines'] = int(num_lines_str)
            if not (values['num_lines'] > 0):
                self._show_message("error","Input Error", "Number of Lines must be a positive integer.")
                return False, {}

            # Crossline Lead-in/out
            bisect_lead_str = self.bisect_lead_entry.text().strip()
            if not bisect_lead_str:
                self._show_message("error","Input Error", "Crossline Lead-in/out cannot be empty.")
                return False, {}
            values['bisect_lead'] = float(bisect_lead_str)
            if not (values['bisect_lead'] >= 0):
                self._show_message("error","Input Error", "Crossline Lead-in/out cannot be negative.")
                return False, {}

            # Survey Speed
            survey_speed_str = self.survey_speed_entry.text().strip()
            if not survey_speed_str:
                self._show_message("error","Input Error", "Survey Speed cannot be empty.")
                return False, {}
            values['survey_speed'] = float(survey_speed_str)
            if not (values['survey_speed'] > 0):
                self._show_message("error","Input Error", "Survey Speed must be positive.")
                return False, {}

            # Export Name
            values['export_name'] = self.export_name_entry.text().strip()
            if not values['export_name']:
                self._show_message("error","Input Error", "Export Name cannot be empty.")
                return False, {}

            values['offset_direction'] = self.offset_direction_var
            values['line_length_multiplier'] = self.line_length_multiplier  # Get slider value
            values['dist_between_lines_multiplier'] = self.dist_between_lines_multiplier  # Get new slider value

        except ValueError as e:
            self._show_message("error","Input Error", f"Invalid numeric input: {str(e)}")
            return False, {}
        except Exception as e:
            self._show_message("error","Error", f"An unexpected error occurred during input validation: {str(e)}")
            return False, {}

        return True, values

    def _on_parameter_changed(self):
        """Handle parameter changes by restarting the auto-regenerate timer."""
        # Stop any existing timer and start a new one
        self.auto_regenerate_timer.stop()
        self.auto_regenerate_timer.start(800)  # 800ms debounce delay
    
    def _auto_regenerate_survey_plan(self):
        """Auto-regenerate survey plan when timer expires, if survey plan already exists."""
        # Only regenerate if a survey plan already exists (has survey lines)
        if hasattr(self, 'survey_lines_data') and len(self.survey_lines_data) > 0:
            try:
                # Validate inputs before regenerating
                is_valid, values = self._validate_inputs()
                if is_valid:
                    self._generate_and_plot(show_success_dialog=False)
            except Exception as e:
                # Silently fail - don't show error dialog for auto-regeneration
                print(f"Auto-regenerate failed: {e}")

    def _generate_and_plot(self, show_success_dialog=True):
        # Clear the elevation profile immediately when generating a new plan
        self._draw_crossline_profile()
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot generate plot.")
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
            self._show_message("error","Projection Error", error_msg)
        except Exception as e:
            error_msg = f"Failed to generate survey lines: {str(e)}"
            print(f"ERROR: {error_msg}")
            print("Full traceback:")
            traceback.print_exc()
            self._show_message("error","Calculation Error", error_msg)

        # In _generate_and_plot and _load_geotiff, after plotting/plan generation, call:
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

            # Plot GeoTIFF if loaded
            if self.geotiff_data_array is not None and self.geotiff_extent is not None:
                if self.hillshade_vs_slope_viz_mode == "hillshade":
                    # Use multidirectional hillshade: average of 4 azimuths
                    masked_data = np.ma.array(self.geotiff_data_array, mask=np.isnan(self.geotiff_data_array))
                    data_for_hillshade = np.nan_to_num(masked_data.filled(np.nanmin(self.geotiff_data_array)))
                    azimuths = [45, 135, 225, 315]
                    altitude = 45
                    hillshades = []
                    for az in azimuths:
                        ls = LightSource(azdeg=az, altdeg=altitude)
                        hillshade = ls.hillshade(data_for_hillshade,vert_exag=0.1)
                        hillshades.append(hillshade)
                    # Average the hillshades
                    hillshade_array = np.mean(hillshades, axis=0)
                    self.geotiff_hillshade_plot = self.ax.imshow(hillshade_array, extent=tuple(self.geotiff_extent),
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
                                                                       orientation='vertical', label='Elevation (m)')

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
                        ls = LightSource(azdeg=315, altdeg=45)
                        hillshade = ls.hillshade(temp_data)
                        self.ax.imshow(hillshade, extent=tuple(self.geotiff_extent), cmap='gray', origin='upper', alpha=0.5, zorder=-2)

                        cmap = 'plasma'
                        self.geotiff_image_plot = self.ax.imshow(display_data, extent=tuple(self.geotiff_extent),
                                                                 cmap=cmap, origin='upper', alpha=0.5,
                                                                 vmin=0, vmax=max_slope_for_cmap,
                                                                 zorder=-1)

                        self.slope_colorbar = self.figure.colorbar(self.geotiff_image_plot, ax=self.ax,
                                                                   orientation='vertical', label='Slope (degrees)')

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
                                                               orientation='vertical', label='Slope (degrees)')

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
                                                                 zorder=0)  # Above background but below other features
                except Exception as e:
                    # Silently fail if contour plotting fails (e.g., invalid data, etc.)
                    if hasattr(self, 'contour_plot') and self.contour_plot is not None:
                        try:
                            for collection in self.contour_plot.collections:
                                collection.remove()
                        except:
                            pass
                        self.contour_plot = None
            else:
                # Remove contours if checkbox is disabled or no GeoTIFF loaded
                if hasattr(self, 'contour_plot') and self.contour_plot is not None:
                    try:
                        for collection in self.contour_plot.collections:
                            collection.remove()
                    except:
                        pass
                    self.contour_plot = None

            # Plot main survey lines
            for i, line in enumerate(self.survey_lines_data):
                latitudes = [p[0] for p in line]
                longitudes = [p[1] for p in line]
                label = "Reference Line" if i == 0 else "_nolegend_"
                self.ax.plot(longitudes, latitudes, color='blue', linewidth=1.5, label=label)
                
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
                self.ax.plot(longitudes, latitudes, color='red', linestyle='--', linewidth=1.5,
                             label='Crossline')
                
                # Add labels for crossline points
                self.ax.annotate('CLS', (longitudes[0], latitudes[0]), 
                                xytext=(5, 5), textcoords='offset points', 
                                fontsize=8, color='red', weight='bold',
                                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
                self.ax.annotate('CLE', (longitudes[1], latitudes[1]), 
                                xytext=(5, 5), textcoords='offset points', 
                                fontsize=8, color='red', weight='bold',
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
                    self.ax.plot(longitudes, latitudes, color='orange', linewidth=2, linestyle='--', marker='x', label=label)
                    
                    # Add labels for heading line points
                    prefix = 'H1' if i == 0 else 'H2'
                    self.ax.annotate(f'{prefix}S', (longitudes[0], latitudes[0]), 
                                    xytext=(5, 5), textcoords='offset points', 
                                    fontsize=8, color='orange', weight='bold',
                                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
                    self.ax.annotate(f'{prefix}E', (longitudes[1], latitudes[1]), 
                                    xytext=(5, 5), textcoords='offset points', 
                                    fontsize=8, color='orange', weight='bold',
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

            self.canvas.draw_idle()

        except Exception as e:
            print(f"Error in _plot_survey_plan: {str(e)}")
            traceback.print_exc()
            self._show_message("error","Plot Error", f"Failed to generate plot: {str(e)}")
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
            self.geotiff_image_plot = None
            self.geotiff_hillshade_plot = None
            # Robustly remove colorbars
            self._remove_colorbar('slope_colorbar')
            self._remove_colorbar('elevation_colorbar')

            # Restored: geotiff_display_mode reset
            self.geotiff_display_mode = "elevation"
            # Update combo box if it exists
            if hasattr(self, 'elevation_slope_combo'):
                self._update_display_mode_combo()

            self.hillshade_vs_slope_viz_mode = "hillshade"
            # Remove all .config() calls for self.slope_viz_btn
            # Remove: self.slope_viz_btn.config(text="Display Slope Visualization")
            # Remove: self.slope_viz_btn.config(text="Display Shaded Relief")

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
            
            # Convert plot limits to dataset coordinates
            if self.geotiff_dataset_original.crs != "EPSG:4326":
                transformer = pyproj.Transformer.from_crs("EPSG:4326", self.geotiff_dataset_original.crs, always_xy=True)
                plot_left, plot_bottom = transformer.transform(xlim[0], ylim[0])
                plot_right, plot_top = transformer.transform(xlim[1], ylim[1])
            else:
                plot_left, plot_bottom = xlim[0], ylim[0]
                plot_right, plot_top = xlim[1], ylim[1]
            
            # Calculate pixel indices for the visible region
            transform = self.geotiff_dataset_original.transform
            row_min, col_min = rowcol(transform, plot_left, plot_top)
            row_max, col_max = rowcol(transform, plot_right, plot_bottom)
            
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
            
            # Add padding around the visible region, but ensure we don't exceed dataset bounds
            padding = max(50, min(100, (row_max - row_min) // 4))  # Adaptive padding
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
            
            # Filter Z-values
            data[data < -11000] = np.nan
            data[data > 0] = np.nan
            
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

    def _load_geotiff(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot load GeoTIFF.")
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
                self.wgs84_to_geotiff_transformer = None
                self.geotiff_to_wgs84_transformer = None

            # Filter Z-values: values < -11000 or > 0 are set to NaN
            progress_label.setText("Processing elevation data...")
            QApplication.processEvents()
            
            self.geotiff_data_array[self.geotiff_data_array < -11000] = np.nan
            self.geotiff_data_array[self.geotiff_data_array > 0] = np.nan

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
            self._zoom_to_geotiff()

        except RasterioIOError as e:
            progress_window.destroy()
            self._show_message("error","GeoTIFF Error", f"Could not open GeoTIFF file: {e}")
            self._clear_plot(full_clear=True)
        except Exception as e:
            progress_window.destroy()
            self._show_message("error","GeoTIFF Error", f"An unexpected error occurred while loading GeoTIFF: {e}")
            self._clear_plot(full_clear=True)

        # Draw crossline profile
        self._draw_crossline_profile()

    def _on_geotiff_display_mode_changed(self, mode_text):
        """Handle display mode selection from combo box."""
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot change display mode.")
            # Reset to previous mode
            if hasattr(self, 'elevation_slope_combo'):
                self._update_display_mode_combo()
            return

        if self.geotiff_data_array is None:
            self._show_message("warning","No GeoTIFF", "Load a GeoTIFF first to change display mode.")
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
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot show contours.")
            self.show_contours_var = False
            if hasattr(self, 'show_contours_checkbox'):
                self.show_contours_checkbox.setChecked(False)
            return

        if self.geotiff_data_array is None:
            self._show_message("warning","No GeoTIFF", "Load a GeoTIFF first to show contours.")
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
            except:
                pass
            finally:
                self._syncing_contour_interval = False
        except:
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
            self._show_message("warning","Disabled Feature",
                                   "Geospatial libraries not loaded. Cannot toggle slope visualization.")
            return

        if self.geotiff_data_array is None:
            self._show_message("warning","No GeoTIFF", "Load a GeoTIFF first to toggle slope visualization.")
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

    def _toggle_pick_pitch_line_mode(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot pick pitch line.")
            return
        if self.geotiff_dataset_original is None:
            self._show_message("warning","No GeoTIFF", "Load a GeoTIFF first to pick a pitch line.")
            return
        self.pick_pitch_line_mode = not self.pick_pitch_line_mode
        if self.pick_pitch_line_mode:
            self.pitch_line_points = []
            self.pick_pitch_line_btn.setText("Drawing Pitch Line: Click Start Point")
            if hasattr(self, 'calibration_frame'):
                self.calibration_frame.setFocus()
            self.canvas_widget.setCursor(Qt.CursorShape.CrossCursor)
            self.set_cal_info_text("Click the start point, and then the end point of the pitch line on the plot.")
            # Add orange message to Activity Log
            if hasattr(self, 'activity_log_text'):
                self.activity_log_text.setReadOnly(False)
                cursor = self.activity_log_text.textCursor()
                cursor.movePosition(QTextCursor.MoveOperation.End)
                # Set orange color
                char_format = QTextCharFormat()
                char_format.setForeground(QColor(255, 165, 0))  # Orange color
                cursor.setCharFormat(char_format)
                cursor.insertText("[Calibration] Left click the starting and ending points of the pitch line on the plot.\n")
                # Reset to default color
                char_format.setForeground(QColor(0, 0, 0))  # Black color
                cursor.setCharFormat(char_format)
                self.activity_log_text.setTextCursor(cursor)
                self.activity_log_text.setReadOnly(True)
            # Remove info text if present
            if hasattr(self, 'pitch_line_info_text') and self.pitch_line_info_text is not None:
                self.pitch_line_info_text.set_visible(False)
                self.pitch_line_info_text = None
                self.canvas.draw_idle()
        else:
            self.pick_pitch_line_btn.setText("Draw a Pitch Line")
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
            # Remove info text if present
            if hasattr(self, 'pitch_line_info_text') and self.pitch_line_info_text is not None:
                self.pitch_line_info_text.set_visible(False)
                self.pitch_line_info_text = None
                self.canvas.draw_idle()
        self._update_cal_line_times()

    def _toggle_edit_pitch_line_mode(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot edit pitch line.")
            return
        if self.geotiff_dataset_original is None:
            self._show_message("warning","No GeoTIFF", "Load a GeoTIFF first to edit a pitch line.")
            return
        if len(self.pitch_line_points) < 2:
            self._show_message("warning","No Pitch Line", "Draw a pitch line first before editing it.")
            return
        
        self.edit_pitch_line_mode = not self.edit_pitch_line_mode
        if self.edit_pitch_line_mode:
            # Clear any existing heading lines as requested
            self.heading_lines = []
            self._plot_survey_plan()  # Redraw to remove heading lines
            
            # Create draggable handles for start and end points
            start_lat, start_lon = self.pitch_line_points[0]
            end_lat, end_lon = self.pitch_line_points[1]
            
            # Create handles as scatter points
            self.pitch_line_start_handle = self.ax.scatter([start_lon], [start_lat], 
                                                         color='red', s=100, marker='o', 
                                                         edgecolors='black', linewidth=2, 
                                                         zorder=10, picker=5)
            self.pitch_line_end_handle = self.ax.scatter([end_lon], [end_lat], 
                                                       color='blue', s=100, marker='o', 
                                                       edgecolors='black', linewidth=2, 
                                                       zorder=10, picker=5)
            
            # Connect pick events
            self.pitch_line_pick_cid = self.canvas.mpl_connect('pick_event', self._on_pitch_line_handle_pick)
            self.pitch_line_motion_cid = self.canvas.mpl_connect('motion_notify_event', self._on_pitch_line_handle_motion)
            self.pitch_line_release_cid = self.canvas.mpl_connect('button_release_event', self._on_pitch_line_handle_release)
            
            self.edit_pitch_line_btn.setText("Click to Stop Editing Pitch Line")
            self.canvas_widget.setCursor(Qt.CursorShape.SizeAllCursor)
            self.set_cal_info_text("You can drag the red (start) and blue (end) points to edit the pitch line. Heading lines have been cleared.")
        else:
            # Remove handles by clearing the axes and redrawing
            if self.pitch_line_start_handle or self.pitch_line_end_handle:
                # Clear the handles by redrawing the plot
                self._plot_survey_plan()
                self.pitch_line_start_handle = None
                self.pitch_line_end_handle = None
            
            # Disconnect events
            if hasattr(self, 'pitch_line_pick_cid'):
                self.canvas.mpl_disconnect(self.pitch_line_pick_cid)
            if hasattr(self, 'pitch_line_motion_cid'):
                self.canvas.mpl_disconnect(self.pitch_line_motion_cid)
            if hasattr(self, 'pitch_line_release_cid'):
                self.canvas.mpl_disconnect(self.pitch_line_release_cid)
            
            self.edit_pitch_line_btn.setText("Edit Pitch Line")
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
            self.dragging_pitch_line_handle = None
            
            # Redraw to update the pitch line
            self._plot_survey_plan()
            self._update_cal_line_offset_from_pitch_line()
            self._draw_pitch_line_profile()
            self._update_cal_export_name_from_pitch_line()
            self._update_cal_line_times()
            
            self.set_cal_info_text("Pitch line editing completed.")

    def _toggle_pick_roll_line_mode(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot pick roll line.")
            return
        if self.geotiff_dataset_original is None:
            self._show_message("warning","No GeoTIFF", "Load a GeoTIFF first to pick a roll line.")
            return
        self.pick_roll_line_mode = not self.pick_roll_line_mode
        if self.pick_roll_line_mode:
            self.roll_line_points = []
            self.pick_roll_line_btn.setText("Drawing Roll Line: Click Start Point")
            if hasattr(self, 'calibration_frame'):
                self.calibration_frame.setFocus()
            self.canvas_widget.setCursor(Qt.CursorShape.CrossCursor)
            self.set_cal_info_text("Click the start point, and then the end point of the roll line on the plot.")
            # Add orange message to Activity Log
            if hasattr(self, 'activity_log_text'):
                self.activity_log_text.setReadOnly(False)
                cursor = self.activity_log_text.textCursor()
                cursor.movePosition(QTextCursor.MoveOperation.End)
                # Set orange color
                char_format = QTextCharFormat()
                char_format.setForeground(QColor(255, 165, 0))  # Orange color
                cursor.setCharFormat(char_format)
                cursor.insertText("[Calibration] Left click the starting and ending points of the roll line on the plot.\n")
                # Reset to default color
                char_format.setForeground(QColor(0, 0, 0))  # Black color
                cursor.setCharFormat(char_format)
                self.activity_log_text.setTextCursor(cursor)
                self.activity_log_text.setReadOnly(True)
            # self._show_message("info","Draw a Roll Line", "Click the start point, then the end point of the roll line on the plot.")
        else:
            self.pick_roll_line_btn.setText("Draw a Roll Line")
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
        self._update_cal_line_times()

    def _toggle_edit_roll_line_mode(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot edit roll line.")
            return
        if self.geotiff_dataset_original is None:
            self._show_message("warning","No GeoTIFF", "Load a GeoTIFF first to edit a roll line.")
            return
        if len(self.roll_line_points) < 2:
            self._show_message("warning","No Roll Line", "Draw a roll line first before editing it.")
            return
        
        self.edit_roll_line_mode = not self.edit_roll_line_mode
        if self.edit_roll_line_mode:
            # Create draggable handles for the roll line endpoints
            (lat1, lon1), (lat2, lon2) = self.roll_line_points
            self.roll_line_start_handle = self.ax.scatter([lon1], [lat1], c='red', s=100, zorder=10, picker=True)
            self.roll_line_end_handle = self.ax.scatter([lon2], [lat2], c='blue', s=100, zorder=10, picker=True)
            
            # Connect pick events for the handles
            self.roll_line_pick_cid = self.canvas.mpl_connect('pick_event', self._on_roll_line_handle_pick)
            self.roll_line_motion_cid = self.canvas.mpl_connect('motion_notify_event', self._on_roll_line_handle_motion)
            self.roll_line_release_cid = self.canvas.mpl_connect('button_release_event', self._on_roll_line_handle_release)
            
            self.edit_roll_line_btn.setText("Click to Stop Editing Roll Line")
            self.canvas_widget.setCursor(Qt.CursorShape.SizeAllCursor)
            self.set_cal_info_text("You can drag the red (start) and blue (end) points to edit the roll line.")
        else:
            # Remove handles by clearing the axes and redrawing
            if self.roll_line_start_handle or self.roll_line_end_handle:
                # Clear the handles by redrawing the plot
                self._plot_survey_plan()
                self.roll_line_start_handle = None
                self.roll_line_end_handle = None
            
            # Disconnect events
            if hasattr(self, 'roll_line_pick_cid'):
                self.canvas.mpl_disconnect(self.roll_line_pick_cid)
            if hasattr(self, 'roll_line_motion_cid'):
                self.canvas.mpl_disconnect(self.roll_line_motion_cid)
            if hasattr(self, 'roll_line_release_cid'):
                self.canvas.mpl_disconnect(self.roll_line_release_cid)
            
            self.edit_roll_line_btn.setText("Edit Roll Line")
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
            self.dragging_roll_line_handle = None
            
            self._update_cal_line_times()
            
            self.set_cal_info_text("Roll line editing completed.")

    def _on_plot_click(self, event):
        # Ignore middle button (button 2) - let it be handled by pan handlers
        if event.button == 2:
            return
        # --- Debugging prints added ---
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

            self.ax.set_xlim(min_lon - buffer_lon, max_lon + buffer_lon)
            self.ax.set_ylim(min_lat - buffer_lat, max_lat + buffer_lat)
            self.canvas.draw_idle()
        else:
            # If nothing to zoom to, just reset view
            self.ax.autoscale_view()
            self.canvas.draw_idle()

    def _remove_geotiff(self):
        """Remove the loaded GeoTIFF and clear it from the plot."""
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded.")
            return

        if self.geotiff_dataset_original is None:
            self._show_message("warning","No GeoTIFF", "No GeoTIFF is currently loaded.")
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
        """Zooms the plot to the GeoTIFF bounds without clearing existing lines."""
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot zoom to GeoTIFF.")
            return

        if self.geotiff_extent is None:
            self._show_message("warning","No GeoTIFF", "No GeoTIFF underlay is loaded to zoom to.")
            return

        # Calculate consistent plot limits based on GeoTIFF extent
        self.current_xlim, self.current_ylim = self._calculate_consistent_plot_limits()
        
        # Store these as the user's intended view limits
        self._last_user_xlim = self.current_xlim
        self._last_user_ylim = self.current_ylim
        
        # Set the axis limits to zoom to GeoTIFF bounds
        self.ax.set_xlim(self.current_xlim)
        self.ax.set_ylim(self.current_ylim)
        
        # Update zoom level if dynamic resolution is enabled
        if self.dynamic_resolution_enabled and hasattr(self, 'geotiff_extent') and self.geotiff_extent is not None:
            full_width = self.geotiff_extent[1] - self.geotiff_extent[0]
            full_height = self.geotiff_extent[3] - self.geotiff_extent[2]
            current_width = self.current_xlim[1] - self.current_xlim[0]
            current_height = self.current_ylim[1] - self.current_ylim[0]
            
            width_ratio = current_width / full_width if full_width > 0 else 1.0
            height_ratio = current_height / full_height if full_height > 0 else 1.0
            new_zoom_level = min(width_ratio, height_ratio)
            new_zoom_level = max(0.01, min(10.0, new_zoom_level))
            self.geotiff_zoom_level = new_zoom_level
        
        # Reload GeoTIFF data to cover the new extent if needed
        if self.geotiff_dataset_original is not None:
            success = self._load_geotiff_at_resolution()
            if not success:
                self._show_message("warning","GeoTIFF Reload", "Could not reload GeoTIFF data for the new view.")
        
        # Just redraw the canvas without replotting everything
        # This preserves all existing lines from other tabs
        self.canvas.draw_idle()

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

    def _export_survey_data(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot export data.")
            return

        is_valid, values = self._validate_inputs()
        if not is_valid:
            return

        export_name = values['export_name']

        if not self.survey_lines_data and not self.cross_line_data:
            self._show_message("warning","No Data", "No survey lines to export. Generate them first.")
            return

        export_dir = QFileDialog.getExistingDirectory(self, "Select Export Directory", self.last_export_dir)
        if not export_dir:
            return
        self.last_export_dir = export_dir
        self._save_last_export_dir()

        try:
            # --- Export to CSV ---
            # CSV files for decimal degrees always use _DD suffix
            csv_file_path = os.path.join(export_dir, f"{export_name}_DD.csv")
            with open(csv_file_path, 'w', newline='') as csvfile:
                csv_writer = csv.writer(csvfile)
                csv_writer.writerow(['Line Number', 'Point Label', 'Latitude', 'Longitude'])
                # Main survey lines
                for i, line in enumerate(self.survey_lines_data):
                    if i % 2 == 0:  # Even lines - normal order
                        start, end = line[0], line[1]
                    else:  # Odd lines - flipped order
                        start, end = line[1], line[0]
                    start_label = f'L{i+1}S'
                    end_label = f'L{i+1}E'
                    csv_writer.writerow([i + 1, start_label, start[0], start[1]])
                    csv_writer.writerow([i + 1, end_label, end[0], end[1]])
                # Crossline
                if self.cross_line_data:
                    csv_writer.writerow([0, 'CLS', self.cross_line_data[0][0], self.cross_line_data[0][1]])
                    csv_writer.writerow([0, 'CLE', self.cross_line_data[1][0], self.cross_line_data[1][1]])
            
            # --- Export to DDM format (Decimal Minutes) ---
            ddm_file_path = os.path.join(export_dir, f"{export_name}_DM.csv")
            with open(ddm_file_path, 'w', newline='', encoding='utf-8') as ddmfile:
                ddm_writer = csv.writer(ddmfile)
                ddm_writer.writerow(['Line Number', 'Point Label', 'Latitude with Decimal Minutes', 'Longitude with Decimal Minutes'])
                # Main survey lines
                for i, line in enumerate(self.survey_lines_data):
                    if i % 2 == 0:  # Even lines - normal order
                        start, end = line[0], line[1]
                    else:  # Odd lines - flipped order
                        start, end = line[1], line[0]
                    start_label = f'L{i+1}S'
                    end_label = f'L{i+1}E'
                    start_lat_ddm = self._decimal_degrees_to_ddm(start[0], is_latitude=True)
                    start_lon_ddm = self._decimal_degrees_to_ddm(start[1], is_latitude=False)
                    end_lat_ddm = self._decimal_degrees_to_ddm(end[0], is_latitude=True)
                    end_lon_ddm = self._decimal_degrees_to_ddm(end[1], is_latitude=False)
                    ddm_writer.writerow([i + 1, start_label, start_lat_ddm, start_lon_ddm])
                    ddm_writer.writerow([i + 1, end_label, end_lat_ddm, end_lon_ddm])
                # Crossline
                if self.cross_line_data:
                    cls_lat_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[0][0], is_latitude=True)
                    cls_lon_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[0][1], is_latitude=False)
                    cle_lat_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[1][0], is_latitude=True)
                    cle_lon_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[1][1], is_latitude=False)
                    ddm_writer.writerow([0, 'CLS', cls_lat_ddm, cls_lon_ddm])
                    ddm_writer.writerow([0, 'CLE', cle_lat_ddm, cle_lon_ddm])
            
            # --- Export to DDM text format (Decimal Minutes) ---
            ddm_txt_file_path = os.path.join(export_dir, f"{export_name}_DM.txt")
            with open(ddm_txt_file_path, 'w', encoding='utf-8') as ddm_txt_file:
                # Main survey lines
                for i, line in enumerate(self.survey_lines_data):
                    if i % 2 == 0:  # Even lines - normal order
                        start, end = line[0], line[1]
                    else:  # Odd lines - flipped order
                        start, end = line[1], line[0]
                    start_label = f'L{i+1}S'
                    end_label = f'L{i+1}E'
                    start_lat_ddm = self._decimal_degrees_to_ddm(start[0], is_latitude=True)
                    start_lon_ddm = self._decimal_degrees_to_ddm(start[1], is_latitude=False)
                    end_lat_ddm = self._decimal_degrees_to_ddm(end[0], is_latitude=True)
                    end_lon_ddm = self._decimal_degrees_to_ddm(end[1], is_latitude=False)
                    ddm_txt_file.write(f"{start_label}, {start_lat_ddm}, {start_lon_ddm}\n")
                    ddm_txt_file.write(f"{end_label}, {end_lat_ddm}, {end_lon_ddm}\n")
                # Crossline
                if self.cross_line_data:
                    cls_lat_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[0][0], is_latitude=True)
                    cls_lon_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[0][1], is_latitude=False)
                    cle_lat_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[1][0], is_latitude=True)
                    cle_lon_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[1][1], is_latitude=False)
                    ddm_txt_file.write(f"CLS, {cls_lat_ddm}, {cls_lon_ddm}\n")
                    ddm_txt_file.write(f"CLE, {cle_lat_ddm}, {cle_lon_ddm}\n")

            # --- Export Input Parameters to TXT ---
            params_file_path = os.path.join(export_dir, f"{export_name}_params_DD.txt")
            with open(params_file_path, 'w') as f:
                f.write("Survey Plan Input Parameters:\n\n")
                for key, value in values.items():
                    f.write(f"{key.replace('_', ' ').title()}: {value}\n")

            # --- Export to ESRI Shapefile (.shp) ---
            schema = {
                'geometry': 'LineString',
                'properties': {'line_num': 'int'},
            }
            crs_epsg = 'EPSG:4326'  # WGS 84

            features = []
            # Add main survey lines
            for i, line_coords in enumerate(self.survey_lines_data):
                # Shapely expects (lon, lat) order
                shapely_line = LineString([(p[1], p[0]) for p in line_coords])
                features.append({
                    'geometry': shapely_line.__geo_interface__,
                    'properties': {'line_num': i + 1},
                })

            # Add bisecting line
            if self.cross_line_data:
                shapely_cross_line = LineString([(p[1], p[0]) for p in self.cross_line_data])
                features.append({
                    'geometry': shapely_cross_line.__geo_interface__,
                    'properties': {'line_num': 0},  # Using 0 for crossline
                })

            shapefile_path = os.path.join(export_dir, f"{export_name}.shp")
            # Fiona needs a directory, not a specific file name for the collection
            # Fiona creates multiple files (.shp, .shx, .dbf, .prj)
            with fiona.open(shapefile_path, 'w', driver='ESRI Shapefile', crs=crs_epsg, schema=schema) as collection:
                collection.writerecords(features)

            self.set_ref_info_text(
                f"Data exported successfully to:\n"
                f"- {os.path.basename(csv_file_path)}\n"
                f"- {os.path.basename(params_file_path)}\n"
                f"- {os.path.basename(shapefile_path)} (and associated files)\n"
                f"in directory: {export_dir}", append=False)

        except Exception as e:
            self._show_message("error","Export Error", f"Failed to export data: {e}")

    def _quit_app(self):
        """Safely quit the application."""
        if self._ask_ok_cancel("Quit", "Do you want to quit the application?"):
            # Close any open GeoTIFF dataset
            if self.geotiff_dataset_original:
                self.geotiff_dataset_original.close()
            # Clean up matplotlib figure
            if hasattr(self, 'figure'):
                plt.close(self.figure)
            # Destroy the main window
            self.destroy()

    def _save_survey_parameters(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot save parameters.")
            return

        is_valid, values = self._validate_inputs()
        if not is_valid:
            return

        save_dir = QFileDialog.getExistingDirectory(self, "Select Directory to Save Survey Parameters", self.last_survey_params_dir)
        if not save_dir:
            return
        self.last_survey_params_dir = save_dir
        self._save_last_survey_params_dir()

        # Use the export name from the form, with .json extension
        export_name = values['export_name']
        if not export_name:
            # Fallback to default naming if export name is empty
            dist_between_lines = int(values['dist_between_lines'])
            heading_int = int(values['heading'])
            export_name = f"Reference_{dist_between_lines}m_{heading_int}deg"
        filename = f"{export_name}_params.json"
        file_path = os.path.join(save_dir, filename)

        try:
            # Create parameters dictionary
            params = {
                'central_lat': float(self.central_lat_entry.text()),
                'central_lon': float(self.central_lon_entry.text()),
                'line_length': float(self.line_length_entry.text()),
                'heading': float(self.heading_entry.text()),
                'dist_between_lines': float(self.dist_between_lines_entry.text()),
                'num_lines': int(self.num_lines_entry.text()),
                'bisect_lead': float(self.bisect_lead_entry.text()),
                'survey_speed': float(self.survey_speed_entry.text()),
                'export_name': self.export_name_entry.text().strip(),
                'offset_direction': self.offset_direction_var,
                'line_length_multiplier': self.line_length_multiplier,
                    'dist_between_lines_multiplier': self.dist_between_lines_multiplier
            }

            # Save to JSON file
            with open(file_path, 'w') as f:
                json.dump(params, f, indent=4)

            self.set_ref_info_text(f"Survey parameters saved successfully to:\n{file_path}", append=False)

        except Exception as e:
            self._show_message("error","Save Error", f"Failed to save survey parameters: {e}")

    def _load_survey_parameters(self, file_path):
        """Load survey parameters from a JSON file."""
        try:
            with open(file_path, 'r') as f:
                params = json.load(f)

            # Update all input fields
            self.central_lat_entry.clear()
            self.central_lat_entry.setText(str(params['central_lat']))

            self.central_lon_entry.clear()
            self.central_lon_entry.setText(str(params['central_lon']))

            self.line_length_entry.setText(str(params['line_length']))

            self.heading_entry.setText(str(params['heading']))

            self.dist_between_lines_entry.setText(str(params['dist_between_lines']))

            self.num_lines_entry.setText(str(params['num_lines']))

            self.bisect_lead_entry.setText(str(params['bisect_lead']))

            self.survey_speed_entry.setText(str(params.get('survey_speed', '')))

            self.export_name_entry.clear()
            self.export_name_entry.setText(params['export_name'])

            self.offset_direction_var.set(params['offset_direction'])

            self.line_length_multiplier.set(params['line_length_multiplier'])
            self._update_multiplier_label_len(params['line_length_multiplier'])

            self.dist_between_lines_multiplier.set(params['dist_between_lines_multiplier'])
            self._update_multiplier_label_dist(params['dist_between_lines_multiplier'])

            # Regenerate the plot with loaded parameters
            self._generate_and_plot()

            self.set_ref_info_text(f"Survey parameters loaded from: {file_path}", append=False)

        except Exception as e:
            self._show_message("error","Load Error", f"Failed to load survey parameters: {e}")

    def _load_survey_parameters_dialog(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot load parameters.")
            return

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Load Survey Parameters",
            self.last_survey_params_dir,
            "JSON files (*.json);;All files (*.*)"
        )
        if file_path:
            self.last_survey_params_dir = os.path.dirname(file_path)
            self._save_last_survey_params_dir()
            self._load_survey_parameters(file_path)

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

    def _on_line_length_or_speed_change(self, *args):
        pass

    def _draw_crossline_profile(self):
        if not (
            self.geotiff_data_array is not None and
            self.geotiff_extent is not None and
            self.cross_line_data and
            len(self.cross_line_data) == 2
        ):
            # No valid data: clear the profile plot completely
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
        # Clear any existing twin axes (slope axes)
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
        left, right, bottom, top = tuple(self.geotiff_extent)
        nrows, ncols = self.geotiff_data_array.shape
        rows = ((top - lats) / (top - bottom) * (nrows - 1)).clip(0, nrows - 1)
        cols = ((lons - left) / (right - left) * (ncols - 1)).clip(0, ncols - 1)
        elevations = []
        slopes = []
        for r, c in zip(rows, cols):
            ir, ic = int(round(r)), int(round(c))
            elevations.append(self.geotiff_data_array[ir, ic])
            # Slope calculation (in degrees)
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
                if window.shape == (3,3) and not np.all(np.isnan(window)):
                    dz_dy, dz_dx = np.gradient(window, dy_m, dx_m)
                    slope_rad = np.arctan(np.sqrt(dz_dx[1,1]**2 + dz_dy[1,1]**2))
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
        self.profile_ax.plot(dists, elevations, color='purple', lw=1, label='Elevation')
        if self.show_slope_profile_var:
            slope_ax = self.profile_ax.twinx()
            slope_ax.plot(dists, slopes, color='blue', lw=1, linestyle='--', label='Slope (deg)')
            slope_ax.set_ylabel('Slope (deg)', fontsize=8)
            slope_ax.tick_params(axis='y', labelsize=7)
            slope_ax.grid(False)
        self.profile_ax.set_xlim(dists[0], dists[-1])
        if np.any(~np.isnan(elevations)):
            self.profile_ax.set_ylim(np.nanmin(elevations), np.nanmax(elevations))
        if slope_ax and np.any(~np.isnan(slopes)):
            slope_ax.set_ylim(0, np.nanmax(slopes[~np.isnan(slopes)])*1.1)
        
        # Add legend to profile plot (include slope axis if present)
        handles, labels = self.profile_ax.get_legend_handles_labels()
        if slope_ax:
            slope_handles, slope_labels = slope_ax.get_legend_handles_labels()
            handles.extend(slope_handles)
            labels.extend(slope_labels)
        if handles:
            self.profile_ax.legend(handles, labels, fontsize=10, loc='upper right')
        
        self.profile_fig.tight_layout(pad=1.0)
        self.profile_canvas.draw_idle()

        # Show pitch line elevation stats in Calibration Planning info/error area
        if (
            hasattr(self, 'cal_info_text') and
            self.geotiff_data_array is not None and
            hasattr(self, 'pitch_line_points') and
            len(self.pitch_line_points) == 2
        ):
            # Use the elevations array from above if available
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
            # Append median depth usage note
            self.set_cal_info_text("Using Median Depth for Line Offset", append=False)

    def _draw_pitch_line_profile(self):
        # Clear any existing twin axes (slope axes)
        for ax in self.profile_fig.get_axes():
            if ax != self.profile_ax:
                ax.remove()
        self.profile_ax.clear()
        self.profile_ax.set_title("Pitch Line Elevation Profile", fontsize=8)
        self.profile_ax.set_xlabel("Distance (m)", fontsize=8)
        self.profile_ax.set_ylabel("Elevation (m)", fontsize=8)
        self.profile_ax.tick_params(axis='both', which='major', labelsize=7)
        slope_ax = None
        if (
            self.geotiff_data_array is not None and
            self.geotiff_extent is not None and
            hasattr(self, 'pitch_line_points') and
            len(self.pitch_line_points) == 2
        ):
            (lat1, lon1), (lat2, lon2) = self.pitch_line_points
            lats = np.linspace(lat1, lat2, 100)
            lons = np.linspace(lon1, lon2, 100)
            left, right, bottom, top = tuple(self.geotiff_extent)
            nrows, ncols = self.geotiff_data_array.shape
            rows = ((top - lats) / (top - bottom) * (nrows - 1)).clip(0, nrows - 1)
            cols = ((lons - left) / (right - left) * (ncols - 1)).clip(0, ncols - 1)
            elevations = []
            slopes = []
            for r, c in zip(rows, cols):
                ir, ic = int(round(r)), int(round(c))
                elevations.append(self.geotiff_data_array[ir, ic])
                # Slope calculation (in degrees)
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
                    if window.shape == (3,3) and not np.all(np.isnan(window)):
                        dz_dy, dz_dx = np.gradient(window, dy_m, dx_m)
                        slope_rad = np.arctan(np.sqrt(dz_dx[1,1]**2 + dz_dy[1,1]**2))
                        slope = np.degrees(slope_rad)
                slopes.append(slope if slope is not None else np.nan)
            elevations = np.array(elevations)
            slopes = np.array(slopes)
            # Calculate true geodetic distance along the pitch line
            try:
                import pyproj
                geod = pyproj.Geod(ellps="WGS84")
                dists = [0.0]
                for i in range(1, len(lats)):
                    _, _, dist = geod.inv(lons[i-1], lats[i-1], lons[i], lats[i])
                    dists.append(dists[-1] + dist)
                dists = np.array(dists)
            except Exception:
                dists = np.linspace(0, 1, 100)
            self.profile_ax.plot(dists, elevations, color='orange', lw=1, label='Elevation')
            if self.show_slope_profile_var:
                slope_ax = self.profile_ax.twinx()
                slope_ax.plot(dists, slopes, color='blue', lw=1, linestyle='--', label='Slope (deg)')
                slope_ax.set_ylabel('Slope (deg)', fontsize=8)
                slope_ax.tick_params(axis='y', labelsize=7)
                slope_ax.grid(False)
            self.profile_ax.set_xlim(dists[0], dists[-1])
            if np.any(~np.isnan(elevations)):
                self.profile_ax.set_ylim(np.nanmin(elevations), np.nanmax(elevations))
            if slope_ax and np.any(~np.isnan(slopes)):
                slope_ax.set_ylim(0, np.nanmax(slopes[~np.isnan(slopes)])*1.1)
        
        # Add legend to profile plot (include slope axis if present)
        handles, labels = self.profile_ax.get_legend_handles_labels()
        if slope_ax:
            slope_handles, slope_labels = slope_ax.get_legend_handles_labels()
            handles.extend(slope_handles)
            labels.extend(slope_labels)
        if handles:
            self.profile_ax.legend(handles, labels, fontsize=10, loc='upper right')
        
        self.profile_fig.tight_layout(pad=1.0)
        self.profile_canvas.draw_idle()

        # Show pitch line elevation stats in Calibration Planning info/error area
        if (
            hasattr(self, 'cal_info_text') and
            self.geotiff_data_array is not None and
            hasattr(self, 'pitch_line_points') and
            len(self.pitch_line_points) == 2
        ):
            # Use the elevations array from above if available
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
            # Append median depth usage note
            self.set_cal_info_text("Using Median Depth for Line Offset", append=False)

    def _export_survey_files(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot export survey files.")
            return

        is_valid, values = self._validate_inputs()
        if not is_valid:
            return

        export_name = values['export_name']

        if not self.survey_lines_data and not self.cross_line_data:
            self._show_message("warning","No Data", "No survey lines to export. Generate them first.")
            return

        export_dir = QFileDialog.getExistingDirectory(self, "Select Export Directory", self.last_export_dir)
        if not export_dir:
            return
        self.last_export_dir = export_dir
        self._save_last_export_dir()

        try:
            import json
            import fiona
            from shapely.geometry import LineString, mapping

            # --- Export to CSV ---
            # CSV files for decimal degrees always use _DD suffix
            csv_file_path = os.path.join(export_dir, f"{export_name}_DD.csv")
            with open(csv_file_path, 'w', newline='') as csvfile:
                csv_writer = csv.writer(csvfile)
                csv_writer.writerow(['Line Number', 'Point Label', 'Latitude', 'Longitude'])
                # Main survey lines
                for i, line in enumerate(self.survey_lines_data):
                    if i % 2 == 0:  # Even lines - normal order
                        start, end = line[0], line[1]
                    else:  # Odd lines - flipped order
                        start, end = line[1], line[0]
                    start_label = f'L{i+1}S'
                    end_label = f'L{i+1}E'
                    csv_writer.writerow([i + 1, start_label, start[0], start[1]])
                    csv_writer.writerow([i + 1, end_label, end[0], end[1]])
                # Crossline
                if self.cross_line_data:
                    csv_writer.writerow([0, 'CLS', self.cross_line_data[0][0], self.cross_line_data[0][1]])
                    csv_writer.writerow([0, 'CLE', self.cross_line_data[1][0], self.cross_line_data[1][1]])
            
            # --- Export to DDM format (Decimal Minutes) ---
            ddm_file_path = os.path.join(export_dir, f"{export_name}_DM.csv")
            with open(ddm_file_path, 'w', newline='', encoding='utf-8') as ddmfile:
                ddm_writer = csv.writer(ddmfile)
                ddm_writer.writerow(['Line Number', 'Point Label', 'Latitude with Decimal Minutes', 'Longitude with Decimal Minutes'])
                # Main survey lines
                for i, line in enumerate(self.survey_lines_data):
                    if i % 2 == 0:  # Even lines - normal order
                        start, end = line[0], line[1]
                    else:  # Odd lines - flipped order
                        start, end = line[1], line[0]
                    start_label = f'L{i+1}S'
                    end_label = f'L{i+1}E'
                    start_lat_ddm = self._decimal_degrees_to_ddm(start[0], is_latitude=True)
                    start_lon_ddm = self._decimal_degrees_to_ddm(start[1], is_latitude=False)
                    end_lat_ddm = self._decimal_degrees_to_ddm(end[0], is_latitude=True)
                    end_lon_ddm = self._decimal_degrees_to_ddm(end[1], is_latitude=False)
                    ddm_writer.writerow([i + 1, start_label, start_lat_ddm, start_lon_ddm])
                    ddm_writer.writerow([i + 1, end_label, end_lat_ddm, end_lon_ddm])
                # Crossline
                if self.cross_line_data:
                    cls_lat_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[0][0], is_latitude=True)
                    cls_lon_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[0][1], is_latitude=False)
                    cle_lat_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[1][0], is_latitude=True)
                    cle_lon_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[1][1], is_latitude=False)
                    ddm_writer.writerow([0, 'CLS', cls_lat_ddm, cls_lon_ddm])
                    ddm_writer.writerow([0, 'CLE', cle_lat_ddm, cle_lon_ddm])
            
            # --- Export to DDM text format (Decimal Minutes) ---
            ddm_txt_file_path = os.path.join(export_dir, f"{export_name}_DM.txt")
            with open(ddm_txt_file_path, 'w', encoding='utf-8') as ddm_txt_file:
                # Main survey lines
                for i, line in enumerate(self.survey_lines_data):
                    if i % 2 == 0:  # Even lines - normal order
                        start, end = line[0], line[1]
                    else:  # Odd lines - flipped order
                        start, end = line[1], line[0]
                    start_label = f'L{i+1}S'
                    end_label = f'L{i+1}E'
                    start_lat_ddm = self._decimal_degrees_to_ddm(start[0], is_latitude=True)
                    start_lon_ddm = self._decimal_degrees_to_ddm(start[1], is_latitude=False)
                    end_lat_ddm = self._decimal_degrees_to_ddm(end[0], is_latitude=True)
                    end_lon_ddm = self._decimal_degrees_to_ddm(end[1], is_latitude=False)
                    ddm_txt_file.write(f"{start_label}, {start_lat_ddm}, {start_lon_ddm}\n")
                    ddm_txt_file.write(f"{end_label}, {end_lat_ddm}, {end_lon_ddm}\n")
                # Crossline
                if self.cross_line_data:
                    cls_lat_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[0][0], is_latitude=True)
                    cls_lon_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[0][1], is_latitude=False)
                    cle_lat_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[1][0], is_latitude=True)
                    cle_lon_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[1][1], is_latitude=False)
                    ddm_txt_file.write(f"CLS, {cls_lat_ddm}, {cls_lon_ddm}\n")
                    ddm_txt_file.write(f"CLE, {cle_lat_ddm}, {cle_lon_ddm}\n")

            # --- Export to ESRI Shapefile (.shp) ---
            schema = {
                'geometry': 'LineString',
                'properties': {'line_num': 'int'},
            }
            crs_epsg = 'EPSG:4326'  # WGS 84
            features = []
            # Add main survey lines
            for i, line_coords in enumerate(self.survey_lines_data):
                shapely_line = LineString([(p[1], p[0]) for p in line_coords])
                features.append({
                    'geometry': mapping(shapely_line),
                    'properties': {'line_num': i + 1},
                })
            # Add crossline
            if self.cross_line_data:
                shapely_cross_line = LineString([(p[1], p[0]) for p in self.cross_line_data])
                features.append({
                    'geometry': mapping(shapely_cross_line),
                    'properties': {'line_num': 0},
                })
            shapefile_path = os.path.join(export_dir, f"{export_name}.shp")
            with fiona.open(shapefile_path, 'w', driver='ESRI Shapefile', crs=crs_epsg, schema=schema) as collection:
                collection.writerecords(features)

            # --- Export to GeoJSON ---
            geojson_file_path = os.path.join(export_dir, f"{export_name}.geojson")
            geojson_features = []
            # Main survey lines
            for i, line in enumerate(self.survey_lines_data):
                geojson_features.append({
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": [[line[0][1], line[0][0]], [line[1][1], line[1][0]]]
                    },
                    "properties": {
                        "line_num": i + 1,
                        "points": [
                            {"point_num": 1, "lat": line[0][0], "lon": line[0][1]},
                            {"point_num": 2, "lat": line[1][0], "lon": line[1][1]}
                        ]
                    }
                })
            # Crossline
            if self.cross_line_data:
                geojson_features.append({
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": [[self.cross_line_data[0][1], self.cross_line_data[0][0]], [self.cross_line_data[1][1], self.cross_line_data[1][0]]]
                    },
                    "properties": {
                        "line_num": 0,
                        "points": [
                            {"point_num": 1, "lat": self.cross_line_data[0][0], "lon": self.cross_line_data[0][1]},
                            {"point_num": 2, "lat": self.cross_line_data[1][0], "lon": self.cross_line_data[1][1]}
                        ]
                    }
                })
            geojson_collection = {
                "type": "FeatureCollection",
                "features": geojson_features
            }
            with open(geojson_file_path, 'w') as f:
                json.dump(geojson_collection, f, indent=2)

            self.set_ref_info_text(
                f"Survey exported successfully to:\n"
                f"- {os.path.basename(csv_file_path)}\n"
                f"- {os.path.basename(shapefile_path)} (and associated files)\n"
                f"- {os.path.basename(geojson_file_path)}\n"
                f"in directory: {export_dir}", append=False)

            # --- Export to Hypack LNW format ---
            lnw_file_path = os.path.join(export_dir, f"{export_name}.lnw")
            with open(lnw_file_path, 'w') as f:
                f.write("LNW 1.0\n")
                # Main survey lines
                for i, line in enumerate(self.survey_lines_data):
                    line_num = i + 1
                    # Start point
                    lat1, lon1 = line[0]
                    lat2, lon2 = line[1]
                    # Get depth from GeoTIFF if available
                    depth1 = self._get_depth_at_point(lat1, lon1) if self.geotiff_data_array is not None else 0.0
                    depth2 = self._get_depth_at_point(lat2, lon2) if self.geotiff_data_array is not None else 0.0
                    # Get speed from survey speed entry
                    try:
                        speed = float(self.survey_speed_entry.text()) if self.survey_speed_entry.text() else 8.0
                    except:
                        speed = 8.0
                    f.write(f"LINE{line_num:03d}_001,{lat1:.6f},{lon1:.6f},{abs(depth1):.1f},{speed:.1f},50.0,{line_num},1\n")
                    f.write(f"LINE{line_num:03d}_002,{lat2:.6f},{lon2:.6f},{abs(depth2):.1f},{speed:.1f},50.0,{line_num},1\n")
                # Crossline
                if self.cross_line_data:
                    lat1, lon1 = self.cross_line_data[0]
                    lat2, lon2 = self.cross_line_data[1]
                    depth1 = self._get_depth_at_point(lat1, lon1) if self.geotiff_data_array is not None else 0.0
                    depth2 = self._get_depth_at_point(lat2, lon2) if self.geotiff_data_array is not None else 0.0
                    try:
                        speed = float(self.survey_speed_entry.text()) if self.survey_speed_entry.text() else 8.0
                    except:
                        speed = 8.0
                    f.write(f"CROSS_001,{lat1:.6f},{lon1:.6f},{abs(depth1):.1f},{speed:.1f},50.0,0,2\n")
                    f.write(f"CROSS_002,{lat2:.6f},{lon2:.6f},{abs(depth2):.1f},{speed:.1f},50.0,0,2\n")

            self.set_ref_info_text(
                f"Survey exported successfully to:\n"
                f"- {os.path.basename(csv_file_path)}\n"
                f"- {os.path.basename(shapefile_path)} (and associated files)\n"
                f"- {os.path.basename(geojson_file_path)}\n"
                f"- {os.path.basename(lnw_file_path)}\n"
                f"in directory: {export_dir}", append=False)

            # --- Export to SIS ASCII Plan format ---
            sis_file_path = os.path.join(export_dir, f"{export_name}.asciiplan")
            with open(sis_file_path, 'w') as f:
                f.write("SIS ASCII Plan\n")
                # Main survey lines
                for i, line in enumerate(self.survey_lines_data):
                    line_num = i + 1
                    # Start point
                    lat1, lon1 = line[0]
                    lat2, lon2 = line[1]
                    # Get depth from GeoTIFF if available
                    depth1 = self._get_depth_at_point(lat1, lon1) if self.geotiff_data_array is not None else 0.0
                    depth2 = self._get_depth_at_point(lat2, lon2) if self.geotiff_data_array is not None else 0.0
                    # Get speed from survey speed entry
                    try:
                        speed = float(self.survey_speed_entry.text()) if self.survey_speed_entry.text() else 8.0
                    except:
                        speed = 8.0
                    f.write(f"LINE{line_num:03d}_001, {lat1:.6f}, {lon1:.6f}, {abs(depth1):.1f}, {speed:.1f}, {line_num}, {line_num}\n")
                    f.write(f"LINE{line_num:03d}_002, {lat2:.6f}, {lon2:.6f}, {abs(depth2):.1f}, {speed:.1f}, {line_num}, {line_num}\n")
                
                # Crossline if present
                if self.cross_line_data:
                    lat1, lon1 = self.cross_line_data[0]
                    lat2, lon2 = self.cross_line_data[1]
                    depth1 = self._get_depth_at_point(lat1, lon1) if self.geotiff_data_array is not None else 0.0
                    depth2 = self._get_depth_at_point(lat2, lon2) if self.geotiff_data_array is not None else 0.0
                    try:
                        speed = float(self.survey_speed_entry.text()) if self.survey_speed_entry.text() else 8.0
                    except:
                        speed = 8.0
                    crossline_num = len(self.survey_lines_data) + 1
                    f.write(f"CROSSLINE_001, {lat1:.6f}, {lon1:.6f}, {abs(depth1):.1f}, {speed:.1f}, {crossline_num}, {crossline_num}\n")
                    f.write(f"CROSSLINE_002, {lat2:.6f}, {lon2:.6f}, {abs(depth2):.1f}, {speed:.1f}, {crossline_num}, {crossline_num}\n")

            # Update success message to include SIS
            self.set_ref_info_text(
                f"Survey exported successfully to:\n"
                f"- {os.path.basename(csv_file_path)}\n"
                f"- {os.path.basename(shapefile_path)} (and associated files)\n"
                f"- {os.path.basename(geojson_file_path)}\n"
                f"- {os.path.basename(lnw_file_path)}\n"
                f"- {os.path.basename(sis_file_path)}\n"
                f"in directory: {export_dir}", append=False)

            # --- Export Comprehensive Survey Statistics ---
            stats_file_path = os.path.join(export_dir, f"{export_name}_stats.txt")
            total_survey_time = self._calculate_total_survey_time()
            
            with open(stats_file_path, 'w') as f:
                f.write("COMPREHENSIVE SURVEY STATISTICS\n")
                f.write("=" * 50 + "\n\n")
                f.write(f"Survey Plan: {export_name}\n")
                f.write(f"Export Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                f.write("SURVEY INPUT PARAMETERS\n")
                f.write("-" * 30 + "\n")
                for key, value in values.items():
                    f.write(f"{key.replace('_', ' ').title()}: {value}\n")
                f.write("\n")
                
                f.write("SURVEY DISTANCE BREAKDOWN\n")
                f.write("-" * 30 + "\n")
                # Calculate single line length and heading
                num_main_lines = len(self.survey_lines_data) if self.survey_lines_data else 0
                if num_main_lines > 0:
                    # Calculate heading for the first main line
                    try:
                        import pyproj
                        geod = pyproj.Geod(ellps="WGS84")
                        first_line = self.survey_lines_data[0]
                        lat1, lon1 = first_line[0]
                        lat2, lon2 = first_line[1]
                        fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)
                        
                        heading = fwd_az % 360
                        reciprocal_heading = back_az % 360
                        
                        f.write(f"Heading: {heading:.1f}°\n")
                        f.write(f"Reciprocal Heading: {reciprocal_heading:.1f}°\n")
                    except Exception:
                        f.write("Heading: Unable to calculate\n")
                    
                    single_line_length_m = total_survey_time['main_lines_distance_m'] / num_main_lines
                    single_line_length_km = single_line_length_m / 1000
                    single_line_length_nm = single_line_length_m / 1852
                    f.write(f"Length of Main Line: {single_line_length_m:.1f} m ({single_line_length_km:.3f} km, {single_line_length_nm:.3f} nm)\n")
                f.write(f"Main Lines Total Distance: {total_survey_time['main_lines_distance_m']:.1f} m\n")
                f.write(f"Main Lines Total Distance: {total_survey_time['main_lines_distance_km']:.3f} km\n")
                f.write(f"Main Lines Total Distance: {total_survey_time['main_lines_distance_nm']:.3f} nm\n")
                f.write(f"Travel Between Lines: {total_survey_time['travel_between_lines_distance_m']:.1f} m\n")
                f.write(f"Travel Between Lines: {total_survey_time['travel_between_lines_distance_km']:.3f} km\n")
                f.write(f"Travel Between Lines: {total_survey_time['travel_between_lines_distance_nm']:.3f} nm\n")
                f.write(f"Travel to Crossline: {total_survey_time['travel_to_crossline_distance_m']:.1f} m\n")
                f.write(f"Travel to Crossline: {total_survey_time['travel_to_crossline_distance_km']:.3f} km\n")
                f.write(f"Travel to Crossline: {total_survey_time['travel_to_crossline_distance_nm']:.3f} nm\n")
                f.write(f"Crossline (per pass): {total_survey_time['crossline_single_pass_distance_m']:.1f} m\n")
                f.write(f"Crossline (per pass): {total_survey_time['crossline_single_pass_distance_km']:.3f} km\n")
                f.write(f"Crossline (per pass): {total_survey_time['crossline_single_pass_distance_nm']:.3f} nm\n")
                f.write(f"Crossline Passes: {total_survey_time['num_crossline_passes']}\n")
                
                # Add crossline heading information
                if self.cross_line_data and len(self.cross_line_data) >= 2:
                    try:
                        if pyproj is not None:
                            geod = pyproj.Geod(ellps="WGS84")
                            lat1, lon1 = self.cross_line_data[0]
                            lat2, lon2 = self.cross_line_data[1]
                            fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)
                            
                            crossline_heading = fwd_az % 360
                            crossline_reciprocal_heading = back_az % 360
                            
                            f.write(f"Crossline Heading: {crossline_heading:.1f}°\n")
                            f.write(f"Crossline Reciprocal Heading: {crossline_reciprocal_heading:.1f}°\n")
                        else:
                            f.write("Crossline Heading: pyproj not available\n")
                    except Exception:
                        f.write("Crossline Heading: Unable to calculate\n")
                
                f.write(f"Crossline Total Distance: {total_survey_time['crossline_total_distance_m']:.1f} m\n")
                f.write(f"Crossline Total Distance: {total_survey_time['crossline_total_distance_km']:.3f} km\n")
                f.write(f"Crossline Total Distance: {total_survey_time['crossline_total_distance_nm']:.3f} nm\n\n")
                
                f.write("TOTAL SURVEY DISTANCE\n")
                f.write("-" * 25 + "\n")
                f.write(f"Total Survey Distance: {total_survey_time['total_distance_m']:.1f} m\n")
                f.write(f"Total Survey Distance: {total_survey_time['total_distance_km']:.3f} km\n")
                f.write(f"Total Survey Distance: {total_survey_time['total_distance_nm']:.3f} nm\n\n")
                
                f.write("SURVEY TIME BREAKDOWN\n")
                f.write("-" * 25 + "\n")
                f.write(f"Main Lines Survey: {total_survey_time['main_lines_minutes']:.1f} min\n")
                f.write(f"Crossline Survey: {total_survey_time['crossline_minutes']:.1f} min\n")
                f.write(f"Travel Between Lines: {total_survey_time['travel_minutes']:.1f} min\n")
                f.write(f"Travel to Crossline: {total_survey_time['travel_to_crossline_minutes']:.1f} min\n\n")
                
                f.write("TOTAL SURVEY TIME\n")
                f.write("-" * 20 + "\n")
                f.write(f"Total Survey Time: {total_survey_time['total_minutes']:.1f} min\n")
                f.write(f"Total Survey Time: {total_survey_time['total_hours']:.2f} hr\n\n")
                
                f.write("SURVEY LINE DETAILS\n")
                f.write("-" * 20 + "\n")
                f.write(f"Number of Main Survey Lines: {len(self.survey_lines_data)}\n")
                if self.cross_line_data:
                    f.write("Crossline: Present\n")
                else:
                    f.write("Crossline: Not present\n")
                f.write(f"Survey Speed: {values.get('survey_speed', '8.0')} knots\n")
                f.write(f"Crossline Passes: {total_survey_time['num_crossline_passes']}\n\n")
                
                f.write("SURVEY PATTERN\n")
                f.write("-" * 15 + "\n")
                f.write("Main lines are surveyed in a zigzag pattern:\n")
                f.write("- Even-numbered lines (2, 4, 6...): Survey from start to end\n")
                f.write("- Odd-numbered lines (1, 3, 5...): Survey from end to start (flipped)\n")
                f.write("- This minimizes travel distance between lines\n")
                f.write("- After completing all main lines, travel to crossline start\n")
                f.write("- Complete crossline survey with specified number of passes\n\n")
                
                # Waypoints with labels
                f.write("WAYPOINTS\n")
                f.write("-" * 10 + "\n")
                if self.survey_lines_data:
                    for i, line in enumerate(self.survey_lines_data):
                        # Determine start and end based on zigzag pattern
                        if i % 2 == 0:  # Even lines - normal order
                            start, end = line[0], line[1]
                        else:  # Odd lines - flipped order
                            start, end = line[1], line[0]
                        start_label = f'L{i+1}S'
                        end_label = f'L{i+1}E'
                        start_lat_ddm = self._decimal_degrees_to_ddm(start[0], is_latitude=True)
                        start_lon_ddm = self._decimal_degrees_to_ddm(start[1], is_latitude=False)
                        end_lat_ddm = self._decimal_degrees_to_ddm(end[0], is_latitude=True)
                        end_lon_ddm = self._decimal_degrees_to_ddm(end[1], is_latitude=False)
                        f.write(f"{start_label}: {start_lat_ddm}, {start_lon_ddm} ({start[0]:.6f}, {start[1]:.6f})\n")
                        f.write(f"{end_label}: {end_lat_ddm}, {end_lon_ddm} ({end[0]:.6f}, {end[1]:.6f})\n")
                if self.cross_line_data:
                    cls_lat_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[0][0], is_latitude=True)
                    cls_lon_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[0][1], is_latitude=False)
                    cle_lat_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[1][0], is_latitude=True)
                    cle_lon_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[1][1], is_latitude=False)
                    f.write(f"CLS: {cls_lat_ddm}, {cls_lon_ddm} ({self.cross_line_data[0][0]:.6f}, {self.cross_line_data[0][1]:.6f})\n")
                    f.write(f"CLE: {cle_lat_ddm}, {cle_lon_ddm} ({self.cross_line_data[1][0]:.6f}, {self.cross_line_data[1][1]:.6f})\n")
                f.write("\n")
                
                f.write("NOTES\n")
                f.write("-" * 5 + "\n")
                f.write("- All distances calculated using geodetic (spherical) geometry\n")
                f.write("- Times calculated based on specified survey speed\n")
                f.write("- Crossline lead-in/out distance extends beyond main survey area\n")
                f.write("- Survey pattern optimized for efficiency\n")

            # --- Export to Text format (point label, latitude, longitude) ---
            txt_file_path = os.path.join(export_dir, f"{export_name}_DD.txt")
            with open(txt_file_path, 'w') as f:
                # Main survey lines
                for i, line in enumerate(self.survey_lines_data):
                    if i % 2 == 0:  # Even lines - normal order
                        start, end = line[0], line[1]
                    else:  # Odd lines - flipped order
                        start, end = line[1], line[0]
                    start_label = f'L{i+1}S'
                    end_label = f'L{i+1}E'
                    f.write(f"{start_label} {start[0]:.6f} {start[1]:.6f}\n")
                    f.write(f"{end_label} {end[0]:.6f} {end[1]:.6f}\n")
                # Crossline
                if self.cross_line_data:
                    f.write(f"CLS {self.cross_line_data[0][0]:.6f} {self.cross_line_data[0][1]:.6f}\n")
                    f.write(f"CLE {self.cross_line_data[1][0]:.6f} {self.cross_line_data[1][1]:.6f}\n")

            # --- Export Survey Plan as PNG ---
            survey_plan_png_path = os.path.join(export_dir, f"{export_name}_survey_plan.png")
            self.figure.savefig(survey_plan_png_path, dpi=300, bbox_inches='tight', facecolor='white')
            
            # --- Export Profile Plot as PNG (if it exists) ---
            profile_png_path = os.path.join(export_dir, f"{export_name}_profile.png")
            if hasattr(self, 'profile_fig') and self.profile_fig is not None:
                self.profile_fig.savefig(profile_png_path, dpi=300, bbox_inches='tight', facecolor='white')

            # --- Export parameters metadata as JSON ---
            json_metadata_path = os.path.join(export_dir, f"{export_name}_params.json")
            try:
                # Get current parameter values (even if not all fields are filled, use what we have)
                params = {}
                try:
                    params['central_lat'] = float(self.central_lat_entry.text())
                except:
                    # Calculate from survey lines if not available
                    if self.survey_lines_data:
                        all_lats = [p[0] for line in self.survey_lines_data for p in line]
                        params['central_lat'] = (min(all_lats) + max(all_lats)) / 2.0
                    else:
                        params['central_lat'] = None
                
                try:
                    params['central_lon'] = float(self.central_lon_entry.text())
                except:
                    if self.survey_lines_data:
                        all_lons = [p[1] for line in self.survey_lines_data for p in line]
                        params['central_lon'] = (min(all_lons) + max(all_lons)) / 2.0
                    else:
                        params['central_lon'] = None
                
                try:
                    params['line_length'] = float(self.line_length_entry.text())
                except:
                    params['line_length'] = None
                
                try:
                    params['heading'] = float(self.heading_entry.get())
                except:
                    params['heading'] = None
                
                try:
                    params['dist_between_lines'] = float(self.dist_between_lines_entry.get())
                except:
                    params['dist_between_lines'] = None
                
                try:
                    params['num_lines'] = int(self.num_lines_entry.get())
                except:
                    params['num_lines'] = len(self.survey_lines_data) if self.survey_lines_data else None
                
                try:
                    params['bisect_lead'] = float(self.bisect_lead_entry.get())
                except:
                    params['bisect_lead'] = None
                
                try:
                    params['survey_speed'] = float(self.survey_speed_entry.get())
                except:
                    params['survey_speed'] = 8.0  # Default
                
                try:
                    params['crossline_passes'] = int(self.crossline_passes_entry.get())
                except:
                    params['crossline_passes'] = 2  # Default
                
                try:
                    params['export_name'] = self.export_name_entry.get().strip()
                except:
                    params['export_name'] = export_name
                
                try:
                    params['offset_direction'] = self.offset_direction_var
                except:
                    params['offset_direction'] = 'North'  # Default
                
                try:
                    params['line_length_multiplier'] = self.line_length_multiplier
                except:
                    params['line_length_multiplier'] = 8.0  # Default
                
                try:
                    params['dist_between_lines_multiplier'] = self.dist_between_lines_multiplier
                except:
                    params['dist_between_lines_multiplier'] = 1.0  # Default
                
                # Save metadata
                with open(json_metadata_path, 'w', encoding='utf-8') as f:
                    json.dump(params, f, indent=2)
            except Exception as e:
                # If metadata export fails, continue without it
                print(f"Warning: Could not export metadata: {e}")
            
            # Update success message to include stats file, text file, PNG files, and metadata
            success_files = [
                f"- {os.path.basename(csv_file_path)}",
                f"- {os.path.basename(shapefile_path)} (and associated files)",
                f"- {os.path.basename(geojson_file_path)}",
                f"- {os.path.basename(lnw_file_path)}",
                f"- {os.path.basename(sis_file_path)}",
                f"- {os.path.basename(txt_file_path)}",
                f"- {os.path.basename(ddm_file_path)}",
                f"- {os.path.basename(ddm_txt_file_path)}",
                f"- {os.path.basename(stats_file_path)}",
                f"- {os.path.basename(survey_plan_png_path)}"
            ]
            
            # Add metadata JSON if it was created
            try:
                if os.path.exists(json_metadata_path):
                    success_files.append(f"- {os.path.basename(json_metadata_path)}")
            except:
                pass
            
            # Add profile PNG if it was created
            if hasattr(self, 'profile_fig') and self.profile_fig is not None:
                success_files.append(f"- {os.path.basename(profile_png_path)}")
            
            self.set_ref_info_text(
                f"Survey exported successfully to:\n" + "\n".join(success_files) + 
                f"\nin directory: {export_dir}", append=False)

        except Exception as e:
            self._show_message("error","Export Error", f"Failed to export survey files: {e}")

    def _import_survey_files(self):
        """Import reference planning survey lines from CSV or GeoJSON file."""
        import csv
        import json
        
        # Open file dialog to select import file
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Survey File to Import",
            self.last_ref_import_dir,
            "Decimal Degree CSV files (*_DD.csv);;CSV files (*.csv);;GeoJSON files (*.geojson);;JSON files (*.json);;All files (*.*)"
        )
        
        if not file_path:
            return
        
        # Save the directory for next time
        import_dir = os.path.dirname(file_path)
        if import_dir and os.path.isdir(import_dir):
            self.last_ref_import_dir = import_dir
            self._save_last_ref_import_dir()
        
        try:
            # Clear existing lines
            self.survey_lines_data = []
            self.cross_line_data = []
            
            # Determine file type and import accordingly
            file_ext = os.path.splitext(file_path)[1].lower()
            
            if file_ext == '.csv':
                # Import from CSV
                with open(file_path, 'r', encoding='utf-8') as csvfile:
                    csv_reader = csv.DictReader(csvfile)
                    
                    # Dictionary to collect points for each line
                    lines_data = {}
                    crossline_points = {}
                    
                    for row in csv_reader:
                        try:
                            line_num = int(row.get('Line Number', -1))
                            point_label = row.get('Point Label', '').strip()
                            lat = float(row.get('Latitude', 0))
                            lon = float(row.get('Longitude', 0))
                        except (ValueError, TypeError):
                            continue
                        
                        # Handle crossline (line_num == 0)
                        if line_num == 0:
                            if point_label == 'CLS':
                                crossline_points['start'] = (lat, lon)
                            elif point_label == 'CLE':
                                crossline_points['end'] = (lat, lon)
                        else:
                            # Main survey line
                            if line_num not in lines_data:
                                lines_data[line_num] = {}
                            lines_data[line_num][point_label] = (lat, lon)
                    
                    # Process main survey lines
                    for line_num in sorted(lines_data.keys()):
                        points = lines_data[line_num]
                        # Look for L{n}S and L{n}E labels
                        start_label = f'L{line_num}S'
                        end_label = f'L{line_num}E'
                        if start_label in points and end_label in points:
                            self.survey_lines_data.append([points[start_label], points[end_label]])
                    
                    # Process crossline
                    if 'start' in crossline_points and 'end' in crossline_points:
                        self.cross_line_data = [crossline_points['start'], crossline_points['end']]
            
            elif file_ext in ['.geojson', '.json']:
                # Import from GeoJSON
                with open(file_path, 'r', encoding='utf-8') as f:
                    geojson_data = json.load(f)
                
                # Ensure geojson_data is a dictionary
                if not isinstance(geojson_data, dict):
                    self._show_message("error","Import Error", "Invalid GeoJSON format: root must be an object")
                    return
                
                # Handle FeatureCollection
                if geojson_data.get('type') == 'FeatureCollection':
                    features = geojson_data.get('features', [])
                    if not isinstance(features, list):
                        features = []
                elif geojson_data.get('type') == 'Feature':
                    features = [geojson_data]
                else:
                    features = []
                
                # Sort features by line_num to maintain order
                features_with_num = []
                for feature in features:
                    if not isinstance(feature, dict) or feature.get('type') != 'Feature':
                        continue
                    geometry = feature.get('geometry', {})
                    if not isinstance(geometry, dict) or geometry.get('type') != 'LineString':
                        continue
                    properties = feature.get('properties', {})
                    if not isinstance(properties, dict):
                        properties = {}
                    line_num = properties.get('line_num', 0)
                    features_with_num.append((line_num, feature))
                
                # Sort by line_num
                features_with_num.sort(key=lambda x: x[0])
                
                for line_num, feature in features_with_num:
                    if not isinstance(feature, dict):
                        continue
                    geometry = feature.get('geometry', {})
                    if not isinstance(geometry, dict):
                        continue
                    coordinates = geometry.get('coordinates', [])
                    if not isinstance(coordinates, list) or len(coordinates) < 2:
                        continue
                    
                    # GeoJSON uses [lon, lat] format, we need [lat, lon]
                    point1 = (coordinates[0][1], coordinates[0][0])
                    point2 = (coordinates[-1][1], coordinates[-1][0])
                    
                    if line_num == 0:
                        # Crossline
                        self.cross_line_data = [point1, point2]
                    else:
                        # Main survey line
                        self.survey_lines_data.append([point1, point2])
            
            else:
                self._show_message("error","Import Error", f"Unsupported file format: {file_ext}")
                return
            
            if not self.survey_lines_data and not self.cross_line_data:
                self._show_message("warning","Import Warning", "No valid survey lines found in the selected file.")
                return
            
            # Try to load parameters metadata JSON file
            base_name = os.path.splitext(os.path.basename(file_path))[0]
            dir_name = os.path.dirname(file_path)
            metadata_path = os.path.join(dir_name, f"{base_name}_params.json")
            params = None
            
            if os.path.exists(metadata_path):
                try:
                    with open(metadata_path, 'r', encoding='utf-8') as f:
                        params = json.load(f)
                except Exception as e:
                    print(f"Warning: Could not load metadata file: {e}")
            
            # Populate parameter fields from metadata or calculate from lines
            try:
                if params:
                    # Use metadata if available
                    if params.get('central_lat') is not None:
                        self.central_lat_entry.clear()
                        self.central_lat_entry.setText(str(params['central_lat']))
                        # Reset "Pick Center from GeoTIFF" button to normal style after importing survey
                        if hasattr(self, 'pick_center_btn'):
                            self.pick_center_btn.setStyleSheet("")
                    
                    if params.get('central_lon') is not None:
                        self.central_lon_entry.clear()
                        self.central_lon_entry.setText(str(params['central_lon']))
                    
                    if params.get('line_length') is not None:
                        self.line_length_entry.setText(str(params['line_length']))
                    
                    if params.get('heading') is not None:
                        self.heading_entry.setText(str(params['heading']))
                    
                    if params.get('dist_between_lines') is not None:
                        self.dist_between_lines_entry.setText(str(params['dist_between_lines']))
                    
                    if params.get('num_lines') is not None:
                        self.num_lines_entry.setText(str(params['num_lines']))
                    
                    if params.get('bisect_lead') is not None:
                        self.bisect_lead_entry.setText(str(params['bisect_lead']))
                    
                    if params.get('survey_speed') is not None:
                        self.survey_speed_entry.setText(str(params['survey_speed']))
                    
                    if params.get('crossline_passes') is not None:
                        self.crossline_passes_entry.setText(str(params['crossline_passes']))
                    
                    if params.get('export_name'):
                        self.export_name_entry.clear()
                        self.export_name_entry.setText(params['export_name'])
                    
                    if params.get('offset_direction'):
                        self.offset_direction_var.set(params['offset_direction'])
                    
                    if params.get('line_length_multiplier') is not None:
                        self.line_length_multiplier.set(params['line_length_multiplier'])
                        self._update_multiplier_label_len(params['line_length_multiplier'])
                    
                    if params.get('dist_between_lines_multiplier') is not None:
                        self.dist_between_lines_multiplier.set(params['dist_between_lines_multiplier'])
                        self._update_multiplier_label_dist(params['dist_between_lines_multiplier'])
                else:
                    # Calculate parameters from imported lines
                    import pyproj
                    
                    # Calculate central lat/lon from all points
                    all_points = []
                    for line in self.survey_lines_data:
                        all_points.extend(line)
                    if self.cross_line_data:
                        all_points.extend(self.cross_line_data)
                    
                    if all_points:
                        all_lats = [p[0] for p in all_points]
                        all_lons = [p[1] for p in all_points]
                        central_lat = (min(all_lats) + max(all_lats)) / 2.0
                        central_lon = (min(all_lons) + max(all_lons)) / 2.0
                        
                        self.central_lat_entry.clear()
                        self.central_lat_entry.setText(f"{central_lat:.6f}")
                        # Reset "Pick Center from GeoTIFF" button to normal style after calculating from imported lines
                        if hasattr(self, 'pick_center_btn'):
                            self.pick_center_btn.setStyleSheet("")
                        
                        self.central_lon_entry.clear()
                        self.central_lon_entry.setText(f"{central_lon:.6f}")
                    
                    # Calculate heading from first line
                    if len(self.survey_lines_data) > 0:
                        first_line = self.survey_lines_data[0]
                        try:
                            geod = pyproj.Geod(ellps="WGS84")
                            lat1, lon1 = first_line[0]
                            lat2, lon2 = first_line[1]
                            fwd_az, back_az, dist = geod.inv(lon1, lat1, lon2, lat2)
                            heading = fwd_az % 360
                            
                            self.heading_entry.clear()
                            self.heading_entry.setText(f"{heading:.1f}")
                            
                            # Calculate line length
                            line_length = dist
                            self.line_length_entry.clear()
                            self.line_length_entry.setText(f"{line_length:.1f}")
                        except Exception:
                            pass
                    
                    # Count number of lines
                    if self.survey_lines_data:
                        self.num_lines_entry.clear()
                        self.num_lines_entry.setText(str(len(self.survey_lines_data)))
                    
                    # Calculate distance between lines (approximate)
                    if len(self.survey_lines_data) > 1:
                        try:
                            geod = pyproj.Geod(ellps="WGS84")
                            # Get perpendicular distance between first two lines
                            line1_mid = ((self.survey_lines_data[0][0][0] + self.survey_lines_data[0][1][0]) / 2,
                                       (self.survey_lines_data[0][0][1] + self.survey_lines_data[0][1][1]) / 2)
                            line2_mid = ((self.survey_lines_data[1][0][0] + self.survey_lines_data[1][1][0]) / 2,
                                       (self.survey_lines_data[1][0][1] + self.survey_lines_data[1][1][1]) / 2)
                            _, _, dist = geod.inv(line1_mid[1], line1_mid[0], line2_mid[1], line2_mid[0])
                            self.dist_between_lines_entry.clear()
                            self.dist_between_lines_entry.setText(f"{dist:.1f}")
                        except Exception:
                            pass
                    
                    # Update export name if not set
                    if not self.export_name_entry.text().strip():
                        self.export_name_entry.clear()
                        self.export_name_entry.setText(base_name)
            except Exception as e:
                print(f"Warning: Error populating parameter fields: {e}")
            
            # Update plot
            self._plot_survey_plan(preserve_view_limits=True)
            
            # Show success message
            imported_items = []
            if self.survey_lines_data:
                imported_items.append(f"{len(self.survey_lines_data)} survey line(s)")
            if self.cross_line_data:
                imported_items.append("crossline")
            
            if imported_items:
                msg = f"Successfully imported: {', '.join(imported_items)}"
                if params:
                    msg += " (parameters loaded from metadata)"
                else:
                    msg += " (parameters calculated from lines)"
                self.set_ref_info_text(msg)
            else:
                self._show_message("warning","Import Warning", "No valid survey lines found in the selected file.")
        
        except Exception as e:
            self._show_message("error","Import Error", f"Failed to import survey files: {e}")
            import traceback
            traceback.print_exc()

    def _update_export_name(self):
        try:
            dist = int(float(self.dist_between_lines_entry.text()))
            heading = int(float(self.heading_entry.text()))
            export_name = f"Reference_{dist}m_{heading}deg"
            self.export_name_entry.clear()
            self.export_name_entry.setText(export_name)
        except Exception:
            pass

    def _update_cal_line_offset_from_pitch_line(self):
        # Only run if we have a loaded GeoTIFF and a valid pitch line
        if (self.geotiff_data_array is None or self.geotiff_extent is None or
            not hasattr(self, 'pitch_line_points') or len(self.pitch_line_points) != 2):
            self.cal_line_offset_entry.clear()
            return
        (lat1, lon1), (lat2, lon2) = self.pitch_line_points
        lats = np.linspace(lat1, lat2, 100)
        lons = np.linspace(lon1, lon2, 100)
        left, right, bottom, top = tuple(self.geotiff_extent)
        nrows, ncols = self.geotiff_data_array.shape
        rows = ((top - lats) / (top - bottom) * (nrows - 1)).clip(0, nrows - 1)
        cols = ((lons - left) / (right - left) * (ncols - 1)).clip(0, ncols - 1)
        elevations = []
        for r, c in zip(rows, cols):
            ir, ic = int(round(r)), int(round(c))
            elevations.append(self.geotiff_data_array[ir, ic])
        elevations = np.array(elevations)
        # Set to the median value below the pitch line (in absolute value)
        if np.any(~np.isnan(elevations)):
            median_val = np.nanmedian(elevations)
            mean_val = np.nanmean(elevations)
            min_val = np.nanmin(elevations)
            max_val = np.nanmax(elevations)
            valid_count = np.sum(~np.isnan(elevations))
            offset_value = abs(median_val)
            self.cal_line_offset_entry.setText(f"{offset_value:.2f}")
            # Log calculation details to activity log
            if hasattr(self, 'set_cal_info_text'):
                self.set_cal_info_text(
                    f"Heading Line Offset calculated: {offset_value:.2f} m (median of {valid_count} points along pitch line). "
                    f"Depth range: {abs(max_val):.2f} m (shallowest) to {abs(min_val):.2f} m (deepest). "
                    f"Mean depth: {abs(mean_val):.2f} m, Median depth: {abs(median_val):.2f} m.",
                    append=False
                )
        else:
            self.cal_line_offset_entry.setText("-")
            if hasattr(self, 'set_cal_info_text'):
                self.set_cal_info_text(
                    "Heading Line Offset: Could not calculate (no valid elevation data along pitch line).",
                    append=False
                )

    def _add_heading_lines_from_pitch_line(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot add heading lines.")
            return
        if self.geotiff_dataset_original is None:
            self._show_message("warning","No GeoTIFF", "Load a GeoTIFF first.")
            return
        if not hasattr(self, 'pitch_line_points') or len(self.pitch_line_points) != 2:
            self._show_message("warning","No Pitch Line", "Pick a pitch line first.")
            return
        try:
            offset_val = float(self.cal_line_offset_entry.text())
            if offset_val <= 0:
                raise ValueError
        except Exception:
            self._show_message("warning","Invalid Offset", "Line Offset must be a positive number.")
            return
        (lat1, lon1), (lat2, lon2) = self.pitch_line_points
        import pyproj
        geod = pyproj.Geod(ellps="WGS84")
        # Calculate azimuth of pitch line
        az12, az21, dist = geod.inv(lon1, lat1, lon2, lat2)
        # Perpendicular azimuths (north and south)
        perp_az_north = (az12 + 90) % 360
        perp_az_south = (az12 - 90) % 360
        # Offset both endpoints north
        n1_lon, n1_lat, _ = geod.fwd(lon1, lat1, perp_az_north, offset_val)
        n2_lon, n2_lat, _ = geod.fwd(lon2, lat2, perp_az_north, offset_val)
        # Offset both endpoints south
        s1_lon, s1_lat, _ = geod.fwd(lon1, lat1, perp_az_south, offset_val)
        s2_lon, s2_lat, _ = geod.fwd(lon2, lat2, perp_az_south, offset_val)
        # Store heading lines
        self.heading_lines = [
            [(n1_lat, n1_lon), (n2_lat, n2_lon)],
            [(s1_lat, s1_lon), (s2_lat, s2_lon)]
        ]
        self._plot_survey_plan(preserve_view_limits=True)
        # Make "Add Heading Lines" button normal after clicking
        if hasattr(self, 'add_heading_lines_btn'):
            self.add_heading_lines_btn.setStyleSheet("")  # Reset to normal
        # Show info in the info/error box instead of a dialog
        if hasattr(self, 'set_cal_info_text'):
            self.set_cal_info_text("Heading lines have been added north and south of the pitch line.")

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
        
        self.ax.set_xlim(min_lon - buffer_lon, max_lon + buffer_lon)
        self.ax.set_ylim(min_lat - buffer_lat, max_lat + buffer_lat)
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
        self.ax.set_xlim(min_lon - buffer_lon, max_lon + buffer_lon)
        self.ax.set_ylim(min_lat - buffer_lat, max_lat + buffer_lat)
        self.canvas.draw_idle()


    def _update_cal_export_name_from_pitch_line(self):
        # Only run if we have a valid pitch line and line offset
        if (not hasattr(self, 'pitch_line_points') or len(self.pitch_line_points) != 2):
            return
        try:
            offset_val = float(self.cal_line_offset_entry.text())
        except Exception:
            return
        (lat1, lon1), (lat2, lon2) = self.pitch_line_points
        try:
            import pyproj
            geod = pyproj.Geod(ellps="WGS84")
            az12, az21, dist = geod.inv(lon1, lat1, lon2, lat2)
            heading = int(round(az12)) % 360
        except Exception:
            heading = 0
        export_name = f"Cal_{int(round(offset_val))}m_{heading}deg"
        self.cal_export_name_entry.setText(export_name)

    def _export_cal_survey_files(self):
        import csv, os, json
        try:
            import fiona
            from shapely.geometry import LineString, mapping
        except ImportError:
            self._show_message("warning","Missing Libraries", "fiona and shapely are required for shapefile export.")
            return
        # Gather lines: Pitch, Heading1, Heading2, Roll (if present)
        lines = []
        line_num = 1
        # Pitch line
        if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2:
            lines.append((line_num, 'Pitch', self.pitch_line_points))
            line_num += 1
        # Heading lines
        if hasattr(self, 'heading_lines') and len(self.heading_lines) == 2:
            lines.append((line_num, 'Heading1', self.heading_lines[0]))
            line_num += 1
            lines.append((line_num, 'Heading2', self.heading_lines[1]))
            line_num += 1
        # Roll line
        if hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2:
            lines.append((line_num, 'Roll', self.roll_line_points))
            line_num += 1
        if not lines:
            self._show_message("warning","No Data", "No calibration lines to export. Define at least one line.")
            return
        export_name = self.cal_export_name_entry.text().strip() or "calibration_survey"
        export_dir = QFileDialog.getExistingDirectory(self, "Select Export Directory", self.last_export_dir)
        if not export_dir:
            return
        self.last_export_dir = export_dir
        self._save_last_export_dir()
        try:
            # --- Export to CSV ---
            # CSV files for decimal degrees always use _DD suffix
            csv_file_path = os.path.join(export_dir, f"{export_name}_DD.csv")
            with open(csv_file_path, 'w', newline='', encoding='utf-8') as csvfile:
                csv_writer = csv.writer(csvfile)
                csv_writer.writerow(['Line Number', 'Point Label', 'Line Name', 'Latitude', 'Longitude'])
                for num, name, pts in lines:
                    # Determine point labels based on line name
                    if name.lower().startswith('pitch'):
                        start_label, end_label = 'PLS', 'PLE'
                    elif name.lower().startswith('roll'):
                        start_label, end_label = 'RLS', 'RLE'
                    elif name.lower().startswith('heading'):
                        # Extract heading number if present
                        import re
                        m = re.search(r'(\d+)', name)
                        n = m.group(1) if m else ''
                        start_label, end_label = f'H{n}S', f'H{n}E'
                    else:
                        start_label, end_label = 'START', 'END'
                    csv_writer.writerow([num, start_label, name, pts[0][0], pts[0][1]])
                    csv_writer.writerow([num, end_label, name, pts[1][0], pts[1][1]])
            
            # --- Export to DDM format (Decimal Minutes) ---
            ddm_file_path = os.path.join(export_dir, f"{export_name}_DM.csv")
            with open(ddm_file_path, 'w', newline='', encoding='utf-8') as ddmfile:
                ddm_writer = csv.writer(ddmfile)
                ddm_writer.writerow(['Line Number', 'Point Label', 'Latitude with Decimal Minutes', 'Longitude with Decimal Minutes'])
                for num, name, pts in lines:
                    # Determine point labels based on line name (same logic as CSV export)
                    if name.lower().startswith('pitch'):
                        start_label, end_label = 'PLS', 'PLE'
                    elif name.lower().startswith('roll'):
                        start_label, end_label = 'RLS', 'RLE'
                    elif name.lower().startswith('heading'):
                        # Extract heading number if present
                        import re
                        m = re.search(r'(\d+)', name)
                        n = m.group(1) if m else ''
                        start_label, end_label = f'H{n}S', f'H{n}E'
                    else:
                        start_label, end_label = 'START', 'END'
                    # Convert to DDM format
                    start_lat_ddm = self._decimal_degrees_to_ddm(pts[0][0], is_latitude=True)
                    start_lon_ddm = self._decimal_degrees_to_ddm(pts[0][1], is_latitude=False)
                    end_lat_ddm = self._decimal_degrees_to_ddm(pts[1][0], is_latitude=True)
                    end_lon_ddm = self._decimal_degrees_to_ddm(pts[1][1], is_latitude=False)
                    ddm_writer.writerow([num, start_label, start_lat_ddm, start_lon_ddm])
                    ddm_writer.writerow([num, end_label, end_lat_ddm, end_lon_ddm])
            
            # --- Export to DDM text format (Decimal Minutes) ---
            ddm_txt_file_path = os.path.join(export_dir, f"{export_name}_DM.txt")
            with open(ddm_txt_file_path, 'w', encoding='utf-8') as ddm_txt_file:
                for num, name, pts in lines:
                    # Determine point labels based on line name (same logic as CSV export)
                    if name.lower().startswith('pitch'):
                        start_label, end_label = 'PLS', 'PLE'
                    elif name.lower().startswith('roll'):
                        start_label, end_label = 'RLS', 'RLE'
                    elif name.lower().startswith('heading'):
                        # Extract heading number if present
                        import re
                        m = re.search(r'(\d+)', name)
                        n = m.group(1) if m else ''
                        start_label, end_label = f'H{n}S', f'H{n}E'
                    else:
                        start_label, end_label = 'START', 'END'
                    # Convert to DDM format
                    start_lat_ddm = self._decimal_degrees_to_ddm(pts[0][0], is_latitude=True)
                    start_lon_ddm = self._decimal_degrees_to_ddm(pts[0][1], is_latitude=False)
                    end_lat_ddm = self._decimal_degrees_to_ddm(pts[1][0], is_latitude=True)
                    end_lon_ddm = self._decimal_degrees_to_ddm(pts[1][1], is_latitude=False)
                    ddm_txt_file.write(f"{start_label}, {start_lat_ddm}, {start_lon_ddm}\n")
                    ddm_txt_file.write(f"{end_label}, {end_lat_ddm}, {end_lon_ddm}\n")
            
            # --- Export to ESRI Shapefile (.shp) ---
            schema = {
                'geometry': 'LineString',
                'properties': {'line_num': 'int', 'line_name': 'str'},
            }
            crs_epsg = 'EPSG:4326'  # WGS 84
            features = []
            for num, name, pts in lines:
                shapely_line = LineString([(p[1], p[0]) for p in pts])
                features.append({
                    'geometry': mapping(shapely_line),
                    'properties': {'line_num': num, 'line_name': name},
                })
            shapefile_path = os.path.join(export_dir, f"{export_name}.shp")
            with fiona.open(shapefile_path, 'w', driver='ESRI Shapefile', crs=crs_epsg, schema=schema) as collection:
                collection.writerecords(features)
            # --- Export to GeoJSON ---
            geojson_file_path = os.path.join(export_dir, f"{export_name}.geojson")
            geojson_features = []
            for num, name, pts in lines:
                geojson_features.append({
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": [[pts[0][1], pts[0][0]], [pts[1][1], pts[1][0]]]
                    },
                    "properties": {
                        "line_num": num,
                        "line_name": name,
                        "points": [
                            {"point_num": 1, "lat": pts[0][0], "lon": pts[0][1]},
                            {"point_num": 2, "lat": pts[1][0], "lon": pts[1][1]}
                        ]
                    }
                })
            geojson_collection = {
                "type": "FeatureCollection",
                "features": geojson_features
            }
            with open(geojson_file_path, 'w', encoding='utf-8') as f:
                json.dump(geojson_collection, f, indent=2)
            # Will be updated below with all files

            # --- Export to Hypack LNW format ---
            lnw_file_path = os.path.join(export_dir, f"{export_name}.lnw")
            with open(lnw_file_path, 'w') as f:
                f.write("LNW 1.0\n")
                # Export each calibration line as LNW waypoints
                for num, name, pts in lines:
                    # Get speed from cal_survey_speed_entry if available, else default to 8 knots
                    try:
                        speed_knots = float(self.cal_survey_speed_entry.text()) if self.cal_survey_speed_entry.text() else 8.0
                    except:
                        speed_knots = 8.0
                    
                    # Get depth at each point if GeoTIFF is loaded
                    depth1 = self._get_depth_at_point(pts[0][0], pts[0][1])
                    depth2 = self._get_depth_at_point(pts[1][0], pts[1][1])
                    
                    # Start point
                    f.write(f"{name}_001, {pts[0][0]:.6f}, {pts[0][1]:.6f}, {depth1:.1f}, {speed_knots:.1f}, 50.0, {num}, {num}\n")
                    # End point
                    f.write(f"{name}_002, {pts[1][0]:.6f}, {pts[1][1]:.6f}, {depth2:.1f}, {speed_knots:.1f}, 50.0, {num}, {num}\n")

            # Will be updated below with all files

            # --- Export to SIS ASCII Plan format ---
            sis_file_path = os.path.join(export_dir, f"{export_name}.asciiplan")
            with open(sis_file_path, 'w') as f:
                f.write("SIS ASCII Plan\n")
                # Export each calibration line as SIS waypoints
                for num, name, pts in lines:
                    # Get speed from cal_survey_speed_entry if available, else default to 8 knots
                    try:
                        speed_knots = float(self.cal_survey_speed_entry.text()) if self.cal_survey_speed_entry.text() else 8.0
                    except:
                        speed_knots = 8.0
                    
                    # Get depth at each point if GeoTIFF is loaded
                    depth1 = self._get_depth_at_point(pts[0][0], pts[0][1])
                    depth2 = self._get_depth_at_point(pts[1][0], pts[1][1])
                    
                    # Start point
                    f.write(f"{name}_001, {pts[0][0]:.6f}, {pts[0][1]:.6f}, {depth1:.1f}, {speed_knots:.1f}, {num}, {num}\n")
                    # End point
                    f.write(f"{name}_002, {pts[1][0]:.6f}, {pts[1][1]:.6f}, {depth2:.1f}, {speed_knots:.1f}, {num}, {num}\n")

            # --- Export to Text format (point label, latitude, longitude) ---
            txt_file_path = os.path.join(export_dir, f"{export_name}_DD.txt")
            with open(txt_file_path, 'w', encoding='utf-8') as f:
                for num, name, pts in lines:
                    # Determine point labels based on line name
                    if name.lower().startswith('pitch'):
                        start_label, end_label = 'PLS', 'PLE'
                    elif name.lower().startswith('roll'):
                        start_label, end_label = 'RLS', 'RLE'
                    elif name.lower().startswith('heading'):
                        # Extract heading number if present
                        import re
                        m = re.search(r'(\d+)', name)
                        n = m.group(1) if m else ''
                        start_label, end_label = f'H{n}S', f'H{n}E'
                    else:
                        start_label, end_label = 'START', 'END'
                    f.write(f"{start_label} {pts[0][0]:.6f} {pts[0][1]:.6f}\n")
                    f.write(f"{end_label} {pts[1][0]:.6f} {pts[1][1]:.6f}\n")

            # Update success message to include SIS and text file
            self.set_cal_info_text(
                f"Survey exported successfully to:\n"
                f"- {os.path.basename(csv_file_path)}\n"
                f"- {os.path.basename(shapefile_path)} (and associated files)\n"
                f"- {os.path.basename(geojson_file_path)}\n"
                f"- {os.path.basename(lnw_file_path)}\n"
                f"- {os.path.basename(sis_file_path)}\n"
                f"- {os.path.basename(txt_file_path)}\n"
                f"in directory: {export_dir}", append=False)
            
            # --- Export comprehensive statistics file ---
            stats = self._calculate_calibration_survey_statistics()
            stats_file_path = None
            if stats:
                stats_file_path = os.path.join(export_dir, f"{export_name}_stats.txt")
                with open(stats_file_path, 'w', encoding='utf-8') as f:
                    f.write("COMPREHENSIVE CALIBRATION SURVEY STATISTICS\n")
                    f.write("=" * 50 + "\n\n")
                    f.write(f"Calibration Survey: {export_name}\n")
                    f.write(f"Export Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                    
                    # Survey parameters
                    try:
                        speed_knots = float(self.cal_survey_speed_entry.text()) if self.cal_survey_speed_entry.text() else 8.0
                        f.write(f"Survey Speed: {speed_knots} knots\n\n")
                    except:
                        f.write("Survey Speed: 8.0 knots (default)\n\n")
                    
                    # Line distances and times
                    f.write("SURVEY LINE DISTANCES AND TIMES\n")
                    f.write("-" * 35 + "\n")
                    
                    if stats['pitch_line_distance_m'] > 0:
                        f.write(f"Pitch Line (2 passes):\n")
                        # Add pitch line heading information
                        if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2:
                            try:
                                if pyproj is not None:
                                    geod = pyproj.Geod(ellps="WGS84")
                                    (lat1, lon1), (lat2, lon2) = self.pitch_line_points
                                    fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)
                                    
                                    pitch_heading = fwd_az % 360
                                    pitch_reciprocal_heading = back_az % 360
                                    
                                    f.write(f"  Heading: {pitch_heading:.1f}°\n")
                                    f.write(f"  Reciprocal Heading: {pitch_reciprocal_heading:.1f}°\n")
                                else:
                                    f.write(f"  Heading: pyproj not available\n")
                            except Exception as e:
                                f.write(f"  Heading: Unable to calculate\n")
                        f.write(f"  Distance (per pass): {stats['pitch_line_distance_m']:.1f} m ({stats['pitch_line_distance_km']:.3f} km, {stats['pitch_line_distance_nm']:.3f} nm)\n")
                        f.write(f"  Time (per line): {stats['pitch_line_time_min']/2:.1f} min\n")
                        f.write(f"  Time (total): {stats['pitch_line_time_min']:.1f} min\n\n")
                    
                    if stats['roll_line_distance_m'] > 0:
                        f.write(f"Roll Line (2 passes):\n")
                        # Add roll line heading information
                        if hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2:
                            try:
                                if pyproj is not None:
                                    geod = pyproj.Geod(ellps="WGS84")
                                    (lat1, lon1), (lat2, lon2) = self.roll_line_points
                                    fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)
                                    
                                    roll_heading = fwd_az % 360
                                    roll_reciprocal_heading = back_az % 360
                                    
                                    f.write(f"  Heading: {roll_heading:.1f}°\n")
                                    f.write(f"  Reciprocal Heading: {roll_reciprocal_heading:.1f}°\n")
                                else:
                                    f.write(f"  Heading: pyproj not available\n")
                            except Exception as e:
                                f.write(f"  Heading: Unable to calculate\n")
                        f.write(f"  Distance (per pass): {stats['roll_line_distance_m']:.1f} m ({stats['roll_line_distance_km']:.3f} km, {stats['roll_line_distance_nm']:.3f} nm)\n")
                        f.write(f"  Time (per line): {stats['roll_line_time_min']/2:.1f} min\n")
                        f.write(f"  Time (total): {stats['roll_line_time_min']:.1f} min\n\n")
                    
                    if stats['heading1_distance_m'] > 0:
                        f.write(f"Heading Line 1 (1 pass):\n")
                        # Add heading line 1 heading information
                        if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 1:
                            try:
                                if pyproj is not None:
                                    geod = pyproj.Geod(ellps="WGS84")
                                    (lat1, lon1), (lat2, lon2) = self.heading_lines[0]
                                    fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)
                                    
                                    heading1_heading = fwd_az % 360
                                    heading1_reciprocal_heading = back_az % 360
                                    
                                    f.write(f"  Heading: {heading1_heading:.1f}°\n")
                                    f.write(f"  Reciprocal Heading: {heading1_reciprocal_heading:.1f}°\n")
                                else:
                                    f.write(f"  Heading: pyproj not available\n")
                            except Exception as e:
                                f.write(f"  Heading: Unable to calculate\n")
                        f.write(f"  Distance: {stats['heading1_distance_m']:.1f} m ({stats['heading1_distance_km']:.3f} km, {stats['heading1_distance_nm']:.3f} nm)\n")
                        f.write(f"  Time: {stats['heading1_time_min']:.1f} min\n\n")
                    
                    if stats['heading2_distance_m'] > 0:
                        f.write(f"Heading Line 2 (1 pass):\n")
                        # Add heading line 2 heading information
                        if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 2:
                            try:
                                if pyproj is not None:
                                    geod = pyproj.Geod(ellps="WGS84")
                                    (lat1, lon1), (lat2, lon2) = self.heading_lines[1]
                                    fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)
                                    
                                    heading2_heading = fwd_az % 360
                                    heading2_reciprocal_heading = back_az % 360
                                    
                                    f.write(f"  Heading: {heading2_heading:.1f}°\n")
                                    f.write(f"  Reciprocal Heading: {heading2_reciprocal_heading:.1f}°\n")
                                else:
                                    f.write(f"  Heading: pyproj not available\n")
                            except Exception as e:
                                f.write(f"  Heading: Unable to calculate\n")
                        f.write(f"  Distance: {stats['heading2_distance_m']:.1f} m ({stats['heading2_distance_km']:.3f} km, {stats['heading2_distance_nm']:.3f} nm)\n")
                        f.write(f"  Time: {stats['heading2_time_min']:.1f} min\n\n")
                    
                    # Travel distances and times
                    f.write("TRAVEL DISTANCES AND TIMES\n")
                    f.write("-" * 30 + "\n")
                    
                    # Combined line execution times and distances
                    if stats['pitch_line_distance_m'] > 0:
                        f.write(f"Pitch Lines (2 passes):\n")
                        f.write(f"  Total Distance: {stats['pitch_line_distance_m'] * 2:.1f} m ({stats['pitch_line_distance_km'] * 2:.3f} km, {stats['pitch_line_distance_nm'] * 2:.3f} nm)\n")
                        f.write(f"  Total Time: {stats['pitch_line_time_min']:.1f} min\n\n")
                    
                    if stats['roll_line_distance_m'] > 0:
                        f.write(f"Roll Lines (2 passes):\n")
                        f.write(f"  Total Distance: {stats['roll_line_distance_m'] * 2:.1f} m ({stats['roll_line_distance_km'] * 2:.3f} km, {stats['roll_line_distance_nm'] * 2:.3f} nm)\n")
                        f.write(f"  Total Time: {stats['roll_line_time_min']:.1f} min\n\n")
                    
                    if stats['heading1_distance_m'] > 0 or stats['heading2_distance_m'] > 0:
                        total_heading_distance = stats['heading1_distance_m'] + stats['heading2_distance_m']
                        total_heading_time = stats['heading1_time_min'] + stats['heading2_time_min']
                        total_heading_distance_km = total_heading_distance / 1000
                        total_heading_distance_nm = total_heading_distance / 1852
                        
                        f.write(f"Heading Lines (2 lines):\n")
                        f.write(f"  Total Distance: {total_heading_distance:.1f} m ({total_heading_distance_km:.3f} km, {total_heading_distance_nm:.3f} nm)\n")
                        f.write(f"  Total Time: {total_heading_time:.1f} min\n\n")
                    
                    if stats['travel_pitch_to_roll_m'] > 0:
                        f.write(f"Pitch Line Start → Roll Line Start:\n")
                        f.write(f"  Distance: {stats['travel_pitch_to_roll_m']:.1f} m ({stats['travel_pitch_to_roll_km']:.3f} km, {stats['travel_pitch_to_roll_nm']:.3f} nm)\n")
                        f.write(f"  Time: {stats['travel_pitch_to_roll_min']:.1f} min\n\n")
                    
                    if stats['travel_roll_to_heading1_m'] > 0:
                        f.write(f"Roll Line Start → Heading Line 1 Start:\n")
                        f.write(f"  Distance: {stats['travel_roll_to_heading1_m']:.1f} m ({stats['travel_roll_to_heading1_km']:.3f} km, {stats['travel_roll_to_heading1_nm']:.3f} nm)\n")
                        f.write(f"  Time: {stats['travel_roll_to_heading1_min']:.1f} min\n\n")
                    
                    if stats['travel_heading1_to_heading2_m'] > 0:
                        f.write(f"Heading Line 1 End → Heading Line 2 Start:\n")
                        f.write(f"  Distance: {stats['travel_heading1_to_heading2_m']:.1f} m ({stats['travel_heading1_to_heading2_km']:.3f} km, {stats['travel_heading1_to_heading2_nm']:.3f} nm)\n")
                        f.write(f"  Time: {stats['travel_heading1_to_heading2_min']:.1f} min\n\n")
                    
                    # Totals
                    f.write("TOTAL SURVEY SUMMARY\n")
                    f.write("-" * 25 + "\n")
                    f.write(f"Total Distance: {stats['total_distance_m']:.1f} m ({stats['total_distance_km']:.3f} km, {stats['total_distance_nm']:.3f} nm)\n")
                    f.write(f"Total Time: {stats['total_time_min']:.1f} min ({stats['total_time_hours']:.2f} hr)\n\n")
                    
                    # Survey pattern explanation
                    f.write("SURVEY PATTERN\n")
                    f.write("-" * 15 + "\n")
                    f.write("1. Run pitch line (2 passes)\n")
                    f.write("2. Travel from pitch line start to roll line start\n")
                    f.write("3. Run roll line (2 passes)\n")
                    f.write("4. Travel from roll line start to heading line 1 start\n")
                    f.write("5. Run heading line 1 (1 pass)\n")
                    f.write("6. Travel from heading line 1 end to heading line 2 start\n")
                    f.write("7. Run heading line 2 (1 pass)\n\n")
                    
                    # Line details
                    f.write("SURVEY LINE DETAILS\n")
                    f.write("-" * 20 + "\n")
                    f.write(f"Number of Calibration Lines: {len(lines)}\n")
                    for num, name, pts in lines:
                        f.write(f"Line {num}: {name} - Start: ({pts[0][0]:.6f}, {pts[0][1]:.6f}), End: ({pts[1][0]:.6f}, {pts[1][1]:.6f})\n")
                    f.write("\n")
                    
                    # Waypoints with labels
                    f.write("WAYPOINTS\n")
                    f.write("-" * 10 + "\n")
                    if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2:
                        pitch_start = self.pitch_line_points[0]
                        pitch_end = self.pitch_line_points[1]
                        pitch_start_lat_ddm = self._decimal_degrees_to_ddm(pitch_start[0], is_latitude=True)
                        pitch_start_lon_ddm = self._decimal_degrees_to_ddm(pitch_start[1], is_latitude=False)
                        pitch_end_lat_ddm = self._decimal_degrees_to_ddm(pitch_end[0], is_latitude=True)
                        pitch_end_lon_ddm = self._decimal_degrees_to_ddm(pitch_end[1], is_latitude=False)
                        f.write(f"Pitch Start: {pitch_start_lat_ddm}, {pitch_start_lon_ddm} ({pitch_start[0]:.6f}, {pitch_start[1]:.6f})\n")
                        f.write(f"Pitch End: {pitch_end_lat_ddm}, {pitch_end_lon_ddm} ({pitch_end[0]:.6f}, {pitch_end[1]:.6f})\n")
                    if hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2:
                        roll_start = self.roll_line_points[0]
                        roll_end = self.roll_line_points[1]
                        roll_start_lat_ddm = self._decimal_degrees_to_ddm(roll_start[0], is_latitude=True)
                        roll_start_lon_ddm = self._decimal_degrees_to_ddm(roll_start[1], is_latitude=False)
                        roll_end_lat_ddm = self._decimal_degrees_to_ddm(roll_end[0], is_latitude=True)
                        roll_end_lon_ddm = self._decimal_degrees_to_ddm(roll_end[1], is_latitude=False)
                        f.write(f"Roll Start: {roll_start_lat_ddm}, {roll_start_lon_ddm} ({roll_start[0]:.6f}, {roll_start[1]:.6f})\n")
                        f.write(f"Roll End: {roll_end_lat_ddm}, {roll_end_lon_ddm} ({roll_end[0]:.6f}, {roll_end[1]:.6f})\n")
                    if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 1 and len(self.heading_lines[0]) == 2:
                        heading1_start = self.heading_lines[0][0]
                        heading1_end = self.heading_lines[0][1]
                        heading1_start_lat_ddm = self._decimal_degrees_to_ddm(heading1_start[0], is_latitude=True)
                        heading1_start_lon_ddm = self._decimal_degrees_to_ddm(heading1_start[1], is_latitude=False)
                        heading1_end_lat_ddm = self._decimal_degrees_to_ddm(heading1_end[0], is_latitude=True)
                        heading1_end_lon_ddm = self._decimal_degrees_to_ddm(heading1_end[1], is_latitude=False)
                        f.write(f"Heading 1 Start: {heading1_start_lat_ddm}, {heading1_start_lon_ddm} ({heading1_start[0]:.6f}, {heading1_start[1]:.6f})\n")
                        f.write(f"Heading 1 End: {heading1_end_lat_ddm}, {heading1_end_lon_ddm} ({heading1_end[0]:.6f}, {heading1_end[1]:.6f})\n")
                    if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 2 and len(self.heading_lines[1]) == 2:
                        heading2_start = self.heading_lines[1][0]
                        heading2_end = self.heading_lines[1][1]
                        heading2_start_lat_ddm = self._decimal_degrees_to_ddm(heading2_start[0], is_latitude=True)
                        heading2_start_lon_ddm = self._decimal_degrees_to_ddm(heading2_start[1], is_latitude=False)
                        heading2_end_lat_ddm = self._decimal_degrees_to_ddm(heading2_end[0], is_latitude=True)
                        heading2_end_lon_ddm = self._decimal_degrees_to_ddm(heading2_end[1], is_latitude=False)
                        f.write(f"Heading 2 Start: {heading2_start_lat_ddm}, {heading2_start_lon_ddm} ({heading2_start[0]:.6f}, {heading2_start[1]:.6f})\n")
                        f.write(f"Heading 2 End: {heading2_end_lat_ddm}, {heading2_end_lon_ddm} ({heading2_end[0]:.6f}, {heading2_end[1]:.6f})\n")
                
                # Will build complete message at end with all files
            
            # --- Export parameters metadata as JSON ---
            json_metadata_path = os.path.join(export_dir, f"{export_name}_params.json")
            try:
                # Get current parameter values
                params = {}
                try:
                    params['line_offset'] = float(self.cal_line_offset_entry.get())
                except:
                    params['line_offset'] = None
                
                try:
                    params['export_name'] = self.cal_export_name_entry.get().strip()
                except:
                    params['export_name'] = export_name
                
                try:
                    params['survey_speed'] = float(self.cal_survey_speed_entry.get())
                except:
                    params['survey_speed'] = 8.0  # Default
                
                # Save metadata
                with open(json_metadata_path, 'w', encoding='utf-8') as f:
                    json.dump(params, f, indent=2)
            except Exception as e:
                # If metadata export fails, continue without it
                print(f"Warning: Could not export metadata: {e}")
                json_metadata_path = None
            
            # --- Export PNG images ---
            # Export main map
            map_png_path = os.path.join(export_dir, f"{export_name}_map.png")
            self.figure.savefig(map_png_path, dpi=300, bbox_inches='tight', facecolor='white')
            
            # Export profile if it exists
            profile_png_path = None
            if hasattr(self, 'profile_fig') and self.profile_fig is not None:
                profile_png_path = os.path.join(export_dir, f"{export_name}_profile.png")
                self.profile_fig.savefig(profile_png_path, dpi=300, bbox_inches='tight', facecolor='white')
            
            # Build complete success message with all exported files
            success_msg = f"Survey exported successfully to:\n"
            success_msg += f"- {os.path.basename(csv_file_path)}\n"
            success_msg += f"- {os.path.basename(shapefile_path)} (and associated files)\n"
            success_msg += f"- {os.path.basename(geojson_file_path)}\n"
            success_msg += f"- {os.path.basename(lnw_file_path)}\n"
            success_msg += f"- {os.path.basename(sis_file_path)}\n"
            if stats_file_path:
                success_msg += f"- {os.path.basename(stats_file_path)}\n"
            if json_metadata_path:
                success_msg += f"- {os.path.basename(json_metadata_path)}\n"
            success_msg += f"- {os.path.basename(map_png_path)}\n"
            if profile_png_path:
                success_msg += f"- {os.path.basename(profile_png_path)}\n"
            success_msg += f"in directory: {export_dir}"
            self.set_cal_info_text(success_msg, append=False)
        except Exception as e:
            self._show_message("error","Export Error", f"Failed to export calibration survey files: {e}")

    def _import_cal_survey_files(self):
        """Import calibration survey lines from CSV or GeoJSON file."""
        import csv
        import json
        
        # Open file dialog to select import file
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Survey File to Import",
            self.last_cal_import_dir,
            "Decimal Degree CSV files (*_DD.csv);;CSV files (*.csv);;GeoJSON files (*.geojson);;JSON files (*.json);;All files (*.*)"
        )
        
        if not file_path:
            return
        
        # Save the directory for next time
        import_dir = os.path.dirname(file_path)
        if import_dir and os.path.isdir(import_dir):
            self.last_cal_import_dir = import_dir
            self._save_last_cal_import_dir()
        
        try:
            # Clear existing lines
            self.pitch_line_points = []
            self.heading_lines = []
            self.roll_line_points = []
            
            # Determine file type and import accordingly
            file_ext = os.path.splitext(file_path)[1].lower()
            
            if file_ext == '.csv':
                # Import from CSV
                with open(file_path, 'r', encoding='utf-8') as csvfile:
                    csv_reader = csv.DictReader(csvfile)
                    
                    # Dictionary to collect points for each line
                    lines_data = {}
                    
                    for row in csv_reader:
                        line_name = row.get('Line Name', '').strip()
                        point_label = row.get('Point Label', '').strip()
                        try:
                            lat = float(row.get('Latitude', 0))
                            lon = float(row.get('Longitude', 0))
                        except (ValueError, TypeError):
                            continue
                        
                        if not line_name:
                            continue
                        
                        # Initialize line if not exists
                        if line_name not in lines_data:
                            lines_data[line_name] = {}
                        
                        # Store point by label
                        lines_data[line_name][point_label] = (lat, lon)
                    
                    # Process collected lines
                    for line_name, points in lines_data.items():
                        line_name_lower = line_name.lower()
                        
                        if line_name_lower.startswith('pitch'):
                            # Pitch line: expect PLS and PLE
                            if 'PLS' in points and 'PLE' in points:
                                self.pitch_line_points = [points['PLS'], points['PLE']]
                        
                        elif line_name_lower.startswith('roll'):
                            # Roll line: expect RLS and RLE
                            if 'RLS' in points and 'RLE' in points:
                                self.roll_line_points = [points['RLS'], points['RLE']]
                        
                        elif line_name_lower.startswith('heading'):
                            # Heading lines: expect H1S/H1E, H2S/H2E, etc.
                            import re
                            match = re.search(r'heading\s*(\d+)', line_name_lower)
                            if match:
                                heading_num = int(match.group(1))
                                start_label = f'H{heading_num}S'
                                end_label = f'H{heading_num}E'
                            else:
                                # Try to extract from point labels
                                start_labels = [p for p in points.keys() if p.endswith('S')]
                                end_labels = [p for p in points.keys() if p.endswith('E')]
                                if start_labels and end_labels:
                                    start_label = start_labels[0]
                                    end_label = end_labels[0]
                                else:
                                    continue
                            
                            if start_label in points and end_label in points:
                                # Ensure heading_lines list is large enough
                                while len(self.heading_lines) < 2:
                                    self.heading_lines.append([])
                                # Heading lines are indexed 0 and 1
                                heading_idx = min(heading_num - 1, 1) if match else 0
                                if heading_idx < 2:
                                    self.heading_lines[heading_idx] = [points[start_label], points[end_label]]
            
            elif file_ext in ['.geojson', '.json']:
                # Import from GeoJSON
                with open(file_path, 'r', encoding='utf-8') as f:
                    geojson_data = json.load(f)
                
                # Ensure geojson_data is a dictionary
                if not isinstance(geojson_data, dict):
                    self._show_message("error","Import Error", "Invalid GeoJSON format: root must be an object")
                    return
                
                # Handle FeatureCollection
                if geojson_data.get('type') == 'FeatureCollection':
                    features = geojson_data.get('features', [])
                    if not isinstance(features, list):
                        features = []
                elif geojson_data.get('type') == 'Feature':
                    features = [geojson_data]
                else:
                    features = []
                
                for feature in features:
                    # Ensure feature is a dictionary
                    if not isinstance(feature, dict):
                        continue
                    if feature.get('type') != 'Feature':
                        continue
                    
                    geometry = feature.get('geometry', {})
                    if not isinstance(geometry, dict):
                        continue
                    if geometry.get('type') != 'LineString':
                        continue
                    
                    coordinates = geometry.get('coordinates', [])
                    if not isinstance(coordinates, list) or len(coordinates) != 2:
                        continue
                    
                    # GeoJSON uses [lon, lat] format, we need [lat, lon]
                    point1 = (coordinates[0][1], coordinates[0][0])
                    point2 = (coordinates[1][1], coordinates[1][0])
                    
                    properties = feature.get('properties', {})
                    if not isinstance(properties, dict):
                        properties = {}
                    line_name = str(properties.get('line_name', '')).strip()
                    line_name_lower = line_name.lower()
                    
                    if line_name_lower.startswith('pitch'):
                        self.pitch_line_points = [point1, point2]
                    elif line_name_lower.startswith('roll'):
                        self.roll_line_points = [point1, point2]
                    elif line_name_lower.startswith('heading'):
                        import re
                        match = re.search(r'heading\s*(\d+)', line_name_lower)
                        heading_idx = 0
                        if match:
                            heading_idx = min(int(match.group(1)) - 1, 1)
                        
                        # Ensure heading_lines list is large enough
                        while len(self.heading_lines) < 2:
                            self.heading_lines.append([])
                        if heading_idx < 2:
                            self.heading_lines[heading_idx] = [point1, point2]
            
            else:
                self._show_message("error","Import Error", f"Unsupported file format: {file_ext}")
                return
            
            # Try to load metadata file (params.json) if it exists
            params = None
            base_path = os.path.splitext(file_path)[0]
            # Try multiple possible locations for the metadata file
            metadata_paths = [
                base_path + "_params.json",  # Same directory, same base name
                os.path.join(os.path.dirname(file_path), os.path.splitext(os.path.basename(file_path))[0] + "_params.json"),
            ]
            for metadata_path in metadata_paths:
                if os.path.exists(metadata_path):
                    try:
                        with open(metadata_path, 'r', encoding='utf-8') as f:
                            params = json.load(f)
                        break
                    except Exception as e:
                        print(f"Warning: Could not read metadata file {metadata_path}: {e}")
            
            # Populate parameter fields from metadata if available
            if params and isinstance(params, dict):
                try:
                    if params.get('line_offset') is not None:
                        self.cal_line_offset_entry.setText(str(params['line_offset']))
                    
                    if params.get('export_name'):
                        self.cal_export_name_entry.setText(params['export_name'])
                    
                    if params.get('survey_speed') is not None:
                        self.cal_survey_speed_entry.setText(str(params['survey_speed']))
                except Exception as e:
                    print(f"Warning: Error populating parameter fields from metadata: {e}")
            
            # If line offset wasn't set from params, try to calculate it from pitch line
            if not params or not isinstance(params, dict) or params.get('line_offset') is None:
                self._update_cal_line_offset_from_pitch_line()
            
            # If export name wasn't set from params, generate it from pitch line
            if not params or not isinstance(params, dict) or not params.get('export_name'):
                self._update_cal_export_name_from_pitch_line()
            
            # Update UI and plot
            self._update_pitch_line_button_states()
            self._update_roll_line_button_states()
            # Enable heading line buttons if heading lines were imported
            if hasattr(self, 'add_heading_lines_btn'):
                has_heading_lines = any(len(line) == 2 for line in self.heading_lines)
                self.add_heading_lines_btn.setEnabled(has_heading_lines)
            # Reset button styles to normal
            if hasattr(self, 'pick_pitch_line_btn'):
                self.pick_pitch_line_btn.setStyleSheet("")
            if hasattr(self, 'add_heading_lines_btn'):
                self.add_heading_lines_btn.setStyleSheet("")
            if hasattr(self, 'pick_roll_line_btn'):
                self.pick_roll_line_btn.setStyleSheet("")
            self._update_cal_line_times()
            self._plot_survey_plan(preserve_view_limits=True)
            # Update profile window
            if hasattr(self, '_draw_current_profile'):
                self._draw_current_profile()
            elif hasattr(self, '_draw_pitch_line_profile'):
                self._draw_pitch_line_profile()
            
            # Show success message
            imported_lines = []
            if len(self.pitch_line_points) == 2:
                imported_lines.append("Pitch line")
            if len(self.roll_line_points) == 2:
                imported_lines.append("Roll line")
            for i, heading_line in enumerate(self.heading_lines):
                if len(heading_line) == 2:
                    imported_lines.append(f"Heading line {i+1}")
            
            if imported_lines:
                self.set_cal_info_text(f"Successfully imported: {', '.join(imported_lines)}")
            else:
                self._show_message("warning","Import Warning", "No valid calibration lines found in the selected file.")
        
        except Exception as e:
            self._show_message("error","Import Error", f"Failed to import calibration survey files: {e}")
            import traceback
            traceback.print_exc()

    def _on_draw_event_update_colormap(self, event=None):
        # Only update if elevation overlay is active and geotiff is loaded
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
            vmax_elev = np.nanmax(display_data) if not np.all(np.isnan(display_data)) else None
        # Update the image clim values
        if vmin_elev is not None and vmax_elev is not None and vmin_elev != vmax_elev:
            self.geotiff_image_plot.set_clim(vmin=vmin_elev, vmax=vmax_elev)
        # Update colorbar if present - just redraw it
        if hasattr(self, 'elevation_colorbar') and self.elevation_colorbar is not None:
            pass  # Remove the draw() call - colorbar will update with canvas.draw_idle()
        # Don't call draw_idle() here - it causes excessive redraws and performance issues
        # The canvas will redraw automatically when needed

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

    def _on_temp_line_motion(self, event):
        # Draw a temporary line from the start point to the current mouse position
        if not hasattr(self, '_temp_line') or not hasattr(self, '_temp_line_start'):
            return
        if event.inaxes != self.ax:
            return
        x0, y0 = self._temp_line_start[1], self._temp_line_start[0]
        x1, y1 = event.xdata, event.ydata
        if x1 is None or y1 is None:
            return
        self._temp_line.set_data([x0, x1], [y0, y1])
        # Calculate and display length, azimuth, and time
        try:
            import pyproj
            geod = pyproj.Geod(ellps="WGS84")
            az12, az21, dist = geod.inv(x0, y0, x1, y1)
        except Exception:
            import numpy as np
            dist = np.sqrt((x1 - x0)**2 + (y1 - y0)**2) * 111320  # rough meters
            az12 = np.degrees(np.arctan2(x1 - x0, y1 - y0)) % 360
        # Determine which speed entry to use
        speed_knots = None
        try:
            if hasattr(self, 'cal_survey_speed_entry') and self.cal_survey_speed_entry.winfo_ismapped():
                speed_knots = float(self.cal_survey_speed_entry.get())
            elif hasattr(self, 'survey_speed_entry'):
                speed_knots = float(self.survey_speed_entry.get())
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

    def _on_pitch_line_handle_pick(self, event):
        """Handle picking of pitch line handles for dragging."""
        if not self.edit_pitch_line_mode:
            return
        
        if event.mouseevent.button != 1:  # Only left mouse button
            return
        
        if event.artist == self.pitch_line_start_handle:
            self.dragging_pitch_line_handle = 'start'
        elif event.artist == self.pitch_line_end_handle:
            self.dragging_pitch_line_handle = 'end'

    def _on_pitch_line_handle_motion(self, event):
        """Handle dragging of pitch line handles."""
        if not self.edit_pitch_line_mode or not self.dragging_pitch_line_handle:
            return
        
        if event.inaxes != self.ax:
            return
        
        if event.xdata is None or event.ydata is None:
            return
        
        # Update the handle position
        if self.dragging_pitch_line_handle == 'start' and self.pitch_line_start_handle is not None:
            self.pitch_line_start_handle.set_offsets([[event.xdata, event.ydata]])
            self.pitch_line_points[0] = (event.ydata, event.xdata)
        elif self.dragging_pitch_line_handle == 'end' and self.pitch_line_end_handle is not None:
            self.pitch_line_end_handle.set_offsets([[event.xdata, event.ydata]])
            self.pitch_line_points[1] = (event.ydata, event.xdata)
        
        # Update the pitch line itself (draw a line between the two handles)
        if hasattr(self, 'pitch_line_edit_line') and self.pitch_line_edit_line is not None:
            self.pitch_line_edit_line.set_data(
                [self.pitch_line_points[0][1], self.pitch_line_points[1][1]],
                [self.pitch_line_points[0][0], self.pitch_line_points[1][0]]
            )
        else:
            (lat1, lon1), (lat2, lon2) = self.pitch_line_points
            self.pitch_line_edit_line, = self.ax.plot([lon1, lon2], [lat1, lat2], color='orange', linewidth=2, zorder=9)
        
        # Update tooltip if in pitch line mode
        if hasattr(self, 'pitch_line_tooltip') and self.pitch_line_tooltip is not None:
            depth = self._get_depth_at_point(event.ydata, event.xdata)
            slope = self._calculate_slope_at_point(event.ydata, event.xdata)[1]
            if depth is not None and slope is not None:
                tooltip_text = f"Depth: {depth:.1f} m\nSlope: {slope:.1f}°"
                if self.pitch_line_tooltip_text is not None:
                    self.pitch_line_tooltip_text.set_text(tooltip_text)
        
        self.canvas.draw_idle()

    def _on_pitch_line_handle_release(self, event):
        """Handle release of pitch line handles."""
        if not self.edit_pitch_line_mode:
            return
        
        if self.dragging_pitch_line_handle:
            # Remove the temp edit line and redraw the full plot to finalize
            if hasattr(self, 'pitch_line_edit_line') and self.pitch_line_edit_line is not None:
                self.pitch_line_edit_line.remove()
                self.pitch_line_edit_line = None
            self._plot_survey_plan(preserve_view_limits=True)
            
            # Recreate the handles after the plot redraw
            start_lat, start_lon = self.pitch_line_points[0]
            end_lat, end_lon = self.pitch_line_points[1]
            
            # Create handles as scatter points
            self.pitch_line_start_handle = self.ax.scatter([start_lon], [start_lat], 
                                                         color='red', s=100, marker='o', 
                                                         edgecolors='black', linewidth=2, 
                                                         zorder=10, picker=5)
            self.pitch_line_end_handle = self.ax.scatter([end_lon], [end_lat], 
                                                       color='blue', s=100, marker='o', 
                                                       edgecolors='black', linewidth=2, 
                                                       zorder=10, picker=5)
            
            # Update all dependent calculations
            self._update_cal_line_offset_from_pitch_line()
            self._draw_pitch_line_profile()
            self._update_cal_export_name_from_pitch_line()
            self._update_cal_line_times()
            
            # Report updated pitch line summary
            try:
                import pyproj
                geod = pyproj.Geod(ellps="WGS84")
                (lat1, lon1), (lat2, lon2) = self.pitch_line_points
                az12, az21, dist_m = geod.inv(lon1, lat1, lon2, lat2)
                dist_nm = dist_m / 1852.0
                speed_knots = float(self.cal_survey_speed_entry.get()) if self.cal_survey_speed_entry.get() else 8.0
                speed_m_per_h = speed_knots * 1852
                time_hours = dist_m / speed_m_per_h if speed_m_per_h > 0 else 0
                time_minutes = time_hours * 60
                summary = (
                    f"Updated Pitch Line:\n"
                    f"Length: {dist_m:.1f} m ({dist_nm:.3f} nm)\n"
                    f"Azimuth: {az12:.1f}°\n"
                    f"Survey Time: {time_minutes:.1f} min"
                )
                self.set_cal_info_text(summary)
            except Exception as e:
                self.set_cal_info_text(f"Error calculating updated pitch line summary: {e}")
            
            self.dragging_pitch_line_handle = None

    def _on_roll_line_handle_pick(self, event):
        """Handle picking of roll line handles for dragging."""
        if not self.edit_roll_line_mode:
            return
        
        if event.mouseevent.button != 1:  # Only left mouse button
            return
        
        if event.artist == self.roll_line_start_handle:
            self.dragging_roll_line_handle = 'start'
        elif event.artist == self.roll_line_end_handle:
            self.dragging_roll_line_handle = 'end'

    def _on_roll_line_handle_motion(self, event):
        """Handle dragging of roll line handles."""
        if not self.edit_roll_line_mode or not self.dragging_roll_line_handle:
            return
        
        if event.inaxes != self.ax:
            return
        
        if event.xdata is None or event.ydata is None:
            return
        
        # Update the handle position
        if self.dragging_roll_line_handle == 'start' and self.roll_line_start_handle is not None:
            self.roll_line_start_handle.set_offsets([[event.xdata, event.ydata]])
            self.roll_line_points[0] = (event.ydata, event.xdata)
        elif self.dragging_roll_line_handle == 'end' and self.roll_line_end_handle is not None:
            self.roll_line_end_handle.set_offsets([[event.xdata, event.ydata]])
            self.roll_line_points[1] = (event.ydata, event.xdata)
        
        # Update the roll line itself (draw a line between the two handles)
        if hasattr(self, 'roll_line_edit_line') and self.roll_line_edit_line is not None:
            self.roll_line_edit_line.set_data(
                [self.roll_line_points[0][1], self.roll_line_points[1][1]],
                [self.roll_line_points[0][0], self.roll_line_points[1][0]]
            )
        else:
            (lat1, lon1), (lat2, lon2) = self.roll_line_points
            self.roll_line_edit_line, = self.ax.plot([lon1, lon2], [lat1, lat2], color='purple', linewidth=2, zorder=9)
        
        self.canvas.draw_idle()

    def _on_roll_line_handle_release(self, event):
        """Handle release of roll line handles."""
        if not self.edit_roll_line_mode:
            return
        
        if self.dragging_roll_line_handle:
            # Remove the temp edit line and redraw the full plot to finalize
            if hasattr(self, 'roll_line_edit_line') and self.roll_line_edit_line is not None:
                self.roll_line_edit_line.remove()
                self.roll_line_edit_line = None
            self._plot_survey_plan(preserve_view_limits=True)
            
            # Recreate the handles after the plot redraw
            start_lat, start_lon = self.roll_line_points[0]
            end_lat, end_lon = self.roll_line_points[1]
            
            # Create handles as scatter points
            self.roll_line_start_handle = self.ax.scatter([start_lon], [start_lat], 
                                                        color='red', s=100, marker='o', 
                                                        edgecolors='black', linewidth=2, 
                                                        zorder=10, picker=5)
            self.roll_line_end_handle = self.ax.scatter([end_lon], [end_lat], 
                                                      color='blue', s=100, marker='o', 
                                                      edgecolors='black', linewidth=2, 
                                                      zorder=10, picker=5)
            
            # Update all dependent calculations
            self._update_cal_line_times()
            
            # Report updated roll line summary
            try:
                import pyproj
                geod = pyproj.Geod(ellps="WGS84")
                (lat1, lon1), (lat2, lon2) = self.roll_line_points
                az12, az21, dist_m = geod.inv(lon1, lat1, lon2, lat2)
                dist_nm = dist_m / 1852.0
                speed_knots = float(self.cal_survey_speed_entry.get()) if self.cal_survey_speed_entry.get() else 8.0
                speed_m_per_h = speed_knots * 1852
                time_hours = dist_m / speed_m_per_h if speed_m_per_h > 0 else 0
                time_minutes = time_hours * 60
                summary = (
                    f"Updated Roll Line:\n"
                    f"Length: {dist_m:.1f} m ({dist_nm:.3f} nm)\n"
                    f"Azimuth: {az12:.1f}°\n"
                    f"Survey Time: {time_minutes:.1f} min"
                )
                self.set_cal_info_text(summary)
            except Exception as e:
                self.set_cal_info_text(f"Error calculating updated roll line summary: {e}")
            
            self.dragging_roll_line_handle = None

    def _calculate_slope_at_point(self, lat, lon):
        """Calculate slope at a given lat/lon point."""
        if (self.geotiff_data_array is None or self.geotiff_extent is None):
            return None, None
        
        try:
            # Convert lat/lon to pixel coordinates
            left, right, bottom, top = tuple(self.geotiff_extent)
            nrows, ncols = self.geotiff_data_array.shape
            
            # Calculate pixel coordinates
            col = int((lon - left) / (right - left) * (ncols - 1))
            row = int((top - lat) / (top - bottom) * (nrows - 1))
            
            # Ensure we're within bounds
            if 0 <= row < nrows and 0 <= col < ncols:
                # Get elevation at this point
                elevation = self.geotiff_data_array[row, col]
                
                # Calculate slope using gradient
                center_lat_geotiff = (self.geotiff_extent[2] + self.geotiff_extent[3]) / 2
                m_per_deg_lat = 111320.0
                m_per_deg_lon = 111320.0 * np.cos(np.radians(center_lat_geotiff))
                
                res_lat_deg = (self.geotiff_extent[3] - self.geotiff_extent[2]) / self.geotiff_data_array.shape[0]
                res_lon_deg = (self.geotiff_extent[1] - self.geotiff_extent[0]) / self.geotiff_data_array.shape[1]
                
                dx_m = res_lon_deg * m_per_deg_lon
                dy_m = res_lat_deg * m_per_deg_lat
                
                # Get a small window around the point for gradient calculation
                window_size = 3
                r_start = max(0, row - window_size // 2)
                r_end = min(nrows, row + window_size // 2 + 1)
                c_start = max(0, col - window_size // 2)
                c_end = min(ncols, col + window_size // 2 + 1)
                
                window_data = self.geotiff_data_array[r_start:r_end, c_start:c_end]
                
                if window_data.size > 1 and not np.all(np.isnan(window_data)):
                    # Calculate gradient
                    dz_dy, dz_dx = np.gradient(window_data, dy_m, dx_m)
                    
                    # Calculate slope magnitude
                    slope_rad = np.arctan(np.sqrt(dz_dx**2 + dz_dy**2))
                    slope_degrees = np.degrees(slope_rad)
                    
                    # Return the center value
                    center_idx = window_size // 2
                    if center_idx < slope_degrees.shape[0] and center_idx < slope_degrees.shape[1]:
                        return elevation, slope_degrees[center_idx, center_idx]
                    else:
                        return elevation, slope_degrees[0, 0]
                else:
                    return elevation, None
            else:
                return None, None
        except Exception as e:
            print(f"Error calculating slope: {e}")
            return None, None

    def _get_depth_at_point(self, lat, lon):
        """Get depth/elevation at a given lat/lon point."""
        if (self.geotiff_data_array is None or self.geotiff_extent is None):
            return 0.0
        
        try:
            # Convert lat/lon to pixel coordinates
            left, right, bottom, top = tuple(self.geotiff_extent)
            nrows, ncols = self.geotiff_data_array.shape
            
            # Calculate pixel coordinates
            col = int((lon - left) / (right - left) * (ncols - 1))
            row = int((top - lat) / (top - bottom) * (nrows - 1))
            
            # Ensure we're within bounds
            if 0 <= row < nrows and 0 <= col < ncols:
                elevation = self.geotiff_data_array[row, col]
                return elevation if not np.isnan(elevation) else 0.0
            else:
                return 0.0
        except Exception as e:
            print(f"Error getting depth: {e}")
            return 0.0

    def _on_mouse_motion(self, event):
        """Handle all mouse motion events including hover info, line planning, pitch line, and pick center."""
        
        if event.inaxes != self.ax:
            # Remove info text if present
            if hasattr(self, 'mouse_hover_info_text') and self.mouse_hover_info_text is not None:
                self.mouse_hover_info_text.set_visible(False)
                self.mouse_hover_info_text = None
                self.canvas.draw_idle()
            return
        
        mouse_lon = event.xdata
        mouse_lat = event.ydata
        
        if mouse_lon is None or mouse_lat is None:
            # Remove info text if present
            if hasattr(self, 'mouse_hover_info_text') and self.mouse_hover_info_text is not None:
                self.mouse_hover_info_text.set_visible(False)
                self.mouse_hover_info_text = None
                self.canvas.draw_idle()
            return
        
        # Handle line planning motion
        if hasattr(self, 'line_planning_mode') and self.line_planning_mode:
            # Show temporary line from last point to current mouse position
            if len(self.line_planning_points) >= 1:
                lats = [p[0] for p in self.line_planning_points] + [mouse_lat]
                lons = [p[1] for p in self.line_planning_points] + [mouse_lon]
                if hasattr(self, 'line_planning_temp_line') and self.line_planning_temp_line is not None:
                    self.line_planning_temp_line.set_data(lons, lats)
                else:
                    self.line_planning_temp_line, = self.ax.plot(lons, lats, color='orange', linewidth=1, linestyle='--', alpha=0.7)
                
                # Calculate and display line information
                if hasattr(self, 'geotiff_data_array') and self.geotiff_data_array is not None:
                    # Calculate total line length
                    import pyproj
                    geod = pyproj.Geod(ellps="WGS84")
                    total_length_m = 0.0
                    for i in range(1, len(lats)):
                        _, _, seg_length = geod.inv(lons[i-1], lats[i-1], lons[i], lats[i])
                        total_length_m += seg_length
                    total_length_nm = total_length_m / 1852.0
                    
                    # Get survey speed
                    try:
                        speed_knots = float(self.line_survey_speed_entry.get()) if hasattr(self, 'line_survey_speed_entry') else 8.0
                    except:
                        speed_knots = 8.0
                    
                    # Calculate survey time
                    speed_m_per_h = speed_knots * 1852
                    time_hours = total_length_m / speed_m_per_h if speed_m_per_h > 0 else 0
                    time_minutes = time_hours * 60
                    
                    # Get depth and slope at current mouse position
                    elevation, slope = self._calculate_slope_at_point(mouse_lat, mouse_lon)
        
                    # Format the information
                    try:
                        depth_str = f"{abs(elevation):.1f} m" if elevation is not None and not np.isnan(elevation) else "-"
                        slope_str = f"{slope:.1f}°" if slope is not None and not np.isnan(slope) else "-"
                    except Exception:
                        depth_str = "-"
                        slope_str = "-"
                    
                    info_str = f"Depth: {depth_str}\nSlope: {slope_str}\nLine Length: {total_length_m:.1f} m\nLine Length: {total_length_nm:.3f} nm\nSurvey Time: {time_minutes:.1f} min"
                    
                    # Update or create the info text
                    if hasattr(self, 'mouse_hover_info_text') and self.mouse_hover_info_text is not None:
                        self.mouse_hover_info_text.set_text(info_str)
                        self.mouse_hover_info_text.set_visible(True)
                    else:
                        self.mouse_hover_info_text = self.ax.text(0.02, 0.98, info_str, transform=self.ax.transAxes, 
                                                                 fontsize=9, va='top', ha='left', 
                                                                 bbox=dict(boxstyle='round', facecolor='orange', alpha=0.8), 
                                                                 zorder=10)
                self.canvas.draw_idle()
            return
        
        # Handle pitch line motion
        if hasattr(self, 'pick_pitch_line_mode') and self.pick_pitch_line_mode:
            print(f"DEBUG: Pitch line mode active, points: {len(self.pitch_line_points)}")
            # Before first click: show depth and slope in upper left corner
            if len(self.pitch_line_points) == 0:
                if hasattr(self, 'geotiff_data_array') and self.geotiff_data_array is not None:
                    elevation, slope = self._calculate_slope_at_point(mouse_lat, mouse_lon)
                    try:
                        depth_str = f"{abs(elevation):.1f} m" if elevation is not None and not np.isnan(elevation) else "-"
                        slope_str = f"{slope:.1f}°" if slope is not None and not np.isnan(slope) else "-"
                    except Exception:
                        depth_str = "-"
                        slope_str = "-"
                    info_str = f"Depth: {depth_str}\nSlope: {slope_str}"
                    print(f"DEBUG: Creating tooltip: {info_str}")
                    # Update or create the info text in upper left corner
                    if hasattr(self, 'mouse_hover_info_text') and self.mouse_hover_info_text is not None:
                        self.mouse_hover_info_text.set_text(info_str)
                        self.mouse_hover_info_text.set_visible(True)
                        print("DEBUG: Updated existing tooltip")
                    else:
                        self.mouse_hover_info_text = self.ax.text(0.02, 0.98, info_str, transform=self.ax.transAxes, 
                                                                 fontsize=9, va='top', ha='left', 
                                                                 bbox=dict(boxstyle='round', facecolor='plum', alpha=0.8), 
                                                                 zorder=10)
                        print("DEBUG: Created new tooltip")
                    self.canvas.draw_idle()
                else:
                    print("DEBUG: No GeoTIFF data available for tooltip")
                return
            # After first click: show line info
            if len(self.pitch_line_points) == 1:
                print(f"DEBUG: After first click, showing line info")
                # Show temporary line from first point to current mouse position
                lats = [self.pitch_line_points[0][0], mouse_lat]
                lons = [self.pitch_line_points[0][1], mouse_lon]
                if hasattr(self, 'pitch_line_temp_line') and self.pitch_line_temp_line is not None:
                    self.pitch_line_temp_line.set_data(lons, lats)
                else:
                    self.pitch_line_temp_line, = self.ax.plot(lons, lats, color='red', linewidth=2, linestyle='--', alpha=0.7)
                
                # Calculate and display line information
                if hasattr(self, 'geotiff_data_array') and self.geotiff_data_array is not None:
                    # Calculate line length
                    import pyproj
                    geod = pyproj.Geod(ellps="WGS84")
                    _, _, line_length_m = geod.inv(lons[0], lats[0], lons[1], lats[1])
                    line_length_nm = line_length_m / 1852.0
                    
                    # Get survey speed
                    try:
                        speed_knots = float(self.cal_survey_speed_entry.text()) if hasattr(self, 'cal_survey_speed_entry') else 8.0
                    except:
                        speed_knots = 8.0
                    
                    # Calculate survey time
                    speed_m_per_h = speed_knots * 1852
                    time_hours = line_length_m / speed_m_per_h if speed_m_per_h > 0 else 0
                    time_minutes = time_hours * 60
                    
                    # Get depth and slope at current mouse position
                    elevation, slope = self._calculate_slope_at_point(mouse_lat, mouse_lon)
                    
                    # Format the information
                    try:
                        depth_str = f"{abs(elevation):.1f} m" if elevation is not None and not np.isnan(elevation) else "-"
                        slope_str = f"{slope:.1f}°" if slope is not None and not np.isnan(slope) else "-"
                    except Exception:
                        depth_str = "-"
                        slope_str = "-"
                    
                    info_str = f"Depth: {depth_str}\nSlope: {slope_str}\nLine Length: {line_length_m:.1f} m\nLine Length: {line_length_nm:.3f} nm\nSurvey Time: {time_minutes:.1f} min"
                    
                    # Update or create the info text in upper left corner
                    if hasattr(self, 'mouse_hover_info_text') and self.mouse_hover_info_text is not None:
                        self.mouse_hover_info_text.set_text(info_str)
                        self.mouse_hover_info_text.set_visible(True)
                    else:
                        self.mouse_hover_info_text = self.ax.text(0.02, 0.98, info_str, transform=self.ax.transAxes, 
                                                                 fontsize=9, va='top', ha='left', 
                                                                 bbox=dict(boxstyle='round', facecolor='plum', alpha=0.8), 
                                                                 zorder=10)
                self.canvas.draw_idle()
            return
        
        # Handle pick center motion
        if hasattr(self, 'pick_center_mode') and self.pick_center_mode:
            # Show crosshair at mouse position
            if hasattr(self, 'pick_center_crosshair') and self.pick_center_crosshair is not None:
                self.pick_center_crosshair.remove()
            self.pick_center_crosshair, = self.ax.plot(mouse_lon, mouse_lat, '+', color='red', markersize=15, markeredgewidth=2)
            self.canvas.draw_idle()
            return
        
        # Handle roll line motion
        if hasattr(self, 'pick_roll_line_mode') and self.pick_roll_line_mode:
            if len(self.roll_line_points) == 1:
                # Show temporary line from first point to current mouse position
                lats = [self.roll_line_points[0][0], mouse_lat]
                lons = [self.roll_line_points[0][1], mouse_lon]
                if hasattr(self, 'roll_line_temp_line') and self.roll_line_temp_line is not None:
                    self.roll_line_temp_line.set_data(lons, lats)
                else:
                    self.roll_line_temp_line, = self.ax.plot(lons, lats, color='purple', linewidth=2, linestyle='--', alpha=0.7)
                
                # Calculate and display line information
                if hasattr(self, 'geotiff_data_array') and self.geotiff_data_array is not None:
                    # Calculate line length
                    import pyproj
                    geod = pyproj.Geod(ellps="WGS84")
                    _, _, line_length_m = geod.inv(lons[0], lats[0], lons[1], lats[1])
                    line_length_nm = line_length_m / 1852.0
                    
                    # Get survey speed
                    try:
                        speed_knots = float(self.cal_survey_speed_entry.text()) if hasattr(self, 'cal_survey_speed_entry') else 8.0
                    except:
                        speed_knots = 8.0
                    
                    # Calculate survey time
                    speed_m_per_h = speed_knots * 1852
                    time_hours = line_length_m / speed_m_per_h if speed_m_per_h > 0 else 0
                    time_minutes = time_hours * 60
                    
                    # Get depth and slope at current mouse position
                    elevation, slope = self._calculate_slope_at_point(mouse_lat, mouse_lon)
                    
                    # Format the information
                    try:
                        depth_str = f"{abs(elevation):.1f} m" if elevation is not None and not np.isnan(elevation) else "-"
                        slope_str = f"{slope:.1f}°" if slope is not None and not np.isnan(slope) else "-"
                    except Exception:
                        depth_str = "-"
                        slope_str = "-"
                    
                    info_str = f"Depth: {depth_str}\nSlope: {slope_str}\nLine Length: {line_length_m:.1f} m\nLine Length: {line_length_nm:.3f} nm\nSurvey Time: {time_minutes:.1f} min"
                    
                    # Update or create the info text
                    if hasattr(self, 'mouse_hover_info_text') and self.mouse_hover_info_text is not None:
                        self.mouse_hover_info_text.set_text(info_str)
                        self.mouse_hover_info_text.set_visible(True)
                    else:
                        self.mouse_hover_info_text = self.ax.text(0.02, 0.98, info_str, transform=self.ax.transAxes, 
                                                                 fontsize=9, va='top', ha='left', 
                                                                 bbox=dict(boxstyle='round', facecolor='plum', alpha=0.8), 
                                                                 zorder=10)
                self.canvas.draw_idle()
            return
        
        # Handle general hover info (only if no special modes are active)
        if not (hasattr(self, 'line_planning_mode') and self.line_planning_mode or
                hasattr(self, 'pick_pitch_line_mode') and self.pick_pitch_line_mode or
                hasattr(self, 'pick_center_mode') and self.pick_center_mode or
                hasattr(self, 'pick_roll_line_mode') and self.pick_roll_line_mode):
            
            # Only show info if we have a GeoTIFF loaded
            if not (hasattr(self, 'geotiff_data_array') and self.geotiff_data_array is not None):
                # Remove info text if present
                if hasattr(self, 'mouse_hover_info_text') and self.mouse_hover_info_text is not None:
                    self.mouse_hover_info_text.set_visible(False)
                    self.mouse_hover_info_text = None
                    self.canvas.draw_idle()
                return
            
            # Calculate elevation and slope at mouse position
            elevation, slope = self._calculate_slope_at_point(mouse_lat, mouse_lon)
            
            # Format the information
            try:
                elev_str = f"{elevation:.1f} m" if elevation is not None and not np.isnan(elevation) else "-"
                slope_str = f"{slope:.1f}°" if slope is not None and not np.isnan(slope) else "-"
            except Exception:
                elev_str = "-"
                slope_str = "-"
            
            info_str = f"Lat: {mouse_lat:.6f}\nLon: {mouse_lon:.6f}\nElevation: {elev_str}\nSlope: {slope_str}"
            
            # Update or create the info text
            if hasattr(self, 'mouse_hover_info_text') and self.mouse_hover_info_text is not None:
                self.mouse_hover_info_text.set_text(info_str)
                self.mouse_hover_info_text.set_visible(True)
            else:
                self.mouse_hover_info_text = self.ax.text(0.02, 0.98, info_str, transform=self.ax.transAxes, 
                                                         fontsize=9, va='top', ha='left', 
                                                         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8), 
                                                         zorder=10)
            self.canvas.draw_idle()

    def _load_last_used_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_used_dir' in config and os.path.isdir(config['last_used_dir']):
                    self.last_used_dir = config['last_used_dir']
        except Exception:
            pass

    def _save_last_used_dir(self):
        try:
            config = {'last_used_dir': self.last_used_dir}
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _load_last_geotiff_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_geotiff_dir' in config and os.path.isdir(config['last_geotiff_dir']):
                    self.last_geotiff_dir = config['last_geotiff_dir']
        except Exception:
            pass

    def _save_last_geotiff_dir(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
            config['last_geotiff_dir'] = self.last_geotiff_dir
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _load_last_survey_params_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_survey_params_dir' in config and os.path.isdir(config['last_survey_params_dir']):
                    self.last_survey_params_dir = config['last_survey_params_dir']
        except Exception:
            pass

    def _save_last_survey_params_dir(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
            config['last_survey_params_dir'] = self.last_survey_params_dir
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _load_last_export_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_export_dir' in config and os.path.isdir(config['last_export_dir']):
                    self.last_export_dir = config['last_export_dir']
        except Exception:
            pass

    def _save_last_export_dir(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
            config['last_export_dir'] = self.last_export_dir
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _load_last_cal_import_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_cal_import_dir' in config and os.path.isdir(config['last_cal_import_dir']):
                    self.last_cal_import_dir = config['last_cal_import_dir']
        except Exception:
            pass

    def _save_last_cal_import_dir(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
            config['last_cal_import_dir'] = self.last_cal_import_dir
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _load_last_ref_import_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_ref_import_dir' in config and os.path.isdir(config['last_ref_import_dir']):
                    self.last_ref_import_dir = config['last_ref_import_dir']
        except Exception:
            pass

    def _save_last_ref_import_dir(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
            config['last_ref_import_dir'] = self.last_ref_import_dir
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _load_last_line_import_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_line_import_dir' in config and os.path.isdir(config['last_line_import_dir']):
                    self.last_line_import_dir = config['last_line_import_dir']
        except Exception:
            pass

    def _save_last_line_import_dir(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
            config['last_line_import_dir'] = self.last_line_import_dir
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _clear_calibration_lines(self):
        """Clear all calibration lines (pitch, heading, roll) and update the plot."""
        self.pitch_line_points = []
        self.heading_lines = []
        self.roll_line_points = []
        self._plot_survey_plan(preserve_view_limits=True)
        self._draw_pitch_line_profile()
        self._update_cal_line_times()
        self._update_cal_export_name_from_pitch_line()
        self._update_pitch_line_button_states()
        self._update_roll_line_button_states()
        self._update_line_planning_button_states()
        # Reset button styles - if GeoTIFF is loaded, make "Draw a Pitch Line" bold and orange
        if hasattr(self, 'geotiff_dataset_original') and self.geotiff_dataset_original is not None:
            if hasattr(self, 'pick_pitch_line_btn'):
                self.pick_pitch_line_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")
        else:
            if hasattr(self, 'pick_pitch_line_btn'):
                self.pick_pitch_line_btn.setStyleSheet("")
        if hasattr(self, 'add_heading_lines_btn'):
            self.add_heading_lines_btn.setStyleSheet("")
        if hasattr(self, 'pick_roll_line_btn'):
            self.pick_roll_line_btn.setStyleSheet("")

    def _update_pitch_line_button_states(self):
        """Update the state of buttons that depend on having a pitch line."""
        has_pitch_line = hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2
        
        # Update button states
        if hasattr(self, 'edit_pitch_line_btn'):
            self.edit_pitch_line_btn.setEnabled(has_pitch_line)
        if hasattr(self, 'add_heading_lines_btn'):
            self.add_heading_lines_btn.setEnabled(has_pitch_line)

    def _update_roll_line_button_states(self):
        """Update the state of buttons that depend on having a roll line."""
        has_roll_line = hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2
        
        # Update button states
        if hasattr(self, 'edit_roll_line_btn'):
            self.edit_roll_line_btn.setEnabled(has_roll_line)

    def _update_line_planning_button_states(self):
        """Update the state of buttons that depend on having a line planning route."""
        has_line_planning = hasattr(self, 'line_planning_points') and len(self.line_planning_points) >= 2
        
        # Update button states
        if hasattr(self, 'line_edit_btn'):
            self.line_edit_btn.setEnabled(has_line_planning)

    def _toggle_edit_line_planning_mode(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot edit line planning.")
            return
        if self.geotiff_dataset_original is None:
            self._show_message("warning","No GeoTIFF", "Load a GeoTIFF first to edit line planning.")
            return
        if len(self.line_planning_points) < 2:
            self._show_message("warning","No Line Planning", "Draw a line planning route first before editing it.")
            return
        
        self.edit_line_planning_mode = not self.edit_line_planning_mode
        if self.edit_line_planning_mode:
            # Create draggable handles for all points
            self.line_planning_handles = []
            for i, (lat, lon) in enumerate(self.line_planning_points):
                color = 'red' if i == 0 else 'blue' if i == len(self.line_planning_points) - 1 else 'green'
                handle = self.ax.scatter([lon], [lat], 
                                       color=color, s=100, marker='o', 
                                       edgecolors='black', linewidth=2, 
                                       zorder=10, picker=5)
                self.line_planning_handles.append(handle)
            
            # Connect pick events
            self.line_planning_pick_cid = self.canvas.mpl_connect('pick_event', self._on_line_planning_handle_pick)
            self.line_planning_motion_cid = self.canvas.mpl_connect('motion_notify_event', self._on_line_planning_handle_motion)
            self.line_planning_release_cid = self.canvas.mpl_connect('button_release_event', self._on_line_planning_handle_release)
            
            self.line_edit_btn.setText("Click to Stop Editing Line Planning")
            self.canvas_widget.setCursor(Qt.CursorShape.SizeAllCursor)
            self.set_line_info_text("You can drag the red (start), green (waypoints), and blue (end) points to edit the line planning route.")
        else:
            # Remove handles by clearing the axes and redrawing
            if self.line_planning_handles:
                # Clear the handles by redrawing the plot
                self._plot_survey_plan(preserve_view_limits=True)
                self.line_planning_handles = []
            
            # Disconnect events
            if hasattr(self, 'line_planning_pick_cid'):
                self.canvas.mpl_disconnect(self.line_planning_pick_cid)
            if hasattr(self, 'line_planning_motion_cid'):
                self.canvas.mpl_disconnect(self.line_planning_motion_cid)
            if hasattr(self, 'line_planning_release_cid'):
                self.canvas.mpl_disconnect(self.line_planning_release_cid)
            
            self.line_edit_btn.setText("Edit Line Planning")
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
            self.dragging_line_planning_handle = None
            
            # Redraw to update the line planning
            self._plot_survey_plan(preserve_view_limits=True)
            self._draw_line_planning_profile()
            
            # Report updated line planning summary
            try:
                import pyproj
                geod = pyproj.Geod(ellps="WGS84")
                num_segments = len(self.line_planning_points) - 1
                total_length_m = 0.0
                for i in range(1, len(self.line_planning_points)):
                    lat1, lon1 = self.line_planning_points[i-1]
                    lat2, lon2 = self.line_planning_points[i]
                    _, _, d = geod.inv(lon1, lat1, lon2, lat2)
                    total_length_m += d
                total_length_km = total_length_m / 1000.0
                total_length_nm = total_length_m / 1852.0
                try:
                    speed_knots = float(self.line_survey_speed_entry.get())
                except Exception:
                    speed_knots = 8.0
                speed_m_per_h = speed_knots * 1852
                time_hours = total_length_m / speed_m_per_h if speed_m_per_h > 0 else 0
                time_minutes = time_hours * 60
                summary = (
                    f"Updated Line Planning Summary:\n"
                    f"Segments: {num_segments}\n"
                    f"Total Length: {total_length_m:.1f} m\n"
                    f"Total Length: {total_length_km:.3f} km\n"
                    f"Total Length: {total_length_nm:.3f} nautical miles\n"
                    f"Time to Run: {time_minutes:.1f} min\n"
                    f"Time to Run: {time_hours:.2f} hr"
                )
                self.set_line_info_text(summary)
            except Exception as e:
                self.set_line_info_text(f"Error calculating updated line planning summary: {e}")

    def _on_line_planning_handle_pick(self, event):
        """Handle picking of line planning handles for dragging."""
        if not self.edit_line_planning_mode:
            return
        
        if event.mouseevent.button != 1:  # Only left mouse button
            return
        
        # Find which handle was picked
        for i, handle in enumerate(self.line_planning_handles):
            if event.artist == handle:
                self.dragging_line_planning_handle = i
                break

    def _on_line_planning_handle_motion(self, event):
        """Handle dragging of line planning handles."""
        if not self.edit_line_planning_mode or self.dragging_line_planning_handle is None:
            return
        
        if event.inaxes != self.ax:
            return
        
        if event.xdata is None or event.ydata is None:
            return
        
        # Update the handle position
        handle_idx = self.dragging_line_planning_handle
        if handle_idx < len(self.line_planning_handles) and self.line_planning_handles[handle_idx] is not None:
            self.line_planning_handles[handle_idx].set_offsets([[event.xdata, event.ydata]])
            self.line_planning_points[handle_idx] = (event.ydata, event.xdata)
        
        # Update the line planning itself
        if hasattr(self, 'line_planning_edit_line') and self.line_planning_edit_line is not None:
            lats = [p[0] for p in self.line_planning_points]
            lons = [p[1] for p in self.line_planning_points]
            self.line_planning_edit_line.set_data(lons, lats)
        else:
            lats = [p[0] for p in self.line_planning_points]
            lons = [p[1] for p in self.line_planning_points]
            self.line_planning_edit_line, = self.ax.plot(lons, lats, color='orange', linewidth=2, zorder=9)
        
        self.canvas.draw_idle()

    def _on_line_planning_handle_release(self, event):
        """Handle release of line planning handles."""
        if not self.edit_line_planning_mode:
            return
        
        if self.dragging_line_planning_handle is not None:
            # Remove the temp edit line and redraw the full plot to finalize
            if hasattr(self, 'line_planning_edit_line') and self.line_planning_edit_line is not None:
                self.line_planning_edit_line.remove()
                self.line_planning_edit_line = None
            self._plot_survey_plan(preserve_view_limits=True)
            
            # Recreate the handles after the plot redraw
            self.line_planning_handles = []
            for i, (lat, lon) in enumerate(self.line_planning_points):
                color = 'red' if i == 0 else 'blue' if i == len(self.line_planning_points) - 1 else 'green'
                handle = self.ax.scatter([lon], [lat], 
                                       color=color, s=100, marker='o', 
                                       edgecolors='black', linewidth=2, 
                                       zorder=10, picker=5)
                self.line_planning_handles.append(handle)
            
            self.dragging_line_planning_handle = None

    def _update_cal_line_times(self):
        """Update the Pitch and Roll line time displays in the Calibration Planning tab."""
        # Pitch Line
        try:
            if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2:
                (lat1, lon1), (lat2, lon2) = self.pitch_line_points
                import pyproj
                geod = pyproj.Geod(ellps="WGS84")
                _, _, dist = geod.inv(lon1, lat1, lon2, lat2)
                speed_knots = float(self.cal_survey_speed_entry.get()) if self.cal_survey_speed_entry.get() else 8.0
                speed_m_per_h = speed_knots * 1852
                time_hours = dist / speed_m_per_h if speed_m_per_h > 0 else 0
                time_minutes = time_hours * 60
                # self.pitch_time_min_value.config(text=f"{time_minutes:.1f}")  # Removed
            else:
                # self.pitch_time_min_value.config(text="-")  # Removed
                pass
        except Exception:
            # self.pitch_time_min_value.config(text="-")  # Removed
            pass
        # Roll Line
        try:
            if hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2:
                (lat1, lon1), (lat2, lon2) = self.roll_line_points
                import pyproj
                geod = pyproj.Geod(ellps="WGS84")
                _, _, dist = geod.inv(lon1, lat1, lon2, lat2)
                speed_knots = float(self.cal_survey_speed_entry.get()) if self.cal_survey_speed_entry.get() else 8.0
                speed_m_per_h = speed_knots * 1852
                time_hours = dist / speed_m_per_h if speed_m_per_h > 0 else 0
                time_minutes = time_hours * 60
                # self.roll_time_value.config(text=f"{time_hours:.2f}")  # Removed
                # self.roll_time_min_value.config(text=f"{time_minutes:.1f}")  # Removed
            else:
                # self.roll_time_value.config(text="-")  # Removed
                # self.roll_time_min_value.config(text="-")  # Removed
                pass
        except Exception:
            # self.roll_time_value.config(text="-")  # Removed
            # self.roll_time_min_value.config(text="-")  # Removed
            pass

    def _draw_current_profile(self):
        # Helper to redraw the correct profile based on the current tab
        if self.param_notebook.currentIndex() == 1:
            self._draw_crossline_profile()
        else:
            self._draw_pitch_line_profile()


    def set_ref_info_text(self, message, append=False):
        """Add a message to the Activity Log with [Reference] prefix. Maintains up to 200 lines of history."""
        if not hasattr(self, 'activity_log_text'):
            return
        prefixed_message = f"[Reference] {message}"
        self.activity_log_text.setReadOnly(False)
        if append:
            # Prepend new message (newest at top) - only when explicitly requested
            current_text = self.activity_log_text.toPlainText().rstrip('\n')
            new_text = prefixed_message + "\n" + current_text if current_text else prefixed_message + "\n"
            # Split into lines and limit to 200
            lines = new_text.splitlines()
            if len(lines) > 200:
                lines = lines[:200]
            self.activity_log_text.setPlainText("\n".join(lines) + "\n")
            # Move cursor to top
            cursor = self.activity_log_text.textCursor()
            cursor.movePosition(QTextCursor.MoveOperation.Start)
            self.activity_log_text.setTextCursor(cursor)
        else:
            # Append to end (default behavior - newest at bottom)
            self.activity_log_text.append(prefixed_message)
            # Maintain up to 200 lines by removing oldest from top
            num_lines = self.activity_log_text.document().blockCount()
            if num_lines > 200:
                cursor = self.activity_log_text.textCursor()
                cursor.movePosition(QTextCursor.MoveOperation.Start)
                cursor.movePosition(QTextCursor.MoveOperation.Down, QTextCursor.MoveMode.MoveAnchor, num_lines - 200)
                cursor.movePosition(QTextCursor.MoveOperation.Start, QTextCursor.MoveMode.KeepAnchor)
                cursor.removeSelectedText()
            # Move cursor to end and scroll to show newest message
            cursor = self.activity_log_text.textCursor()
            cursor.movePosition(QTextCursor.MoveOperation.End)
            self.activity_log_text.setTextCursor(cursor)
        self.activity_log_text.ensureCursorVisible()
        self.activity_log_text.setReadOnly(True)

    def _reset_reference_tab(self):
        """Reset Reference tab data and entry fields to defaults."""
        # Clear line data
        self.survey_lines_data = []
        self.cross_line_data = []
        self.central_point_coords = (None, None)
        
        # Reset entry fields to defaults
        if hasattr(self, 'central_lat_entry'):
            self.central_lat_entry.clear()
        if hasattr(self, 'central_lon_entry'):
            self.central_lon_entry.clear()
        if hasattr(self, 'line_length_entry'):
            self.line_length_entry.clear()
        if hasattr(self, 'heading_entry'):
            self.heading_entry.setText("0")
        if hasattr(self, 'dist_between_lines_entry'):
            self.dist_between_lines_entry.clear()
        if hasattr(self, 'num_lines_entry'):
            self.num_lines_entry.setText("5")
        if hasattr(self, 'bisect_lead_entry'):
            self.bisect_lead_entry.setText("100")
        if hasattr(self, 'survey_speed_entry'):
            self.survey_speed_entry.setText("8")
        if hasattr(self, 'crossline_passes_entry'):
            self.crossline_passes_entry.setText("2")
        if hasattr(self, 'export_name_entry'):
            self.export_name_entry.clear()
            try:
                dist = int(float(self.dist_between_lines_entry.get())) if self.dist_between_lines_entry.get() else 0
                heading = int(float(self.heading_entry.get())) if self.heading_entry.get() else 0
                export_name = f"Reference_{dist}m_{heading}deg"
            except Exception:
                export_name = "Reference_0m_0deg"
            self.export_name_entry.setText(export_name)
        if hasattr(self, 'offset_direction_var'):
            self.offset_direction_var.set("North")
            if hasattr(self, 'offset_direction_combo'):
                self.offset_direction_combo.setCurrentIndex(0)
        if hasattr(self, 'line_length_multiplier'):
            self.line_length_multiplier.set(8.0)
        if hasattr(self, 'dist_between_lines_multiplier'):
            self.dist_between_lines_multiplier.set(1.0)
        
        # Redraw plot (without clearing GeoTIFF)
        self._plot_survey_plan(preserve_view_limits=True)
    
    def _reset_calibration_tab(self):
        """Reset Calibration tab data and entry fields to defaults."""
        # Clear line data
        self.pitch_line_points = []
        self.heading_lines = []
        self.roll_line_points = []
        
        # Reset mode flags
        self.pick_pitch_line_mode = False
        self.edit_pitch_line_mode = False
        self.pick_roll_line_mode = False
        self.edit_roll_line_mode = False
        
        # Clear visual handles/annotations
        self.pitch_line_start_handle = None
        self.pitch_line_end_handle = None
        self.dragging_pitch_line_handle = None
        self.pitch_line_temp_line = None
        self.pitch_line_annotation = None
        self.pitch_line_tooltip = None
        self.pitch_line_tooltip_text = None
        self.pitch_line_edit_line = None
        self.roll_line_start_handle = None
        self.roll_line_end_handle = None
        self.dragging_roll_line_handle = None
        self.roll_line_edit_line = None
        
        # Reset entry fields to defaults
        if hasattr(self, 'cal_survey_speed_entry'):
            self.cal_survey_speed_entry.setText("8")
        if hasattr(self, 'cal_line_offset_entry'):
            self.cal_line_offset_entry.clear()
        if hasattr(self, 'cal_export_name_entry'):
            self.cal_export_name_entry.clear()
        
        # Update button states
        self._update_pitch_line_button_states()
        self._update_roll_line_button_states()
        
        # Reset button styles - if GeoTIFF is loaded, make "Draw a Pitch Line" bold and orange
        if hasattr(self, 'geotiff_dataset_original') and self.geotiff_dataset_original is not None:
            if hasattr(self, 'pick_pitch_line_btn'):
                self.pick_pitch_line_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")
        else:
            if hasattr(self, 'pick_pitch_line_btn'):
                self.pick_pitch_line_btn.setStyleSheet("")
        if hasattr(self, 'add_heading_lines_btn'):
            self.add_heading_lines_btn.setStyleSheet("")
        if hasattr(self, 'pick_roll_line_btn'):
            self.pick_roll_line_btn.setStyleSheet("")
        
        # Redraw plot and profiles
        self._plot_survey_plan(preserve_view_limits=True)
        if hasattr(self, '_draw_pitch_line_profile'):
            self._draw_pitch_line_profile()
        if hasattr(self, '_update_cal_line_times'):
            self._update_cal_line_times()
    
    def _reset_line_planning_tab(self):
        """Reset Line Planning tab data and entry fields to defaults."""
        # Clear line data
        self.line_planning_points = []
        self.line_planning_mode = False
        self.edit_line_planning_mode = False
        
        # Clear visual elements
        if hasattr(self, 'line_planning_line') and self.line_planning_line is not None:
            self.line_planning_line.remove()
            self.line_planning_line = None
        if hasattr(self, 'line_planning_temp_line') and self.line_planning_temp_line is not None:
            self.line_planning_temp_line.remove()
            self.line_planning_temp_line = None
        if hasattr(self, 'line_planning_edit_line') and self.line_planning_edit_line is not None:
            self.line_planning_edit_line.remove()
            self.line_planning_edit_line = None
        if hasattr(self, 'line_planning_info_text') and self.line_planning_info_text is not None:
            self.line_planning_info_text.set_visible(False)
            self.line_planning_info_text = None
        
        # Clear handles
        self.line_planning_handles = []
        self.dragging_line_planning_handle = None
        
        # Reset entry fields to defaults
        if hasattr(self, 'line_survey_speed_entry'):
            self.line_survey_speed_entry.clear()
            self.line_survey_speed_entry.setText("8")
        
        # Update button states
        if hasattr(self, 'line_start_draw_btn'):
            self.line_start_draw_btn.setText("Start Drawing Line")
        if hasattr(self, 'line_edit_btn'):
            self.line_edit_btn.setEnabled(False)
        if hasattr(self, 'canvas_widget'):
            self.canvas_widget.setCursor(Qt.CursorShape.ArrowCursor)
        
        # Update button states
        self._update_line_planning_button_states()
        
        # Redraw plot
        self._plot_survey_plan(preserve_view_limits=True)
    
    def _on_tab_changed(self, event=None):
        """Handle tab change event - clear previous tab's data and settings."""
        if not hasattr(self, 'param_notebook'):
            return
        
        try:
            current_tab = self.param_notebook.currentIndex()
            
            # Only clear if this is not the first tab change (i.e., user has actually switched tabs)
            if self.tab_switch_initialized:
                # Clear the previous tab's data (not the current one)
                if self.previous_tab_index == 0:  # Calibration tab
                    self._reset_calibration_tab()
                elif self.previous_tab_index == 1:  # Reference tab
                    self._reset_reference_tab()
                elif self.previous_tab_index == 2:  # Line Planning tab
                    self._reset_line_planning_tab()
            
            # Mark as initialized and update previous tab index
            self.tab_switch_initialized = True
            self.previous_tab_index = current_tab
        except Exception as e:
            print(f"Error in tab change handler: {e}")

    def _create_widgets(self):
        print("Creating widgets...")
        # --- 1. Tabbed Parameters Section ---
        # Use QScrollArea for scrolling
        self.param_scroll = QScrollArea()
        self.param_scroll.setWidgetResizable(True)
        self.param_scroll.setMaximumWidth(420)
        self.param_scroll.setMinimumWidth(300)  # Minimum width for visibility
        self.param_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)  # Disable horizontal scrollbar
        
        param_widget = QWidget()
        param_widget.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Preferred)
        param_layout = QVBoxLayout(param_widget)
        
        # --- GeoTIFF Control GroupBox (above tabs) ---
        geotiff_groupbox = QGroupBox("GeoTIFF Control")
        geotiff_layout = QVBoxLayout(geotiff_groupbox)
        
        # Load GeoTIFF and Remove GeoTIFF buttons
        geotiff_button_frame = QWidget()
        geotiff_button_layout = QHBoxLayout(geotiff_button_frame)
        geotiff_button_layout.setContentsMargins(0, 0, 0, 0)
        self.load_geotiff_btn = QPushButton("Load GeoTIFF")
        self.load_geotiff_btn.clicked.connect(self._load_geotiff)
        # Set button to orange and bold initially
        self.load_geotiff_btn.setStyleSheet("QPushButton { color: rgb(255, 165, 0); font-weight: bold; }")
        geotiff_button_layout.addWidget(self.load_geotiff_btn)
        self.remove_geotiff_btn = QPushButton("Remove GeoTIFF")
        self.remove_geotiff_btn.clicked.connect(self._remove_geotiff)
        self.remove_geotiff_btn.setEnabled(False)  # Disabled initially (no GeoTIFF loaded)
        geotiff_button_layout.addWidget(self.remove_geotiff_btn)
        geotiff_layout.addWidget(geotiff_button_frame)
        
        # Display mode dropdown
        geotiff_layout.addWidget(QLabel("Display Mode:"))
        self.elevation_slope_combo = QComboBox()
        self.elevation_slope_combo.addItems(["Shaded Relief", "Shaded Slope", "Hillshade", "Slope"])
        self.elevation_slope_combo.setCurrentText("Shaded Relief")  # Set default
        self.elevation_slope_combo.currentTextChanged.connect(self._on_geotiff_display_mode_changed)
        geotiff_layout.addWidget(self.elevation_slope_combo)
        
        # Dynamic Resolution button
        self.dynamic_resolution_btn = QPushButton("Dynamic Resolution: ON")
        self.dynamic_resolution_btn.clicked.connect(self._toggle_dynamic_resolution)
        geotiff_layout.addWidget(self.dynamic_resolution_btn)
        
        # Show Contours checkbox and Interval on same row
        contours_interval_frame = QWidget()
        contours_interval_layout = QHBoxLayout(contours_interval_frame)
        contours_interval_layout.setContentsMargins(0, 0, 0, 0)
        self.show_contours_checkbox = QCheckBox("Contours")
        self.show_contours_checkbox.setChecked(self.show_contours_var)
        self.show_contours_checkbox.stateChanged.connect(self._on_contour_checkbox_changed)
        contours_interval_layout.addWidget(self.show_contours_checkbox)
        contours_interval_layout.addStretch()  # Add stretch to push interval to the right
        contours_interval_layout.addWidget(QLabel("Interval (m):"))
        self.contour_interval_entry = QLineEdit("200")
        self.contour_interval_entry.textChanged.connect(self._on_contour_interval_changed)
        self.contour_interval_entry.setMaximumWidth(80)  # Limit width of entry field
        contours_interval_layout.addWidget(self.contour_interval_entry)
        geotiff_layout.addWidget(contours_interval_frame)
        
        param_layout.addWidget(geotiff_groupbox)
        
        # --- Test Planning GroupBox (contains tabs) ---
        test_planning_groupbox = QGroupBox("Test Planning")
        test_planning_layout = QVBoxLayout(test_planning_groupbox)
        
        self.param_notebook = QTabWidget()
        
        # Make the selected tab text bold
        self.param_notebook.setStyleSheet("""
            QTabBar::tab:selected {
                font-weight: bold;
            }
        """)
        
        self.reference_frame = QWidget()
        self.calibration_frame = QWidget()
        self.line_planning_frame = QWidget()
        self.param_notebook.addTab(self.calibration_frame, "Calibration")
        self.param_notebook.addTab(self.reference_frame, "Reference")
        self.param_notebook.addTab(self.line_planning_frame, "Line")
        
        test_planning_layout.addWidget(self.param_notebook)
        param_layout.addWidget(test_planning_groupbox)
        
        # --- Activity Log GroupBox ---
        activity_log_groupbox = QGroupBox("Activity Log")
        activity_log_layout = QVBoxLayout(activity_log_groupbox)
        
        self.activity_log_text = QTextEdit()
        self.activity_log_text.setReadOnly(True)
        # Set light yellow background and black text color for visibility
        self.activity_log_text.setStyleSheet("background-color: #ffffe0; color: black;")
        # Let it expand to fill the groupbox - no height constraint
        activity_log_layout.addWidget(self.activity_log_text, 1)  # Stretch factor 1 to fill available space
        
        param_layout.addWidget(activity_log_groupbox)
        
        self.param_scroll.setWidget(param_widget)
        
        # Connect tab change signal
        self.param_notebook.currentChanged.connect(self._on_tab_changed)

        # --- Reference Planning Tab (was input_frame) ---
        self.input_frame = self.reference_frame  # For compatibility with rest of code
        ref_layout = QGridLayout(self.input_frame)
        ref_layout.setColumnStretch(0, 1)  # Labels
        ref_layout.setColumnStretch(1, 2)  # Entries
        ref_layout.setSpacing(3)  # Reduce spacing between groupboxes
        
        row = 0
        
        # --- Reference Line Parameters GroupBox ---
        ref_line_params_groupbox = QGroupBox("Reference Line Parameters")
        ref_line_params_groupbox.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        ref_line_params_layout = QGridLayout(ref_line_params_groupbox)
        ref_line_params_layout.setSpacing(3)
        ref_line_params_layout.setContentsMargins(9, 9, 9, 9)
        ref_line_params_layout.setColumnStretch(0, 1)
        ref_line_params_layout.setColumnStretch(1, 2)
        
        ref_line_row = 0
        self.pick_center_btn = QPushButton("Pick Center from GeoTIFF")
        self.pick_center_btn.clicked.connect(self._toggle_pick_center_mode)
        self.pick_center_btn.setEnabled(False)  # Disabled initially (no GeoTIFF loaded)
        ref_line_params_layout.addWidget(self.pick_center_btn, ref_line_row, 0, 1, 2)
        ref_line_row += 1
        
        ref_line_params_layout.addWidget(QLabel("Central Latitude:"), ref_line_row, 0)
        self.central_lat_entry = QLineEdit()
        ref_line_params_layout.addWidget(self.central_lat_entry, ref_line_row, 1)
        ref_line_row += 1
        
        ref_line_params_layout.addWidget(QLabel("Central Longitude:"), ref_line_row, 0)
        self.central_lon_entry = QLineEdit()
        ref_line_params_layout.addWidget(self.central_lon_entry, ref_line_row, 1)
        ref_line_row += 1
        
        ref_line_params_layout.addWidget(QLabel("Number of Lines:"), ref_line_row, 0)
        self.num_lines_entry = QLineEdit("5")
        ref_line_params_layout.addWidget(self.num_lines_entry, ref_line_row, 1)
        ref_line_row += 1
        
        ref_line_params_layout.addWidget(QLabel("Offset Direction:"), ref_line_row, 0)
        self.offset_direction_combo = QComboBox()
        self.offset_direction_combo.addItems(["North", "South"])
        self.offset_direction_combo.setCurrentIndex(0)
        self.offset_direction_var = "North"
        self.offset_direction_combo.currentTextChanged.connect(lambda text: setattr(self, 'offset_direction_var', text))
        ref_line_params_layout.addWidget(self.offset_direction_combo, ref_line_row, 1)
        ref_line_row += 1
        
        ref_line_params_layout.addWidget(QLabel("Line Length (m):"), ref_line_row, 0)
        self.line_length_entry = QLineEdit()
        ref_line_params_layout.addWidget(self.line_length_entry, ref_line_row, 1)
        ref_line_row += 1
        
        ref_line_params_layout.addWidget(QLabel("Heading (deg, 0-360):"), ref_line_row, 0)
        self.heading_entry = QLineEdit("0")
        ref_line_params_layout.addWidget(self.heading_entry, ref_line_row, 1)
        ref_line_row += 1
        
        ref_line_params_layout.addWidget(QLabel("Distance Between Lines (m):"), ref_line_row, 0)
        self.dist_between_lines_entry = QLineEdit()
        ref_line_params_layout.addWidget(self.dist_between_lines_entry, ref_line_row, 1)
        ref_line_row += 1
        
        ref_line_params_layout.addWidget(QLabel("Crossline Lead-in/out (m):"), ref_line_row, 0)
        self.bisect_lead_entry = QLineEdit("100")
        ref_line_params_layout.addWidget(self.bisect_lead_entry, ref_line_row, 1)
        ref_line_row += 1
        
        ref_line_params_layout.addWidget(QLabel("Line Length Multiplier:"), ref_line_row, 0)
        slider_frame_len = QWidget()
        slider_layout_len = QHBoxLayout(slider_frame_len)
        slider_layout_len.setContentsMargins(0, 0, 0, 0)
        self.multiplier_slider_len = QSlider(Qt.Orientation.Horizontal)
        self.multiplier_slider_len.setMinimum(10)  # 1.0 * 10
        self.multiplier_slider_len.setMaximum(100)  # 10.0 * 10
        self.multiplier_slider_len.setValue(int(self.line_length_multiplier * 10))
        self.multiplier_slider_len.valueChanged.connect(lambda val: [setattr(self, 'line_length_multiplier', val/10.0), self._update_multiplier_label_len(val/10.0), self._on_parameter_changed()])
        slider_layout_len.addWidget(self.multiplier_slider_len)
        self.multiplier_label_len = QLabel(f"{self.line_length_multiplier:.1f}")
        slider_layout_len.addWidget(self.multiplier_label_len)
        ref_line_params_layout.addWidget(slider_frame_len, ref_line_row, 1)
        ref_line_row += 1
        
        ref_line_params_layout.addWidget(QLabel("Separation Multiplier:"), ref_line_row, 0)
        slider_frame_dist = QWidget()
        slider_layout_dist = QHBoxLayout(slider_frame_dist)
        slider_layout_dist.setContentsMargins(0, 0, 0, 0)
        self.multiplier_slider_dist = QSlider(Qt.Orientation.Horizontal)
        self.multiplier_slider_dist.setMinimum(0)  # 0.0 * 10
        self.multiplier_slider_dist.setMaximum(20)  # 2.0 * 10
        self.multiplier_slider_dist.setValue(int(self.dist_between_lines_multiplier * 10))
        self.multiplier_slider_dist.valueChanged.connect(lambda val: [setattr(self, 'dist_between_lines_multiplier', val/10.0), self._update_multiplier_label_dist(val/10.0), self._on_parameter_changed()])
        slider_layout_dist.addWidget(self.multiplier_slider_dist)
        self.multiplier_label_dist = QLabel(f"{self.dist_between_lines_multiplier:.1f}")
        slider_layout_dist.addWidget(self.multiplier_label_dist)
        ref_line_params_layout.addWidget(slider_frame_dist, ref_line_row, 1)
        ref_line_row += 1
        
        self.generate_plot_btn = QPushButton("Generate Survey Plan")
        self.generate_plot_btn.clicked.connect(self._generate_and_plot)
        ref_line_params_layout.addWidget(self.generate_plot_btn, ref_line_row, 0, 1, 2)
        
        # Connect parameter changes to auto-regenerate (after all widgets are created)
        self.central_lat_entry.textChanged.connect(self._on_parameter_changed)
        self.central_lon_entry.textChanged.connect(self._on_parameter_changed)
        self.line_length_entry.textChanged.connect(self._on_parameter_changed)
        self.heading_entry.textChanged.connect(self._on_parameter_changed)
        self.dist_between_lines_entry.textChanged.connect(self._on_parameter_changed)
        self.num_lines_entry.textChanged.connect(self._on_parameter_changed)
        self.bisect_lead_entry.textChanged.connect(self._on_parameter_changed)
        self.offset_direction_combo.currentTextChanged.connect(self._on_parameter_changed)
        
        ref_layout.addWidget(ref_line_params_groupbox, row, 0, 1, 2)
        ref_layout.setRowStretch(row, 0)
        row += 1
        
        # --- Plot Control GroupBox ---
        ref_plot_control_groupbox = QGroupBox("Plot Control")
        ref_plot_control_groupbox.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        ref_plot_control_layout = QVBoxLayout(ref_plot_control_groupbox)
        ref_plot_control_layout.setSpacing(0)
        ref_plot_control_layout.setContentsMargins(9, 9, 9, 9)
        
        self.zoom_to_plan_btn = QPushButton("Zoom to Plan")
        self.zoom_to_plan_btn.clicked.connect(self._zoom_to_plan)
        ref_plot_control_layout.addWidget(self.zoom_to_plan_btn)
        ref_plot_control_layout.addSpacing(3)
        
        self.zoom_to_geotiff_btn_ref = QPushButton("Zoom to GeoTIFF")
        self.zoom_to_geotiff_btn_ref.clicked.connect(self._zoom_to_geotiff)
        ref_plot_control_layout.addWidget(self.zoom_to_geotiff_btn_ref)
        ref_plot_control_layout.addSpacing(3)
        
        self.clear_plot_btn = QPushButton("Clear Plot")
        self.clear_plot_btn.clicked.connect(self._clear_plot)
        ref_plot_control_layout.addWidget(self.clear_plot_btn)
        
        ref_layout.addWidget(ref_plot_control_groupbox, row, 0, 1, 2)
        ref_layout.setRowStretch(row, 0)
        row += 1
        
        # --- Test Plan Info GroupBox ---
        ref_test_plan_info_groupbox = QGroupBox("Test Plan Info")
        ref_test_plan_info_groupbox.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        ref_test_plan_info_layout = QGridLayout(ref_test_plan_info_groupbox)
        ref_test_plan_info_layout.setSpacing(3)
        ref_test_plan_info_layout.setContentsMargins(9, 9, 9, 9)
        ref_test_plan_info_layout.setColumnStretch(0, 1)
        ref_test_plan_info_layout.setColumnStretch(1, 2)
        
        ref_test_plan_row = 0
        ref_test_plan_info_layout.addWidget(QLabel("Survey Speed (knots):"), ref_test_plan_row, 0)
        self.survey_speed_entry = QLineEdit("8")
        ref_test_plan_info_layout.addWidget(self.survey_speed_entry, ref_test_plan_row, 1)
        ref_test_plan_row += 1
        
        ref_test_plan_info_layout.addWidget(QLabel("Number of Crossline Passes:"), ref_test_plan_row, 0)
        self.crossline_passes_entry = QLineEdit("2")
        ref_test_plan_info_layout.addWidget(self.crossline_passes_entry, ref_test_plan_row, 1)
        ref_test_plan_row += 1
        
        self.ref_show_info_btn = QPushButton("Show Reference Test Info")
        self.ref_show_info_btn.clicked.connect(self._show_reference_planning_info)
        ref_test_plan_info_layout.addWidget(self.ref_show_info_btn, ref_test_plan_row, 0, 1, 2)
        
        ref_layout.addWidget(ref_test_plan_info_groupbox, row, 0, 1, 2)
        ref_layout.setRowStretch(row, 0)
        row += 1
        
        # --- Import/Export GroupBox ---
        ref_import_export_groupbox = QGroupBox("Import/Export")
        ref_import_export_groupbox.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        ref_import_export_layout = QVBoxLayout(ref_import_export_groupbox)
        ref_import_export_layout.setSpacing(0)
        ref_import_export_layout.setContentsMargins(9, 9, 9, 9)
        
        # Import/Export buttons on same row
        ref_import_export_button_frame = QWidget()
        ref_import_export_button_layout = QHBoxLayout(ref_import_export_button_frame)
        ref_import_export_button_layout.setSpacing(3)
        ref_import_export_button_layout.setContentsMargins(0, 0, 0, 0)
        
        self.import_survey_btn = QPushButton("Import Survey")
        self.import_survey_btn.clicked.connect(self._import_survey_files)
        ref_import_export_button_layout.addWidget(self.import_survey_btn)
        self.export_survey_btn = QPushButton("Export Survey")
        self.export_survey_btn.clicked.connect(self._export_survey_files)
        ref_import_export_button_layout.addWidget(self.export_survey_btn)
        
        ref_import_export_layout.addWidget(ref_import_export_button_frame)
        ref_import_export_layout.addSpacing(3)
        
        # Export Name at the bottom
        ref_export_name_frame = QWidget()
        ref_export_name_layout = QGridLayout(ref_export_name_frame)
        ref_export_name_layout.setSpacing(0)
        ref_export_name_layout.setContentsMargins(0, 0, 0, 0)
        ref_export_name_layout.setColumnStretch(0, 1)
        ref_export_name_layout.setColumnStretch(1, 2)
        ref_export_name_layout.addWidget(QLabel("Export Name:"), 0, 0)
        self.export_name_entry = QLineEdit()
        try:
            dist = int(float(self.dist_between_lines_entry.text() or "0"))
            heading = int(float(self.heading_entry.text() or "0"))
            export_name = f"Reference_{dist}m_{heading}deg"
        except Exception:
            export_name = "Reference_0m_0deg"
        self.export_name_entry.setText(export_name)
        ref_export_name_layout.addWidget(self.export_name_entry, 0, 1)
        ref_import_export_layout.addWidget(ref_export_name_frame)
        
        ref_layout.addWidget(ref_import_export_groupbox, row, 0, 1, 2)
        ref_layout.setRowStretch(row, 0)
        row += 1
        
        # Add stretch at the bottom to push all groupboxes to the top
        ref_layout.setRowStretch(row, 1)
        
        self.dist_between_lines_entry.textChanged.connect(self._update_export_name)
        self.heading_entry.textChanged.connect(self._update_export_name)

        # --- Calibration Planning Tab ---
        cal_layout = QGridLayout(self.calibration_frame)
        cal_layout.setColumnStretch(0, 1)
        cal_layout.setColumnStretch(1, 2)
        cal_layout.setSpacing(3)  # Reduce spacing between groupboxes
        
        cal_row = 0
        # GeoTIFF controls moved to groupbox above tabs - removed from here
        
        # --- Calibration Line Parameters GroupBox ---
        cal_line_params_groupbox = QGroupBox("Calibration Line Parameters")
        cal_line_params_groupbox.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        cal_line_params_layout = QGridLayout(cal_line_params_groupbox)
        cal_line_params_layout.setSpacing(3)
        cal_line_params_layout.setContentsMargins(9, 9, 9, 9)
        cal_line_params_layout.setColumnStretch(0, 1)
        cal_line_params_layout.setColumnStretch(1, 2)
        
        cal_line_row = 0
        self.pick_pitch_line_btn = QPushButton("Draw a Pitch Line")
        self.pick_pitch_line_btn.clicked.connect(self._toggle_pick_pitch_line_mode)
        cal_line_params_layout.addWidget(self.pick_pitch_line_btn, cal_line_row, 0, 1, 2)
        cal_line_row += 1
        
        self.edit_pitch_line_btn = QPushButton("Edit Pitch Line")
        self.edit_pitch_line_btn.setEnabled(False)
        self.edit_pitch_line_btn.clicked.connect(self._toggle_edit_pitch_line_mode)
        cal_line_params_layout.addWidget(self.edit_pitch_line_btn, cal_line_row, 0, 1, 2)
        cal_line_row += 1
        
        cal_line_params_layout.addWidget(QLabel("Heading Line Offset (m):"), cal_line_row, 0)
        self.cal_line_offset_entry = QLineEdit()
        self.cal_line_offset_entry.textChanged.connect(self._update_cal_export_name_from_pitch_line)
        cal_line_params_layout.addWidget(self.cal_line_offset_entry, cal_line_row, 1)
        cal_line_row += 1
        
        self.add_heading_lines_btn = QPushButton("Add Heading Lines")
        self.add_heading_lines_btn.setEnabled(False)
        self.add_heading_lines_btn.clicked.connect(self._add_heading_lines_from_pitch_line)
        cal_line_params_layout.addWidget(self.add_heading_lines_btn, cal_line_row, 0, 1, 2)
        cal_line_row += 1
        
        self.pick_roll_line_btn = QPushButton("Draw a Roll Line")
        self.pick_roll_line_btn.clicked.connect(self._toggle_pick_roll_line_mode)
        cal_line_params_layout.addWidget(self.pick_roll_line_btn, cal_line_row, 0, 1, 2)
        cal_line_row += 1
        
        self.edit_roll_line_btn = QPushButton("Edit Roll Line")
        self.edit_roll_line_btn.setEnabled(False)
        self.edit_roll_line_btn.clicked.connect(self._toggle_edit_roll_line_mode)
        cal_line_params_layout.addWidget(self.edit_roll_line_btn, cal_line_row, 0, 1, 2)
        cal_line_row += 1
        
        self.clear_cal_lines_btn = QPushButton("Clear Lines")
        self.clear_cal_lines_btn.clicked.connect(self._clear_calibration_lines)
        cal_line_params_layout.addWidget(self.clear_cal_lines_btn, cal_line_row, 0, 1, 2)
        
        cal_layout.addWidget(cal_line_params_groupbox, cal_row, 0, 1, 2)
        cal_layout.setRowStretch(cal_row, 0)
        cal_row += 1
        
        # --- Plot Control GroupBox ---
        cal_plot_control_groupbox = QGroupBox("Plot Control")
        cal_plot_control_groupbox.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        cal_plot_control_layout = QVBoxLayout(cal_plot_control_groupbox)
        cal_plot_control_layout.setSpacing(0)
        cal_plot_control_layout.setContentsMargins(9, 9, 9, 9)
        
        self.zoom_to_all_lines_btn = QPushButton("Zoom to Calibration Lines")
        self.zoom_to_all_lines_btn.clicked.connect(self._zoom_to_any_lines)
        cal_plot_control_layout.addWidget(self.zoom_to_all_lines_btn)
        cal_plot_control_layout.addSpacing(3)
        
        self.zoom_to_geotiff_btn_cal = QPushButton("Zoom to GeoTIFF")
        self.zoom_to_geotiff_btn_cal.clicked.connect(self._zoom_to_geotiff)
        cal_plot_control_layout.addWidget(self.zoom_to_geotiff_btn_cal)
        cal_plot_control_layout.addSpacing(3)
        
        self.clear_plot_btn_cal = QPushButton("Clear Plot")
        self.clear_plot_btn_cal.clicked.connect(self._clear_plot)
        cal_plot_control_layout.addWidget(self.clear_plot_btn_cal)
        
        cal_layout.addWidget(cal_plot_control_groupbox, cal_row, 0, 1, 2)
        cal_layout.setRowStretch(cal_row, 0)
        cal_row += 1
        
        # --- Test Plan Info GroupBox ---
        cal_test_plan_info_groupbox = QGroupBox("Test Plan Info")
        cal_test_plan_info_groupbox.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        cal_test_plan_info_layout = QGridLayout(cal_test_plan_info_groupbox)
        cal_test_plan_info_layout.setSpacing(3)
        cal_test_plan_info_layout.setContentsMargins(9, 9, 9, 9)
        cal_test_plan_info_layout.setColumnStretch(0, 1)
        cal_test_plan_info_layout.setColumnStretch(1, 2)
        
        cal_test_plan_row = 0
        cal_test_plan_info_layout.addWidget(QLabel("Survey Speed (knots):"), cal_test_plan_row, 0)
        self.cal_survey_speed_entry = QLineEdit("8")
        cal_test_plan_info_layout.addWidget(self.cal_survey_speed_entry, cal_test_plan_row, 1)
        cal_test_plan_row += 1
        
        self.cal_show_stats_btn = QPushButton("Show Calibration Test Info")
        self.cal_show_stats_btn.clicked.connect(self._show_calibration_statistics)
        cal_test_plan_info_layout.addWidget(self.cal_show_stats_btn, cal_test_plan_row, 0, 1, 2)
        
        cal_layout.addWidget(cal_test_plan_info_groupbox, cal_row, 0, 1, 2)
        cal_layout.setRowStretch(cal_row, 0)
        cal_row += 1
        
        # --- Import/Export GroupBox ---
        cal_import_export_groupbox = QGroupBox("Import/Export")
        cal_import_export_groupbox.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        cal_import_export_layout = QVBoxLayout(cal_import_export_groupbox)
        cal_import_export_layout.setSpacing(0)
        cal_import_export_layout.setContentsMargins(9, 9, 9, 9)
        
        # Import/Export buttons on same row
        import_export_button_frame = QWidget()
        import_export_button_layout = QHBoxLayout(import_export_button_frame)
        import_export_button_layout.setSpacing(3)
        import_export_button_layout.setContentsMargins(0, 0, 0, 0)
        
        self.cal_import_survey_btn = QPushButton("Import Survey")
        self.cal_import_survey_btn.clicked.connect(self._import_cal_survey_files)
        import_export_button_layout.addWidget(self.cal_import_survey_btn)
        self.cal_export_survey_btn = QPushButton("Export Survey")
        self.cal_export_survey_btn.clicked.connect(self._export_cal_survey_files)
        import_export_button_layout.addWidget(self.cal_export_survey_btn)
        
        cal_import_export_layout.addWidget(import_export_button_frame)
        cal_import_export_layout.addSpacing(3)
        
        # Export Name at the bottom
        export_name_frame = QWidget()
        export_name_layout = QGridLayout(export_name_frame)
        export_name_layout.setSpacing(0)
        export_name_layout.setContentsMargins(0, 0, 0, 0)
        export_name_layout.setColumnStretch(0, 1)
        export_name_layout.setColumnStretch(1, 2)
        export_name_layout.addWidget(QLabel("Export Name:"), 0, 0)
        self.cal_export_name_entry = QLineEdit()
        export_name_layout.addWidget(self.cal_export_name_entry, 0, 1)
        cal_import_export_layout.addWidget(export_name_frame)
        
        cal_layout.addWidget(cal_import_export_groupbox, cal_row, 0, 1, 2)
        cal_layout.setRowStretch(cal_row, 0)
        cal_row += 1
        
        # Add stretch at the bottom to push all groupboxes to the top
        cal_layout.setRowStretch(cal_row, 1)

        # --- Line Planning Tab ---
        line_layout = QGridLayout(self.line_planning_frame)
        line_layout.setColumnStretch(0, 1)
        line_layout.setColumnStretch(1, 2)
        
        line_row = 0
        # GeoTIFF controls moved to groupbox above tabs - removed from here
        
        # --- Line Planning GroupBox ---
        line_planning_groupbox = QGroupBox("Line Planning")
        line_planning_groupbox_layout = QVBoxLayout(line_planning_groupbox)
        
        self.line_start_draw_btn = QPushButton("Start Drawing Line")
        self.line_start_draw_btn.clicked.connect(self._toggle_line_planning_mode)
        line_planning_groupbox_layout.addWidget(self.line_start_draw_btn)
        line_planning_groupbox_layout.addSpacing(3)
        
        self.line_edit_btn = QPushButton("Edit Line Planning")
        self.line_edit_btn.clicked.connect(self._toggle_edit_line_planning_mode)
        line_planning_groupbox_layout.addWidget(self.line_edit_btn)
        line_planning_groupbox_layout.addSpacing(3)
        
        self.line_clear_btn = QPushButton("Clear Line")
        self.line_clear_btn.clicked.connect(self._clear_line_planning)
        line_planning_groupbox_layout.addWidget(self.line_clear_btn)
        
        line_layout.addWidget(line_planning_groupbox, line_row, 0, 1, 2)
        line_layout.setRowStretch(line_row, 0)
        line_row += 1
        
        # --- Plot Control GroupBox ---
        line_plot_control_groupbox = QGroupBox("Plot Control")
        line_plot_control_groupbox.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        line_plot_control_layout = QVBoxLayout(line_plot_control_groupbox)
        line_plot_control_layout.setSpacing(0)
        line_plot_control_layout.setContentsMargins(9, 9, 9, 9)
        
        self.zoom_to_line_btn = QPushButton("Zoom to Line")
        self.zoom_to_line_btn.clicked.connect(self._zoom_to_line)
        line_plot_control_layout.addWidget(self.zoom_to_line_btn)
        line_plot_control_layout.addSpacing(3)
        
        self.zoom_to_geotiff_btn_line = QPushButton("Zoom to GeoTIFF")
        self.zoom_to_geotiff_btn_line.clicked.connect(self._zoom_to_geotiff)
        line_plot_control_layout.addWidget(self.zoom_to_geotiff_btn_line)
        
        line_layout.addWidget(line_plot_control_groupbox, line_row, 0, 1, 2)
        line_layout.setRowStretch(line_row, 0)
        line_row += 1
        
        # --- Test Plan Info GroupBox ---
        line_test_plan_info_groupbox = QGroupBox("Test Plan Info")
        line_test_plan_info_groupbox.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        line_test_plan_info_layout = QGridLayout(line_test_plan_info_groupbox)
        line_test_plan_info_layout.setSpacing(3)
        line_test_plan_info_layout.setContentsMargins(9, 9, 9, 9)
        line_test_plan_info_layout.setColumnStretch(0, 1)
        line_test_plan_info_layout.setColumnStretch(1, 2)
        
        line_test_plan_row = 0
        line_test_plan_info_layout.addWidget(QLabel("Survey Speed (knots):"), line_test_plan_row, 0)
        self.line_survey_speed_entry = QLineEdit("8")
        line_test_plan_info_layout.addWidget(self.line_survey_speed_entry, line_test_plan_row, 1)
        line_test_plan_row += 1
        
        self.line_show_info_btn = QPushButton("Show Survey Info")
        self.line_show_info_btn.clicked.connect(self._show_line_information)
        line_test_plan_info_layout.addWidget(self.line_show_info_btn, line_test_plan_row, 0, 1, 2)
        
        line_layout.addWidget(line_test_plan_info_groupbox, line_row, 0, 1, 2)
        line_layout.setRowStretch(line_row, 0)
        line_row += 1
        
        # --- Import/Export GroupBox ---
        line_import_export_groupbox = QGroupBox("Import/Export")
        line_import_export_groupbox.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        line_import_export_layout = QVBoxLayout(line_import_export_groupbox)
        line_import_export_layout.setSpacing(0)
        line_import_export_layout.setContentsMargins(9, 9, 9, 9)
        
        # Import/Export buttons on same row
        line_import_export_button_frame = QWidget()
        line_import_export_button_layout = QHBoxLayout(line_import_export_button_frame)
        line_import_export_button_layout.setSpacing(3)
        line_import_export_button_layout.setContentsMargins(0, 0, 0, 0)
        
        self.line_import_btn = QPushButton("Import Survey")
        self.line_import_btn.clicked.connect(self._import_drawn_line)
        line_import_export_button_layout.addWidget(self.line_import_btn)
        self.line_export_btn = QPushButton("Export Survey")
        self.line_export_btn.clicked.connect(self._export_drawn_line)
        line_import_export_button_layout.addWidget(self.line_export_btn)
        
        line_import_export_layout.addWidget(line_import_export_button_frame)
        line_import_export_layout.addSpacing(3)
        
        # Export Name at the bottom
        line_export_name_frame = QWidget()
        line_export_name_layout = QGridLayout(line_export_name_frame)
        line_export_name_layout.setSpacing(0)
        line_export_name_layout.setContentsMargins(0, 0, 0, 0)
        line_export_name_layout.setColumnStretch(0, 1)
        line_export_name_layout.setColumnStretch(1, 2)
        line_export_name_layout.addWidget(QLabel("Export Name:"), 0, 0)
        self.line_export_name_entry = QLineEdit()
        default_export_name = f"Line_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
        self.line_export_name_entry.setText(default_export_name)
        line_export_name_layout.addWidget(self.line_export_name_entry, 0, 1)
        line_import_export_layout.addWidget(line_export_name_frame)
        
        line_layout.addWidget(line_import_export_groupbox, line_row, 0, 1, 2)
        line_layout.setRowStretch(line_row, 0)
        line_row += 1
        
        # Add stretch at the bottom to push all groupboxes to the top
        line_layout.setRowStretch(line_row, 1)

        # --- 2. Main Plot Area (right side) ---
        self.plot_frame = QWidget()
        self.plot_frame.setMinimumSize(400, 300)  # Ensure minimum size
        plot_layout = QVBoxLayout(self.plot_frame)
        plot_layout.setContentsMargins(0, 0, 0, 0)
        self.canvas = FigureCanvas(self.figure)
        plot_layout.addWidget(self.canvas)
        self.canvas_widget = self.canvas  # For compatibility
        
        print(f"Widgets created - param_scroll: {self.param_scroll is not None}, plot_frame: {self.plot_frame is not None}, canvas: {self.canvas is not None}")

        # Add stub for _on_parameter_change if missing
        if not hasattr(self, '_on_parameter_change'):
            def _on_parameter_change(event=None):
                pass
            self._on_parameter_change = _on_parameter_change

    # Add a method to draw the elevation profile for the drawn line
    def _draw_line_planning_profile(self):
        self.profile_ax.clear()
        # Clear any existing twin axes (slope axes)
        for ax in self.profile_fig.get_axes():
            if ax != self.profile_ax:
                ax.remove()
        self.profile_ax.set_title("Line Planning Elevation Profile", fontsize=8)
        self.profile_ax.set_xlabel("Distance (m)", fontsize=8)
        self.profile_ax.set_ylabel("Elevation (m)", fontsize=8)
        self.profile_ax.tick_params(axis='both', which='major', labelsize=7)
        slope_ax = None
        if (
            hasattr(self, 'geotiff_data_array') and self.geotiff_data_array is not None and
            hasattr(self, 'geotiff_extent') and self.geotiff_extent is not None and
            hasattr(self, 'line_planning_points') and len(self.line_planning_points) >= 2
        ):
            # Interpolate points along the line
            import numpy as np
            import pyproj
            geod = pyproj.Geod(ellps="WGS84")
            lats = []
            lons = []
            dists = [0.0]
            total_dist = 0.0
            # For mapping waypoints to profile distances
            waypoint_distances = [0.0]
            for i in range(1, len(self.line_planning_points)):
                lat1, lon1 = self.line_planning_points[i-1]
                lat2, lon2 = self.line_planning_points[i]
                # Interpolate 50 points per segment
                seg_lats = np.linspace(lat1, lat2, 50)
                seg_lons = np.linspace(lon1, lon2, 50)
                if i == 1:
                    lats.extend(seg_lats)
                    lons.extend(seg_lons)
                else:
                    lats.extend(seg_lats[1:])
                    lons.extend(seg_lons[1:])
                # Compute cumulative distance
                for j in range(1, 50):
                    _, _, d = geod.inv(seg_lons[j-1], seg_lats[j-1], seg_lons[j], seg_lats[j])
                    if not np.isnan(d) and not np.isinf(d):
                        total_dist += d
                    dists.append(total_dist)
                waypoint_distances.append(total_dist)
            lats = np.array(lats)
            lons = np.array(lons)
            dists = np.array(dists)
            # Sample elevations and slopes
            left, right, bottom, top = tuple(self.geotiff_extent)
            nrows, ncols = self.geotiff_data_array.shape
            rows = ((top - lats) / (top - bottom) * (nrows - 1)).clip(0, nrows - 1)
            cols = ((lons - left) / (right - left) * (ncols - 1)).clip(0, ncols - 1)
            elevations = []
            slopes = []
            for idx, (r, c) in enumerate(zip(rows, cols)):
                ir, ic = int(round(r)), int(round(c))
                elev = self.geotiff_data_array[ir, ic]
                elevations.append(elev)
                # Slope calculation (in degrees)
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
                    if window.shape == (3,3) and not np.all(np.isnan(window)):
                        dz_dy, dz_dx = np.gradient(window, dy_m, dx_m)
                        slope_rad = np.arctan(np.sqrt(dz_dx[1,1]**2 + dz_dy[1,1]**2))
                        slope = np.degrees(slope_rad)
                slopes.append(slope if slope is not None else np.nan)
            elevations = np.array(elevations)
            slopes = np.array(slopes)
            self.profile_ax.plot(dists, elevations, color='orange', lw=1, label='Elevation')
            # Plot waypoints as circles on the elevation profile
            waypoint_elevations = []
            for i, (lat, lon) in enumerate(self.line_planning_points):
                # Find the closest point in the interpolated profile
                _, idx = min((abs(lat - lats[j]) + abs(lon - lons[j]), j) for j in range(len(lats)))
                waypoint_elevations.append(elevations[idx])
            self.profile_ax.plot(waypoint_distances, waypoint_elevations, 'o', color='red', markersize=6, label='Waypoints')
            if np.any(~np.isnan(slopes)):
                slope_ax = self.profile_ax.twinx()
                slope_ax.plot(dists, slopes, color='blue', lw=1, linestyle='--', label='Slope (deg)')
                slope_ax.set_ylabel('Slope (deg)', fontsize=8)
                slope_ax.tick_params(axis='y', labelsize=7)
                slope_ax.grid(False)
            # Check for valid limits before setting
            if len(dists) > 0 and not (np.isnan(dists[0]) or np.isinf(dists[0]) or np.isnan(dists[-1]) or np.isinf(dists[-1])):
                self.profile_ax.set_xlim(dists[0], dists[-1])
            if np.any(~np.isnan(elevations)):
                self.profile_ax.set_ylim(np.nanmin(elevations), np.nanmax(elevations))
            if slope_ax and np.any(~np.isnan(slopes)):
                slope_ax.set_ylim(0, np.nanmax(slopes[~np.isnan(slopes)])*1.1)
        
        # Add legend to profile plot (include slope axis if present)
        handles, labels = self.profile_ax.get_legend_handles_labels()
        if slope_ax:
            slope_handles, slope_labels = slope_ax.get_legend_handles_labels()
            handles.extend(slope_handles)
            labels.extend(slope_labels)
        if handles:
            self.profile_ax.legend(handles, labels, fontsize=10, loc='upper right')
        
        self.profile_fig.tight_layout(pad=1.0)
        self.profile_canvas.draw_idle()

    def set_line_info_text(self, message, append=False):
        """Add a message to the Activity Log with [Line] prefix. Maintains up to 200 lines of history."""
        if not hasattr(self, 'activity_log_text'):
            return
        prefixed_message = f"[Line] {message}"
        self.activity_log_text.setReadOnly(False)
        if append:
            # Prepend new message (newest at top) - only when explicitly requested
            current_text = self.activity_log_text.toPlainText().rstrip('\n')
            new_text = prefixed_message + "\n" + current_text if current_text else prefixed_message + "\n"
            # Split into lines and limit to 200
            lines = new_text.splitlines()
            if len(lines) > 200:
                lines = lines[:200]
            self.activity_log_text.setPlainText("\n".join(lines) + "\n")
            # Move cursor to top
            cursor = self.activity_log_text.textCursor()
            cursor.movePosition(QTextCursor.MoveOperation.Start)
            self.activity_log_text.setTextCursor(cursor)
        else:
            # Append to end (default behavior - newest at bottom)
            self.activity_log_text.append(prefixed_message)
            # Maintain up to 200 lines by removing oldest from top
            num_lines = self.activity_log_text.document().blockCount()
            if num_lines > 200:
                cursor = self.activity_log_text.textCursor()
                cursor.movePosition(QTextCursor.MoveOperation.Start)
                cursor.movePosition(QTextCursor.MoveOperation.Down, QTextCursor.MoveMode.MoveAnchor, num_lines - 200)
                cursor.movePosition(QTextCursor.MoveOperation.Start, QTextCursor.MoveMode.KeepAnchor)
                cursor.removeSelectedText()
            # Move cursor to end and scroll to show newest message
            cursor = self.activity_log_text.textCursor()
            cursor.movePosition(QTextCursor.MoveOperation.End)
            self.activity_log_text.setTextCursor(cursor)
        self.activity_log_text.ensureCursorVisible()
        self.activity_log_text.setReadOnly(True)

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
            # We want: (width * geographic_aspect) / new_height = figure_aspect
            # So: new_height = (width * geographic_aspect) / figure_aspect
            new_height = (base_width * geographic_aspect) / figure_aspect
            new_width = base_width  # Keep width at GeoTIFF size
        else:
            # Display would be taller than figure - expand width to fill, keep height at GeoTIFF size
            # We want: (new_width * geographic_aspect) / height = figure_aspect
            # So: new_width = (height * figure_aspect) / geographic_aspect
            new_width = (base_height * figure_aspect) / geographic_aspect
            new_height = base_height  # Keep height at GeoTIFF size
        
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

    def _calculate_total_survey_time(self):
        """Calculate total survey time including travel between lines and crossline."""
        try:
            geod = pyproj.Geod(ellps="WGS84")
            speed_knots = float(self.survey_speed_entry.text()) if self.survey_speed_entry.text() else 8.0
            speed_m_per_h = speed_knots * 1852
            
            # Initialize distance calculations
            main_lines_total_distance = 0
            travel_between_lines_total_distance = 0
            travel_to_crossline_distance = 0
            crossline_single_pass_distance = 0
            
            # Calculate main lines survey time and distance (all main lines)
            main_lines_survey_time_minutes = 0
            for line in self.survey_lines_data:
                lat1, lon1 = line[0]
                lat2, lon2 = line[1]
                _, _, line_length = geod.inv(lon1, lat1, lon2, lat2)
                main_lines_total_distance += line_length
                line_time_hours = line_length / speed_m_per_h if speed_m_per_h > 0 else 0
                main_lines_survey_time_minutes += line_time_hours * 60
            
            # Calculate travel time and distance between main lines (zigzag pattern)
            travel_between_lines_minutes = 0
            if len(self.survey_lines_data) > 1:
                for i in range(1, len(self.survey_lines_data)):
                    # For zigzag pattern: from end of line i-1 to start of line i
                    # Even lines (0,2,4...) use normal order, odd lines (1,3,5...) use flipped order
                    if (i-1) % 2 == 0:  # Even line - end point
                        end_lat, end_lon = self.survey_lines_data[i-1][1]
                    else:  # Odd line - start point (flipped)
                        end_lat, end_lon = self.survey_lines_data[i-1][0]
                    
                    if i % 2 == 0:  # Even line - start point
                        start_lat, start_lon = self.survey_lines_data[i][0]
                    else:  # Odd line - end point (flipped)
                        start_lat, start_lon = self.survey_lines_data[i][1]
                    
                    _, _, travel_distance = geod.inv(end_lon, end_lat, start_lon, start_lat)
                    travel_between_lines_total_distance += travel_distance
                    travel_time_hours = travel_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                    travel_between_lines_minutes += travel_time_hours * 60
            
            # Calculate travel time and distance from last main line to crossline
            travel_to_crossline_minutes = 0
            if self.cross_line_data and self.survey_lines_data:
                # Find the end point of the last main line (according to zigzag order)
                last_line_index = len(self.survey_lines_data) - 1
                if last_line_index % 2 == 0:  # Even line - end point
                    last_end_lat, last_end_lon = self.survey_lines_data[last_line_index][1]
                else:  # Odd line - start point (flipped)
                    last_end_lat, last_end_lon = self.survey_lines_data[last_line_index][0]
                
                # Find the crossline endpoint closest to the last main line endpoint
                crossline_start_lat, crossline_start_lon = self.cross_line_data[0]
                crossline_end_lat, crossline_end_lon = self.cross_line_data[1]
                
                # Calculate distance to both crossline endpoints
                _, _, dist_to_crossline_start = geod.inv(last_end_lon, last_end_lat, crossline_start_lon, crossline_start_lat)
                _, _, dist_to_crossline_end = geod.inv(last_end_lon, last_end_lat, crossline_end_lon, crossline_end_lat)
                
                # Use the closer crossline endpoint
                if dist_to_crossline_start <= dist_to_crossline_end:
                    crossline_target_lat, crossline_target_lon = crossline_start_lat, crossline_start_lon
                    travel_distance = dist_to_crossline_start
                else:
                    crossline_target_lat, crossline_target_lon = crossline_end_lat, crossline_end_lon
                    travel_distance = dist_to_crossline_end
                
                travel_to_crossline_distance = travel_distance
                travel_time_hours = travel_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                travel_to_crossline_minutes = travel_time_hours * 60
                
                # Debug output for travel to crossline calculation
                print(f"DEBUG: Travel to Crossline - Last line end: ({last_end_lat:.6f}, {last_end_lon:.6f})")
                print(f"DEBUG: Travel to Crossline - Crossline target: ({crossline_target_lat:.6f}, {crossline_target_lon:.6f})")
                print(f"DEBUG: Travel to Crossline - Distance: {travel_distance:.1f} m")
            
            # Calculate crossline survey time and distance
            crossline_survey_minutes = 0
            num_passes = 1
            if self.cross_line_data:
                lat1, lon1 = self.cross_line_data[0]
                lat2, lon2 = self.cross_line_data[1]
                _, _, crossline_length = geod.inv(lon1, lat1, lon2, lat2)
                crossline_single_pass_distance = crossline_length
                crossline_time_hours = crossline_length / speed_m_per_h if speed_m_per_h > 0 else 0
                crossline_survey_minutes = crossline_time_hours * 60
                
                # Multiply by number of crossline passes
                try:
                    num_passes = int(self.crossline_passes_entry.text()) if self.crossline_passes_entry.text() else 2
                    crossline_survey_minutes *= num_passes
                except (ValueError, AttributeError):
                    # If the entry doesn't exist or has invalid value, default to 2 passes
                    num_passes = 2
                    crossline_survey_minutes *= num_passes
            
            # Calculate total time and distance
            total_minutes = (main_lines_survey_time_minutes + 
                           travel_between_lines_minutes + 
                           travel_to_crossline_minutes + 
                           crossline_survey_minutes)
            total_hours = total_minutes / 60
            
            # Calculate total distance
            total_distance = (main_lines_total_distance + 
                            travel_between_lines_total_distance + 
                            travel_to_crossline_distance + 
                            (crossline_single_pass_distance * num_passes))
            
            return {
                'total_minutes': total_minutes,
                'total_hours': total_hours,
                'main_lines_minutes': main_lines_survey_time_minutes,
                'crossline_minutes': crossline_survey_minutes,
                'travel_minutes': travel_between_lines_minutes,
                'travel_to_crossline_minutes': travel_to_crossline_minutes,
                'total_distance_m': total_distance,
                'total_distance_km': total_distance / 1000.0,
                'total_distance_nm': total_distance / 1852.0,
                'main_lines_distance_m': main_lines_total_distance,
                'main_lines_distance_km': main_lines_total_distance / 1000.0,
                'main_lines_distance_nm': main_lines_total_distance / 1852.0,
                'travel_between_lines_distance_m': travel_between_lines_total_distance,
                'travel_between_lines_distance_km': travel_between_lines_total_distance / 1000.0,
                'travel_between_lines_distance_nm': travel_between_lines_total_distance / 1852.0,
                'travel_to_crossline_distance_m': travel_to_crossline_distance,
                'travel_to_crossline_distance_km': travel_to_crossline_distance / 1000.0,
                'travel_to_crossline_distance_nm': travel_to_crossline_distance / 1852.0,
                'crossline_single_pass_distance_m': crossline_single_pass_distance,
                'crossline_single_pass_distance_km': crossline_single_pass_distance / 1000.0,
                'crossline_single_pass_distance_nm': crossline_single_pass_distance / 1852.0,
                'crossline_total_distance_m': crossline_single_pass_distance * num_passes,
                'crossline_total_distance_km': (crossline_single_pass_distance * num_passes) / 1000.0,
                'crossline_total_distance_nm': (crossline_single_pass_distance * num_passes) / 1852.0,
                'num_crossline_passes': num_passes
            }
        except Exception as e:
            print(f"Error calculating total survey time: {e}")
            return {
                'total_minutes': 0,
                'total_hours': 0,
                'main_lines_minutes': 0,
                'crossline_minutes': 0,
                'travel_minutes': 0,
                'travel_to_crossline_minutes': 0,
                'total_distance_m': 0,
                'total_distance_km': 0,
                'total_distance_nm': 0,
                'main_lines_distance_m': 0,
                'main_lines_distance_km': 0,
                'main_lines_distance_nm': 0,
                'travel_between_lines_distance_m': 0,
                'travel_between_lines_distance_km': 0,
                'travel_between_lines_distance_nm': 0,
                'travel_to_crossline_distance_m': 0,
                'travel_to_crossline_distance_km': 0,
                'travel_to_crossline_distance_nm': 0,
                'crossline_single_pass_distance_m': 0,
                'crossline_single_pass_distance_km': 0,
                'crossline_single_pass_distance_nm': 0,
                'crossline_total_distance_m': 0,
                'crossline_total_distance_km': 0,
                'crossline_total_distance_nm': 0,
                'num_crossline_passes': 1
            }

    def _calculate_calibration_survey_statistics(self):
        """Calculate comprehensive calibration survey statistics including pitch, roll, heading lines and travel distances."""
        try:
            import pyproj
            geod = pyproj.Geod(ellps="WGS84")
            speed_knots = float(self.cal_survey_speed_entry.text()) if self.cal_survey_speed_entry.text() else 8.0
            speed_m_per_h = speed_knots * 1852
            
            # Initialize statistics
            stats = {
                'pitch_line_distance_m': 0.0,
                'pitch_line_time_min': 0.0,
                'roll_line_distance_m': 0.0,
                'roll_line_time_min': 0.0,
                'heading1_distance_m': 0.0,
                'heading1_time_min': 0.0,
                'heading2_distance_m': 0.0,
                'heading2_time_min': 0.0,
                'travel_pitch_to_roll_m': 0.0,
                'travel_pitch_to_roll_min': 0.0,
                'travel_roll_to_heading1_m': 0.0,
                'travel_roll_to_heading1_min': 0.0,
                'travel_heading1_to_heading2_m': 0.0,
                'travel_heading1_to_heading2_min': 0.0,
                'total_distance_m': 0.0,
                'total_time_min': 0.0,
                'total_time_hours': 0.0
            }
            
            # Calculate pitch line statistics (2 passes)
            if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2:
                (lat1, lon1), (lat2, lon2) = self.pitch_line_points
                _, _, pitch_distance = geod.inv(lon1, lat1, lon2, lat2)
                stats['pitch_line_distance_m'] = pitch_distance
                pitch_time_hours = pitch_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                stats['pitch_line_time_min'] = pitch_time_hours * 60 * 2  # 2 passes
            
            # Calculate roll line statistics (2 passes)
            if hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2:
                (lat1, lon1), (lat2, lon2) = self.roll_line_points
                _, _, roll_distance = geod.inv(lon1, lat1, lon2, lat2)
                stats['roll_line_distance_m'] = roll_distance
                roll_time_hours = roll_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                stats['roll_line_time_min'] = roll_time_hours * 60 * 2  # 2 passes
            
            # Calculate heading line 1 statistics (1 pass)
            if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 1:
                (lat1, lon1), (lat2, lon2) = self.heading_lines[0]
                _, _, heading1_distance = geod.inv(lon1, lat1, lon2, lat2)
                stats['heading1_distance_m'] = heading1_distance
                heading1_time_hours = heading1_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                stats['heading1_time_min'] = heading1_time_hours * 60
            
            # Calculate heading line 2 statistics (1 pass)
            if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 2:
                (lat1, lon1), (lat2, lon2) = self.heading_lines[1]
                _, _, heading2_distance = geod.inv(lon1, lat1, lon2, lat2)
                stats['heading2_distance_m'] = heading2_distance
                heading2_time_hours = heading2_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                stats['heading2_time_min'] = heading2_time_hours * 60
            
            # Calculate travel distances and times
            # Travel from pitch line start to roll line start
            if (hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2 and
                hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2):
                pitch_start_lat, pitch_start_lon = self.pitch_line_points[0]
                roll_start_lat, roll_start_lon = self.roll_line_points[0]
                _, _, travel_distance = geod.inv(pitch_start_lon, pitch_start_lat, roll_start_lon, roll_start_lat)
                stats['travel_pitch_to_roll_m'] = travel_distance
                travel_time_hours = travel_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                stats['travel_pitch_to_roll_min'] = travel_time_hours * 60
            
            # Travel from roll line start to heading line 1 start
            if (hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2 and
                hasattr(self, 'heading_lines') and len(self.heading_lines) >= 1):
                roll_start_lat, roll_start_lon = self.roll_line_points[0]
                heading1_start_lat, heading1_start_lon = self.heading_lines[0][0]
                _, _, travel_distance = geod.inv(roll_start_lon, roll_start_lat, heading1_start_lon, heading1_start_lat)
                stats['travel_roll_to_heading1_m'] = travel_distance
                travel_time_hours = travel_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                stats['travel_roll_to_heading1_min'] = travel_time_hours * 60
            
            # Travel from heading line 1 end to heading line 2 start
            if (hasattr(self, 'heading_lines') and len(self.heading_lines) >= 2):
                heading1_end_lat, heading1_end_lon = self.heading_lines[0][1]
                heading2_start_lat, heading2_start_lon = self.heading_lines[1][0]
                _, _, travel_distance = geod.inv(heading1_end_lon, heading1_end_lat, heading2_start_lon, heading2_start_lat)
                stats['travel_heading1_to_heading2_m'] = travel_distance
                travel_time_hours = travel_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                stats['travel_heading1_to_heading2_min'] = travel_time_hours * 60
            
            # Calculate totals
            stats['total_distance_m'] = (
                stats['pitch_line_distance_m'] * 2 +  # 2 passes
                stats['roll_line_distance_m'] * 2 +  # 2 passes
                stats['heading1_distance_m'] +  # 1 pass
                stats['heading2_distance_m'] +  # 1 pass
                stats['travel_pitch_to_roll_m'] +
                stats['travel_roll_to_heading1_m'] +
                stats['travel_heading1_to_heading2_m']
            )
            
            stats['total_time_min'] = (
                stats['pitch_line_time_min'] +
                stats['roll_line_time_min'] +
                stats['heading1_time_min'] +
                stats['heading2_time_min'] +
                stats['travel_pitch_to_roll_min'] +
                stats['travel_roll_to_heading1_min'] +
                stats['travel_heading1_to_heading2_min']
            )
            
            stats['total_time_hours'] = stats['total_time_min'] / 60.0
            
            # Add derived units
            stats['total_distance_km'] = stats['total_distance_m'] / 1000.0
            stats['total_distance_nm'] = stats['total_distance_m'] / 1852.0
            stats['pitch_line_distance_km'] = stats['pitch_line_distance_m'] / 1000.0
            stats['pitch_line_distance_nm'] = stats['pitch_line_distance_m'] / 1852.0
            stats['roll_line_distance_km'] = stats['roll_line_distance_m'] / 1000.0
            stats['roll_line_distance_nm'] = stats['roll_line_distance_m'] / 1852.0
            stats['heading1_distance_km'] = stats['heading1_distance_m'] / 1000.0
            stats['heading1_distance_nm'] = stats['heading1_distance_m'] / 1852.0
            stats['heading2_distance_km'] = stats['heading2_distance_m'] / 1000.0
            stats['heading2_distance_nm'] = stats['heading2_distance_m'] / 1852.0
            stats['travel_pitch_to_roll_km'] = stats['travel_pitch_to_roll_m'] / 1000.0
            stats['travel_pitch_to_roll_nm'] = stats['travel_pitch_to_roll_m'] / 1852.0
            stats['travel_roll_to_heading1_km'] = stats['travel_roll_to_heading1_m'] / 1000.0
            stats['travel_roll_to_heading1_nm'] = stats['travel_roll_to_heading1_m'] / 1852.0
            stats['travel_heading1_to_heading2_km'] = stats['travel_heading1_to_heading2_m'] / 1000.0
            stats['travel_heading1_to_heading2_nm'] = stats['travel_heading1_to_heading2_m'] / 1852.0
            
            return stats
            
        except Exception as e:
            print(f"Error calculating calibration survey statistics: {e}")
            return None

    def _show_calibration_statistics(self):
        """Display comprehensive calibration survey statistics in a message box."""
        stats = self._calculate_calibration_survey_statistics()
        if not stats:
            self._show_message("warning","Statistics Error", "Unable to calculate calibration survey statistics. Please ensure all required lines are drawn.")
            return
        
        # Format the statistics for display
        stats_text = "COMPREHENSIVE CALIBRATION SURVEY STATISTICS\n"
        stats_text += "=" * 50 + "\n\n"
        
        # Survey parameters
        try:
            speed_knots = float(self.cal_survey_speed_entry.text()) if self.cal_survey_speed_entry.text() else 8.0
            stats_text += f"Survey Speed: {speed_knots} knots\n\n"
        except:
            stats_text += "Survey Speed: 8.0 knots (default)\n\n"
        
        # Line distances and times
        stats_text += "SURVEY LINE DISTANCES AND TIMES\n"
        stats_text += "-" * 35 + "\n"
        
        if stats['pitch_line_distance_m'] > 0:
            stats_text += f"Pitch Line (2 passes):\n"
            # Add pitch line heading information
            if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2:
                try:
                    if pyproj is not None:
                        geod = pyproj.Geod(ellps="WGS84")
                        (lat1, lon1), (lat2, lon2) = self.pitch_line_points
                        fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)
                        
                        pitch_heading = fwd_az % 360
                        pitch_reciprocal_heading = back_az % 360
                        
                        stats_text += f"  Heading: {pitch_heading:.1f}°\n"
                        stats_text += f"  Reciprocal Heading: {pitch_reciprocal_heading:.1f}°\n"
                    else:
                        stats_text += f"  Heading: pyproj not available\n"
                except Exception as e:
                    stats_text += f"  Heading: Unable to calculate\n"
            stats_text += f"  Distance (per pass): {stats['pitch_line_distance_m']:.1f} m ({stats['pitch_line_distance_km']:.3f} km, {stats['pitch_line_distance_nm']:.3f} nm)\n"
            stats_text += f"  Time (per line): {stats['pitch_line_time_min']/2:.1f} min\n"
            stats_text += f"  Time (total): {stats['pitch_line_time_min']:.1f} min\n\n"
        
        if stats['roll_line_distance_m'] > 0:
            stats_text += f"Roll Line (2 passes):\n"
            # Add roll line heading information
            if hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2:
                try:
                    if pyproj is not None:
                        geod = pyproj.Geod(ellps="WGS84")
                        (lat1, lon1), (lat2, lon2) = self.roll_line_points
                        fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)
                        
                        roll_heading = fwd_az % 360
                        roll_reciprocal_heading = back_az % 360
                        
                        stats_text += f"  Heading: {roll_heading:.1f}°\n"
                        stats_text += f"  Reciprocal Heading: {roll_reciprocal_heading:.1f}°\n"
                    else:
                        stats_text += f"  Heading: pyproj not available\n"
                except Exception as e:
                    stats_text += f"  Heading: Unable to calculate\n"
            stats_text += f"  Distance (per pass): {stats['roll_line_distance_m']:.1f} m ({stats['roll_line_distance_km']:.3f} km, {stats['roll_line_distance_nm']:.3f} nm)\n"
            stats_text += f"  Time (per line): {stats['roll_line_time_min']/2:.1f} min\n"
            stats_text += f"  Time (total): {stats['roll_line_time_min']:.1f} min\n\n"
        
        if stats['heading1_distance_m'] > 0:
            stats_text += f"Heading Line 1 (1 pass):\n"
            # Add heading line 1 heading information
            if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 1:
                try:
                    if pyproj is not None:
                        geod = pyproj.Geod(ellps="WGS84")
                        (lat1, lon1), (lat2, lon2) = self.heading_lines[0]
                        fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)
                        
                        heading1_heading = fwd_az % 360
                        heading1_reciprocal_heading = back_az % 360
                        
                        stats_text += f"  Heading: {heading1_heading:.1f}°\n"
                        stats_text += f"  Reciprocal Heading: {heading1_reciprocal_heading:.1f}°\n"
                    else:
                        stats_text += f"  Heading: pyproj not available\n"
                except Exception as e:
                    stats_text += f"  Heading: Unable to calculate\n"
            stats_text += f"  Distance: {stats['heading1_distance_m']:.1f} m ({stats['heading1_distance_km']:.3f} km, {stats['heading1_distance_nm']:.3f} nm)\n"
            stats_text += f"  Time: {stats['heading1_time_min']:.1f} min\n\n"
        
        if stats['heading2_distance_m'] > 0:
            stats_text += f"Heading Line 2 (1 pass):\n"
            # Add heading line 2 heading information
            if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 2:
                try:
                    if pyproj is not None:
                        geod = pyproj.Geod(ellps="WGS84")
                        (lat1, lon1), (lat2, lon2) = self.heading_lines[1]
                        fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)
                        
                        heading2_heading = fwd_az % 360
                        heading2_reciprocal_heading = back_az % 360
                        
                        stats_text += f"  Heading: {heading2_heading:.1f}°\n"
                        stats_text += f"  Reciprocal Heading: {heading2_reciprocal_heading:.1f}°\n"
                    else:
                        stats_text += f"  Heading: pyproj not available\n"
                except Exception as e:
                    stats_text += f"  Heading: Unable to calculate\n"
            stats_text += f"  Distance: {stats['heading2_distance_m']:.1f} m ({stats['heading2_distance_km']:.3f} km, {stats['heading2_distance_nm']:.3f} nm)\n"
            stats_text += f"  Time: {stats['heading2_time_min']:.1f} min\n\n"
        
        # Travel distances and times
        stats_text += "TRAVEL DISTANCES AND TIMES\n"
        stats_text += "-" * 30 + "\n"
        
        # Combined line execution times and distances
        if stats['pitch_line_distance_m'] > 0:
            stats_text += f"Pitch Lines (2 passes):\n"
            stats_text += f"  Total Distance: {stats['pitch_line_distance_m'] * 2:.1f} m ({stats['pitch_line_distance_km'] * 2:.3f} km, {stats['pitch_line_distance_nm'] * 2:.3f} nm)\n"
            stats_text += f"  Total Time: {stats['pitch_line_time_min']:.1f} min\n\n"
        
        if stats['roll_line_distance_m'] > 0:
            stats_text += f"Roll Lines (2 passes):\n"
            stats_text += f"  Total Distance: {stats['roll_line_distance_m'] * 2:.1f} m ({stats['roll_line_distance_km'] * 2:.3f} km, {stats['roll_line_distance_nm'] * 2:.3f} nm)\n"
            stats_text += f"  Total Time: {stats['roll_line_time_min']:.1f} min\n\n"
        
        if stats['heading1_distance_m'] > 0 or stats['heading2_distance_m'] > 0:
            total_heading_distance = stats['heading1_distance_m'] + stats['heading2_distance_m']
            total_heading_time = stats['heading1_time_min'] + stats['heading2_time_min']
            total_heading_distance_km = total_heading_distance / 1000
            total_heading_distance_nm = total_heading_distance / 1852
            
            stats_text += f"Heading Lines (2 lines):\n"
            stats_text += f"  Total Distance: {total_heading_distance:.1f} m ({total_heading_distance_km:.3f} km, {total_heading_distance_nm:.3f} nm)\n"
            stats_text += f"  Total Time: {total_heading_time:.1f} min\n\n"
        
        if stats['travel_pitch_to_roll_m'] > 0:
            stats_text += f"Pitch Line Start → Roll Line Start:\n"
            stats_text += f"  Distance: {stats['travel_pitch_to_roll_m']:.1f} m ({stats['travel_pitch_to_roll_km']:.3f} km, {stats['travel_pitch_to_roll_nm']:.3f} nm)\n"
            stats_text += f"  Time: {stats['travel_pitch_to_roll_min']:.1f} min\n\n"
        
        if stats['travel_roll_to_heading1_m'] > 0:
            stats_text += f"Roll Line Start → Heading Line 1 Start:\n"
            stats_text += f"  Distance: {stats['travel_roll_to_heading1_m']:.1f} m ({stats['travel_roll_to_heading1_km']:.3f} km, {stats['travel_roll_to_heading1_nm']:.3f} nm)\n"
            stats_text += f"  Time: {stats['travel_roll_to_heading1_min']:.1f} min\n\n"
        
        if stats['travel_heading1_to_heading2_m'] > 0:
            stats_text += f"Heading Line 1 End → Heading Line 2 Start:\n"
            stats_text += f"  Distance: {stats['travel_heading1_to_heading2_m']:.1f} m ({stats['travel_heading1_to_heading2_km']:.3f} km, {stats['travel_heading1_to_heading2_nm']:.3f} nm)\n"
            stats_text += f"  Time: {stats['travel_heading1_to_heading2_min']:.1f} min\n\n"
        
        # Totals
        stats_text += "TOTAL SURVEY SUMMARY\n"
        stats_text += "-" * 25 + "\n"
        stats_text += f"Total Distance: {stats['total_distance_m']:.1f} m ({stats['total_distance_km']:.3f} km, {stats['total_distance_nm']:.3f} nm)\n"
        stats_text += f"Total Time: {stats['total_time_min']:.1f} min ({stats['total_time_hours']:.2f} hr)\n\n"
        
        # Survey pattern explanation
        stats_text += "SURVEY PATTERN\n"
        stats_text += "-" * 15 + "\n"
        stats_text += "1. Run pitch line (2 passes)\n"
        stats_text += "2. Travel from pitch line start to roll line start\n"
        stats_text += "3. Run roll line (2 passes)\n"
        stats_text += "4. Travel from roll line start to heading line 1 start\n"
        stats_text += "5. Run heading line 1 (1 pass)\n"
        stats_text += "6. Travel from heading line 1 end to heading line 2 start\n"
        stats_text += "7. Run heading line 2 (1 pass)\n"
        stats_text += "\n"
        
        # Waypoints with labels
        stats_text += "WAYPOINTS\n"
        stats_text += "-" * 10 + "\n"
        if hasattr(self, 'pitch_line_points') and len(self.pitch_line_points) == 2:
            pitch_start = self.pitch_line_points[0]
            pitch_end = self.pitch_line_points[1]
            pitch_start_lat_ddm = self._decimal_degrees_to_ddm(pitch_start[0], is_latitude=True)
            pitch_start_lon_ddm = self._decimal_degrees_to_ddm(pitch_start[1], is_latitude=False)
            pitch_end_lat_ddm = self._decimal_degrees_to_ddm(pitch_end[0], is_latitude=True)
            pitch_end_lon_ddm = self._decimal_degrees_to_ddm(pitch_end[1], is_latitude=False)
            stats_text += f"Pitch Start: {pitch_start_lat_ddm}, {pitch_start_lon_ddm} ({pitch_start[0]:.6f}, {pitch_start[1]:.6f})\n"
            stats_text += f"Pitch End: {pitch_end_lat_ddm}, {pitch_end_lon_ddm} ({pitch_end[0]:.6f}, {pitch_end[1]:.6f})\n"
        if hasattr(self, 'roll_line_points') and len(self.roll_line_points) == 2:
            roll_start = self.roll_line_points[0]
            roll_end = self.roll_line_points[1]
            roll_start_lat_ddm = self._decimal_degrees_to_ddm(roll_start[0], is_latitude=True)
            roll_start_lon_ddm = self._decimal_degrees_to_ddm(roll_start[1], is_latitude=False)
            roll_end_lat_ddm = self._decimal_degrees_to_ddm(roll_end[0], is_latitude=True)
            roll_end_lon_ddm = self._decimal_degrees_to_ddm(roll_end[1], is_latitude=False)
            stats_text += f"Roll Start: {roll_start_lat_ddm}, {roll_start_lon_ddm} ({roll_start[0]:.6f}, {roll_start[1]:.6f})\n"
            stats_text += f"Roll End: {roll_end_lat_ddm}, {roll_end_lon_ddm} ({roll_end[0]:.6f}, {roll_end[1]:.6f})\n"
        if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 1:
            heading1_start = self.heading_lines[0][0]
            heading1_end = self.heading_lines[0][1]
            heading1_start_lat_ddm = self._decimal_degrees_to_ddm(heading1_start[0], is_latitude=True)
            heading1_start_lon_ddm = self._decimal_degrees_to_ddm(heading1_start[1], is_latitude=False)
            heading1_end_lat_ddm = self._decimal_degrees_to_ddm(heading1_end[0], is_latitude=True)
            heading1_end_lon_ddm = self._decimal_degrees_to_ddm(heading1_end[1], is_latitude=False)
            stats_text += f"Heading 1 Start: {heading1_start_lat_ddm}, {heading1_start_lon_ddm} ({heading1_start[0]:.6f}, {heading1_start[1]:.6f})\n"
            stats_text += f"Heading 1 End: {heading1_end_lat_ddm}, {heading1_end_lon_ddm} ({heading1_end[0]:.6f}, {heading1_end[1]:.6f})\n"
        if hasattr(self, 'heading_lines') and len(self.heading_lines) >= 2:
            heading2_start = self.heading_lines[1][0]
            heading2_end = self.heading_lines[1][1]
            heading2_start_lat_ddm = self._decimal_degrees_to_ddm(heading2_start[0], is_latitude=True)
            heading2_start_lon_ddm = self._decimal_degrees_to_ddm(heading2_start[1], is_latitude=False)
            heading2_end_lat_ddm = self._decimal_degrees_to_ddm(heading2_end[0], is_latitude=True)
            heading2_end_lon_ddm = self._decimal_degrees_to_ddm(heading2_end[1], is_latitude=False)
            stats_text += f"Heading 2 Start: {heading2_start_lat_ddm}, {heading2_start_lon_ddm} ({heading2_start[0]:.6f}, {heading2_start[1]:.6f})\n"
            stats_text += f"Heading 2 End: {heading2_end_lat_ddm}, {heading2_end_lon_ddm} ({heading2_end[0]:.6f}, {heading2_end[1]:.6f})\n"
        
        # Create custom dialog window with copy functionality
        self._show_statistics_dialog("Calibration Planning Statistics", stats_text)

    def _calculate_line_planning_statistics(self):
        """Calculate statistics for the drawn line in line planning mode."""
        if not self.line_planning_points or len(self.line_planning_points) < 2:
            return None
        
        try:
            if pyproj is not None:
                geod = pyproj.Geod(ellps="WGS84")
                
                # Calculate total distance
                total_distance = 0.0
                segment_distances = []
                
                for i in range(len(self.line_planning_points) - 1):
                    lat1, lon1 = self.line_planning_points[i]
                    lat2, lon2 = self.line_planning_points[i + 1]
                    fwd_az, back_az, dist = geod.inv(lon1, lat1, lon2, lat2)
                    segment_distances.append(dist)
                    total_distance += dist
                
                # Calculate heading for the first segment
                if len(self.line_planning_points) >= 2:
                    lat1, lon1 = self.line_planning_points[0]
                    lat2, lon2 = self.line_planning_points[1]
                    fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)
                    heading = fwd_az % 360
                    reciprocal_heading = back_az % 360
                else:
                    heading = 0
                    reciprocal_heading = 0
                
                # Calculate time based on survey speed
                try:
                    speed_knots = float(self.line_survey_speed_entry.get()) if self.line_survey_speed_entry.get() else 8.0
                except:
                    speed_knots = 8.0
                
                speed_ms = speed_knots * 0.514444  # Convert knots to m/s
                total_time_seconds = total_distance / speed_ms
                total_time_minutes = total_time_seconds / 60
                total_time_hours = total_time_minutes / 60
                
                # Get depth information if GeoTIFF is available
                depth_info = ""
                if self.geotiff_data_array is not None:
                    depths = []
                    for lat, lon in self.line_planning_points:
                        depth = self._get_depth_at_point(lat, lon)
                        depths.append(depth)
                    
                    if depths:
                        min_depth = min(depths)
                        max_depth = max(depths)
                        avg_depth = sum(depths) / len(depths)
                        depth_info = f"Depth Range: {min_depth:.1f} to {max_depth:.1f} m (avg: {avg_depth:.1f} m)\n"
                
                return {
                    'num_points': len(self.line_planning_points),
                    'total_distance_m': total_distance,
                    'total_distance_km': total_distance / 1000,
                    'total_distance_nm': total_distance / 1852,
                    'heading': heading,
                    'reciprocal_heading': reciprocal_heading,
                    'speed_knots': speed_knots,
                    'total_time_minutes': total_time_minutes,
                    'total_time_hours': total_time_hours,
                    'depth_info': depth_info,
                    'segment_distances': segment_distances
                }
            else:
                return None
        except Exception as e:
            print(f"Error calculating line planning statistics: {e}")
            return None

    def _show_line_information(self):
        """Display line planning statistics in a dialog."""
        stats = self._calculate_line_planning_statistics()
        if not stats:
            self._show_message("warning","No Line Data", "No line has been drawn. Please draw a line first.")
            return
        
        # Format the statistics for display
        stats_text = "LINE PLANNING STATISTICS\n"
        stats_text += "=" * 30 + "\n\n"
        
        stats_text += f"Number of Points: {stats['num_points']}\n"
        stats_text += f"Total Distance: {stats['total_distance_m']:.1f} m ({stats['total_distance_km']:.3f} km, {stats['total_distance_nm']:.3f} nm)\n"
        stats_text += f"Survey Speed: {stats['speed_knots']:.1f} knots\n"
        stats_text += f"Estimated Time: {stats['total_time_minutes']:.1f} min ({stats['total_time_hours']:.2f} hr)\n\n"
        
        if stats['depth_info']:
            stats_text += stats['depth_info']
        
        # Segment information
        if len(stats['segment_distances']) > 1:
            stats_text += f"\nSEGMENT DISTANCE AND HEADINGS\n"
            stats_text += "-" * 30 + "\n"
            for i in range(len(self.line_planning_points) - 1):
                lat1, lon1 = self.line_planning_points[i]
                lat2, lon2 = self.line_planning_points[i + 1]
                if pyproj is not None:
                    geod = pyproj.Geod(ellps="WGS84")
                    fwd_az, back_az, dist = geod.inv(lon1, lat1, lon2, lat2)
                    heading = fwd_az % 360
                    stats_text += f"Segment {i+1}: {dist:.1f} m ({dist/1000:.3f} km, {dist/1852:.3f} nm) - Heading: {heading:.1f}°\n"
                else:
                    dist = stats['segment_distances'][i]
                    stats_text += f"Segment {i+1}: {dist:.1f} m ({dist/1000:.3f} km, {dist/1852:.3f} nm) - Heading: Unable to calculate\n"
        
        # Waypoints with labels
        if hasattr(self, 'line_planning_points') and self.line_planning_points:
            stats_text += f"\nWAYPOINTS\n"
            stats_text += "-" * 10 + "\n"
            for i, point in enumerate(self.line_planning_points):
                lat, lon = point
                lat_ddm = self._decimal_degrees_to_ddm(lat, is_latitude=True)
                lon_ddm = self._decimal_degrees_to_ddm(lon, is_latitude=False)
                stats_text += f"WP{i+1}: {lat_ddm}, {lon_ddm} ({lat:.6f}, {lon:.6f})\n"
        
        # Create custom dialog window with copy functionality
        self._show_statistics_dialog("Line Planning Statistics", stats_text)


    def _decimal_degrees_to_ddm(self, decimal_deg, is_latitude=True):
        """Convert decimal degrees to degrees and decimal minutes format.
        
        Args:
            decimal_deg: Decimal degrees value
            is_latitude: True for latitude, False for longitude
        
        Returns:
            Formatted string like "DD°MM.mmm'N" or "DDD°MM.mmm'E"
        """
        try:
            # Determine sign and direction
            if is_latitude:
                direction = 'N' if decimal_deg >= 0 else 'S'
                degrees_width = 2
            else:
                direction = 'E' if decimal_deg >= 0 else 'W'
                degrees_width = 3
            
            # Work with absolute value
            abs_deg = abs(decimal_deg)
            
            # Extract degrees and minutes
            degrees = int(abs_deg)
            decimal_minutes = (abs_deg - degrees) * 60.0
            
            # Format with leading zeros for degrees, 3 decimal places for minutes
            return f"{degrees:0{degrees_width}d}°{decimal_minutes:06.3f}'{direction}"
        except:
            return f"{decimal_deg:.6f}"

    def _calculate_reference_survey_statistics(self):
        """Calculate comprehensive reference planning survey statistics using the existing calculation function."""
        try:
            # Use the existing _calculate_total_survey_time function which provides all the correct statistics
            stats = self._calculate_total_survey_time()
            if not stats:
                return None
            
            # Add number of main lines to the stats
            stats['num_main_lines'] = len(self.survey_lines_data)
            
            return stats
            
        except Exception as e:
            print(f"Error calculating reference survey statistics: {e}")
            return None

    def _show_reference_planning_info(self):
        """Display comprehensive reference planning survey statistics in a custom dialog with copy functionality."""
        stats = self._calculate_reference_survey_statistics()
        if not stats:
            self._show_message("warning","Statistics Error", "Unable to calculate reference planning statistics. Please ensure survey lines are generated.")
            return
        
        # Format the statistics for display
        stats_text = "COMPREHENSIVE REFERENCE PLANNING STATISTICS\n"
        stats_text += "=" * 50 + "\n\n"
        
        # Survey parameters
        try:
            speed_knots = float(self.survey_speed_entry.text()) if self.survey_speed_entry.text() else 8.0
            stats_text += f"Survey Speed: {speed_knots} knots\n\n"
        except:
            stats_text += "Survey Speed: 8.0 knots (default)\n\n"
        
        # Main survey lines
        stats_text += "MAIN SURVEY LINES\n"
        stats_text += "-" * 20 + "\n"
        stats_text += f"Number of Main Lines: {stats['num_main_lines']}\n"
        
        # Add heading information for the first main line
        if self.survey_lines_data and len(self.survey_lines_data) > 0:
            try:
                if pyproj is not None:
                    geod = pyproj.Geod(ellps="WGS84")
                    first_line = self.survey_lines_data[0]
                    lat1, lon1 = first_line[0]
                    lat2, lon2 = first_line[1]
                    fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)
                    
                    # Convert azimuth to heading (0-360 degrees)
                    heading = fwd_az % 360
                    reciprocal_heading = back_az % 360
                    
                    stats_text += f"Heading: {heading:.1f}°\n"
                    stats_text += f"Reciprocal Heading: {reciprocal_heading:.1f}°\n"
                else:
                    stats_text += f"Heading: pyproj not available\n"
            except Exception as e:
                stats_text += f"Heading: Unable to calculate\n"
        
        if stats['num_main_lines'] > 0:
            single_line_length_m = stats['main_lines_distance_m'] / stats['num_main_lines']
            single_line_length_km = single_line_length_m / 1000
            single_line_length_nm = single_line_length_m / 1852
            stats_text += f"Length of Main Line: {single_line_length_m:.1f} m ({single_line_length_km:.3f} km, {single_line_length_nm:.3f} nm)\n"
        stats_text += f"Total Distance: {stats['main_lines_distance_m']:.1f} m ({stats['main_lines_distance_km']:.3f} km, {stats['main_lines_distance_nm']:.3f} nm)\n"
        stats_text += f"Survey Time: {stats['main_lines_minutes']:.1f} min\n\n"
        
        # Travel between lines
        if stats['travel_between_lines_distance_m'] > 0:
            stats_text += "TRAVEL BETWEEN LINES\n"
            stats_text += "-" * 25 + "\n"
            stats_text += f"Total Travel Distance: {stats['travel_between_lines_distance_m']:.1f} m ({stats['travel_between_lines_distance_km']:.3f} km, {stats['travel_between_lines_distance_nm']:.3f} nm)\n"
            stats_text += f"Travel Time: {stats['travel_minutes']:.1f} min\n\n"
        
        # Crossline information
        if stats['crossline_single_pass_distance_m'] > 0:
            stats_text += "CROSSLINE\n"
            stats_text += "-" * 10 + "\n"
            stats_text += f"Crossline Passes: {stats['num_crossline_passes']}\n"
            
            # Add crossline heading information
            if self.cross_line_data and len(self.cross_line_data) >= 2:
                try:
                    if pyproj is not None:
                        geod = pyproj.Geod(ellps="WGS84")
                        lat1, lon1 = self.cross_line_data[0]
                        lat2, lon2 = self.cross_line_data[1]
                        fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)
                        
                        # Convert azimuth to heading (0-360 degrees)
                        crossline_heading = fwd_az % 360
                        crossline_reciprocal_heading = back_az % 360
                        
                        stats_text += f"Crossline Heading: {crossline_heading:.1f}°\n"
                        stats_text += f"Crossline Reciprocal Heading: {crossline_reciprocal_heading:.1f}°\n"
                    else:
                        stats_text += f"Crossline Heading: pyproj not available\n"
                except Exception as e:
                    stats_text += f"Crossline Heading: Unable to calculate\n"
            
            stats_text += f"Crossline (per pass): {stats['crossline_single_pass_distance_m']:.1f} m ({stats['crossline_single_pass_distance_km']:.3f} km, {stats['crossline_single_pass_distance_nm']:.3f} nm)\n"
            stats_text += f"Crossline Total Distance: {stats['crossline_total_distance_m']:.1f} m ({stats['crossline_total_distance_km']:.3f} km, {stats['crossline_total_distance_nm']:.3f} nm)\n"
            stats_text += f"Crossline Time: {stats['crossline_minutes']:.1f} min\n"
            if stats['travel_to_crossline_distance_m'] > 0:
                stats_text += f"Travel to Crossline: {stats['travel_to_crossline_distance_m']:.1f} m ({stats['travel_to_crossline_distance_km']:.3f} km, {stats['travel_to_crossline_distance_nm']:.3f} nm)\n"
                stats_text += f"Travel to Crossline Time: {stats['travel_to_crossline_minutes']:.1f} min\n"
            stats_text += "\n"
        
        # Total summary
        stats_text += "TOTAL SURVEY SUMMARY\n"
        stats_text += "-" * 25 + "\n"
        stats_text += f"Total Distance: {stats['total_distance_m']:.1f} m ({stats['total_distance_km']:.3f} km, {stats['total_distance_nm']:.3f} nm)\n"
        stats_text += f"Total Time: {stats['total_minutes']:.1f} min ({stats['total_hours']:.2f} hr)\n\n"
        
        # Survey pattern explanation
        stats_text += "SURVEY PATTERN\n"
        stats_text += "-" * 15 + "\n"
        stats_text += "1. Run main survey lines (zigzag pattern)\n"
        if stats['crossline_single_pass_distance_m'] > 0:
            stats_text += "2. Travel from last main line to crossline start\n"
            stats_text += f"3. Run crossline ({stats['num_crossline_passes']} passes)\n"
        stats_text += "\n"
        
        # Waypoints with labels
        stats_text += "WAYPOINTS\n"
        stats_text += "-" * 10 + "\n"
        if self.survey_lines_data:
            for i, line in enumerate(self.survey_lines_data):
                # Determine start and end based on zigzag pattern
                if i % 2 == 0:  # Even lines - normal order
                    start, end = line[0], line[1]
                else:  # Odd lines - flipped order
                    start, end = line[1], line[0]
                start_label = f'L{i+1}S'
                end_label = f'L{i+1}E'
                start_lat_ddm = self._decimal_degrees_to_ddm(start[0], is_latitude=True)
                start_lon_ddm = self._decimal_degrees_to_ddm(start[1], is_latitude=False)
                end_lat_ddm = self._decimal_degrees_to_ddm(end[0], is_latitude=True)
                end_lon_ddm = self._decimal_degrees_to_ddm(end[1], is_latitude=False)
                stats_text += f"{start_label}: {start_lat_ddm}, {start_lon_ddm} ({start[0]:.6f}, {start[1]:.6f})\n"
                stats_text += f"{end_label}: {end_lat_ddm}, {end_lon_ddm} ({end[0]:.6f}, {end[1]:.6f})\n"
        if self.cross_line_data:
            cls_lat_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[0][0], is_latitude=True)
            cls_lon_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[0][1], is_latitude=False)
            cle_lat_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[1][0], is_latitude=True)
            cle_lon_ddm = self._decimal_degrees_to_ddm(self.cross_line_data[1][1], is_latitude=False)
            stats_text += f"CLS: {cls_lat_ddm}, {cls_lon_ddm} ({self.cross_line_data[0][0]:.6f}, {self.cross_line_data[0][1]:.6f})\n"
            stats_text += f"CLE: {cle_lat_ddm}, {cle_lon_ddm} ({self.cross_line_data[1][0]:.6f}, {self.cross_line_data[1][1]:.6f})\n"
        
        # Create custom dialog window with copy functionality
        self._show_statistics_dialog("Reference Planning Statistics", stats_text)

    def _show_statistics_dialog(self, title, text):
        """Create a custom dialog window with copy functionality for statistics display."""
        dialog = QDialog(self)
        dialog.setWindowTitle(title)
        dialog.setMinimumSize(600, 500)
        dialog.setModal(True)
        
        # Create main layout
        main_layout = QVBoxLayout(dialog)
        main_layout.setContentsMargins(10, 10, 10, 10)
        
        # Create text widget with scrollbar
        text_widget = QTextEdit()
        text_widget.setReadOnly(True)
        text_widget.setPlainText(text)
        text_widget.setFont(QApplication.font())
        main_layout.addWidget(text_widget)
        
        # Create button frame
        button_frame = QWidget()
        button_layout = QHBoxLayout(button_frame)
        button_layout.setContentsMargins(0, 0, 0, 0)
        
        # Copy button
        copy_btn = QPushButton("Copy to Clipboard")
        def copy_to_clipboard():
            try:
                # Get selected text or all text
                cursor = text_widget.textCursor()
                if cursor.hasSelection():
                    selected_text = cursor.selectedText()
                else:
                    selected_text = text_widget.toPlainText()
                
                # Copy to clipboard
                clipboard = QApplication.clipboard()
                clipboard.setText(selected_text)
                
                # Show brief confirmation
                copy_btn.setText("Copied!")
                QTimer.singleShot(1000, lambda: copy_btn.setText("Copy to Clipboard"))
                
            except Exception as e:
                self._show_message("error", "Copy Error", f"Failed to copy to clipboard: {e}")
        
        copy_btn.clicked.connect(copy_to_clipboard)
        button_layout.addWidget(copy_btn)
        button_layout.addStretch()
        
        # Close button
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.accept)
        button_layout.addWidget(close_btn)
        
        main_layout.addWidget(button_frame)
        
        dialog.exec()

if __name__ == "__main__":
    app = QApplication([])
    window = SurveyPlanApp()
    window.show()
    app.exec()
    
    
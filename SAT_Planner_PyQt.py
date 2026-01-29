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

# Constants, geospatial libs, and geo utils from sat_planner package
from sat_planner import (
    __version__,
    CONFIG_FILENAME,
    GEOSPATIAL_LIBS_AVAILABLE,
    decimal_degrees_to_ddm,
)
from sat_planner.constants import (
    rasterio,
    transform,
    RasterioIOError,
    Window,
    window_transform,
    window_bounds,
    rowcol,
    reproject,
    Resampling,
    pyproj,
    CRSError,
    LineString,
    fiona,
)
from sat_planner.mixins.geotiff_mixin import GeoTIFFMixin
from sat_planner.mixins.plotting_mixin import PlottingMixin
from sat_planner.mixins.reference_mixin import ReferenceMixin
from sat_planner.mixins.calibration_mixin import CalibrationMixin
from sat_planner.mixins.line_planning_mixin import LinePlanningMixin

"""
UNH/CCOM-JHC Shipboard Acceptance Testing (SAT) and Quality Assurance Testing (QAT) Planner
A Python application to help plan and visualize shipboard acceptance testing and quality assurance testing.

![Example Plot](media/SAT_Planner.jpg)

Program by Paul Johnson, pjohnson@ccom.unh.edu
Date: 2025-09-12

Center for Coastal and Ocean Mapping/Joint Hydrographic Center, University of New Hampshire

This program was developed at the University of New Hampshire, Center for Coastal and Ocean Mapping - Joint Hydrographic Center (UNH/CCOM-JHC) under the grant NA20NOS4000196 from the National Oceanic and Atmospheric Administration (NOAA).

This software is released for general use under the BSD 3-Clause License.

"""

class SurveyPlanApp(GeoTIFFMixin, PlottingMixin, ReferenceMixin, CalibrationMixin, LinePlanningMixin, QMainWindow):
    CONFIG_FILENAME = CONFIG_FILENAME  # from sat_planner

    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"UNH/CCOM-JHC - SAT Planner - v{__version__} - pjohnson@ccom.unh.edu")
        
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

        # Pitch/center/hover motion (line planning in LinePlanningMixin)
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

        # Patch pitch/center/hover motion methods into self
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

    def _get_profile_data_from_geotiff(self, lats, lons):
        """Get elevation and slope data from GeoTIFF for profile calculations.
        Always uses the original full dataset to ensure consistent results regardless of zoom level.
        Returns (elevations, slopes, extent) or (None, None, None) if data unavailable."""
        if self.geotiff_dataset_original is None:
            return None, None, None
        
        # Use original full extent for profile calculations
        if hasattr(self, 'geotiff_original_extent') and self.geotiff_original_extent is not None:
            extent = self.geotiff_original_extent
        elif self.geotiff_extent is not None:
            extent = self.geotiff_extent
        else:
            return None, None, None
        
        # Read data directly from original dataset for full extent
        try:
            from rasterio.warp import reproject, Resampling, calculate_default_transform
            from rasterio.transform import from_bounds
            import pyproj
            
            # Get bounds in original CRS
            if self.geotiff_dataset_original.crs != "EPSG:4326":
                transformer = pyproj.Transformer.from_crs("EPSG:4326", self.geotiff_dataset_original.crs, always_xy=True)
                # Transform extent to original CRS
                left, bottom = transformer.transform(extent[0], extent[2])
                right, top = transformer.transform(extent[1], extent[3])
            else:
                left, right, bottom, top = extent
            
            # Read full dataset
            data = self.geotiff_dataset_original.read(1).astype(float)
            
            # If not WGS84, we need to reproject the data
            if self.geotiff_dataset_original.crs != "EPSG:4326":
                # Reproject to WGS84 for profile calculations
                src_crs = self.geotiff_dataset_original.crs
                dst_crs = "EPSG:4326"
                src_bounds = self.geotiff_dataset_original.bounds
                src_width = self.geotiff_dataset_original.width
                src_height = self.geotiff_dataset_original.height
                
                # Calculate transform for reprojection
                dst_transform, dst_width, dst_height = calculate_default_transform(
                    src_crs, dst_crs, src_width, src_height,
                    left=src_bounds.left, bottom=src_bounds.bottom,
                    right=src_bounds.right, top=src_bounds.top
                )
                
                # Reproject
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
                
                # Get reprojected bounds
                from rasterio.transform import array_bounds
                dst_bounds = array_bounds(dst_height, dst_width, dst_transform)
                extent = [dst_bounds.left, dst_bounds.right, dst_bounds.bottom, dst_bounds.top]
            
            # Filter Z-values
            data[data < -11000] = np.nan
            data[data > 0] = np.nan
            
            # Sample elevations and slopes along the line
            left, right, bottom, top = extent
            nrows, ncols = data.shape
            rows = ((top - lats) / (top - bottom) * (nrows - 1)).clip(0, nrows - 1)
            cols = ((lons - left) / (right - left) * (ncols - 1)).clip(0, ncols - 1)
            
            elevations = []
            slopes = []
            for r, c in zip(rows, cols):
                ir, ic = int(round(r)), int(round(c))
                elevations.append(data[ir, ic])
                
                # Slope calculation
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
                    if window.shape == (3,3) and not np.all(np.isnan(window)):
                        dz_dy, dz_dx = np.gradient(window, dy_m, dx_m)
                        slope_rad = np.arctan(np.sqrt(dz_dx[1,1]**2 + dz_dy[1,1]**2))
                        slope = np.degrees(slope_rad)
                slopes.append(slope if slope is not None else np.nan)
            
            return np.array(elevations), np.array(slopes), extent
            
        except Exception as e:
            print(f"Error reading profile data from GeoTIFF: {e}")
            return None, None, None

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
        
        # Use helper function to get data from original full dataset
        elevations, slopes, profile_extent = self._get_profile_data_from_geotiff(lats, lons)
        if elevations is None:
            # Fallback to current data if helper fails
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
        
        # Always ensure canvas is updated, even if no data
        if not (
            hasattr(self, 'pitch_line_points') and
            len(self.pitch_line_points) == 2
        ):
            # No data - just clear and update
            self.profile_fig.tight_layout(pad=1.0)
            if hasattr(self, 'profile_canvas'):
                self.profile_canvas.draw_idle()
            return
        
        if (
            hasattr(self, 'pitch_line_points') and
            len(self.pitch_line_points) == 2
        ):
            (lat1, lon1), (lat2, lon2) = self.pitch_line_points
            lats = np.linspace(lat1, lat2, 100)
            lons = np.linspace(lon1, lon2, 100)
            
            # Use helper function to get data from original full dataset
            elevations, slopes, profile_extent = self._get_profile_data_from_geotiff(lats, lons)
            if elevations is None:
                # Fallback to current data if helper fails
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
                    start_lat_ddm = decimal_degrees_to_ddm(start[0], is_latitude=True)
                    start_lon_ddm = decimal_degrees_to_ddm(start[1], is_latitude=False)
                    end_lat_ddm = decimal_degrees_to_ddm(end[0], is_latitude=True)
                    end_lon_ddm = decimal_degrees_to_ddm(end[1], is_latitude=False)
                    ddm_writer.writerow([i + 1, start_label, start_lat_ddm, start_lon_ddm])
                    ddm_writer.writerow([i + 1, end_label, end_lat_ddm, end_lon_ddm])
                # Crossline
                if self.cross_line_data:
                    cls_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][0], is_latitude=True)
                    cls_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][1], is_latitude=False)
                    cle_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][0], is_latitude=True)
                    cle_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][1], is_latitude=False)
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
                    start_lat_ddm = decimal_degrees_to_ddm(start[0], is_latitude=True)
                    start_lon_ddm = decimal_degrees_to_ddm(start[1], is_latitude=False)
                    end_lat_ddm = decimal_degrees_to_ddm(end[0], is_latitude=True)
                    end_lon_ddm = decimal_degrees_to_ddm(end[1], is_latitude=False)
                    ddm_txt_file.write(f"{start_label}, {start_lat_ddm}, {start_lon_ddm}\n")
                    ddm_txt_file.write(f"{end_label}, {end_lat_ddm}, {end_lon_ddm}\n")
                # Crossline
                if self.cross_line_data:
                    cls_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][0], is_latitude=True)
                    cls_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][1], is_latitude=False)
                    cle_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][0], is_latitude=True)
                    cle_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][1], is_latitude=False)
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
                        start_lat_ddm = decimal_degrees_to_ddm(start[0], is_latitude=True)
                        start_lon_ddm = decimal_degrees_to_ddm(start[1], is_latitude=False)
                        end_lat_ddm = decimal_degrees_to_ddm(end[0], is_latitude=True)
                        end_lon_ddm = decimal_degrees_to_ddm(end[1], is_latitude=False)
                        f.write(f"{start_label}: {start_lat_ddm}, {start_lon_ddm} ({start[0]:.6f}, {start[1]:.6f})\n")
                        f.write(f"{end_label}: {end_lat_ddm}, {end_lon_ddm} ({end[0]:.6f}, {end[1]:.6f})\n")
                if self.cross_line_data:
                    cls_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][0], is_latitude=True)
                    cls_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][1], is_latitude=False)
                    cle_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][0], is_latitude=True)
                    cle_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][1], is_latitude=False)
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
            
            # Update Pitch Line Info labels
            if hasattr(self, 'pitch_shallowest_depth_label'):
                self.pitch_shallowest_depth_label.setText(f"{abs(max_val):.2f}")
            if hasattr(self, 'pitch_max_depth_label'):
                self.pitch_max_depth_label.setText(f"{abs(min_val):.2f}")
            if hasattr(self, 'pitch_mean_depth_label'):
                self.pitch_mean_depth_label.setText(f"{abs(mean_val):.2f}")
            if hasattr(self, 'pitch_median_depth_label'):
                self.pitch_median_depth_label.setText(f"{abs(median_val):.2f}")
            
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
            
            # Clear Pitch Line Info labels
            if hasattr(self, 'pitch_shallowest_depth_label'):
                self.pitch_shallowest_depth_label.setText("-")
            if hasattr(self, 'pitch_max_depth_label'):
                self.pitch_max_depth_label.setText("-")
            if hasattr(self, 'pitch_mean_depth_label'):
                self.pitch_mean_depth_label.setText("-")
            if hasattr(self, 'pitch_median_depth_label'):
                self.pitch_median_depth_label.setText("-")
            
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
        
        # Handle line planning motion (delegated to LinePlanningMixin)
        if hasattr(self, 'line_planning_mode') and self.line_planning_mode and len(self.line_planning_points) >= 1:
            self._on_line_planning_motion(event)
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

    def _draw_current_profile(self):
        # Helper to redraw the correct profile based on the current tab
        if not hasattr(self, 'param_notebook'):
            return
        
        if not hasattr(self, 'profile_ax') or not hasattr(self, 'profile_canvas'):
            return
        
        current_tab = self.param_notebook.currentIndex()
        if current_tab == 0:  # Calibration tab - show pitch line profile
            self._draw_pitch_line_profile()
        elif current_tab == 1:  # Reference tab - show crossline profile
            self._draw_crossline_profile()
        elif current_tab == 2:  # Line tab - show line planning profile
            self._draw_line_planning_profile()
        else:
            # Unknown tab - clear profile
            self.profile_ax.clear()
            self.profile_ax.set_title("Elevation Profile", fontsize=8)
            self.profile_ax.set_xlabel("Distance (m)", fontsize=8)
            self.profile_ax.set_ylabel("Elevation (m)", fontsize=8)
            self.profile_fig.tight_layout(pad=1.0)
            self.profile_canvas.draw()

    def _on_tab_changed(self, event=None):
        """Handle tab change event - update profile plot for active tab."""
        if not hasattr(self, 'param_notebook'):
            return
        
        try:
            current_tab = self.param_notebook.currentIndex()
            
            # Mark as initialized and update previous tab index
            self.tab_switch_initialized = True
            self.previous_tab_index = current_tab
        except Exception as e:
            print(f"Error in tab change handler: {e}")
        
        # Always update profile plot to show the profile for the active tab
        try:
            self._draw_current_profile()
            # Force immediate canvas update
            if hasattr(self, 'profile_canvas'):
                self.profile_canvas.draw()
        except Exception as e:
            print(f"Error updating profile plot: {e}")

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
        geotiff_layout.setSpacing(0)
        geotiff_layout.setContentsMargins(9, 9, 9, 9)
        
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
        geotiff_layout.addSpacing(3)
        
        # Display mode dropdown - label and combo on same line
        display_frame = QWidget()
        display_layout = QHBoxLayout(display_frame)
        display_layout.setContentsMargins(0, 0, 0, 0)
        display_layout.addWidget(QLabel("Display:"))
        self.elevation_slope_combo = QComboBox()
        self.elevation_slope_combo.addItems(["Shaded Relief", "Shaded Slope", "Hillshade", "Slope"])
        self.elevation_slope_combo.setCurrentText("Shaded Relief")  # Set default
        self.elevation_slope_combo.currentTextChanged.connect(self._on_geotiff_display_mode_changed)
        display_layout.addWidget(self.elevation_slope_combo)
        geotiff_layout.addWidget(display_frame)
        geotiff_layout.addSpacing(3)
        
        # Dynamic Resolution button
        self.dynamic_resolution_btn = QPushButton("Dynamic Resolution: ON")
        self.dynamic_resolution_btn.clicked.connect(self._toggle_dynamic_resolution)
        geotiff_layout.addWidget(self.dynamic_resolution_btn)
        geotiff_layout.addSpacing(3)
        
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
        self.activity_log_text.setStyleSheet("background-color: #e0e0e0; color: black;")
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
        ref_plot_control_groupbox = QGroupBox("Reference Plot Control")
        ref_plot_control_groupbox.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        ref_plot_control_layout = QVBoxLayout(ref_plot_control_groupbox)
        ref_plot_control_layout.setSpacing(0)
        ref_plot_control_layout.setContentsMargins(9, 9, 9, 9)
        
        self.zoom_to_plan_btn = QPushButton("Zoom to Reference Plan")
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
        ref_test_plan_info_groupbox = QGroupBox("Reference Info")
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
        ref_import_export_groupbox = QGroupBox("Reference Import/Export")
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
        
        # --- Pitch Line Info GroupBox ---
        pitch_line_info_groupbox = QGroupBox("Pitch Line Info")
        pitch_line_info_groupbox.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        pitch_line_info_layout = QGridLayout(pitch_line_info_groupbox)
        pitch_line_info_layout.setSpacing(3)
        pitch_line_info_layout.setContentsMargins(9, 9, 9, 9)
        pitch_line_info_layout.setColumnStretch(0, 1)
        pitch_line_info_layout.setColumnStretch(1, 2)
        
        pitch_info_row = 0
        pitch_line_info_layout.addWidget(QLabel("Shallowest Depth (m):"), pitch_info_row, 0)
        self.pitch_shallowest_depth_label = QLabel("-")
        pitch_line_info_layout.addWidget(self.pitch_shallowest_depth_label, pitch_info_row, 1)
        pitch_info_row += 1
        
        pitch_line_info_layout.addWidget(QLabel("Maximum Depth (m):"), pitch_info_row, 0)
        self.pitch_max_depth_label = QLabel("-")
        pitch_line_info_layout.addWidget(self.pitch_max_depth_label, pitch_info_row, 1)
        pitch_info_row += 1
        
        pitch_line_info_layout.addWidget(QLabel("Mean Depth (m):"), pitch_info_row, 0)
        self.pitch_mean_depth_label = QLabel("-")
        pitch_line_info_layout.addWidget(self.pitch_mean_depth_label, pitch_info_row, 1)
        pitch_info_row += 1
        
        pitch_line_info_layout.addWidget(QLabel("Median Depth (m):"), pitch_info_row, 0)
        self.pitch_median_depth_label = QLabel("-")
        pitch_line_info_layout.addWidget(self.pitch_median_depth_label, pitch_info_row, 1)
        
        cal_layout.addWidget(pitch_line_info_groupbox, cal_row, 0, 1, 2)
        cal_layout.setRowStretch(cal_row, 0)
        cal_row += 1
        
        # --- Plot Control GroupBox ---
        cal_plot_control_groupbox = QGroupBox("Calibration Plot Control")
        cal_plot_control_groupbox.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        cal_plot_control_layout = QVBoxLayout(cal_plot_control_groupbox)
        cal_plot_control_layout.setSpacing(0)
        cal_plot_control_layout.setContentsMargins(9, 9, 9, 9)
        
        self.zoom_to_all_lines_btn = QPushButton("Zoom to Calibration Plan")
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
        cal_test_plan_info_groupbox = QGroupBox("Calibration Info")
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
        cal_import_export_groupbox = QGroupBox("Calibration Import/Export")
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
        line_plot_control_groupbox = QGroupBox("Line Plot Control")
        line_plot_control_groupbox.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        line_plot_control_layout = QVBoxLayout(line_plot_control_groupbox)
        line_plot_control_layout.setSpacing(0)
        line_plot_control_layout.setContentsMargins(9, 9, 9, 9)
        
        self.zoom_to_line_btn = QPushButton("Zoom to Line Plan")
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
        line_test_plan_info_groupbox = QGroupBox("Line Info")
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
        line_import_export_groupbox = QGroupBox("Line Import/Export")
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
    
    
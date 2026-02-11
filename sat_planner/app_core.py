"""
SAT Planner application core: SurveyPlanApp main window and core setup.
UNH/CCOM-JHC Shipboard Acceptance Testing (SAT) and Quality Assurance Testing (QAT) Planner.

Copyright (c) 2025, UNH/CCOM-JHC. BSD 3-Clause License (see LICENSE in repo root).
"""
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                             QGridLayout, QLabel, QLineEdit, QPushButton, QCheckBox, QTextEdit,
                             QScrollArea, QTabWidget, QFileDialog, QMessageBox, QDialog,
                             QDialogButtonBox, QSlider, QComboBox, QFrame, QSizePolicy, QProgressBar, QGroupBox)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt6.QtGui import QTextCursor, QColor, QTextCharFormat, QPixmap
import matplotlib
matplotlib.use('QtAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.colors import LightSource
from matplotlib.figure import Figure
import numpy as np
import os
import sys
import csv
import traceback
import json
import platform
import datetime
import re
import threading
import time

from . import __version__, CONFIG_FILENAME, GEOSPATIAL_LIBS_AVAILABLE
from .utils_ui import show_message as _show_message_fn, ask_yes_no as _ask_yes_no_fn, ask_ok_cancel as _ask_ok_cancel_fn
from .utils_geo import decimal_degrees_to_ddm
from .mixins.basemap_mixin import BasemapMixin
from .mixins.geotiff_mixin import GeoTIFFMixin
from .mixins.plotting_mixin import PlottingMixin
from .mixins.reference_mixin import ReferenceMixin
from .mixins.calibration_mixin import CalibrationMixin
from .mixins.line_planning_mixin import LinePlanningMixin
from .mixins.profiles_mixin import ProfilesMixin
from .mixins.map_interaction_mixin import MapInteractionMixin
from .mixins.export_import_mixin import ExportImportMixin
from .mixins.config_mixin import ConfigMixin


class SurveyPlanApp(BasemapMixin, GeoTIFFMixin, PlottingMixin, ReferenceMixin, CalibrationMixin, LinePlanningMixin, ProfilesMixin, MapInteractionMixin, ExportImportMixin, ConfigMixin, QMainWindow):
    CONFIG_FILENAME = CONFIG_FILENAME  # from package

    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"UNH/CCOM-JHC - SAT Planner - v{__version__} - pjohnson@ccom.unh.edu")

        # Set minimum window size
        self.setMinimumSize(1600, 1110)

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

        # Slope overlay display variables
        self.show_slope_overlay_var = False  # Default to off
        self.slope_overlay_min_var = 10.0  # Default min slope in degrees
        self.slope_overlay_max_var = 20.0  # Default max slope in degrees
        self.slope_overlay_image_plot = None  # Store slope overlay plot object for removal/update

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
            lat_str = decimal_degrees_to_ddm(cur_lat, is_latitude=True)
            lon_str = decimal_degrees_to_ddm(cur_lon, is_latitude=False)
            info_str = f"Lat: {lat_str}\nLon: {lon_str}\nLength: {dist:.1f} m\nAzimuth: {az12:.1f}°\nTime: {time_minutes:.1f} min"
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
            lat_str = decimal_degrees_to_ddm(cur_lat, is_latitude=True)
            lon_str = decimal_degrees_to_ddm(cur_lon, is_latitude=False)
            info_str = f"Lat: {lat_str}\nLon: {lon_str}\nDepth: {depth_str}\nSlope: {slope_str}"
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

        # Imagery Basemap checkbox
        if not hasattr(self, 'imagery_basemap_checkbox'):
            self.show_imagery_basemap_var = False
            self.imagery_basemap_checkbox = QCheckBox("Imagery Basemap")
            self.imagery_basemap_checkbox.setChecked(self.show_imagery_basemap_var)
            self.imagery_basemap_checkbox.stateChanged.connect(self._toggle_imagery_basemap)
            if not hasattr(self, 'basemap_image_plot'):
                self.basemap_image_plot = None

        # NOAA ENC Charts checkbox and opacity slider
        if not hasattr(self, 'noaa_charts_checkbox'):
            self.show_noaa_charts_var = False
            self.noaa_charts_checkbox = QCheckBox("NOAA ENC Charts")
            self.noaa_charts_checkbox.setChecked(self.show_noaa_charts_var)
            self.noaa_charts_checkbox.stateChanged.connect(self._toggle_noaa_charts)
            self.noaa_charts_opacity = 50
            self.noaa_charts_opacity_slider = QSlider(Qt.Orientation.Horizontal)
            self.noaa_charts_opacity_slider.setMinimum(0)
            self.noaa_charts_opacity_slider.setMaximum(100)
            self.noaa_charts_opacity_slider.setValue(self.noaa_charts_opacity)
            self.noaa_charts_opacity_slider.setMaximumWidth(100)
            self.noaa_charts_opacity_slider.setEnabled(False)
            self.noaa_charts_opacity_slider.valueChanged.connect(self._update_noaa_charts_opacity)
            self.noaa_charts_opacity_label = QLabel(f"Opacity: {self.noaa_charts_opacity}%")
            self.noaa_charts_opacity_label.setEnabled(False)
            if not hasattr(self, 'noaa_charts_image_plot'):
                self.noaa_charts_image_plot = None

        # About button (only create once)
        if not hasattr(self, 'about_btn'):
            self.about_btn = QPushButton("About This Program")
            self.about_btn.clicked.connect(self._show_about_dialog)

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
        # Connect axis limit change events for dynamic resolution (toolbar zoom/pan support)
        self.cid_xlim_changed = self.ax.callbacks.connect('xlim_changed', self._on_axis_limits_changed)
        self.cid_ylim_changed = self.ax.callbacks.connect('ylim_changed', self._on_axis_limits_changed)

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
            if hasattr(self, 'zoom_to_geotiff_btn'):
                self.zoom_to_geotiff_btn.setEnabled(False)

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
        self.resize(1600, 1110)
        print("Initialization complete, window should be visible")

    # Helper methods for dialog conversions (delegate to utils_ui)
    def _show_message(self, msg_type, title, message):
        return _show_message_fn(self, msg_type, title, message)

    def _ask_yes_no(self, title, message):
        return _ask_yes_no_fn(self, title, message)

    def _ask_ok_cancel(self, title, message):
        return _ask_ok_cancel_fn(self, title, message)

    def _show_about_dialog(self):
        """Show About dialog with program information."""
        dialog = QDialog(self)
        dialog.setWindowTitle("About This Program")
        dialog.setMinimumWidth(500)
        dialog.setMinimumHeight(400)

        layout = QVBoxLayout(dialog)
        layout.setSpacing(15)
        layout.setContentsMargins(20, 20, 20, 20)

        program_name = QLabel(f"UNH/CCOM-JHC SAT PLanner v{__version__}")
        program_name.setStyleSheet("font-size: 16pt; font-weight: bold;")
        program_name.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(program_name)

        compile_date = "Unknown"
        try:
            if getattr(sys, 'frozen', False):
                exe_path = sys.executable
                if os.path.exists(exe_path):
                    mod_time = os.path.getmtime(exe_path)
                    compile_date = datetime.datetime.fromtimestamp(mod_time).strftime("%B %d, %Y")
            else:
                script_path = __file__
                if os.path.exists(script_path):
                    mod_time = os.path.getmtime(script_path)
                    compile_date = datetime.datetime.fromtimestamp(mod_time).strftime("%B %d, %Y")
        except Exception:
            compile_date = "Unknown"

        date_label = QLabel(f"Compiled: {compile_date}")
        date_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        date_label.setStyleSheet("font-size: 10pt; color: gray;")
        layout.addWidget(date_label)

        # Media folder is at project root (parent of sat_planner package)
        _pkg_dir = os.path.dirname(os.path.abspath(__file__))
        _project_root = os.path.dirname(_pkg_dir)
        logo_path = os.path.join(_project_root, "media", "CCOM.png")
        if os.path.exists(logo_path):
            logo_label = QLabel()
            pixmap = QPixmap(logo_path)
            if pixmap.width() > 300:
                pixmap = pixmap.scaledToWidth(300, Qt.TransformationMode.SmoothTransformation)
            logo_label.setPixmap(pixmap)
            logo_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            layout.addWidget(logo_label)

        author_label = QLabel("Paul Johnson")
        author_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        author_label.setStyleSheet("font-size: 12pt; font-weight: bold; margin-top: 10px;")
        layout.addWidget(author_label)

        email_label = QLabel("pjohnson@ccom.unh.edu")
        email_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        email_label.setStyleSheet("font-size: 10pt; margin-top: 3px;")
        layout.addWidget(email_label)

        institution_label = QLabel("Center for Coastal and Ocean Mapping/Joint Hydrographic Center, University of New Hampshire")
        institution_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        institution_label.setWordWrap(True)
        institution_label.setStyleSheet("font-size: 10pt; margin-top: 5px;")
        layout.addWidget(institution_label)

        grant_text = ("This program was developed at the University of New Hampshire, "
                      "Center for Coastal and Ocean Mapping - Joint Hydrographic Center "
                      "(UNH/CCOM-JHC) under the grant NA25NOSX400C0001-T1-01 from the National "
                      "Oceanic and Atmospheric Administration (NOAA).")
        grant_label = QLabel(grant_text)
        grant_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        grant_label.setWordWrap(True)
        grant_label.setStyleSheet("font-size: 9pt; margin-top: 15px;")
        layout.addWidget(grant_label)

        license_label = QLabel("This software is released for general use under the BSD 3-Clause License.")
        license_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        license_label.setWordWrap(True)
        license_label.setStyleSheet("font-size: 9pt; margin-top: 10px;")
        layout.addWidget(license_label)

        layout.addStretch()

        ok_button = QPushButton("OK")
        ok_button.setMaximumWidth(100)
        ok_button.clicked.connect(dialog.accept)
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        button_layout.addWidget(ok_button)
        button_layout.addStretch()
        layout.addLayout(button_layout)

        dialog.exec()

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
                # Horizontal layout for checkbox and About button
                checkbox_button_layout = QHBoxLayout()
                checkbox_button_layout.setContentsMargins(0, 0, 0, 0)
                checkbox_button_layout.addWidget(self.slope_profile_checkbox)
                if hasattr(self, 'imagery_basemap_checkbox'):
                    checkbox_button_layout.addWidget(self.imagery_basemap_checkbox)
                if hasattr(self, 'noaa_charts_checkbox'):
                    checkbox_button_layout.addWidget(self.noaa_charts_checkbox)
                    if hasattr(self, 'noaa_charts_opacity_label'):
                        checkbox_button_layout.addWidget(self.noaa_charts_opacity_label)
                    if hasattr(self, 'noaa_charts_opacity_slider'):
                        self.noaa_charts_opacity_slider.setMaximumWidth(100)
                        checkbox_button_layout.addWidget(self.noaa_charts_opacity_slider)
                checkbox_button_layout.addStretch()
                if hasattr(self, 'about_btn'):
                    checkbox_button_layout.addWidget(self.about_btn)
                profile_layout.addLayout(checkbox_button_layout)
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

        # Zoom to GeoTIFF button - zooms map to original bounds when GeoTIFF was first loaded
        self.zoom_to_geotiff_btn = QPushButton("Zoom to GeoTIFF")
        self.zoom_to_geotiff_btn.clicked.connect(self._zoom_to_geotiff)
        self.zoom_to_geotiff_btn.setEnabled(False)  # Enabled when a GeoTIFF is loaded
        geotiff_layout.addWidget(self.zoom_to_geotiff_btn)
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
        geotiff_layout.addSpacing(3)

        # Show Slope overlay checkbox and Min/Max on same row
        slope_overlay_frame = QWidget()
        slope_overlay_layout = QHBoxLayout(slope_overlay_frame)
        slope_overlay_layout.setContentsMargins(0, 0, 0, 0)
        self.show_slope_overlay_checkbox = QCheckBox("Slopes")
        self.show_slope_overlay_checkbox.setChecked(self.show_slope_overlay_var)
        self.show_slope_overlay_checkbox.stateChanged.connect(self._on_slope_overlay_checkbox_changed)
        slope_overlay_layout.addWidget(self.show_slope_overlay_checkbox)
        slope_overlay_layout.addStretch()  # Add stretch to push min/max to the right
        slope_overlay_layout.addWidget(QLabel("Min:"))
        self.slope_overlay_min_entry = QLineEdit("10")
        self.slope_overlay_min_entry.textChanged.connect(self._on_slope_overlay_min_changed)
        self.slope_overlay_min_entry.setMaximumWidth(60)  # Limit width of entry field
        slope_overlay_layout.addWidget(self.slope_overlay_min_entry)
        slope_overlay_layout.addWidget(QLabel("Max:"))
        self.slope_overlay_max_entry = QLineEdit("20")
        self.slope_overlay_max_entry.textChanged.connect(self._on_slope_overlay_max_changed)
        self.slope_overlay_max_entry.setMaximumWidth(60)  # Limit width of entry field
        slope_overlay_layout.addWidget(self.slope_overlay_max_entry)
        geotiff_layout.addWidget(slope_overlay_frame)

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

        self.clear_plot_btn = QPushButton("Remove Reference Surface Plan")
        self.clear_plot_btn.clicked.connect(self._clear_reference_plan)  # Remove only reference plan lines
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

        ref_test_plan_info_layout.addWidget(QLabel("Turn Time (min):"), ref_test_plan_row, 0)
        self.ref_turn_time_entry = QLineEdit("5")
        ref_test_plan_info_layout.addWidget(self.ref_turn_time_entry, ref_test_plan_row, 1)
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

        self.clear_plot_btn_cal = QPushButton("Remove Calibration Plan")
        self.clear_plot_btn_cal.clicked.connect(self._clear_calibration_plan)  # Remove only calibration lines
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

        cal_test_plan_info_layout.addWidget(QLabel("Turn Time (min):"), cal_test_plan_row, 0)
        self.cal_turn_time_entry = QLineEdit("5")
        cal_test_plan_info_layout.addWidget(self.cal_turn_time_entry, cal_test_plan_row, 1)
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

        self.clear_plot_btn_line = QPushButton("Remove Line Plan")
        self.clear_plot_btn_line.clicked.connect(self._clear_line_planning)  # Clear and remove the line plan
        line_plot_control_layout.addWidget(self.clear_plot_btn_line)

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
        # Add navigation toolbar (zoom/pan/home etc.)
        self.toolbar = NavigationToolbar(self.canvas, self.plot_frame)
        plot_layout.addWidget(self.toolbar)

        print(f"Widgets created - param_scroll: {self.param_scroll is not None}, plot_frame: {self.plot_frame is not None}, canvas: {self.canvas is not None}")

        # Add stub for _on_parameter_change if missing
        if not hasattr(self, '_on_parameter_change'):
            def _on_parameter_change(event=None):
                pass
            self._on_parameter_change = _on_parameter_change

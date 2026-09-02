"""Interactive Download Bathymetry dialog for SAT Planner (vendored from Bathymetry Downloader)."""

import json
import os
from datetime import datetime

from PyQt6.QtCore import QEvent, Qt, QTimer, QUrl, pyqtSignal
from PyQt6.QtGui import QDesktopServices
from PyQt6.QtWidgets import (
    QApplication,
    QDialog,
    QFileDialog,
    QLabel,
    QMessageBox,
)

from .data_sources import (
    DATA_SOURCES,
    DEFAULT_DATA_SOURCE,
    UI_DATA_SOURCES,
    UI_DATA_SOURCE_ORDER,
    bbox_sr_for_data_source,
    download_filename_prefix,
    format_cell_size_degrees,
    is_gmrt_source,
    output_mode_filename_suffix,
    shows_output_data_types,
    uses_meter_cell_size,
    uses_web_mercator_map,
)
from .download_module import BathymetryDownloader
from .geo_utils import (
    bbox_to_degrees,
    bbox_to_meters,
    bboxes_overlap,
    clamp_extent_to_bounds,
    extent_compare_tolerance,
    extents_equal,
    transform_bbox,
)
from .gmrt_module import (
    GMRTDownloader,
    MAX_TILES_PER_DOWNLOAD,
    clamp_gmrt_bbox,
    estimate_gmrt_pixels,
    generate_gmrt_tiles,
    needs_tiling_for_spans,
)
from .map_widget import MapWidget
from .service_loader import ServiceInfoLoader
from .ui_layout import build_main_window_ui


class DownloadBathymetryDialog(QDialog):
    """Full interactive bathymetry download dialog with map preview and AOI selection."""

    geotiff_downloaded = pyqtSignal(str)

    def __init__(self, parent=None, initial_source=None):
        super().__init__(parent)
        self.setModal(False)
        self.data_sources = DATA_SOURCES
        self.ui_source_names = list(UI_DATA_SOURCE_ORDER)
        if initial_source in UI_DATA_SOURCES:
            self.current_data_source = initial_source
        else:
            self.current_data_source = DEFAULT_DATA_SOURCE
        self.base_url = self.data_sources[self.current_data_source]["url"]
        # Use known extent as fallback (will be updated when service info loads)
        self.service_extent = self.data_sources[self.current_data_source]["default_extent"]
        self.pixel_size_x = None  # Pixel size in X direction from service
        self.pixel_size_y = None  # Pixel size in Y direction from service
        self.downloader = None
        self.service_loader = None
        self._updating_coordinates = False  # Flag to prevent recursive updates
        self.output_directory = None  # Store selected output directory
        self.config_file = os.path.join(os.path.expanduser("~"), ".sat_planner_bathy_download.json")
        self._data_source_changing = False  # Flag to track when data source is changing
        self._zoom_history = []
        self._zoom_history_index = -1
        self._zoom_history_suppress = False

        self.init_ui()
        self.load_config()  # Load saved output directory
        self.load_service_info()

    def init_ui(self):
        """Initialize the user interface."""
        build_main_window_ui(self)
        if hasattr(self, "data_source_combo"):
            idx = self.data_source_combo.findText(self.current_data_source)
            if idx >= 0:
                self.data_source_combo.setCurrentIndex(idx)

    def current_source_name(self):
        return self.current_data_source

    def load_service_info(self):
        """Load service information from REST endpoint in background thread."""
        # Abort any previous metadata loader so a late reply can't clobber a newer source
        if getattr(self, "service_loader", None) and self.service_loader.isRunning():
            try:
                self.service_loader.loaded.disconnect()
                self.service_loader.error.disconnect()
                self.service_loader.statusMessage.disconnect()
            except TypeError:
                pass

        # Start with default extent and initialize map immediately
        self.log_message("Initializing map widget with default extent...")
        self.init_map_widget()
        if self.map_widget:
            self.log_message("Map widget created successfully on startup")
        else:
            self.log_message("WARNING: Map widget is None after initial creation attempt")

        # GMRT has no ArcGIS service metadata endpoint
        if self._is_gmrt_source():
            ds = self.data_sources[self.current_data_source]
            west, south, east, north = ds["default_extent"]
            self.log_message("GMRT source: using fixed geographic extent (±180°, ±85°)")
            self.on_service_info_loaded({
                "extent": {"xmin": west, "ymin": south, "xmax": east, "ymax": north},
                "pixel_size_x": None,
                "pixel_size_y": None,
            })
            return
        
        # Try to load actual service info in background
        self.service_loader = ServiceInfoLoader(self.base_url)
        self.service_loader.loaded.connect(self.on_service_info_loaded)
        self.service_loader.error.connect(self.on_service_info_error)
        self.service_loader.statusMessage.connect(self.log_message)
        self.service_loader.start()
        
    def on_service_info_loaded(self, service_data):
        """Handle successful service info load."""
        extent_dict = service_data.get("extent", {})
        ds = self.data_sources.get(self.current_data_source, {})
        self.service_extent = (
            extent_dict["xmin"],
            extent_dict["ymin"],
            extent_dict["xmax"],
            extent_dict["ymax"]
        )
        # Keep WGOM (and other Web Mercator) extents in native EPSG:3857
        self.log_message("Service info loaded successfully")
        self.log_message(f"REST endpoint extent (bathymetry data bounds): {self.service_extent}")
        
        # Don't set default bounds here - wait until map loads so extent is correct
        # Default bounds will be set in on_map_first_loaded after map loads with correct extent
        
        # Update cell size dropdown based on pixel size from service
        pixel_size_x = service_data.get("pixel_size_x")
        pixel_size_y = service_data.get("pixel_size_y")
        self.pixel_size_x = pixel_size_x
        self.pixel_size_y = pixel_size_y
        if ds.get("native_resolution_only"):
            if ds.get("configurable_cell_size_degrees"):
                self._set_cell_size_degrees_default(force=self._data_source_changing)
            else:
                self._set_native_cell_size_only()
        elif self._is_gmrt_source():
            self._set_gmrt_cell_size_options(force_default=self._data_source_changing)
        elif pixel_size_x is not None and pixel_size_y is not None:
            base_cell_size = max(abs(pixel_size_x), abs(pixel_size_y))
            self.update_cell_size_options(base_cell_size, force_highest_resolution=self._data_source_changing)
        else:
            self.log_message("Warning: Pixel size not available from service, using default cell sizes")
            self.update_cell_size_options(4.0, force_highest_resolution=self._data_source_changing)
        
        self._data_source_changing = False
        
        # Ensure map widget is initialized (this will remove loading label)
        if self.map_widget is None:
            self.log_message("Initializing map widget...")
            self.init_map_widget()
        
        # Update map extent (but don't reload if map widget already exists and is loading)
        if self.map_widget:
            # Prefer explicit source-change flag; base_url is often updated before this runs
            force_reload = getattr(self, "_force_map_reload", False)
            self._force_map_reload = False
            url_changed = force_reload or (self.map_widget.base_url != self.base_url)
            
            # Update pixel sizes in map widget from service info (this happens after service loads)
            self.map_widget.pixel_size_x = self.pixel_size_x
            self.map_widget.pixel_size_y = self.pixel_size_y
            
            self.map_widget.base_url = self.base_url
            
            # Update raster functions and display URL from current data source
            new_raster_function = self.data_sources[self.current_data_source]["bathymetry_raster_function"]
            new_hillshade_raster_function = self.data_sources[self.current_data_source]["hillshade_raster_function"]
            self.map_widget.raster_function = new_raster_function
            self.map_widget.hillshade_raster_function = new_hillshade_raster_function
            self.map_widget.display_url = self.data_sources[self.current_data_source].get("display_url")
            self.map_widget.land_display_url = self.data_sources[self.current_data_source].get("land_display_url")
            self.map_widget.bbox_sr = self._bbox_sr_for_data_source(self.current_data_source)
            self.map_widget.preview_url = self.data_sources[self.current_data_source].get("preview_url")
            self.map_widget.gmrt_mask = (
                hasattr(self, "gmrt_mask_checkbox") and self.gmrt_mask_checkbox.isChecked()
            )
            
            # Check if there's a pending selection to preserve
            if hasattr(self, '_pending_selection') and self._pending_selection:
                # Use the pending selection extent instead of full service extent
                selection_extent = self._pending_selection
                self.log_message(f"Preserving selection, will zoom to it: {selection_extent}")
                self.map_widget.extent = selection_extent
                self.map_widget._requested_extent = selection_extent
            else:
                # Always set extent to REST endpoint service extent as a baseline
                # This ensures the map shows exactly the bathymetry data bounds from the REST endpoint
                self.log_message(f"Updating map extent to REST endpoint extent: {self.service_extent}")
                self.map_widget.extent = self.service_extent
                # Also update _requested_extent to ensure coordinate conversion is correct
                # This ensures the map displays exactly the REST endpoint bounds, not a rounded or adjusted version
                self.map_widget._requested_extent = self.service_extent
            
            # Update service extent in map widget
            self.map_widget.service_extent = self.service_extent
            
            # CRITICAL: Always update selected_bbox_world to REST endpoint extent if it matches default extent
            # This ensures the box shows the exact REST endpoint bounds, not the default extent
            # Check if selected_bbox_world is None OR if it matches the default extent (needs update)
            default_extent = self.data_sources[self.current_data_source]["default_extent"]
            needs_update = (
                self.map_widget.selected_bbox_world is None or
                self.map_widget.selected_bbox_world == default_extent
            )
            if needs_update and not (hasattr(self, '_pending_selection') and self._pending_selection):
                self.log_message(f"Updating selected_bbox_world from {self.map_widget.selected_bbox_world} to REST endpoint extent {self.service_extent}")
                self.map_widget.selected_bbox_world = self.service_extent
                self.map_widget.set_selection_validity(True)
                self.selected_bbox = self.service_extent
                # Update coordinate display to show REST endpoint bounds
                self.update_coordinate_display(*self.service_extent, update_map=False)
            
            # CRITICAL: Always reload map with REST endpoint extent to ensure it shows exact bathymetry data bounds
            # Check if the map was loaded with a different extent (e.g., default extent)
            current_extent = self.map_widget.extent
            default_extent = self.data_sources[self.current_data_source]["default_extent"]
            
            # Check if map was loaded with default extent
            _tol = 1.0 if self._uses_web_mercator_map() else 1e-5
            extent_matches_default = (
                abs(current_extent[0] - default_extent[0]) < _tol and
                abs(current_extent[1] - default_extent[1]) < _tol and
                abs(current_extent[2] - default_extent[2]) < _tol and
                abs(current_extent[3] - default_extent[3]) < _tol
            )
            
            # Ensure default bounds are set
            if self.map_widget.selected_bbox_world is None:
                self.map_widget.selected_bbox_world = self.service_extent
                self.map_widget.set_selection_validity(True)
                self.selected_bbox = self.service_extent
                self.map_widget.service_extent = self.service_extent
                self.log_message(f"Set default bounds to REST endpoint extent: {self.service_extent}")
            
            # Check if there's a pending selection to restore
            if hasattr(self, '_pending_selection') and self._pending_selection:
                # Restore the selection - this will zoom to it and reload the map
                self.log_message(f"Restoring pending selection: {self._pending_selection}")
                QTimer.singleShot(300, lambda: self._restore_selection())
            # Force reload if URL changed (data source switch) or if extent differs
            elif url_changed or current_extent != self.service_extent or extent_matches_default:
                if url_changed:
                    self.log_message(f"Data source URL changed, reloading map with new service...")
                else:
                    self.log_message(f"Map extent ({current_extent}) differs from REST endpoint extent ({self.service_extent}), zooming to REST endpoint bounds...")
                if not getattr(self.map_widget, '_loading', False):
                    # CRITICAL: Use zoom_to_selection directly (like when user hits return)
                    # This recalculates the extent with padding and positions the box correctly
                    # Don't reload first - zoom_to_selection will reload with the correct extent
                    # Wait a bit longer to ensure widget is fully sized
                    QTimer.singleShot(300, lambda: self.zoom_to_selection(*self.service_extent, record_history=False))
            elif not self.map_widget.map_loaded and not getattr(self.map_widget, '_loading', False):
                self.log_message("Map not loaded yet, will load and zoom to REST endpoint extent...")
                # Map hasn't loaded yet - use zoom_to_selection which will load the map with correct extent
                # zoom_to_selection will:
                # 1. Calculate extent with padding
                # 2. Set map widget extent
                # 3. Call load_map() to reload basemap and raster layers
                # Wait a bit longer to ensure widget is fully sized
                QTimer.singleShot(300, lambda: self.zoom_to_selection(*self.service_extent, record_history=False))
            else:
                self.log_message("Map already loaded with REST endpoint extent")
        else:
            self.log_message("ERROR: Map widget is None after initialization attempt")
            
    def on_service_info_error(self, error_message):
        """Handle service info load error with helpful message."""
        # Log the error
        self.log_message(error_message)
        
        # Show error message with suggestion to check for updates
        msg = QMessageBox(self)
        msg.setIcon(QMessageBox.Icon.Warning)
        msg.setWindowTitle("Connection Error")
        msg.setText(
            f"Unable to connect to the REST endpoint.\n\n"
            f"{error_message}\n\n"
            f"If this problem persists, please:\n"
            f"1. Check for a new version at: https://github.com/seamapper/GEBCO_Downloader\n"
            f"2. Contact: pjohnson@ccom.unh.edu"
        )
        msg.setStandardButtons(QMessageBox.StandardButton.Ok)
        msg.exec()
        
        # Continue with default extent - map should already be initialized
            
    def init_map_widget(self):
        """Initialize the map widget."""
        if self.service_extent is None:
            self.log_message("ERROR: service_extent is None, cannot initialize map")
            return
            
        # Host map inside the black canvas (falls back to Map group layout)
        if hasattr(self, "map_canvas") and self.map_canvas and self.map_canvas.layout():
            layout = self.map_canvas.layout()
        elif hasattr(self, "map_group") and self.map_group and self.map_group.layout():
            layout = self.map_group.layout()
        else:
            self.log_message("ERROR: map canvas/layout not found")
            return
            
        # Remove loading label if it exists - try multiple approaches
        label_removed = False
        
        # First, try to remove via the stored reference
        if hasattr(self, 'loading_label') and self.loading_label:
            try:
                layout.removeWidget(self.loading_label)
                self.loading_label.hide()
                self.loading_label.setParent(None)
                self.loading_label.deleteLater()
                self.loading_label = None
                label_removed = True
            except:
                pass
        
        # Also search for any QLabel with "Loading" text in the layout
        if not label_removed:
            for i in reversed(range(layout.count())):
                item = layout.itemAt(i)
                if item and item.widget():
                    widget = item.widget()
                    if isinstance(widget, QLabel) and "Loading" in widget.text():
                        try:
                            layout.removeWidget(widget)
                            widget.hide()
                            widget.setParent(None)
                            widget.deleteLater()
                            label_removed = True
                            break
                        except:
                            pass
                    
        # Create map widget if it doesn't exist
        if self.map_widget is None:
            try:
                # Get raster functions from current data source
                raster_function = self.data_sources[self.current_data_source]["bathymetry_raster_function"]
                hillshade_raster_function = self.data_sources[self.current_data_source]["hillshade_raster_function"]
                ds = self.data_sources[self.current_data_source]
                show_basemap = bool(ds.get("show_basemap", False))
                show_hillshade = bool(ds.get("show_hillshade", False))
                # Multiply blend lets hillshade show through the colored bathymetry (CCOM behavior)
                use_blend = show_hillshade
                self.log_message(f"Creating MapWidget with extent: {self.service_extent}, raster function: {raster_function}")
                display_url = ds.get("display_url")
                land_display_url = ds.get("land_display_url")
                self.map_widget = MapWidget(self.base_url, self.service_extent, raster_function=raster_function, show_basemap=show_basemap, show_hillshade=show_hillshade, use_blend=use_blend, hillshade_raster_function=hillshade_raster_function, display_url=display_url, land_display_url=land_display_url)
                self.map_widget.bbox_sr = self._bbox_sr_for_data_source(self.current_data_source)
                self.map_widget.preview_url = ds.get("preview_url")
                self.map_widget.gmrt_mask = (
                    hasattr(self, "gmrt_mask_checkbox") and self.gmrt_mask_checkbox.isChecked()
                )
                self.map_widget.bathymetry_opacity = 1.0  # Full opacity
                # Sync legend visibility with checkbox state
                if hasattr(self, 'legend_checkbox'):
                    self.map_widget.show_legend = self.legend_checkbox.isChecked()
                    self.map_widget.show_aoi = self.aoi_checkbox.isChecked()
                # Store service extent in map widget
                self.map_widget.service_extent = self.service_extent
                # Store pixel sizes for raster function selection
                self.map_widget.pixel_size_x = self.pixel_size_x
                self.map_widget.pixel_size_y = self.pixel_size_y
                self.map_widget.selectionChanged.connect(self.on_selection_changed)
                self.map_widget.selectionCompleted.connect(self.on_selection_completed)
                self.map_widget.mapFirstLoaded.connect(self.on_map_first_loaded)
                self.map_widget.userViewChanged.connect(self._on_user_map_view_changed)
                self.map_widget.statusMessage.connect(self.log_message)  # Connect status messages to log
                layout.addWidget(self.map_widget)
                self.map_widget.show()
                # Force UI update
                if hasattr(self, "map_canvas") and self.map_canvas:
                    self.map_canvas.update()
                self.map_group.update()
                layout.update()
                self.log_message("MapWidget created and added to layout successfully")
                
                # Don't set default bounds here - wait until map loads so extent is correct
                # Default bounds will be set in on_map_first_loaded after map loads
                
                # Trigger map load after a short delay to ensure widget is sized
                self.log_message("Scheduling map load in 200ms...")
                QTimer.singleShot(200, lambda: self.trigger_map_load())
            except Exception as e:
                self.log_message(f"ERROR creating MapWidget: {e}")
                import traceback
                self.log_message(traceback.format_exc())
                self.map_widget = None
                
    def trigger_map_load(self):
        """Trigger map load - called via timer."""
        # Don't load here - let on_service_info_loaded handle it
        # This ensures REST endpoint extent is available and widget is properly sized
        if self.map_widget:
            if self.service_extent:
                self.log_message(f"Widget ready, waiting for service info to trigger map load...")
                self.log_message(f"Widget size: {self.map_widget.width()}x{self.map_widget.height()}")
                self.log_message(f"REST endpoint extent: {self.service_extent}")
            else:
                self.log_message("REST endpoint extent not available yet, waiting for service info to load...")
        else:
            self.log_message("ERROR: map_widget is None when trying to trigger load")
            
    def fit_to_extent(self):
        """Fit map to full service extent."""
        if self.map_widget and self.service_extent:
            self.map_widget.selected_bbox_world = self.service_extent
            self.map_widget.set_selection_validity(True)
            self.selected_bbox = self.service_extent
            self.map_widget.service_extent = self.service_extent
            self.zoom_to_selection(*self.service_extent)
            # Update coordinate display to show the service extent
            self.update_coordinate_display(*self.service_extent, update_map=False)

    def _current_zoom_state(self):
        """Return the current map view extent and AOI as a history snapshot."""
        if not self.map_widget:
            return None
        extent = self.map_widget._requested_extent or self.map_widget.extent
        if extent is None:
            return None
        aoi = self.map_widget.selected_bbox_world
        if aoi is not None:
            aoi = tuple(aoi)
        return (tuple(extent), aoi)

    def _seed_zoom_history(self):
        """Initialize zoom history with the current map state."""
        state = self._current_zoom_state()
        if state is None:
            return
        self._zoom_history = [state]
        self._zoom_history_index = 0
        self._update_zoom_nav_buttons()

    def _append_zoom_history(self):
        """Record the current map state in zoom history."""
        if self._zoom_history_suppress or not self.map_widget:
            return
        state = self._current_zoom_state()
        if state is None:
            return
        if (
            self._zoom_history
            and 0 <= self._zoom_history_index < len(self._zoom_history)
            and self._zoom_history[self._zoom_history_index] == state
        ):
            return
        self._zoom_history = self._zoom_history[: self._zoom_history_index + 1]
        self._zoom_history.append(state)
        self._zoom_history_index = len(self._zoom_history) - 1
        self._update_zoom_nav_buttons()

    def _reset_zoom_history(self):
        """Clear zoom history (e.g. after switching data source)."""
        self._zoom_history = []
        self._zoom_history_index = -1
        self._update_zoom_nav_buttons()

    def _update_zoom_nav_buttons(self):
        """Enable Zoom Prev/Next based on the current position in history."""
        if hasattr(self, "zoom_back_btn"):
            self.zoom_back_btn.setEnabled(self._zoom_history_index > 0)
        if hasattr(self, "zoom_next_btn"):
            self.zoom_next_btn.setEnabled(
                0 <= self._zoom_history_index < len(self._zoom_history) - 1
            )

    def _apply_zoom_state(self, state):
        """Restore a saved map view extent and AOI without recording history."""
        if not self.map_widget:
            return
        extent, aoi = state
        self._zoom_history_suppress = True
        try:
            self.map_widget.extent = extent
            self.map_widget._requested_extent = extent
            self.map_widget.selection_start = None
            self.map_widget.selection_end = None
            self.map_widget.is_selecting = False
            if aoi:
                self.map_widget.selected_bbox_world = tuple(aoi)
                self.map_widget.set_selection_validity(True)
                self.selected_bbox = tuple(aoi)
                self.update_coordinate_display(*aoi, update_map=False)
            else:
                self.map_widget.selected_bbox_world = None
                self.selected_bbox = None
                self._updating_coordinates = True
                try:
                    for field in self._coordinate_fields:
                        field.clear()
                finally:
                    self._updating_coordinates = False
            self.map_widget.load_map()
        finally:
            self._zoom_history_suppress = False
        self.check_and_update_download_button()
        self._update_zoom_nav_buttons()

    def _on_user_map_view_changed(self):
        """Record map wheel/pan changes in zoom history."""
        if self._zoom_history_suppress:
            return
        aoi = self.map_widget.selected_bbox_world if self.map_widget else None
        self.selected_bbox = tuple(aoi) if aoi else None
        if not self._zoom_history:
            self._seed_zoom_history()
        else:
            self._append_zoom_history()

    def zoom_back(self):
        """Restore the previous map view and AOI."""
        if self._zoom_history_index <= 0:
            return
        self._zoom_history_index -= 1
        self._apply_zoom_state(self._zoom_history[self._zoom_history_index])

    def zoom_next(self):
        """Move forward to the next map view and AOI in history."""
        if self._zoom_history_index >= len(self._zoom_history) - 1:
            return
        self._zoom_history_index += 1
        self._apply_zoom_state(self._zoom_history[self._zoom_history_index])
            
    def clear_selection(self):
        """Clear the map selection."""
        if self.map_widget:
            self.map_widget.clear_selection()
        # Remove bold formatting from download button when selection is cleared
        if hasattr(self, 'download_btn'):
            font = self.download_btn.font()
            font.setBold(False)
            self.download_btn.setFont(font)
    
    def refresh_map(self):
        """Refresh the map display for the currently shown area."""
        if self.map_widget:
            self.log_message("Refreshing map display...")
            self.map_widget.load_map()
        else:
            self.log_message("Warning: Map widget not available for refresh")
    
    def export_map_image(self):
        """Export the current map display as a PNG image."""
        if not self.map_widget:
            self.log_message("Warning: Map widget not available for export")
            QMessageBox.warning(self, "Export Error", "Map widget not available.")
            return
        
        if not self.map_widget.map_loaded:
            self.log_message("Warning: Map not loaded yet")
            QMessageBox.warning(self, "Export Error", "Map is not loaded yet. Please wait for the map to load.")
            return
        
        # Generate default filename with timestamp
        from datetime import datetime
        date_time_str = datetime.now().strftime("%Y%m%d_%H%M%S")
        default_filename = f"GEBCO_Map_{date_time_str}.png"
        
        # Determine save location
        if self.output_directory and os.path.isdir(self.output_directory):
            default_path = os.path.join(self.output_directory, default_filename)
        else:
            default_path = default_filename
        
        # Prompt for save location
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Map Image",
            default_path,
            "PNG Files (*.png);;All Files (*)"
        )
        
        if not file_path:
            return  # User cancelled
        
        try:
            # Grab the map widget as a pixmap
            pixmap = self.map_widget.grab()
            
            if pixmap.isNull():
                raise Exception("Failed to capture map widget")
            
            # Save to file
            if not pixmap.save(file_path, "PNG"):
                raise Exception("Failed to save PNG file")
            
            self.log_message(f"✓ Map image exported: {file_path}")
            QMessageBox.information(self, "Success", f"Map image saved to:\n{file_path}")
            
        except Exception as e:
            error_msg = f"Error exporting map image: {str(e)}"
            self.log_message(f"✗ {error_msg}")
            QMessageBox.critical(self, "Export Error", error_msg)
            
    # Raster function is fixed to "DAR - StdDev - BlueGreen" - no handler needed
            
    def on_legend_toggled(self, state):
        """Handle legend checkbox toggle."""
        if self.map_widget:
            show_legend = (state == Qt.CheckState.Checked.value or state == 2)
            self.map_widget.show_legend = show_legend
            # Just update display (no need to reload map)
            self.map_widget.update()
    
    def on_aoi_toggled(self, state):
        """Handle AOI checkbox toggle."""
        show_aoi = (state == Qt.CheckState.Checked.value or state == 2)
        if self.map_widget:
            self.map_widget.show_aoi = show_aoi
            # Just update display (no need to reload map)
            self.map_widget.update()
        if not show_aoi:
            # When AoI is unchecked, remember Legend state and uncheck Legend
            self._legend_was_on_before_aoi_off = self.legend_checkbox.isChecked()
            if self._legend_was_on_before_aoi_off:
                self.legend_checkbox.setChecked(False)
        else:
            # When AoI is re-enabled, restore Legend if it was on before
            if self._legend_was_on_before_aoi_off:
                self.legend_checkbox.setChecked(True)
                
    def check_and_update_download_button(self):
        """Check if selection is valid and within size limits, update download button state."""
        # Check if there's a valid selection
        bbox = None
        if hasattr(self, 'selected_bbox') and self.selected_bbox:
            bbox = self.selected_bbox
        elif self.map_widget:
            bbox = self.map_widget.get_selection_bbox()
        
        if not bbox:
            # No selection - disable button and clear selection validity
            self.download_btn.setEnabled(False)
            # Remove bold formatting when no selection
            font = self.download_btn.font()
            font.setBold(False)
            self.download_btn.setFont(font)
            if self.map_widget:
                self.map_widget.set_selection_validity(True)  # Default to valid (no selection shown)
            return
        
        # Check if selection exceeds maximum size (bbox is in GCS: west, south, east, north)
        try:
            west, south, east, north = bbox
            ds = self.data_sources.get(self.current_data_source, {})
            if self._uses_web_mercator_map():
                # selected_bbox is already in EPSG:3857
                xmin_m, ymin_m, xmax_m, ymax_m = bbox
                cell_size = self._get_selected_cell_size_meters()
                pixels_width = int((xmax_m - xmin_m) / cell_size)
                pixels_height = int((ymax_m - ymin_m) / cell_size)
            elif self._is_gmrt_source():
                pixels_width, pixels_height = estimate_gmrt_pixels(
                    west, south, east, north, self._get_selected_cell_size_meters()
                )
            elif ds.get("native_resolution_only") or ds.get("configurable_cell_size_degrees"):
                deg_per_pixel = self._get_degrees_per_pixel_for_output()
                pixels_width = int((east - west) / deg_per_pixel)
                pixels_height = int((north - south) / abs(deg_per_pixel))
            else:
                cell_size = self._get_selected_cell_size_meters()
                width_m = (east - west) * 111320 * 0.5
                height_m = (north - south) * 110540
                pixels_width = int(width_m / cell_size)
                pixels_height = int(height_m / cell_size)
            
            # No size limit - always enable download button
            # Warning dialog will be shown when downloading large datasets
            is_valid = True
            
            # Update map widget selection color (always valid now)
            if self.map_widget:
                self.map_widget.set_selection_validity(is_valid)
            
            # Enable download button unless a multi-output source has no type selected
            self.download_btn.setEnabled(True)
            if self._shows_output_data_types() and hasattr(self, "check_combined") and not self._collect_output_type_requests():
                self.download_btn.setEnabled(False)
            # Make text bold only if this is a user manual selection (not initial dataset bounds)
            is_initial_bounds = False
            if hasattr(self, 'service_extent') and self.service_extent:
                se = self.service_extent
                tol = self._extent_compare_tolerance()
                if (abs(se[0] - west) < tol and abs(se[1] - south) < tol and
                    abs(se[2] - east) < tol and abs(se[3] - north) < tol):
                    is_initial_bounds = True
            
            # Only make bold if it's NOT the initial dataset bounds
            font = self.download_btn.font()
            font.setBold(not is_initial_bounds)
            self.download_btn.setFont(font)
        except Exception:
            # Error calculating - disable button to be safe
            self.download_btn.setEnabled(False)
            # Remove bold formatting on error
            font = self.download_btn.font()
            font.setBold(False)
            self.download_btn.setFont(font)
            if self.map_widget:
                self.map_widget.set_selection_validity(True)  # Default to valid on error
    
    def _shows_output_data_types(self):
        """Return True when the current data source supports multiple output grid types."""
        return shows_output_data_types(self.data_sources, self.current_data_source)

    def _is_gmrt_source(self, data_source_name=None):
        """Return True for Lamont GMRT GridServer sources."""
        return is_gmrt_source(self.data_sources, data_source_name or self.current_data_source)

    def _uses_web_mercator_map(self, data_source_name=None):
        """Return True when map extents/selection are in EPSG:3857 (WGOM)."""
        return uses_web_mercator_map(self.data_sources, data_source_name or self.current_data_source)

    def _uses_meter_cell_size(self, data_source_name=None):
        """Return True when the UI shows meter cell-size options (WGOM or GMRT)."""
        return uses_meter_cell_size(self.data_sources, data_source_name or self.current_data_source)

    def _transform_bbox(self, bbox, from_crs, to_crs):
        """Transform (xmin, ymin, xmax, ymax) between CRS definitions."""
        return transform_bbox(bbox, from_crs, to_crs)

    def _bbox_to_meters(self, bbox_4326):
        """Convert a GCS selection bbox to Web Mercator meters."""
        return bbox_to_meters(bbox_4326)

    def _bbox_to_degrees(self, bbox_3857):
        """Convert a Web Mercator selection bbox to GCS degrees."""
        return bbox_to_degrees(bbox_3857)

    def _map_bbox_to_display(self, bbox):
        """Convert map-CRS bbox to GCS for West/South/East/North fields."""
        if self._uses_web_mercator_map():
            return self._bbox_to_degrees(bbox)
        return bbox

    def _display_bbox_to_map(self, bbox_4326):
        """Convert GCS field values to the active map CRS."""
        if self._uses_web_mercator_map():
            return self._bbox_to_meters(bbox_4326)
        if self._is_gmrt_source():
            return clamp_gmrt_bbox(*bbox_4326)
        return bbox_4326

    def _extent_compare_tolerance(self):
        """Tolerance for comparing extents in the active map CRS."""
        return extent_compare_tolerance(self._uses_web_mercator_map())

    def _get_selected_cell_size_meters(self):
        """Return the selected meter cell size, or source default if unavailable."""
        ds = self.data_sources.get(self.current_data_source, {})
        default = float(ds.get("default_cell_size_meters", 4.0))
        if hasattr(self, "cell_size_combo") and self.cell_size_combo.count():
            try:
                return float(self.cell_size_combo.currentText())
            except (ValueError, TypeError):
                pass
        return default

    def _set_gmrt_cell_size_options(self, force_default=False):
        """Populate meter cell-size combo with GMRT presets."""
        if not hasattr(self, "cell_size_combo"):
            return
        ds = self.data_sources.get(self.current_data_source, {})
        options = ds.get("cell_size_meters_options", [60, 120, 240, 480, 960])
        default = ds.get("default_cell_size_meters", 480)
        current = None if force_default else (self.cell_size_combo.currentText() if self.cell_size_combo.count() else None)
        if hasattr(self, "cell_size_label"):
            self.cell_size_label.setText("Cell Size (m):")
        self.cell_size_combo.blockSignals(True)
        self.cell_size_combo.clear()
        self.cell_size_combo.addItems([str(int(v)) if float(v).is_integer() else str(v) for v in options])
        target = str(int(default)) if float(default).is_integer() else str(default)
        if current and current in [self.cell_size_combo.itemText(i) for i in range(self.cell_size_combo.count())]:
            self.cell_size_combo.setCurrentText(current)
        else:
            self.cell_size_combo.setCurrentText(target)
        self.cell_size_combo.blockSignals(False)
        if hasattr(self, "selected_bbox") and self.selected_bbox:
            self.update_coordinate_display(*self.selected_bbox, update_map=False)

    def _get_output_type_tid_url(self):
        """Return the TID ImageServer URL for masking bathymetry/land/direct outputs."""
        tid_name = f"{self.current_data_source} TID"
        tid_ds = self.data_sources.get(tid_name, {})
        return tid_ds.get("url")

    def _get_download_filename_prefix(self):
        """Return the filename prefix for downloaded GeoTIFFs."""
        return download_filename_prefix(self.data_sources, self.current_data_source)

    def _collect_output_type_requests(self):
        """Return selected output grid modes for sources with multiple output types."""
        if not self._shows_output_data_types() or not hasattr(self, "check_combined"):
            return []
        requests = []
        if self.check_combined.isChecked():
            requests.append(("combined", None))
        if self.check_bathymetry_only.isChecked():
            requests.append(("bathymetry_only", None))
        if self.check_land_only.isChecked():
            requests.append(("land_only", None))
        if self.check_direct_measurements_only.isChecked():
            requests.append(("direct_measurements_only", None))
        if self.check_direct_unknown_measurements_only.isChecked():
            requests.append(("direct_unknown_measurements_only", None))
        return requests

    def _output_mode_filename_suffix(self, mode):
        """Map an output mode to its filename suffix."""
        return output_mode_filename_suffix(mode)

    def _get_native_pixel_size_degrees(self):
        """Return native cell size in degrees from service info or data source config."""
        ds = self.data_sources.get(self.current_data_source, {})
        if self.pixel_size_x is not None:
            return max(abs(self.pixel_size_x), abs(self.pixel_size_y or self.pixel_size_x))
        return ds.get("native_pixel_size_degrees", 0.004166666666666667)

    def _format_cell_size_degrees(self, value):
        """Format a cell size in degrees for display in the UI."""
        return format_cell_size_degrees(value)

    def _get_output_cell_size_degrees(self):
        """Return the user-selected cell size in degrees, or native if unset/invalid."""
        if hasattr(self, "cell_size_degrees_edit") and self.cell_size_degrees_container.isVisible():
            text = self.cell_size_degrees_edit.text().strip()
            if text:
                try:
                    value = float(text)
                    if value > 0:
                        return value
                except ValueError:
                    pass
        return self._get_native_pixel_size_degrees()

    def _get_degrees_per_pixel_for_output(self):
        """Return degrees-per-pixel used for output grid sizing and downloads."""
        ds = self.data_sources.get(self.current_data_source, {})
        if ds.get("configurable_cell_size_degrees"):
            return self._get_output_cell_size_degrees()
        if ds.get("native_resolution_only"):
            return self._get_native_pixel_size_degrees()
        return ds.get("native_pixel_size_degrees", 0.004166666666666667)

    def _set_cell_size_degrees_default(self, force=False):
        """Set Cell Size (deg) field to the service native resolution."""
        if not hasattr(self, "cell_size_degrees_edit"):
            return
        if force or not self.cell_size_degrees_edit.text().strip():
            self.cell_size_degrees_edit.setText(
                self._format_cell_size_degrees(self._get_native_pixel_size_degrees())
            )
        if hasattr(self, "selected_bbox") and self.selected_bbox:
            xmin, ymin, xmax, ymax = self.selected_bbox
            self.update_coordinate_display(xmin, ymin, xmax, ymax, update_map=False)

    def _on_cell_size_degrees_changed(self):
        """Handle manual entry of cell size in degrees."""
        ds = self.data_sources.get(self.current_data_source, {})
        if not ds.get("configurable_cell_size_degrees"):
            return
        text = self.cell_size_degrees_edit.text().strip()
        if not text:
            self._set_cell_size_degrees_default(force=True)
            return
        try:
            value = float(text)
            if value <= 0:
                raise ValueError("Cell size must be positive")
            self.cell_size_degrees_edit.setText(self._format_cell_size_degrees(value))
        except ValueError:
            QMessageBox.warning(self, "Invalid Cell Size", "Enter a positive numeric cell size in degrees.")
            self._set_cell_size_degrees_default(force=True)
            return
        if hasattr(self, "selected_bbox") and self.selected_bbox:
            xmin, ymin, xmax, ymax = self.selected_bbox
            self.update_coordinate_display(xmin, ymin, xmax, ymax, update_map=False)
        self.check_and_update_download_button()

    def _set_native_cell_size_only(self):
        """Set cell size dropdown to single 'Native' option (for sources with native_resolution_only)."""
        if not hasattr(self, 'cell_size_combo'):
            return
        if hasattr(self, 'cell_size_label'):
            self.cell_size_label.setText("Resolution:")
        self.cell_size_combo.clear()
        self.cell_size_combo.addItems(["Native"])
        if hasattr(self, 'selected_bbox') and self.selected_bbox:
            xmin, ymin, xmax, ymax = self.selected_bbox
            self.update_coordinate_display(xmin, ymin, xmax, ymax, update_map=False)
    
    def update_cell_size_options(self, base_cell_size, force_highest_resolution=False):
        """Update cell size dropdown options based on base cell size from service.
        
        Args:
            base_cell_size: The base cell size (max of pixelSizeX and pixelSizeY)
            force_highest_resolution: If True, always select the highest resolution (smallest cell size)
        """
        if not hasattr(self, 'cell_size_combo'):
            return
        if hasattr(self, 'cell_size_label'):
            self.cell_size_label.setText("Cell Size (m):")
        
        # Calculate the five options: base, 2x, 3x, 4x, 5x
        option1 = base_cell_size  # Highest resolution (smallest cell size)
        option2 = base_cell_size * 2
        option3 = base_cell_size * 3
        option4 = base_cell_size * 4
        option5 = base_cell_size * 5
        
        # Store current selection if exists (only if not forcing highest resolution)
        current_text = self.cell_size_combo.currentText() if not force_highest_resolution else None
        
        # Clear and repopulate dropdown
        self.cell_size_combo.clear()
        self.cell_size_combo.addItems([f"{option1:.1f}", f"{option2:.1f}", f"{option3:.1f}", f"{option4:.1f}", f"{option5:.1f}"])
        
        if force_highest_resolution:
            # Always select the highest resolution (first option, smallest cell size)
            self.cell_size_combo.setCurrentIndex(0)
        else:
            # Try to restore previous selection if it matches one of the new options
            # Otherwise, select the first (smallest) option
            try:
                current_value = float(current_text)
                # Find closest match
                options = [option1, option2, option3, option4, option5]
                closest_idx = min(range(len(options)), key=lambda i: abs(options[i] - current_value))
                self.cell_size_combo.setCurrentIndex(closest_idx)
            except (ValueError, TypeError):
                # If previous selection was invalid, default to first option
                self.cell_size_combo.setCurrentIndex(0)
        
        # Update pixel count if selection exists
        if hasattr(self, 'selected_bbox') and self.selected_bbox:
            xmin, ymin, xmax, ymax = self.selected_bbox
            self.update_coordinate_display(xmin, ymin, xmax, ymax, update_map=False)
    
    def on_cell_size_changed(self, cell_size_text):
        """Handle cell size change - update pixel count if selection exists."""
        if not hasattr(self, 'cell_size_combo'):
            return
        # Update pixel count display if there's a current selection
        if hasattr(self, 'selected_bbox') and self.selected_bbox:
            xmin, ymin, xmax, ymax = self.selected_bbox
            self.update_coordinate_display(xmin, ymin, xmax, ymax, update_map=False)
        # Update download button state
        self.check_and_update_download_button()
            
    def eventFilter(self, watched, event):
        """Commit coordinate edits when focus leaves the selected-area fields."""
        if (
            watched in getattr(self, "_coordinate_fields", ())
            and event.type() == QEvent.Type.FocusOut
        ):
            # Defer so tab order between fields is resolved before checking focus
            QTimer.singleShot(0, self._commit_geographic_if_focus_left_group)
        return super().eventFilter(watched, event)

    def _commit_geographic_if_focus_left_group(self):
        """Apply coordinate field values only after focus leaves all four fields."""
        if self._updating_coordinates:
            return
        focus = QApplication.focusWidget()
        if focus in self._coordinate_fields:
            return
        self.on_geographic_changed()

    def on_geographic_changed(self):
        """Handle committed manual entry in geographic coordinate fields (always GCS)."""
        if self._updating_coordinates:
            return
            
        try:
            # Get values from Geographic fields
            west_text = self.west_edit.text().strip()
            south_text = self.south_edit.text().strip()
            east_text = self.east_edit.text().strip()
            north_text = self.north_edit.text().strip()
            
            # Check if all fields have values
            if not (west_text and south_text and east_text and north_text):
                return
            
            west = float(west_text)
            south = float(south_text)
            east = float(east_text)
            north = float(north_text)
            
            # Validate that min < max
            if west >= east or south >= north:
                QMessageBox.warning(self, "Invalid Coordinates", "West must be less than East and South must be less than North.")
                return
            
            # Convert GCS entry to map CRS, then snap in map units
            map_bbox = self._display_bbox_to_map((west, south, east, north))
            snapped = self._snap_bounds_to_cell_size(*map_bbox)
            snapped_bbox = snapped[:4]
            was_adjusted = snapped[4]
            
            if was_adjusted:
                disp_before = (west, south, east, north)
                disp_after = self._map_bbox_to_display(snapped_bbox)
                self.log_message(
                    f"Selection bounds adjusted to align with cell size grid: "
                    f"({disp_before[0]:.6f}, {disp_before[1]:.6f}, {disp_before[2]:.6f}, {disp_before[3]:.6f}) → "
                    f"({disp_after[0]:.6f}, {disp_after[1]:.6f}, {disp_after[2]:.6f}, {disp_after[3]:.6f})"
                )
            
            self.update_coordinate_display(*snapped_bbox, update_map=True)
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter valid numeric coordinates.")
            self.check_and_update_download_button()  # Disable button on invalid input
            
    def update_coordinate_display(self, xmin, ymin, xmax, ymax, update_map=True):
        """Update coordinate fields from a map-CRS bbox; fields always show GCS."""
        if self._updating_coordinates:
            return
        west, south, east, north = self._map_bbox_to_display((xmin, ymin, xmax, ymax))
        self._updating_coordinates = True
        try:
            self.west_edit.setText(f"{west:.6f}")
            self.south_edit.setText(f"{south:.6f}")
            self.east_edit.setText(f"{east:.6f}")
            self.north_edit.setText(f"{north:.6f}")
            if update_map:
                self.selected_bbox = (xmin, ymin, xmax, ymax)
                if self.map_widget:
                    self.zoom_to_selection(xmin, ymin, xmax, ymax)
            else:
                self.selected_bbox = (xmin, ymin, xmax, ymax)
        finally:
            self._updating_coordinates = False
        
        # Update download button state
        self.check_and_update_download_button()
        
        # Calculate expected number of pixels in map CRS
        try:
            ds = self.data_sources.get(self.current_data_source, {})
            if self._uses_web_mercator_map():
                cell_size = self._get_selected_cell_size_meters()
                pixels_width = int((xmax - xmin) / cell_size)
                pixels_height = int((ymax - ymin) / cell_size)
                cell_size_label = f"{cell_size}m"
            elif self._is_gmrt_source():
                cell_size = self._get_selected_cell_size_meters()
                pixels_width, pixels_height = estimate_gmrt_pixels(xmin, ymin, xmax, ymax, cell_size)
                cell_size_label = f"{cell_size}m"
            elif ds.get("native_resolution_only") or ds.get("configurable_cell_size_degrees"):
                deg_per_pixel = self._get_degrees_per_pixel_for_output()
                pixels_width = int((xmax - xmin) / deg_per_pixel)
                pixels_height = int((ymax - ymin) / abs(deg_per_pixel))
                cell_size_label = "native"
            else:
                width_meters = (east - west) * 111320 * 0.5  # approx at mid-lat
                height_meters = (north - south) * 110540
                cell_size = self._get_selected_cell_size_meters()
                pixels_width = int(width_meters / cell_size)
                pixels_height = int(height_meters / cell_size)
                cell_size_label = f"{cell_size}m"
            
            total_pixels = pixels_width * pixels_height
            pixels_width_str = f"{pixels_width:,}"
            pixels_height_str = f"{pixels_height:,}"
            total_pixels_str = f"{total_pixels:,}"
            large_size_threshold = 10000
            is_large = pixels_width > large_size_threshold or pixels_height > large_size_threshold
            
            if is_large:
                self.pixel_count_label.setText(
                    f"⚠️ Output Grid Pixels : {pixels_width_str} × {pixels_height_str} = {total_pixels_str} "
                    f"(LARGE DATASET!)"
                )
                self.pixel_count_label.setStyleSheet("font-weight: bold; padding: 5px; color: orange;")
            else:
                self.pixel_count_label.setText(
                    f"Output Grid Pixels : {pixels_width_str} × {pixels_height_str} = {total_pixels_str}"
                )
                self.pixel_count_label.setStyleSheet("font-weight: bold; padding: 5px;")
        except Exception:
            self.pixel_count_label.setText("Pixels: --")
            self.pixel_count_label.setStyleSheet("font-weight: bold; padding: 5px;")
            
    def on_selection_changed(self, xmin, ymin, xmax, ymax):
        """Handle selection change from map (during dragging)."""
        if xmin == 0 and ymin == 0 and xmax == 0 and ymax == 0:
            # Selection cleared
            self.west_edit.clear()
            self.south_edit.clear()
            self.east_edit.clear()
            self.north_edit.clear()
            self.pixel_count_label.setText("Pixels: --")
            self.selected_bbox = None
            self.download_btn.setEnabled(False)
        else:
            # Show real-time values while selecting (without updating map to avoid recursion)
            self.update_coordinate_display(xmin, ymin, xmax, ymax, update_map=False)
            # Button state will be updated by update_coordinate_display (which calls check_and_update_download_button)
            
    def _snap_bounds_to_cell_size(self, xmin, ymin, xmax, ymax):
        """Snap bounding box to align with cell size grid (in active map CRS).
        
        Returns:
            tuple: (snapped_xmin, snapped_ymin, snapped_xmax, snapped_ymax, was_adjusted)
        """
        import math
        ds = self.data_sources.get(self.current_data_source, {})
        
        if self._uses_web_mercator_map():
            pixel_size = self._get_selected_cell_size_meters()
        elif self._is_gmrt_source():
            # Snap in degrees using approximate meters-per-degree conversion
            pixel_size = self._get_selected_cell_size_meters() / 111320.0
        elif ds.get("native_resolution_only") or ds.get("configurable_cell_size_degrees"):
            pixel_size = self._get_degrees_per_pixel_for_output()
        else:
            pixel_size = self._get_selected_cell_size_meters() / 111320.0
        
        snapped_xmin = math.floor(xmin / pixel_size) * pixel_size
        snapped_xmax = math.ceil(xmax / pixel_size) * pixel_size
        snapped_ymin = math.floor(ymin / pixel_size) * pixel_size
        snapped_ymax = math.ceil(ymax / pixel_size) * pixel_size
        
        was_adjusted = (
            abs(snapped_xmin - xmin) > 1e-10 or
            abs(snapped_ymin - ymin) > 1e-10 or
            abs(snapped_xmax - xmax) > 1e-10 or
            abs(snapped_ymax - ymax) > 1e-10
        )
        
        return snapped_xmin, snapped_ymin, snapped_xmax, snapped_ymax, was_adjusted
    
    def on_selection_completed(self, xmin, ymin, xmax, ymax):
        """Handle selection completion (when mouse is released) - zoom to selection."""
        if xmin != 0 or ymin != 0 or xmax != 0 or ymax != 0:
            # Snap bounds to cell size grid (map CRS)
            snapped_xmin, snapped_ymin, snapped_xmax, snapped_ymax, was_adjusted = self._snap_bounds_to_cell_size(xmin, ymin, xmax, ymax)
            
            if was_adjusted:
                before = self._map_bbox_to_display((xmin, ymin, xmax, ymax))
                after = self._map_bbox_to_display((snapped_xmin, snapped_ymin, snapped_xmax, snapped_ymax))
                self.log_message(
                    f"Selection bounds adjusted to align with cell size grid: "
                    f"({before[0]:.6f}, {before[1]:.6f}, {before[2]:.6f}, {before[3]:.6f}) → "
                    f"({after[0]:.6f}, {after[1]:.6f}, {after[2]:.6f}, {after[3]:.6f})"
                )
            
            # Store the selected bbox for download
            self.selected_bbox = (snapped_xmin, snapped_ymin, snapped_xmax, snapped_ymax)
            # Temporarily disconnect the selection changed signal to prevent clearing
            self.map_widget.selectionChanged.disconnect()
            self.zoom_to_selection(snapped_xmin, snapped_ymin, snapped_xmax, snapped_ymax)
            # Reconnect the signal
            self.map_widget.selectionChanged.connect(self.on_selection_changed)
            # Set the final bounds in the text fields after zoom
            self.update_coordinate_display(snapped_xmin, snapped_ymin, snapped_xmax, snapped_ymax)
            # Button state will be updated by update_coordinate_display (which calls check_and_update_download_button)
            
    def _get_service_extent(self):
        """Return the active service extent from the map widget or main window."""
        if self.map_widget and getattr(self.map_widget, "service_extent", None):
            return self.map_widget.service_extent
        return self.service_extent

    @staticmethod
    def _extents_equal(extent_a, extent_b, tol=1e-5):
        return extents_equal(extent_a, extent_b, tol)

    @staticmethod
    def _clamp_extent_to_bounds(extent, bounds):
        """Clamp a view extent to valid geographic/service bounds."""
        return clamp_extent_to_bounds(extent, bounds)

    def zoom_to_selection(self, xmin, ymin, xmax, ymax, record_history=True):
        """Zoom map to the selected area."""
        if self.map_widget:
            # Get widget aspect ratio - ensure widget is properly sized
            widget_width = self.map_widget.width()
            widget_height = self.map_widget.height()
            
            # If widget isn't sized yet, wait a bit and try again
            if widget_width <= 0 or widget_height <= 0:
                self.log_message(f"Widget not sized yet ({widget_width}x{widget_height}), retrying zoom_to_selection in 200ms...")
                QTimer.singleShot(
                    200,
                    lambda x=xmin, y=ymin, X=xmax, Y=ymax, rh=record_history: self.zoom_to_selection(
                        x, y, X, Y, record_history=rh
                    ),
                )
                return

            service_extent = self._get_service_extent()

            # Calculate the selected area dimensions
            selection_width = xmax - xmin
            selection_height = ymax - ymin

            # Add 5% padding around the selection
            padding_x = selection_width * 0.05
            padding_y = selection_height * 0.05

            # Start with padded extent
            padded_xmin = xmin - padding_x
            padded_ymin = ymin - padding_y
            padded_xmax = xmax + padding_x
            padded_ymax = ymax + padding_y

            padded_width = padded_xmax - padded_xmin
            padded_height = padded_ymax - padded_ymin

            if widget_width > 0 and widget_height > 0:
                widget_aspect = widget_width / widget_height
                padded_aspect = padded_width / padded_height

                # Calculate center of padded area
                center_x = (padded_xmin + padded_xmax) / 2
                center_y = (padded_ymin + padded_ymax) / 2

                # Adjust extent to match widget aspect ratio while containing the padded selection
                if padded_aspect > widget_aspect:
                    # Padded area is wider than widget - use padded width, adjust height
                    new_width = padded_width
                    new_height = new_width / widget_aspect
                else:
                    # Padded area is taller than widget - use padded height, adjust width
                    new_height = padded_height
                    new_width = new_height * widget_aspect

                # Create new extent centered on the padded selection
                new_extent = (
                    center_x - new_width / 2,
                    center_y - new_height / 2,
                    center_x + new_width / 2,
                    center_y + new_height / 2
                )
            else:
                new_extent = (padded_xmin, padded_ymin, padded_xmax, padded_ymax)

            new_extent = self._clamp_extent_to_bounds(new_extent, service_extent)
            # Set the extent FIRST, then store the selection bbox
            # This ensures the selection bbox is stored with the correct extent context
            self.map_widget.extent = new_extent
            self.map_widget._requested_extent = new_extent  # Also update _requested_extent for accurate coordinate conversion
            # Store the selected bbox in world coordinates for drawing (original selection, no modifications)
            # This is what will be shown in the yellow/green box and used for download
            self.map_widget.selected_bbox_world = (xmin, ymin, xmax, ymax)
            self.selected_bbox = (xmin, ymin, xmax, ymax)
            # Ensure service_extent is preserved
            if not hasattr(self.map_widget, 'service_extent') or self.map_widget.service_extent is None:
                self.map_widget.service_extent = self.service_extent
            # Don't clear selection - keep it visible
            self.map_widget.load_map()

            if not self._zoom_history:
                self._seed_zoom_history()
            elif record_history and not self._zoom_history_suppress:
                self._append_zoom_history()
            else:
                self._update_zoom_nav_buttons()
            
    def start_download(self):
        """Start downloading the selected area."""
        bbox = None
        if hasattr(self, 'selected_bbox') and self.selected_bbox:
            bbox = self.selected_bbox
        elif self.map_widget:
            bbox = self.map_widget.get_selection_bbox()
        else:
            try:
                west_text = self.west_edit.text().strip()
                south_text = self.south_edit.text().strip()
                east_text = self.east_edit.text().strip()
                north_text = self.north_edit.text().strip()
                if west_text and south_text and east_text and north_text:
                    bbox = self._display_bbox_to_map((
                        float(west_text), float(south_text), float(east_text), float(north_text)
                    ))
            except ValueError:
                QMessageBox.warning(self, "Invalid Input", "Please enter valid coordinates.")
                return
        if not bbox:
            QMessageBox.warning(self, "No Selection", "Please select an area on the map.")
            return
        output_crs = "EPSG:4326"
        ds = self.data_sources.get(self.current_data_source, {})
        native_only = ds.get("native_resolution_only", False)
        is_gmrt = self._is_gmrt_source()
        if is_gmrt:
            west, south, east, north = clamp_gmrt_bbox(*bbox)
            bbox = (west, south, east, north)
            cell_size = self._get_selected_cell_size_meters()
            pixels_width, pixels_height = estimate_gmrt_pixels(west, south, east, north, cell_size)
            cell_size_for_filename = int(cell_size) if float(cell_size).is_integer() else cell_size
            lon_span = east - west
            lat_span = north - south
            if self.tile_download_checkbox.isChecked() and needs_tiling_for_spans(lon_span, lat_span, cell_size):
                n_tiles = len(generate_gmrt_tiles(west, east, south, north, cell_size))
                if n_tiles > MAX_TILES_PER_DOWNLOAD:
                    QMessageBox.warning(
                        self,
                        "Too Many GMRT Tiles",
                        f"This area would require {n_tiles} tiles at {cell_size:g} m "
                        f"(limit {MAX_TILES_PER_DOWNLOAD}).\n\n"
                        f"Select a smaller area or a coarser cell size.",
                    )
                    return
        elif native_only:
            # Bbox is (west, south, east, north) in 4326
            bbox_4326 = bbox
            lon_min, lat_min, lon_max, lat_max = bbox
            pixel_size_degrees = self._get_degrees_per_pixel_for_output()
            pixels_width = int((lon_max - lon_min) / pixel_size_degrees)
            pixels_height = int((lat_max - lat_min) / abs(pixel_size_degrees))
            if ds.get("configurable_cell_size_degrees"):
                cell_size_for_filename = self._format_cell_size_degrees(pixel_size_degrees).replace(".", "p")
            else:
                cell_size_for_filename = "native"
        else:
            # Meter-based Web Mercator source; selection is already EPSG:3857
            cell_size = self._get_selected_cell_size_meters()
            xmin, ymin, xmax, ymax = bbox
            width_meters = xmax - xmin
            height_meters = ymax - ymin
            pixels_width = int(width_meters / cell_size)
            pixels_height = int(height_meters / cell_size)
            cell_size_for_filename = int(cell_size) if float(cell_size).is_integer() else cell_size
        
        large_size_threshold = 10000
        if pixels_width > large_size_threshold or pixels_height > large_size_threshold:
            total_pixels = pixels_width * pixels_height
            msg = QMessageBox(self)
            msg.setIcon(QMessageBox.Icon.Warning)
            msg.setWindowTitle("Large Dataset Warning")
            msg.setText(
                f"You are about to download a very large dataset.\n\n"
                f"Requested size: {pixels_width:,} × {pixels_height:,} pixels\n"
                f"Total pixels: {total_pixels:,}\n\n"
                f"This download may take a significant amount of time and disk space.\n"
                f"{'Tiled download is enabled and will break this into multiple requests.' if self.tile_download_checkbox.isChecked() else 'Consider enabling Tile Download for better reliability.'}\n\n"
                f"Do you want to continue?"
            )
            msg.setStandardButtons(QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel)
            msg.setDefaultButton(QMessageBox.StandardButton.Cancel)
            result = msg.exec()
            if result == QMessageBox.StandardButton.Cancel:
                return
        
        current_time = datetime.now()
        date_time_str = current_time.strftime("%Y-%m-%d_%H-%M-%S")
        
        # Build list of requested outputs for multi-type GEBCO sources
        output_requests = self._collect_output_type_requests()
        tid_url = self._get_output_type_tid_url() if output_requests else None
        if self._shows_output_data_types() and not output_requests:
            QMessageBox.warning(
                self,
                "No Output Selected",
                "Select at least one output: Combined Bathymetry && Land, Bathymetry Only, "
                "Land Only, Direct Measurements Only, or Direct && Unknown Measurement Only.",
            )
            return
        
        filename_prefix = self._get_download_filename_prefix()
        # Resolve output path(s) for multi-output sources
        if native_only and self._shows_output_data_types() and output_requests:
            if len(output_requests) > 1:
                if not self.output_directory or not os.path.isdir(self.output_directory):
                    QMessageBox.warning(self, "Output Directory Required", "Select an output directory when saving multiple grids.")
                    return
                out_dir = self.output_directory
                resolved = []
                for mode, _ in output_requests:
                    mode_name = self._output_mode_filename_suffix(mode)
                    fn = f"{filename_prefix}_{mode_name}_{date_time_str}.tif"
                    resolved.append((mode, os.path.join(out_dir, fn)))
                output_requests = resolved
            else:
                mode = output_requests[0][0]
                mode_name = self._output_mode_filename_suffix(mode)
                default_name = f"{filename_prefix}_{mode_name}_{date_time_str}.tif"
                if self.output_directory and os.path.isdir(self.output_directory):
                    output_path = os.path.join(self.output_directory, default_name)
                else:
                    output_path, _ = QFileDialog.getSaveFileName(self, "Save GeoTIFF", default_name, "GeoTIFF Files (*.tif *.tiff);;All Files (*)")
                    if not output_path:
                        return
                output_requests = [(mode, output_path)]
        elif native_only and not self._shows_output_data_types():
            prefix = ds.get("download_filename_prefix", "bathymetry")
            if ds.get("configurable_cell_size_degrees"):
                default_filename = f"{prefix}_{cell_size_for_filename}deg_{date_time_str}.tif"
            else:
                default_filename = f"{prefix}_{date_time_str}.tif"
            if self.output_directory and os.path.isdir(self.output_directory):
                output_path = os.path.join(self.output_directory, default_filename)
            else:
                output_path, _ = QFileDialog.getSaveFileName(self, "Save GeoTIFF", default_filename, "GeoTIFF Files (*.tif *.tiff);;All Files (*)")
                if not output_path:
                    return
            output_requests = [("combined", output_path)]
        elif native_only:
            default_filename = f"{filename_prefix}_{date_time_str}.tif"
            if self.output_directory and os.path.isdir(self.output_directory):
                output_path = os.path.join(self.output_directory, default_filename)
            else:
                output_path, _ = QFileDialog.getSaveFileName(self, "Save GeoTIFF", default_filename, "GeoTIFF Files (*.tif *.tiff);;All Files (*)")
                if not output_path:
                    return
            output_requests = [("combined", output_path)]
        else:
            default_filename = f"{filename_prefix}_{cell_size_for_filename}m_{date_time_str}.tif"
            if self.output_directory and os.path.isdir(self.output_directory):
                output_path = os.path.join(self.output_directory, default_filename)
            else:
                output_path, _ = QFileDialog.getSaveFileName(self, "Save GeoTIFF", default_filename, "GeoTIFF Files (*.tif *.tiff);;All Files (*)")
                if not output_path:
                    return
            output_requests = [("combined", output_path)]
        
        # Disable download button
        self.download_btn.setEnabled(False)
        self.progress_bar.setValue(0)
        self.status_label.setText("Starting download...")
        
        # Get tile download setting
        use_tile_download = self.tile_download_checkbox.isChecked()
        
        max_size = 14000
        ignore_source_nodata = ds.get("ignore_source_nodata", False)
        if is_gmrt:
            self.downloader = GMRTDownloader(
                bbox_4326=bbox,
                output_path=output_path,
                gmrt_layer=ds.get("gmrt_layer", "topo"),
                mresolution=cell_size,
                use_tile_download=use_tile_download,
            )
        elif native_only:
            self.downloader = BathymetryDownloader(
                self.base_url,
                bbox_4326,
                None,  # output_path unused when output_requests provided
                output_crs,
                pixel_size=None,
                max_size=max_size,
                use_tile_download=use_tile_download,
                bbox_in_4326=True,
                pixel_size_degrees=pixel_size_degrees,
                tid_url=tid_url,
                output_requests=output_requests,
                ignore_source_nodata=ignore_source_nodata,
            )
        else:
            self.downloader = BathymetryDownloader(
                self.base_url,
                bbox,
                output_path,
                output_crs,
                pixel_size=cell_size,
                max_size=max_size,
                use_tile_download=use_tile_download,
                ignore_source_nodata=ignore_source_nodata,
            )
        self.downloader.progress.connect(self.progress_bar.setValue)
        self.downloader.status.connect(self.on_status_update)
        self.downloader.finished.connect(self.on_download_finished)
        self.downloader.error.connect(self.on_download_error)
        self.downloader.start()
        
    def on_status_update(self, message):
        """Handle status update from downloader."""
        self.status_label.setText(message)
        self.log_message(message)
        
    def on_download_finished(self, file_path):
        """Handle download completion. file_path may be newline-separated for multiple files."""
        paths = [p.strip() for p in file_path.splitlines() if p.strip()]
        if not paths:
            paths = [file_path]
        display = "\n".join(paths)
        self.status_label.setText(f"Download complete: {paths[0]}" if len(paths) == 1 else f"Download complete: {len(paths)} files")
        self.log_message(f"✓ Download complete: {display}")
        self.download_btn.setEnabled(True)
        font = self.download_btn.font()
        font.setBold(False)
        self.download_btn.setFont(font)

        path_to_load = paths[0]
        split_cb = getattr(self, "split_topo_depths_checkbox", None)
        if split_cb is not None and split_cb.isChecked() and split_cb.isEnabled():
            try:
                from sat_planner.gmrt_split import split_topo_bathy

                result = split_topo_bathy(path_to_load)
                if result.bathy_path and os.path.isfile(result.bathy_path):
                    path_to_load = result.bathy_path
                    try:
                        if os.path.isfile(paths[0]) and os.path.abspath(paths[0]) != os.path.abspath(path_to_load):
                            os.remove(paths[0])
                    except OSError:
                        pass
                    self.log_message(f"Split complete. Bathymetry file: {path_to_load}")
                else:
                    self.log_message("Split produced no bathymetry file; keeping combined grid.")
            except Exception as exc:
                self.log_message(f"Split failed ({exc}); keeping combined grid.")

        QMessageBox.information(self, "Success", f"GeoTIFF saved to:\n{display}")
        self.geotiff_downloaded.emit(path_to_load)
        
    def on_download_error(self, error_message):
        """Handle download error."""
        self.status_label.setText(f"Error: {error_message}")
        self.log_message(f"✗ Error: {error_message}")
        self.download_btn.setEnabled(True)
        # Remove bold formatting after download error
        font = self.download_btn.font()
        font.setBold(False)
        self.download_btn.setFont(font)
        
        # Check if it's a connection error and show helpful message
        if "connection" in error_message.lower() or "timeout" in error_message.lower() or "network" in error_message.lower() or "rest endpoint" in error_message.lower():
            msg = QMessageBox(self)
            msg.setIcon(QMessageBox.Icon.Warning)
            msg.setWindowTitle("Connection Error")
            msg.setText(
                f"Unable to connect to the REST endpoint.\n\n"
                f"{error_message}\n\n"
                f"If this problem persists, please:\n"
                f"1. Check for a new version at: https://github.com/seamapper/GEBCO_Downloader\n"
                f"2. Contact: pjohnson@ccom.unh.edu"
            )
            msg.setStandardButtons(QMessageBox.StandardButton.Ok)
            msg.exec()
        QMessageBox.critical(self, "Download Error", error_message)
        
    def log_message(self, message, bold=False, color=None):
        """Add message to log. If bold is True, the message is shown in bold (HTML). If color is set, wrap in span (e.g. 'orange')."""
        if bold:
            message = f"<b>{message}</b>"
        if color:
            message = f'<span style="color: {color};">{message}</span>'
        self.log_text.append(message)
        # Auto-scroll to bottom
        scrollbar = self.log_text.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())
        
    def resizeEvent(self, event):
        """Handle window resize event - refresh map display."""
        super().resizeEvent(event)
        # Refresh map when window is resized (with a small delay to avoid multiple refreshes)
        if self.map_widget and self.map_widget.map_loaded:
            # Use a timer to debounce rapid resize events
            if not hasattr(self, '_resize_timer'):
                self._resize_timer = QTimer()
                self._resize_timer.setSingleShot(True)
                self._resize_timer.timeout.connect(self._refresh_map_on_resize)
            
            # Restart timer - will trigger refresh 300ms after resize stops
            self._resize_timer.stop()
            self._resize_timer.start(300)
    
    def _refresh_map_on_resize(self):
        """Refresh map display after window resize."""
        if self.map_widget and self.map_widget.map_loaded:
            # Ensure widget size is updated before calculating new extent
            # Process events to ensure Qt has updated the widget size
            QApplication.processEvents()
            
            # Verify widget size is valid before proceeding
            widget_width = self.map_widget.width()
            widget_height = self.map_widget.height()
            
            if widget_width <= 0 or widget_height <= 0:
                # Widget not sized yet, skip this resize
                return
            
            # If there's a selected area, zoom to it to maintain constant visual size
            # This treats the resize as if the user made a new selection with the same bounds
            # The zoom_to_selection function will recalculate the extent based on the new widget size
            # making the selection box appear the same visual size
            if self.map_widget.selected_bbox_world:
                xmin, ymin, xmax, ymax = self.map_widget.selected_bbox_world
                # Zoom to the selection - this will recalculate the extent based on new widget size
                # making the selection box appear the same visual size
                self.zoom_to_selection(xmin, ymin, xmax, ymax, record_history=False)
            else:
                # No selection - just reload the map with current extent
                self.map_widget.load_map()
    
    def closeEvent(self, event):
        """Handle window close event."""
        if self.downloader and self.downloader.isRunning():
            reply = QMessageBox.question(
                self,
                "Download in Progress",
                "A download is in progress. Do you want to cancel it and exit?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
            )
            if reply == QMessageBox.StandardButton.Yes:
                self.downloader.cancel()
                self.downloader.wait(3000)  # Wait up to 3 seconds
            else:
                event.ignore()
                return
        event.accept()
    
    def load_config(self):
        """Load configuration from JSON file."""
        try:
            if os.path.exists(self.config_file):
                with open(self.config_file, 'r') as f:
                    config = json.load(f)
                    self.output_directory = config.get('output_directory')
                    # Update edit field if it exists (it should after init_ui)
                    if hasattr(self, 'output_dir_edit'):
                        if self.output_directory and os.path.isdir(self.output_directory):
                            self.output_dir_edit.setText(self.output_directory)
                        else:
                            self.output_directory = None
                            self.output_dir_edit.clear()
                    elif not (self.output_directory and os.path.isdir(self.output_directory)):
                        self.output_directory = None
        except Exception as e:
            # If config file is corrupted or can't be read, just use defaults
            self.output_directory = None
            if hasattr(self, 'output_dir_edit'):
                self.output_dir_edit.clear()
    
    def save_config(self):
        """Save configuration to JSON file."""
        try:
            config = {
                'output_directory': self.output_directory
            }
            with open(self.config_file, 'w') as f:
                json.dump(config, f, indent=2)
        except Exception as e:
            # If we can't save config, just continue - it's not critical
            pass
    
    def select_output_directory(self):
        """Open dialog to select output directory."""
        # Start with current directory or saved directory
        start_dir = self.output_directory if self.output_directory and os.path.isdir(self.output_directory) else os.getcwd()
        
        directory = QFileDialog.getExistingDirectory(
            self,
            "Select Output Directory",
            start_dir,
            QFileDialog.Option.ShowDirsOnly | QFileDialog.Option.DontResolveSymlinks
        )
        
        if directory:
            self.output_directory = directory
            self.output_dir_edit.setText(directory)
            self.save_config()  # Save to config file
    
    def on_map_first_loaded(self):
        """Handle first successful map load - show instructions and set default bounds."""
        # If no selection exists yet, set default to REST endpoint service extent bounds
        # CRITICAL: Use the REST endpoint extent (service_extent) for the box, NOT the map's displayed extent
        # The map might show a slightly different area due to rounding or basemap coverage,
        # but the box should show the exact bathymetry data bounds from the REST endpoint
        if self.map_widget and self.map_widget.selected_bbox_world is None:
            if self.service_extent:
                # Set default selection to REST endpoint service extent bounds (exact bathymetry data bounds)
                # This is the correct extent from the REST endpoint, not the map's displayed extent
                self.map_widget.selected_bbox_world = self.service_extent
                self.map_widget.set_selection_validity(True)
                self.selected_bbox = self.service_extent
                # Ensure service_extent is stored in map widget (this is the REST endpoint extent)
                self.map_widget.service_extent = self.service_extent
                
                # CRITICAL: Zoom to the REST endpoint bounds using zoom_to_selection
                # This ensures the map extent is recalculated with padding and the box is positioned correctly
                # This mimics what happens when the user hits return in a coordinate field
                QTimer.singleShot(300, lambda: self.zoom_to_selection(*self.service_extent, record_history=False))
                self.log_message("Default selection set to service extent bounds, will zoom to dataset bounds")
        # Show map interaction instructions once when map first loads (in orange)
        self.log_message(
            "Map tips — Pan: middle mouse button + drag. Zoom: mouse wheel. Select Area of Interest: left mouse button + drag.",
            color="orange"
        )
    
    def _zoom_to_service_extent(self):
        """Zoom to service extent bounds - helper method for delayed zoom."""
        if self.map_widget and self.service_extent:
            # Ensure the selected bbox is set to service extent (this is the dataset bounds)
            self.map_widget.selected_bbox_world = self.service_extent
            self.map_widget.set_selection_validity(True)
            self.selected_bbox = self.service_extent
            # Ensure service extent is stored in map widget for color distinction
            self.map_widget.service_extent = self.service_extent
            
            # Zoom to the service extent - this will reload the map with the correct extent
            self.zoom_to_selection(*self.service_extent, record_history=False)
            # After zoom completes, the map will reload and the box should be visible
            # The box will be repainted in paintEvent when the new map loads
        
        instructions = [
            "",
            "=" * 60,
            "Map loaded successfully!",
            "=" * 60,
            "",
            "To select an area:",
            "  1. Click and drag with the left mouse button on the map",
            "  2. The selected area will be shown with a purple dashed box",
            "  3. You can also manually enter coordinates in the West/South/East/North (GCS) fields",
            "",
            "To download the selected area:",
            "  1. Select an area on the map (or enter coordinates)",
            "  2. For GEBCO 2025, choose Combined Bathymetry && Land, Bathymetry Only, Land Only, Direct Measurements Only, or Direct && Unknown Measurement Only if needed",
            "  3. Click 'Download Selected Area' button",
            "  4. Choose a filename and location (defaults to selected directory)",
            "",
            "Map controls:",
            "  - Mouse wheel: Zoom in/out (centered on window)",
            "  - Middle mouse button + drag: Pan the map",
            "  - Left mouse button + drag: Select area",
            "",
            "=" * 60
        ]
        
        for line in instructions:
            self.log_message(line)
    
    def _bboxes_overlap(self, bbox1, bbox2):
        """Check if two bounding boxes overlap."""
        return bboxes_overlap(bbox1, bbox2)
    
    def on_data_source_changed(self, data_source_name):
        """Handle data source selection change."""
        if data_source_name not in self.data_sources:
            return
        
        # Get new data source extent
        new_service_extent = self.data_sources[data_source_name]["default_extent"]
        
        # Preserve selection when staying in the same map CRS (GCS↔GCS or WGOM↔WGOM).
        # Clear when crossing between Web Mercator WGOM and geographic sources.
        saved_selection = None
        switching_to_wgom = self._uses_web_mercator_map(data_source_name)
        switching_from_wgom = self._uses_web_mercator_map(self.current_data_source)
        if switching_to_wgom != switching_from_wgom:
            self.selected_bbox = None
            if self.map_widget:
                self.map_widget.selected_bbox_world = None
        elif hasattr(self, 'selected_bbox') and self.selected_bbox:
            saved_selection = self.selected_bbox
        
        # Update current data source
        self._reset_zoom_history()
        self.current_data_source = data_source_name
        self.base_url = self.data_sources[data_source_name]["url"]
        self.service_extent = new_service_extent
        
        # Set flag to force highest resolution when cell size options are updated
        self._data_source_changing = True
        
        # Update map widget settings if it exists
        if self.map_widget:
            # Drop in-flight preview loaders before swapping source URLs
            self.map_widget._stop_all_loaders()
            # Update raster functions
            new_raster_function = self.data_sources[data_source_name]["bathymetry_raster_function"]
            new_hillshade_raster_function = self.data_sources[data_source_name]["hillshade_raster_function"]
            self.map_widget.raster_function = new_raster_function
            self.map_widget.hillshade_raster_function = new_hillshade_raster_function
            self.map_widget.base_url = self.base_url
            self.map_widget.display_url = self.data_sources[data_source_name].get("display_url")
            self.map_widget.land_display_url = self.data_sources[data_source_name].get("land_display_url")
            self.map_widget.show_basemap = bool(self.data_sources[data_source_name].get("show_basemap", False))
            self.map_widget.show_hillshade = bool(self.data_sources[data_source_name].get("show_hillshade", False))
            # Multiply blend with hillshade (same as CCOM Downloader)
            self.map_widget.use_blend = self.map_widget.show_hillshade
            self.map_widget.bbox_sr = self._bbox_sr_for_data_source(data_source_name)
            self.map_widget.preview_url = self.data_sources[data_source_name].get("preview_url")
            self.map_widget.gmrt_mask = (
                hasattr(self, "gmrt_mask_checkbox") and self.gmrt_mask_checkbox.isChecked()
            )
            # Update service extent in map widget
            self.map_widget.service_extent = self.service_extent
            # Don't update pixel sizes here - they will be updated in on_service_info_loaded after the new service loads
            # This ensures we get the correct pixel sizes for the new data source
        
        # Show/hide output option controls based on data source
        self._update_output_options_visibility()
        if self.data_sources[data_source_name].get("configurable_cell_size_degrees"):
            self._set_cell_size_degrees_default(force=True)
        
        # Update attribution text
        self._update_attribution()
        
        # Must set pending selection BEFORE load_service_info — GMRT loads
        # synchronously and would otherwise reset the AOI to the full extent.
        self._pending_selection = saved_selection
        # Always force a map reload on source change (GMRT sources share the same
        # GridServer/ImageServer URLs, so URL-diff checks would skip the refresh).
        self._force_map_reload = True

        # Observed Only → show hi-res mask; Topo-Bathy → clear it
        if self._is_gmrt_source(data_source_name) and hasattr(self, "gmrt_mask_checkbox"):
            observed_only = (
                self.data_sources[data_source_name].get("gmrt_layer") == "topo-mask"
            )
            self.gmrt_mask_checkbox.blockSignals(True)
            self.gmrt_mask_checkbox.setChecked(observed_only)
            self.gmrt_mask_checkbox.blockSignals(False)
            if self.map_widget:
                self.map_widget.gmrt_mask = observed_only
        
        # Reload service info (this will update extent and reload map)
        self.load_service_info()
    
    def _bbox_sr_for_data_source(self, data_source_name):
        """Return bboxSR for ImageServer export, or None for native service CRS."""
        return bbox_sr_for_data_source(self.data_sources, data_source_name)

    def _update_output_options_visibility(self):
        """Show/hide output controls based on the active data source."""
        ds = self.data_sources.get(self.current_data_source, {})
        if hasattr(self, "output_data_types_group"):
            self.output_data_types_group.setVisible(ds.get("show_output_data_types", True))
        if hasattr(self, "cell_size_meters_container"):
            self.cell_size_meters_container.setVisible(self._uses_meter_cell_size())
        if hasattr(self, "cell_size_degrees_container"):
            self.cell_size_degrees_container.setVisible(ds.get("configurable_cell_size_degrees", False))
        if hasattr(self, "bathymetry_only_notice_label"):
            self.bathymetry_only_notice_label.setVisible(ds.get("show_bathymetry_only_notice", False))
        if hasattr(self, "gmrt_mask_checkbox"):
            self.gmrt_mask_checkbox.setVisible(self._is_gmrt_source())
        if hasattr(self, "split_topo_depths_checkbox"):
            is_topo_bathy = self._is_gmrt_source() or self.current_data_source.startswith("GEBCO")
            self.split_topo_depths_checkbox.setVisible(is_topo_bathy)
            self.split_topo_depths_checkbox.setEnabled(is_topo_bathy)
            if is_topo_bathy:
                self.split_topo_depths_checkbox.setChecked(True)

    def on_gmrt_mask_toggled(self, state):
        """Reload GMRT preview when the hi-res mask checkbox changes."""
        if not self.map_widget or not self._is_gmrt_source():
            return
        checked = (state == Qt.CheckState.Checked.value or state == 2)
        self.map_widget.gmrt_mask = checked
        if not getattr(self.map_widget, "_loading", False):
            self.map_widget.load_map()
    
    def _update_attribution(self):
        """Update the attribution text based on the current data source."""
        if not hasattr(self, 'attribution_label'):
            return
        
        ds = self.data_sources.get(self.current_data_source, {})
        attribution_text = ds.get("attribution", "")
        attribution_url = ds.get("attribution_url", "")
        
        if attribution_text:
            self.attribution_label.setText(attribution_text)
            self.attribution_label.setToolTip(f"Click to open: {attribution_url}")
            self._current_attribution_url = attribution_url
            self.attribution_label.setVisible(True)
        else:
            self.attribution_label.setVisible(False)
            self._current_attribution_url = None
    
    def _open_attribution_url(self):
        """Open the attribution URL in the default web browser."""
        if hasattr(self, '_current_attribution_url') and self._current_attribution_url:
            QDesktopServices.openUrl(QUrl(self._current_attribution_url))
    
    def _reload_map_with_selection(self):
        """Reload map and restore selection if it exists."""
        # If there's a pending selection, let zoom_to_selection handle loading the map
        # Otherwise, load the map with current extent
        if self.map_widget:
            if hasattr(self, '_pending_selection') and self._pending_selection:
                # Wait a moment for map widget to be ready, then restore selection
                # zoom_to_selection will load the map with the correct extent
                QTimer.singleShot(100, lambda: self._restore_selection())
            else:
                # No selection - just reload the map with current extent
                self.map_widget.load_map()
    
    def _restore_selection(self):
        """Restore a previously saved selection and zoom to it."""
        if hasattr(self, '_pending_selection') and self._pending_selection:
            bbox = self._pending_selection
            self.selected_bbox = bbox
            
            # Update map widget's selected_bbox_world so the selection box is drawn
            if self.map_widget:
                # Allow reload even if a previous load is still marked in-progress
                if getattr(self.map_widget, "_loading", False):
                    self.map_widget._stop_all_loaders()
                    self.map_widget._loading = False
                self.map_widget.selected_bbox_world = bbox
                self.map_widget.set_selection_validity(True)
            
            # Zoom to the selection to maintain visual size
            self.zoom_to_selection(bbox[0], bbox[1], bbox[2], bbox[3], record_history=False)
            
            # Update coordinate displays (without updating map since zoom_to_selection already does)
            self.update_coordinate_display(bbox[0], bbox[1], bbox[2], bbox[3], update_map=False)
            
            # Clear pending selection
            self._pending_selection = None


"""
Map widget for displaying bathymetry data and selecting areas of interest.

License: BSD 3-Clause License
Copyright (c) 2025, Center for Coastal and Ocean Mapping, University of New Hampshire
All rights reserved.

See LICENSE file for full license text.
"""
from PyQt6.QtWidgets import QWidget, QLabel, QVBoxLayout, QHBoxLayout, QPushButton
from PyQt6.QtCore import Qt, QRect, QPoint, pyqtSignal, QSize
from PyQt6.QtGui import QPainter, QPen, QBrush, QColor, QPixmap, QImage, QPalette
import math
import pyproj
from urllib.parse import urlencode

from .map_loaders import BasemapLoader, MapServerLoader, MapTileLoader


def _format_rest_url(url, params):
    """Build a full REST request URL for activity log display."""
    return f"{url}?{urlencode(params)}"


class MapWidget(QWidget):
    """Interactive map widget for displaying bathymetry and selecting areas."""
    
    selectionChanged = pyqtSignal(float, float, float, float)  # xmin, ymin, xmax, ymax
    selectionCompleted = pyqtSignal(float, float, float, float)  # xmin, ymin, xmax, ymax - emitted when selection is finished
    mapFirstLoaded = pyqtSignal()  # Emitted when map is successfully loaded for the first time
    statusMessage = pyqtSignal(str)  # Emit status/log messages
    userViewChanged = pyqtSignal()  # Emitted after user pan/wheel changes the map view
    
    def set_selection_validity(self, is_valid):
        """Set whether the current selection is within size limits."""
        if self.selection_is_valid != is_valid:
            self.selection_is_valid = is_valid
            self.update()  # Trigger repaint to update color
    
    def __init__(self, base_url, initial_extent, parent=None, raster_function="Shaded Relief - Haxby - MD Hillshade 2", show_basemap=True, show_hillshade=True, use_blend=False, hillshade_raster_function="Multidirectional Hillshade 3x", display_url=None, land_display_url=None):
        super().__init__(parent)
        self.base_url = base_url
        self.display_url = display_url  # When set, map is drawn from this MapServer (e.g. GEBCO Haxby) instead of ImageServer layers
        self.land_display_url = land_display_url  # Land layer MapServer (e.g. GEBCO Land Grey) shown as basemap
        self.extent = initial_extent  # (west, south, east, north) in GCS (4326)
        self.current_pixmap = QPixmap()
        self.basemap_pixmap = QPixmap()
        self.hillshade_pixmap = QPixmap()
        self.show_basemap = show_basemap
        self.show_hillshade = show_hillshade
        self.show_legend = True  # Legend visibility (on by default)
        self.show_aoi = True  # Area of Interest (selection rectangle) visibility (on by default)
        self.use_blend = use_blend  # Use Multiply blend mode for top layer
        self.bathymetry_opacity = 1.0  # Opacity for bathymetry layer (0.0 to 1.0) - default 100%
        self.selection_start = None
        self.selection_end = None
        self.is_selecting = False
        self.is_panning = False
        self.pan_start = None
        self.pan_origin = None  # Track original pan start position for drawing pan line
        self.pan_end = None  # Track current pan position for drawing pan line
        self.raster_function = raster_function
        self.hillshade_raster_function = hillshade_raster_function  # Raster function for hillshade layer
        self.pixel_size_x = None  # Pixel size in X direction from service (meters)
        self.pixel_size_y = None  # Pixel size in Y direction from service (meters)
        self.map_loaded = False
        self._first_load_complete = False  # Track if first load has completed
        self._loading = False  # Flag to prevent multiple simultaneous loads
        self._active_loaders = []  # Track active loaders
        self._orphan_loaders = []  # Keep abandoned QThreads alive until they finish
        self._load_timer = None  # Timer for debouncing zoom operations
        self.selected_bbox_world = None  # (west, south, east, north) in GCS (4326)
        self.selection_is_valid = True  # Track if selection is within size limits (True = valid/green, False = too large/red)
        self.service_extent = None  # Store service extent (for reference)
        self._extent_locked = False  # Flag to prevent extent changes during resize
        self._original_pixmap_size = None  # Store original pixmap size before scaling for coordinate conversion
        self._scaled_pixmap_size = None  # Store scaled pixmap size (what's actually drawn)
        self._requested_extent = initial_extent  # (west, south, east, north) GCS
        self.bbox_sr = None  # bboxSR for ImageServer export (e.g. "4326")
        self.preview_url = None  # GMRT ImageServer URL when api == gmrt
        self.gmrt_mask = False  # GMRT high-res coverage mask overlay
        print(f"MapWidget initialized with raster function: {self.raster_function}, show_basemap: {self.show_basemap}, show_hillshade: {self.show_hillshade}, use_blend: {self.use_blend}")
        
        # Set a smaller minimum size to allow 60/40 split (60% of 1200 = 720px)
        self.setMinimumSize(600, 400)
        self.setMouseTracking(True)
        self.setAutoFillBackground(True)
        palette = self.palette()
        palette.setColor(self.backgroundRole(), QColor(0, 0, 0))
        palette.setColor(QPalette.ColorRole.Window, QColor(0, 0, 0))
        self.setPalette(palette)
        self.setStyleSheet("background-color: #000000;")
        
        # Don't load map immediately - wait for widget to be shown and sized
        # The load will be triggered by showEvent or when explicitly called
        
    def set_raster_function(self, raster_function):
        """Set the raster function for map display."""
        self.raster_function = raster_function
        self.load_map()
        
    def showEvent(self, event):
        """Handle widget being shown - trigger map load if not already loaded."""
        super().showEvent(event)
        print(f"MapWidget showEvent called, map_loaded={self.map_loaded}, size={self.width()}x{self.height()}")
        # Don't auto-load here - let MainWindow control when to load
        # This ensures the REST endpoint extent is available before loading
        # The map will be loaded explicitly by MainWindow after service info loads
        if not self.map_loaded:
            print("MapWidget showEvent: Waiting for MainWindow to trigger map load after service info loads")
            
    def _on_orphan_loader_finished(self):
        """Drop references to abandoned loaders once their threads exit."""
        sender = self.sender()
        if sender is None:
            return
        try:
            sender.finished.disconnect(self._on_orphan_loader_finished)
        except TypeError:
            pass
        if sender in self._orphan_loaders:
            self._orphan_loaders.remove(sender)

    def _stop_all_loaders(self):
        """Disconnect active loaders without blocking or destroying running QThreads.

        Using terminate()/wait() from the GUI thread can deadlock inside requests/SSL.
        Dropping the last Python reference while a QThread is still running crashes with
        'QThread: Destroyed while thread is still running'. Keep orphans until finished.
        """
        loaders_to_stop = []
        if getattr(self, "loader", None):
            loaders_to_stop.append(self.loader)
        if getattr(self, "basemap_loader", None):
            loaders_to_stop.append(self.basemap_loader)
        if getattr(self, "hillshade_loader", None):
            loaders_to_stop.append(self.hillshade_loader)

        for loader in loaders_to_stop:
            for signal_name in ("tileLoaded", "statusMessage"):
                signal = getattr(loader, signal_name, None)
                if signal is None:
                    continue
                try:
                    signal.disconnect()
                except TypeError:
                    pass
            # Detach completion handlers used by the current load cycle
            try:
                loader.finished.disconnect(self._check_all_loaders_finished)
            except TypeError:
                pass
            try:
                loader.finished.disconnect(self.on_loader_finished)
            except TypeError:
                pass

            if loader.isRunning():
                if loader not in self._orphan_loaders:
                    self._orphan_loaders.append(loader)
                try:
                    loader.finished.connect(
                        self._on_orphan_loader_finished,
                        Qt.ConnectionType.UniqueConnection,
                    )
                except TypeError:
                    pass

        self.loader = None
        self.basemap_loader = None
        self.hillshade_loader = None
        self._active_loaders = []
        self._loading = False

    def _check_all_loaders_finished(self):
        """Check if all loaders are finished and reset loading flag."""
        if not self._loading:
            return

        # Ignore stale finished signals from abandoned loaders
        active = {id(loader) for loader in self._active_loaders}
        sender = self.sender()
        if sender is not None and id(sender) not in active:
            return

        all_finished = True
        if getattr(self, "loader", None) and self.loader.isRunning():
            all_finished = False
        if getattr(self, "basemap_loader", None) and self.basemap_loader.isRunning():
            all_finished = False
        if getattr(self, "hillshade_loader", None) and self.hillshade_loader.isRunning():
            all_finished = False

        if all_finished:
            self._loading = False

    def _start_land_basemap_loader(self, requested_extent, size):
        """Load GEBCO land MapServer layer as an underlay when configured."""
        if not self.land_display_url:
            return
        print("Loading land basemap (GCS)...")
        self.basemap_loader = MapServerLoader(self.land_display_url, requested_extent, size, transparent=True, purpose="land display")
        self.basemap_loader.statusMessage.connect(self.statusMessage.emit)
        self.basemap_loader.tileLoaded.connect(lambda pixmap, *args: self.on_basemap_loaded(pixmap))
        self.basemap_loader.finished.connect(self._check_all_loaders_finished)
        self._active_loaders.append(self.basemap_loader)
        self.basemap_loader.start()

    @staticmethod
    def _wgom_source_pixels(bbox, pixel_size_x, pixel_size_y):
        """Return source pixel dimensions for a WGOM AOI (bbox already in EPSG:3857 meters)."""
        xmin, ymin, xmax, ymax = bbox
        pixels_x = int((xmax - xmin) / abs(pixel_size_x))
        pixels_y = int((ymax - ymin) / abs(pixel_size_y))
        return pixels_x, pixels_y
    
    def load_map(self):
        """Load map for current extent."""
        # Cancel any pending load timer
        if self._load_timer:
            self._load_timer.stop()
            self._load_timer = None

        # Replace any in-progress load (source switches / resize) instead of skipping
        if self._loading:
            print("load_map() interrupting in-progress load...")
        self._stop_all_loaders()
        self._loading = True
        
        # Preserve the current extent - we'll use this for the request
        # and restore it after loading to prevent the selection from moving
        requested_extent = self.extent
        
        print("=" * 50)
        print("load_map() called!")
        print(f"Widget visible: {self.isVisible()}")
        print(f"Widget size: {self.width()}x{self.height()}")
        print(f"Extent: {self.extent}")
        print(f"Base URL: {self.base_url}")
        
        # Ensure widget has a valid size
        widget_width = self.width()
        widget_height = self.height()
        
        # If widget has no size yet, use minimum size
        if widget_width <= 0 or widget_height <= 0:
            widget_width = 800
            widget_height = 600
            print(f"Widget has no size yet, using default: {widget_width}x{widget_height}")
        else:
            print(f"Widget size: {widget_width}x{widget_height}")
        
        # Use widget size to fill the window completely
        size = (widget_width, widget_height)
        
        # Determine raster function based on area of interest pixel dimensions in source data
        # Check if this is the Hi Resolution or Regional service by checking the base_url
        is_hi_resolution = "WGOM_LI_SNE_BTY_4m" in self.base_url
        is_regional = "WGOM_LI_SNE_BTY" in self.base_url and "16m" in self.base_url
        uses_dynamic_raster_function = is_hi_resolution or is_regional
        
        if uses_dynamic_raster_function:
            # Get the area of interest (selected bbox) or use current extent if no selection
            area_bbox = self.selected_bbox_world if self.selected_bbox_world else requested_extent
            
            # Get pixel size from service (default based on service type if not available)
            # IMPORTANT: Use actual pixel sizes from service, not defaults, unless they're None
            if is_hi_resolution:
                default_pixel_size = 4.0
            else:  # Regional
                default_pixel_size = 16.0
            
            # Use actual pixel sizes from service if available, otherwise use service-specific default
            pixel_size_x = self.pixel_size_x if self.pixel_size_x is not None else default_pixel_size
            pixel_size_y = self.pixel_size_y if self.pixel_size_y is not None else default_pixel_size
            
            # AOI is GCS; WGOM pixel sizes are meters — convert before dividing
            pixels_x, pixels_y = self._wgom_source_pixels(area_bbox, pixel_size_x, pixel_size_y)
            
            # Use "StdDev - BlueGreen" for areas > 4000 pixels in either dimension
            # Use "DAR - StdDev - BlueGreen" for areas <= 4000 pixels in both dimensions
            if pixels_x > 4000 or pixels_y > 4000:
                new_raster_function = "StdDev - BlueGreen"
            else:
                new_raster_function = "DAR - StdDev - BlueGreen"
            
            # Update raster function if it changed
            if self.raster_function != new_raster_function:
                msg = f"Updating raster function based on area of interest size ({pixels_x}x{pixels_y} source pixels): {self.raster_function} -> {new_raster_function}"
                print(msg)
                # Format message with green color for raster function info
                green_msg = f'<span style="color: green;">{msg}</span>'
                self.statusMessage.emit(green_msg)
                self.raster_function = new_raster_function
        
        # Always log which raster function is being used
        if uses_dynamic_raster_function and self.selected_bbox_world:
            area_bbox = self.selected_bbox_world
            # Get pixel size from service (default based on service type if not available)
            if is_hi_resolution:
                default_pixel_size = 4.0
            else:  # Regional
                default_pixel_size = 16.0
            pixel_size_x = self.pixel_size_x if self.pixel_size_x is not None else default_pixel_size
            pixel_size_y = self.pixel_size_y if self.pixel_size_y is not None else default_pixel_size
            pixels_x, pixels_y = self._wgom_source_pixels(area_bbox, pixel_size_x, pixel_size_y)
            raster_info_msg = f"Map display using raster function: {self.raster_function} (area of interest: {pixels_x}x{pixels_y} source pixels)"
        else:
            raster_info_msg = f"Map display using raster function: {self.raster_function}"
        print(raster_info_msg)
        # Format message with green color for raster function info
        green_raster_info_msg = f'<span style="color: green;">{raster_info_msg}</span>'
        self.statusMessage.emit(green_raster_info_msg)
        
        print(f"Starting map load with extent: {requested_extent}, size: {size}")
        print(f"Using raster function: {self.raster_function}")
        
        # Store the requested extent so we can restore it after loading
        self._requested_extent = requested_extent

        # GMRT custom ImageServer preview (JPEG; not ArcGIS)
        if getattr(self, "preview_url", None):
            print("Loading GMRT ImageServer preview...")
            from .gmrt_module import GMRTMapLoader
            self.loader = GMRTMapLoader(
                requested_extent, size,
                mask=bool(getattr(self, "gmrt_mask", False)),
                purpose="GMRT preview",
            )
            self.loader.statusMessage.connect(self.statusMessage.emit)
            self.loader.tileLoaded.connect(self.on_tile_loaded)
            self.loader.finished.connect(self.on_loader_finished)
            self.loader.finished.connect(self._check_all_loaders_finished)
            self._active_loaders.append(self.loader)
            self.loader.start()
            self.map_loaded = True
            print("=" * 50)
            return
        
        # When display_url is set (e.g. GEBCO MapServer), use GCS extent: land basemap + display layer
        if self.display_url:
            self._start_land_basemap_loader(requested_extent, size)
            print("Loading display layer (GCS)...")
            # Bathymetry layer: transparent (transparent=True) so land shows through
            self.loader = MapServerLoader(self.display_url, requested_extent, size, transparent=True, purpose="bathymetry display")
            self.loader.statusMessage.connect(self.statusMessage.emit)
            self.loader.tileLoaded.connect(self.on_tile_loaded)
            self.loader.finished.connect(self.on_loader_finished)
            self.loader.finished.connect(self._check_all_loaders_finished)
            self._active_loaders.append(self.loader)
            self.loader.start()
            self.map_loaded = True
            print("=" * 50)
            return
        
        # ImageServer path (e.g. NCEI multibeam): optional land underlay + raster function overlay
        self._start_land_basemap_loader(requested_extent, size)

        # Load basemap if enabled (World Imagery — same CRS/extent as bathymetry overlay)
        if self.show_basemap:
            print("Loading basemap...")
            self.basemap_loader = BasemapLoader(requested_extent, size, bbox_sr=self.bbox_sr)
            self.basemap_loader.statusMessage.connect(self.statusMessage.emit)
            self.basemap_loader.tileLoaded.connect(self.on_basemap_loaded)
            self.basemap_loader.finished.connect(self._check_all_loaders_finished)
            self._active_loaders.append(self.basemap_loader)
            self.basemap_loader.start()
        
        # Load hillshade layer if enabled (as underlay)
        if self.show_hillshade:
            print("Loading hillshade layer...")
            self.hillshade_loader = MapTileLoader(
                self.base_url, requested_extent, size, self.hillshade_raster_function,
                bbox_sr=self.bbox_sr, purpose="hillshade display",
            )
            self.hillshade_loader.statusMessage.connect(self.statusMessage.emit)
            self.hillshade_loader.tileLoaded.connect(self.on_hillshade_loaded)
            self.hillshade_loader.finished.connect(self._check_all_loaders_finished)
            self._active_loaders.append(self.hillshade_loader)
            self.hillshade_loader.start()
        
        # Load bathymetry layer (main layer)
        print("Creating MapTileLoader...")
        self.loader = MapTileLoader(
            self.base_url, requested_extent, size, self.raster_function,
            bbox_sr=self.bbox_sr, purpose="bathymetry display",
        )
        print("Connecting signals...")
        self.loader.statusMessage.connect(self.statusMessage.emit)
        self.loader.tileLoaded.connect(self.on_tile_loaded)
        self.loader.finished.connect(self.on_loader_finished)
        self.loader.finished.connect(self._check_all_loaders_finished)
        self._active_loaders.append(self.loader)
        print("Starting loader thread...")
        self.loader.start()
        print(f"Loader thread started, isRunning: {self.loader.isRunning()}")
        self.map_loaded = True
        print("=" * 50)
        
    def on_basemap_loaded(self, pixmap):
        """Handle basemap tile loaded."""
        sender = self.sender()
        if sender is not None and sender is not getattr(self, "basemap_loader", None):
            print("Ignoring stale basemap from abandoned loader")
            return
        if not pixmap.isNull():
            # Ensure basemap is fully opaque (remove alpha channel if present)
            # Convert to QImage, then to RGB format to remove transparency
            qimage = pixmap.toImage()
            if qimage.hasAlphaChannel():
                # Convert to RGB888 format (removes alpha channel, truly opaque)
                qimage = qimage.convertToFormat(QImage.Format.Format_RGB888)
                pixmap = QPixmap.fromImage(qimage)
            
            widget_size = self.size()
            if widget_size.width() > 0 and widget_size.height() > 0:
                if widget_size.width() == pixmap.width() and widget_size.height() == pixmap.height():
                    self.basemap_pixmap = pixmap
                else:
                    self.basemap_pixmap = pixmap.scaled(
                        widget_size,
                        Qt.AspectRatioMode.KeepAspectRatio,
                        Qt.TransformationMode.SmoothTransformation
                    )
            else:
                self.basemap_pixmap = pixmap
            print(f"Basemap loaded: {self.basemap_pixmap.width()}x{self.basemap_pixmap.height()}")
            self._sync_basemap_to_current_pixmap()
            self.update()  # Trigger repaint
            
    def on_hillshade_loaded(self, pixmap, xmin, ymin, xmax, ymax):
        """Handle hillshade tile loaded."""
        if not pixmap.isNull():
            widget_size = self.size()
            if widget_size.width() > 0 and widget_size.height() > 0:
                if widget_size.width() == pixmap.width() and widget_size.height() == pixmap.height():
                    self.hillshade_pixmap = pixmap
                else:
                    self.hillshade_pixmap = pixmap.scaled(
                        widget_size,
                        Qt.AspectRatioMode.KeepAspectRatio,
                        Qt.TransformationMode.SmoothTransformation
                    )
            else:
                self.hillshade_pixmap = pixmap
            print(f"Hillshade loaded: {self.hillshade_pixmap.width()}x{self.hillshade_pixmap.height()}")
            self.update()  # Trigger repaint
        
    def on_loader_finished(self):
        """Handle loader thread finishing."""
        # If we still have an empty pixmap, the load might have failed
        if self.current_pixmap.isNull():
            print("Warning: Map tile loader finished but no pixmap was loaded")
        
    def on_tile_loaded(self, pixmap, xmin, ymin, xmax, ymax):
        """Handle loaded tile."""
        # Ignore callbacks from abandoned loaders after a source/extent switch
        sender = self.sender()
        if sender is not None and sender is not getattr(self, "loader", None):
            print("Ignoring stale tileLoaded from abandoned loader")
            return
        print(f"on_tile_loaded called! pixmap.isNull: {pixmap.isNull()}, size: {pixmap.width()}x{pixmap.height()}")
        if not pixmap.isNull():
            # Check if pixmap has actual content (not all white/transparent)
            # Sample a few pixels to verify
            sample_image = pixmap.toImage()
            if not sample_image.isNull():
                # Sample a few pixels
                colors = []
                for x in [10, pixmap.width()//2, pixmap.width()-10]:
                    for y in [10, pixmap.height()//2, pixmap.height()-10]:
                        if x < pixmap.width() and y < pixmap.height():
                            color = sample_image.pixelColor(x, y)
                            colors.append((color.red(), color.green(), color.blue()))
                print(f"Sample pixel colors: {colors[:3]}...")  # Print first 3
            
            widget_size = self.size()
            print(f"Widget size in on_tile_loaded: {widget_size.width()}x{widget_size.height()}")
            
            # Store original pixmap size before any scaling
            # This is needed for accurate world-to-screen coordinate conversion
            # The original pixmap size represents the actual pixel dimensions requested from the server
            # which directly correspond to the geographic extent
            self._original_pixmap_size = (pixmap.width(), pixmap.height())
            
            # CRITICAL: For coordinate conversion, we need to use a consistent reference size
            # When the pixmap matches the widget size, we should still use the widget size
            # for coordinate conversion to maintain consistent visual size of the selection box
            # When the pixmap is scaled, we use the scaled size
            
            # Don't scale if sizes match - use pixmap directly
            if widget_size.width() == pixmap.width() and widget_size.height() == pixmap.height():
                print("Pixmap size matches widget, using directly")
                self.current_pixmap = pixmap
                # Use widget size for coordinate conversion to maintain consistent visual size
                # This ensures the selection box doesn't change size when pixmap pixel size changes
                self._scaled_pixmap_size = (widget_size.width(), widget_size.height())
            elif widget_size.width() > 0 and widget_size.height() > 0:
                # Scale to fit widget while maintaining aspect ratio
                scaled_pixmap = pixmap.scaled(
                    widget_size, 
                    Qt.AspectRatioMode.KeepAspectRatio,
                    Qt.TransformationMode.SmoothTransformation
                )
                print(f"Scaled pixmap: {scaled_pixmap.width()}x{scaled_pixmap.height()}")
                self.current_pixmap = scaled_pixmap
                # Use the scaled size (what's actually drawn) for coordinate conversion
                self._scaled_pixmap_size = (scaled_pixmap.width(), scaled_pixmap.height())
            else:
                # Widget not sized yet, use pixmap as-is
                print("Widget not sized, using pixmap as-is")
                self.current_pixmap = pixmap
                self._scaled_pixmap_size = (pixmap.width(), pixmap.height())  # No scaling
                
            # ALWAYS preserve the requested extent instead of using the server's response
            # This prevents the selection from moving when the window is resized
            # The _requested_extent is set in load_map() and should be preserved
            # until the user explicitly changes the extent (via pan/zoom)
            if self._extent_locked:
                # Extent is locked (during resize) - never change it
                pass  # Keep current extent unchanged
            elif hasattr(self, '_requested_extent') and self._requested_extent is not None:
                # Restore the extent we requested, not what the server returned
                # This ensures coordinate conversion uses the correct extent
                self.extent = self._requested_extent
                # Keep _requested_extent so it persists across multiple loads during resize
            else:
                # First load - use server response for accurate coordinate conversion
                # The server's returned extent matches what's actually displayed
                # However, if we have a _requested_extent set (from initial load), use that instead
                # to ensure coordinate conversion matches what we intended to display
                if hasattr(self, '_requested_extent') and self._requested_extent is not None:
                    # Use the requested extent (which should be the service extent)
                    self.extent = self._requested_extent
                    # Keep _requested_extent for consistency
                else:
                    # Fallback to server response if no requested extent
                    server_extent = (xmin, ymin, xmax, ymax)
                    self.extent = server_extent
                    self._requested_extent = server_extent
            print(f"Setting current_pixmap, isNull: {self.current_pixmap.isNull()}, size: {self.current_pixmap.width()}x{self.current_pixmap.height()}")
            self._sync_basemap_to_current_pixmap()
            print(f"Calling update() to repaint widget")
            
            # CRITICAL: Force immediate repaint to ensure selection box is redrawn with new pixmap size
            # This is especially important during window resize when pixmap size changes
            # The selection box will be recalculated in paintEvent using the new pixmap size
            self.update()  # Trigger repaint immediately
            self.repaint()  # Force immediate repaint
            
            # Also schedule a delayed repaint to ensure box is visible after map loads
            # This is especially important after zoom operations
            if self.selected_bbox_world:
                from PyQt6.QtCore import QTimer
                QTimer.singleShot(100, lambda: self.update())
            
            print(f"Map tile loaded successfully: {pixmap.width()}x{pixmap.height()}")
            if self.selected_bbox_world:
                print(f"Selected bbox exists: {self.selected_bbox_world}, will be repainted")
            
            # Emit signal on first successful load
            if not self._first_load_complete:
                self._first_load_complete = True
                # Don't set default bounds here - let MainWindow handle it after extent is confirmed
                # The extent at this point matches what's displayed, so coordinate conversion will be accurate
                self.mapFirstLoaded.emit()
        else:
            print("Error: Received null pixmap from tile loader")
            
    def _sync_basemap_to_current_pixmap(self):
        """Keep land basemap aligned with the bathymetry layer dimensions."""
        if self.basemap_pixmap.isNull() or self.current_pixmap.isNull():
            return
        if self.basemap_pixmap.size() != self.current_pixmap.size():
            self.basemap_pixmap = self.basemap_pixmap.scaled(
                self.current_pixmap.size(),
                Qt.AspectRatioMode.IgnoreAspectRatio,
                Qt.TransformationMode.SmoothTransformation,
            )

    def _extent_for_conversion(self):
        """Return the geographic extent used for screen/world coordinate conversion."""
        if hasattr(self, '_requested_extent') and self._requested_extent is not None:
            return self._requested_extent
        return self.extent

    @staticmethod
    def _lat_to_mercator_y(lat):
        """Web Mercator Y (unitless) used by GMRT ImageServer latitude mapping."""
        lat = max(min(lat, 85.05112878), -85.05112878)
        return math.log(math.tan(math.pi / 4.0 + math.radians(lat) / 2.0))

    @staticmethod
    def _mercator_y_to_lat(merc_y):
        """Inverse of _lat_to_mercator_y."""
        return math.degrees(2.0 * math.atan(math.exp(merc_y)) - math.pi / 2.0)

    def _uses_gmrt_preview(self):
        """True when the active map image is a GMRT ImageServer JPEG."""
        return bool(getattr(self, "preview_url", None))

    def _pixmap_draw_rect(self):
        """Return the on-screen rectangle where the map pixmap is painted (centered)."""
        widget_rect = self.rect()
        if self.current_pixmap.isNull():
            return widget_rect
        pr = self.current_pixmap.rect()
        x = (widget_rect.width() - pr.width()) // 2
        y = (widget_rect.height() - pr.height()) // 2
        return QRect(x, y, pr.width(), pr.height())

    def _geographic_data_rect(self):
        """
        Return the screen rectangle where geographic coordinates map.

        ArcGIS exports letterbox when geographic aspect differs from image aspect
        (e.g. global 2:1 extent in a ~4:3 widget), leaving transparent margins.

        GMRT ImageServer JPEGs use Web Mercator Y and include fixed chrome margins
        (~33px L/R, ~13px T/B at native size), matching GMRT_Downloader.
        """
        map_rect = self._pixmap_draw_rect()
        if map_rect.isEmpty():
            return map_rect

        if self._uses_gmrt_preview():
            # Scale GMRT chrome insets with the displayed pixmap
            ow, oh = self._original_pixmap_size or (map_rect.width(), map_rect.height())
            if ow <= 0 or oh <= 0:
                return map_rect
            sx = map_rect.width() / float(ow)
            sy = map_rect.height() / float(oh)
            left = int(round(33 * sx))
            right = int(round(33 * sx))
            top = int(round(13 * sy))
            bottom = int(round(13 * sy))
            content_w = map_rect.width() - left - right
            content_h = map_rect.height() - top - bottom
            if content_w <= 1 or content_h <= 1:
                return map_rect
            return QRect(map_rect.left() + left, map_rect.top() + top, content_w, content_h)

        extent = self._extent_for_conversion()
        if extent is None:
            return map_rect

        xmin, ymin, xmax, ymax = extent
        geo_w = xmax - xmin
        geo_h = ymax - ymin
        if geo_w <= 0 or geo_h <= 0:
            return map_rect

        geo_aspect = geo_w / geo_h
        pix_aspect = map_rect.width() / map_rect.height()

        if geo_aspect > pix_aspect:
            draw_w = map_rect.width()
            draw_h = int(round(map_rect.width() / geo_aspect))
            x = map_rect.left()
            y = map_rect.top() + (map_rect.height() - draw_h) // 2
        else:
            draw_h = map_rect.height()
            draw_w = int(round(map_rect.height() * geo_aspect))
            x = map_rect.left() + (map_rect.width() - draw_w) // 2
            y = map_rect.top()

        return QRect(x, y, draw_w, draw_h)

    def screen_to_world(self, point):
        """Convert screen coordinates to world coordinates."""
        if self.current_pixmap.isNull():
            return None

        data_rect = self._geographic_data_rect()
        if not data_rect.contains(point):
            return None

        rel_x = (point.x() - data_rect.left()) / data_rect.width()
        rel_y = (point.y() - data_rect.top()) / data_rect.height()

        extent = self._extent_for_conversion()
        if extent is None:
            return None
        xmin, ymin, xmax, ymax = extent
        
        world_x = xmin + rel_x * (xmax - xmin)
        if self._uses_gmrt_preview():
            # GMRT ImageServer maps Y in Web Mercator space
            merc_north = self._lat_to_mercator_y(ymax)
            merc_south = self._lat_to_mercator_y(ymin)
            merc_y = merc_north - rel_y * (merc_north - merc_south)
            world_y = self._mercator_y_to_lat(merc_y)
        else:
            world_y = ymax - rel_y * (ymax - ymin)  # equirectangular (GEBCO/NCEI/WGOM display)
        
        return (world_x, world_y)
        
    def world_to_screen(self, world_x, world_y):
        """Convert world coordinates to screen coordinates."""
        if self.current_pixmap.isNull():
            return None

        target_rect = self._geographic_data_rect()
        extent = self._extent_for_conversion()
        if extent is None:
            return None
        xmin, ymin, xmax, ymax = extent
        
        # Calculate relative position within the extent (0.0 to 1.0)
        rel_x = (world_x - xmin) / (xmax - xmin) if (xmax - xmin) != 0 else 0
        if self._uses_gmrt_preview():
            merc_north = self._lat_to_mercator_y(ymax)
            merc_south = self._lat_to_mercator_y(ymin)
            merc_y = self._lat_to_mercator_y(world_y)
            rel_y = (merc_north - merc_y) / (merc_north - merc_south) if (merc_north - merc_south) != 0 else 0
        else:
            rel_y = (ymax - world_y) / (ymax - ymin) if (ymax - ymin) != 0 else 0  # Y is inverted
        
        screen_x = target_rect.left() + rel_x * target_rect.width()
        screen_y = target_rect.top() + rel_y * target_rect.height()
        
        # Clamp coordinates to widget bounds to prevent drawing outside the widget
        widget_rect = self.rect()
        clamped_x = max(0, min(int(screen_x), widget_rect.width() - 1))
        clamped_y = max(0, min(int(screen_y), widget_rect.height() - 1))
        return QPoint(clamped_x, clamped_y)
        
    def get_selection_bbox(self):
        """Get the bounding box of the current selection in world coordinates."""
        if not self.selection_start or not self.selection_end:
            return None
            
        start_world = self.screen_to_world(self.selection_start)
        end_world = self.screen_to_world(self.selection_end)
        
        if not start_world or not end_world:
            return None
            
        xmin = min(start_world[0], end_world[0])
        xmax = max(start_world[0], end_world[0])
        ymin = min(start_world[1], end_world[1])
        ymax = max(start_world[1], end_world[1])
        
        return (xmin, ymin, xmax, ymax)
        
    def clear_selection(self):
        """Clear the current selection."""
        self.selection_start = None
        self.selection_end = None
        self.is_selecting = False
        self.selected_bbox_world = None  # Also clear persistent selection
        self.update()
        self.selectionChanged.emit(0, 0, 0, 0)
        
    def world_bbox_to_screen_rect(self, bbox_world):
        """Convert world bbox to screen rectangle for drawing."""
        if bbox_world is None:
            return None
        
        # CRITICAL: Ensure pixmap is valid before converting
        if self.current_pixmap.isNull() or not self.map_loaded:
            return None
            
        xmin, ymin, xmax, ymax = bbox_world
        
        # Use _requested_extent if available (matches what's displayed), otherwise use extent
        # This ensures coordinate conversion uses the correct extent
        conversion_extent = self._requested_extent if (hasattr(self, '_requested_extent') and self._requested_extent is not None) else self.extent
        
        # Use conversion_extent directly for coordinate conversion
        # Store original values but don't modify self.extent (world_to_screen will use conversion_extent via parameter)
        original_extent = self.extent
        original_requested = getattr(self, '_requested_extent', None)
        
        # Temporarily set extent for world_to_screen to use
        self.extent = conversion_extent
        if not hasattr(self, '_requested_extent'):
            self._requested_extent = None
        self._requested_extent = conversion_extent
        
        try:
            # Convert corners to screen coordinates
            # world_to_screen will use self._requested_extent which we just set to conversion_extent
            top_left = self.world_to_screen(xmin, ymax)
            bottom_right = self.world_to_screen(xmax, ymin)
            
            if top_left is not None and bottom_right is not None:
                screen_rect = QRect(top_left, bottom_right)
                return screen_rect
        finally:
            # Restore original extent
            self.extent = original_extent
            if hasattr(self, '_requested_extent'):
                self._requested_extent = original_requested
        
        return None
        
    def mousePressEvent(self, event):
        """Handle mouse press for selection or panning."""
        if event.button() == Qt.MouseButton.LeftButton:
            # Left button for selection only
            # Clear previous selection when starting a new one
            self.clear_selection()
            self.selection_start = event.position().toPoint()
            self.selection_end = self.selection_start
            self.is_selecting = True
            self.update()
        elif event.button() == Qt.MouseButton.MiddleButton:
            # Middle button for panning
            self.is_panning = True
            self.pan_start = event.position().toPoint()
            self.pan_origin = event.position().toPoint()  # Store original position for pan line
            self.pan_end = event.position().toPoint()  # Initialize pan_end
            self.update()
            
    def mouseMoveEvent(self, event):
        """Handle mouse move for selection or panning."""
        if self.is_selecting:
            self.selection_end = event.position().toPoint()
            bbox = self.get_selection_bbox()
            if bbox:
                self.selectionChanged.emit(*bbox)
            self.update()
        elif self.is_panning and self.pan_start:
            # Calculate pan delta
            current_pos = event.position().toPoint()
            self.pan_end = current_pos  # Track current position for pan line
            delta = current_pos - self.pan_start
            
            # Convert screen delta to world delta using the geographic data rect
            data_rect = self._geographic_data_rect()
            
            if not data_rect.isNull():
                extent = self._extent_for_conversion()
                if extent is None:
                    return
                xmin, ymin, xmax, ymax = extent
                world_width = xmax - xmin
                world_height = ymax - ymin
                
                rel_delta_x = -delta.x() / data_rect.width() * world_width
                rel_delta_y = delta.y() / data_rect.height() * world_height
                
                # Update extent
                self.extent = (
                    xmin + rel_delta_x,
                    ymin + rel_delta_y,
                    xmax + rel_delta_x,
                    ymax + rel_delta_y
                )
                # Update _requested_extent to match the new extent for accurate coordinate conversion
                self._requested_extent = self.extent
                self.pan_start = current_pos
                self.clear_selection()
                
                # Update display to show pan line
                self.update()
                
                # Debounce: Cancel any pending load and schedule a new one after a delay
                if self._load_timer:
                    self._load_timer.stop()
                
                from PyQt6.QtCore import QTimer
                self._load_timer = QTimer()
                self._load_timer.setSingleShot(True)
                self._load_timer.timeout.connect(self.load_map)
                self._load_timer.start(300)  # Wait 300ms before loading (debounce)
            
    def mouseReleaseEvent(self, event):
        """Handle mouse release for selection or panning."""
        if event.button() == Qt.MouseButton.LeftButton:
            if self.is_selecting:
                self.selection_end = event.position().toPoint()
                self.is_selecting = False
                bbox = self.get_selection_bbox()
                if bbox:
                    # Store the selected bbox in world coordinates for persistent display
                    self.selected_bbox_world = bbox
                    self.selectionChanged.emit(*bbox)
                    # Emit selection completed signal for zooming
                    self.selectionCompleted.emit(*bbox)
                # Clear the active selection rectangle (red dashed line)
                self.selection_start = None
                self.selection_end = None
                self.update()
            elif self.is_panning:
                self.is_panning = False
                self.pan_start = None
                self.pan_origin = None
                self.pan_end = None
                self.update()  # Clear pan line
        elif event.button() == Qt.MouseButton.MiddleButton and self.is_panning:
            self.is_panning = False
            self.pan_start = None
            self.pan_origin = None
            self.pan_end = None
            self.update()  # Clear pan line
            self.userViewChanged.emit()
            
    def wheelEvent(self, event):
        """Handle mouse wheel for zooming."""
        if self.current_pixmap.isNull():
            return
            
        # Get center of widget in world coordinates
        widget_center = QPoint(self.width() // 2, self.height() // 2)
        world_pos = self.screen_to_world(widget_center)
        
        if not world_pos:
            # Fallback: use center of extent if screen_to_world fails
            xmin, ymin, xmax, ymax = self.extent
            center_x = (xmin + xmax) / 2
            center_y = (ymin + ymax) / 2
            world_pos = (center_x, center_y)
            
        # Calculate zoom factor
        zoom_factor = 1.2 if event.angleDelta().y() > 0 else 1 / 1.2
        
        # Calculate new extent centered on window center
        xmin, ymin, xmax, ymax = self.extent
        width = (xmax - xmin) / zoom_factor
        height = (ymax - ymin) / zoom_factor
        
        center_x, center_y = world_pos
        new_xmin = center_x - width / 2
        new_xmax = center_x + width / 2
        new_ymin = center_y - height / 2
        new_ymax = center_y + height / 2
        
        self.extent = (new_xmin, new_ymin, new_xmax, new_ymax)
        # Update _requested_extent to match the new extent for accurate coordinate conversion
        self._requested_extent = self.extent
        self.clear_selection()
        
        # Debounce: Cancel any pending load and schedule a new one after a delay
        if self._load_timer:
            self._load_timer.stop()
        
        from PyQt6.QtCore import QTimer
        self._load_timer = QTimer()
        self._load_timer.setSingleShot(True)
        self._load_timer.timeout.connect(self.load_map)
        self._load_timer.start(300)  # Wait 300ms before loading (debounce)
        self.userViewChanged.emit()
            
    def resizeEvent(self, event):
        """Handle widget resize."""
        super().resizeEvent(event)
        # Don't scale pixmap here - let the map loading handle it properly
        # Scaling here interferes with coordinate conversion and causes selection box size issues
        # The map will be reloaded by _refresh_map_on_resize in main.py
            
    def paintEvent(self, event):
        """Paint the map and selection rectangle."""
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        
        widget_rect = self.rect()
        painter.fillRect(widget_rect, QColor(0, 0, 0))
        
        # Draw basemap first (if available) - bottom layer
        # Draw if show_basemap is True OR if land_display_url is set (land layer for GEBCO)
        should_draw_basemap = (self.show_basemap or self.land_display_url) and not self.basemap_pixmap.isNull()
        if should_draw_basemap:
            basemap_rect = self.basemap_pixmap.rect()
            x = (widget_rect.width() - basemap_rect.width()) // 2
            y = (widget_rect.height() - basemap_rect.height()) // 2
            target_rect = QRect(x, y, basemap_rect.width(), basemap_rect.height())
            painter.drawPixmap(target_rect, self.basemap_pixmap)
        
        # Draw hillshade layer (if available) - middle layer (underlay) at full opacity
        if self.show_hillshade and not self.hillshade_pixmap.isNull():
            hillshade_rect = self.hillshade_pixmap.rect()
            x = (widget_rect.width() - hillshade_rect.width()) // 2
            y = (widget_rect.height() - hillshade_rect.height()) // 2
            target_rect = QRect(x, y, hillshade_rect.width(), hillshade_rect.height())
            # Ensure full opacity for hillshade (opacity only affects top layer)
            painter.setOpacity(1.0)
            painter.drawPixmap(target_rect, self.hillshade_pixmap)
        
        # Draw bathymetry layer on top (with opacity and/or blend mode) - top layer
        # Only this layer uses the opacity setting and blend mode
        if not self.current_pixmap.isNull():
            pixmap_rect = self.current_pixmap.rect()
            x = (widget_rect.width() - pixmap_rect.width()) // 2
            y = (widget_rect.height() - pixmap_rect.height()) // 2
            target_rect = QRect(x, y, pixmap_rect.width(), pixmap_rect.height())
            
            # Set blend mode if enabled (Multiply mode for natural blending with hillshade)
            # Multiply darkens the colors while preserving hillshade detail, better than Overlay for cartography
            if self.use_blend:
                painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Multiply)
            
            # Draw with opacity (only the top layer uses opacity)
            painter.setOpacity(self.bathymetry_opacity)
            painter.drawPixmap(target_rect, self.current_pixmap)
            
            # Reset opacity and composition mode for subsequent drawing
            painter.setOpacity(1.0)
            painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceOver)
        elif not self.show_basemap and not self.land_display_url:
            # Draw placeholder only if no map layers are configured
            status_text = "Loading map..."
            if hasattr(self, 'loader') and self.loader and self.loader.isRunning():
                status_text = "Loading map..."
            else:
                status_text = "No map data available"
            painter.setPen(QColor(255, 255, 255))  # White text on black background
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, status_text)
            
        # Draw selection rectangle (always on top)
        # First draw the persistent selected bbox if it exists
        # CRITICAL: Only draw if pixmap is loaded and valid to avoid drawing with stale data during resize
        # Also ensure pixmap size matches widget size (or is being scaled correctly)
        if self.show_aoi and self.selected_bbox_world and self.map_loaded and not self.current_pixmap.isNull():
            # Only draw if pixmap is valid and has been loaded
            # Check that pixmap size is reasonable (not stale)
            pixmap_width = self.current_pixmap.width()
            pixmap_height = self.current_pixmap.height()
            widget_width = self.width()
            widget_height = self.height()
            
            # Only draw if pixmap dimensions are valid (greater than 0)
            if pixmap_width > 0 and pixmap_height > 0:
                bbox_screen = self.world_bbox_to_screen_rect(self.selected_bbox_world)
                if bbox_screen:
                    # Draw selection rectangle based on validity
                    if self.selection_is_valid:
                        # Valid selection - use green dashed line (no fill)
                        pen = QPen(QColor(0, 255, 0), 2, Qt.PenStyle.DashLine)  # Green dashed line
                    else:
                        # Selection too large - use red dashed line (no fill)
                        pen = QPen(QColor(255, 0, 0), 2, Qt.PenStyle.DashLine)  # Red dashed line
                    painter.setPen(pen)
                    painter.setBrush(Qt.BrushStyle.NoBrush)  # No fill - outline only
                    painter.drawRect(bbox_screen)
        
        # Draw active selection rectangle (while dragging)
        if self.show_aoi and self.selection_start and self.selection_end:
            selection_rect = QRect(self.selection_start, self.selection_end).normalized()
            pen = QPen(QColor(0, 255, 0), 2, Qt.PenStyle.DashLine)  # Green dashed line
            painter.setPen(pen)
            painter.setBrush(Qt.BrushStyle.NoBrush)  # No fill - outline only
            painter.drawRect(selection_rect)
        
        # Draw pan line (red line showing pan direction and distance)
        if self.is_panning and self.pan_origin and self.pan_end:
            pen = QPen(QColor(255, 0, 0), 2, Qt.PenStyle.DashLine)  # Red dashed line
            painter.setPen(pen)
            painter.drawLine(self.pan_origin, self.pan_end)
        
        # Draw legend in upper left corner
        self._draw_legend(painter)
    
    def _draw_legend(self, painter):
        """Draw a legend in the upper left corner showing box color meanings."""
        if not self.map_loaded or not self.show_legend:
            return  # Don't draw legend until map is loaded or if legend is disabled
        
        # Legend configuration
        margin = 10
        padding = 8
        line_height = 20
        line_width = 30
        legend_width = 150
        # Height calculation:
        # - Top padding: padding
        # - Item 1: line_height (text drawn at item_y + line_height - 4)
        # - Bottom padding: padding + 4 (extra space for text)
        # Total: 2*padding + line_height + 4
        legend_height = padding * 2 + line_height + 4
        
        # Position in upper left corner
        x = margin
        y = margin
        
        # Draw semi-transparent background
        legend_rect = QRect(x, y, legend_width, legend_height)
        bg_color = QColor(0, 0, 0, 140)  # Black with 140/255 opacity (~55% opaque)
        painter.fillRect(legend_rect, bg_color)
        
        # Draw border
        border_pen = QPen(QColor(255, 255, 255), 1)
        painter.setPen(border_pen)
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawRect(legend_rect)
        
        # Draw legend items (Area of Interest)
        items = [
            (QColor(0, 255, 0), "Area of Interest")
        ]
        
        start_y = y + padding
        for i, (color, label) in enumerate(items):
            item_y = start_y + i * line_height
            
            # Draw colored line sample (dashed)
            line_x = x + padding
            line_y = item_y + line_height // 2
            pen = QPen(color, 2, Qt.PenStyle.DashLine)
            painter.setPen(pen)
            painter.drawLine(line_x, line_y, line_x + line_width, line_y)
            
            # Draw label (text baseline at item_y + line_height - 4 to account for text height)
            painter.setPen(QColor(255, 255, 255))
            text_x = line_x + line_width + 8
            painter.drawText(text_x, item_y + line_height - 4, label)


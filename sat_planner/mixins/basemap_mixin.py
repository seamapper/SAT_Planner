"""
Basemap and web imagery: ArcGIS World Imagery tiles, NOAA ENC Charts overlay.
Toggle, load, and display; reload on zoom/pan is coordinated from map_interaction_mixin.
"""
import traceback
import json
import numpy as np
from urllib.request import urlopen, Request
from urllib.parse import urlencode
from io import BytesIO

from PIL import Image

from sat_planner.constants import GEOSPATIAL_LIBS_AVAILABLE, pyproj, transform, reproject, Resampling


class BasemapMixin:
    """Mixin for basemap imagery and NOAA ENC Charts overlay (web-served tiles/export)."""

    def _toggle_imagery_basemap(self):
        """Toggle imagery basemap on/off."""
        if hasattr(self, 'imagery_basemap_checkbox'):
            # Store current view limits before enabling basemap
            current_xlim = self.ax.get_xlim()
            current_ylim = self.ax.get_ylim()
            
            self.show_imagery_basemap_var = self.imagery_basemap_checkbox.isChecked()
            if not hasattr(self, '_last_basemap_zoom_level'):
                self._last_basemap_zoom_level = None
            if not hasattr(self, '_last_basemap_extent'):
                self._last_basemap_extent = None
            if not self.show_imagery_basemap_var:
                if hasattr(self, 'basemap_image_plots') and isinstance(self.basemap_image_plots, list):
                    for _plot in self.basemap_image_plots:
                        try:
                            _plot.remove()
                        except Exception:
                            pass
                    self.basemap_image_plots = []
            
            # Ensure view limits are stored for preservation
            self._last_user_xlim = current_xlim
            self._last_user_ylim = current_ylim
            
            self._plot_survey_plan(preserve_view_limits=True)
            
            # Explicitly restore the view limits after plotting to ensure they match what was shown before
            if current_xlim and current_ylim:
                self.ax.set_xlim(current_xlim)
                self.ax.set_ylim(current_ylim)
                self.canvas.draw_idle()

    def _toggle_noaa_charts(self):
        """Toggle NOAA ENC Charts overlay on/off."""
        if hasattr(self, 'noaa_charts_checkbox'):
            self.show_noaa_charts_var = self.noaa_charts_checkbox.isChecked()
            if hasattr(self, 'noaa_charts_opacity_slider'):
                self.noaa_charts_opacity_slider.setEnabled(self.show_noaa_charts_var)
            if hasattr(self, 'noaa_charts_opacity_label'):
                self.noaa_charts_opacity_label.setEnabled(self.show_noaa_charts_var)
            if not self.show_noaa_charts_var:
                if hasattr(self, 'noaa_charts_image_plot') and self.noaa_charts_image_plot is not None:
                    try:
                        self.noaa_charts_image_plot.remove()
                    except Exception:
                        pass
                    self.noaa_charts_image_plot = None
                self.canvas.draw_idle()
            else:
                self._load_and_plot_noaa_charts()

    def _update_noaa_charts_opacity(self, value):
        """Update the opacity of NOAA charts overlay."""
        self.noaa_charts_opacity = value
        if hasattr(self, 'noaa_charts_opacity_label'):
            self.noaa_charts_opacity_label.setText(f"Opacity: {value}%")
        if hasattr(self, 'noaa_charts_image_plot') and self.noaa_charts_image_plot is not None:
            try:
                self.noaa_charts_image_plot.set_alpha(value / 100.0)
                self.canvas.draw_idle()
            except Exception as e:
                print(f"Error updating NOAA charts opacity: {e}")

    def _toggle_eez_layer(self):
        """Toggle World EEZ overlay on/off."""
        if hasattr(self, 'eez_checkbox'):
            self.show_eez_var = self.eez_checkbox.isChecked()
            if hasattr(self, 'eez_opacity_slider'):
                self.eez_opacity_slider.setEnabled(self.show_eez_var)
            if hasattr(self, 'eez_opacity_label'):
                self.eez_opacity_label.setEnabled(self.show_eez_var)
            if not self.show_eez_var:
                if hasattr(self, 'eez_image_plots') and isinstance(self.eez_image_plots, list):
                    for _plot in self.eez_image_plots:
                        try:
                            _plot.remove()
                        except Exception:
                            pass
                    self.eez_image_plots = []
                if hasattr(self, 'eez_image_plot') and self.eez_image_plot is not None:
                    try:
                        self.eez_image_plot.remove()
                    except Exception:
                        pass
                    self.eez_image_plot = None
                self.canvas.draw_idle()
            else:
                self._load_and_plot_eez(force_reload=True)

    def _update_eez_opacity(self, value):
        """Update EEZ layer opacity."""
        self.eez_opacity = value
        if hasattr(self, 'eez_opacity_label'):
            self.eez_opacity_label.setText(f"Opacity: {value}%")
        if hasattr(self, 'eez_image_plots') and isinstance(self.eez_image_plots, list) and self.eez_image_plots:
            try:
                for _plot in self.eez_image_plots:
                    _plot.set_alpha(value / 100.0)
                self.canvas.draw_idle()
            except Exception as e:
                print(f"Error updating EEZ opacity: {e}")
            return
        if hasattr(self, 'eez_image_plot') and self.eez_image_plot is not None:
            try:
                self.eez_image_plot.set_alpha(value / 100.0)
                self.canvas.draw_idle()
            except Exception as e:
                print(f"Error updating EEZ opacity: {e}")

    def _load_and_plot_basemap(self, force_reload=False):
        """Load and display ArcGIS World Imagery basemap tiles."""
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            if not xlim or not ylim or xlim[0] >= xlim[1] or ylim[0] >= ylim[1]:
                return
            extent_width = xlim[1] - xlim[0]
            extent_height = ylim[1] - ylim[0]
            current_extent = (xlim[0], xlim[1], ylim[0], ylim[1])
            if not force_reload:
                if hasattr(self, '_last_basemap_extent') and self._last_basemap_extent is not None:
                    last = self._last_basemap_extent
                    lw, lh = last[1] - last[0], last[3] - last[2]
                    wc = abs(extent_width - lw) / max(extent_width, lw) if max(extent_width, lw) > 0 else 1.0
                    hc = abs(extent_height - lh) / max(extent_height, lh) if max(extent_height, lh) > 0 else 1.0
                    cx = (xlim[0] + xlim[1]) / 2
                    cy = (ylim[0] + ylim[1]) / 2
                    lcx = (last[0] + last[1]) / 2
                    lcy = (last[2] + last[3]) / 2
                    cxc = abs(cx - lcx) / extent_width if extent_width > 0 else 0
                    cyc = abs(cy - lcy) / extent_height if extent_height > 0 else 0
                    if wc < 0.15 and hc < 0.15 and cxc < 0.25 and cyc < 0.25:
                        return
            self._last_basemap_extent = current_extent
            if not GEOSPATIAL_LIBS_AVAILABLE or pyproj is None:
                return
            if hasattr(self, 'basemap_image_plots') and isinstance(self.basemap_image_plots, list):
                for _plot in self.basemap_image_plots:
                    try:
                        _plot.remove()
                    except Exception:
                        pass
                self.basemap_image_plots = []

            # Keep request bounds valid for Web Mercator and ArcGIS tiling.
            # For wrapped/world-scale views, fetch a normalized world window and then
            # draw repeated copies across longitude.
            view_lon_width = float(xlim[1] - xlim[0])
            if view_lon_width >= 360.0:
                req_lon_min, req_lon_max = -180.0, 180.0
            else:
                req_lon_min = float(xlim[0])
                req_lon_max = float(xlim[1])
                while req_lon_min < -180.0 and req_lon_max < -180.0:
                    req_lon_min += 360.0
                    req_lon_max += 360.0
                while req_lon_min > 180.0 and req_lon_max > 180.0:
                    req_lon_min -= 360.0
                    req_lon_max -= 360.0
                if req_lon_min < -180.0 or req_lon_max > 180.0:
                    req_lon_min, req_lon_max = -180.0, 180.0
            req_lat_min = max(-85.0511, min(85.0511, float(ylim[0])))
            req_lat_max = max(-85.0511, min(85.0511, float(ylim[1])))
            if req_lon_min > req_lon_max:
                req_lon_min, req_lon_max = req_lon_max, req_lon_min
            if req_lat_min > req_lat_max:
                req_lat_min, req_lat_max = req_lat_max, req_lat_min
            if abs(req_lon_max - req_lon_min) < 1e-9 or abs(req_lat_max - req_lat_min) < 1e-9:
                return

            wgs84_to_mercator = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
            x_min_merc, y_min_merc = wgs84_to_mercator.transform(req_lon_min, req_lat_min)
            x_max_merc, y_max_merc = wgs84_to_mercator.transform(req_lon_max, req_lat_max)
            extent_width_m = abs(x_max_merc - x_min_merc)
            extent_height_m = abs(y_max_merc - y_min_merc)
            min_extent = min(extent_width_m, extent_height_m)
            target_resolution = min_extent / 512
            zoom_level = int(np.round(np.clip(np.log2(156543.03392 / target_resolution), 0, 19)))
            origin_x, origin_y = -20037508.34, 20037508.34
            tile_size = 256
            res = 156543.03392 / (2 ** zoom_level)
            tile_x_min = int(np.floor((x_min_merc - origin_x) / (tile_size * res)))
            tile_x_max = int(np.ceil((x_max_merc - origin_x) / (tile_size * res)))
            tile_y_min = int(np.floor((origin_y - y_max_merc) / (tile_size * res)))
            tile_y_max = int(np.ceil((origin_y - y_min_merc) / (tile_size * res)))
            max_tile = 2 ** zoom_level
            tile_x_min = max(0, min(tile_x_min, max_tile - 1))
            tile_x_max = max(0, min(tile_x_max, max_tile - 1))
            tile_y_min = max(0, min(tile_y_min, max_tile - 1))
            tile_y_max = max(0, min(tile_y_max, max_tile - 1))
            base_url = "https://services.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}"
            tiles = []
            for ty in range(tile_y_min, tile_y_max + 1):
                for tx in range(tile_x_min, tile_x_max + 1):
                    try:
                        url = base_url.format(z=zoom_level, y=ty, x=tx)
                        req = Request(url, headers={'User-Agent': 'SAT_Planner'})
                        with urlopen(req, timeout=5) as response:
                            img_data = response.read()
                            img = Image.open(BytesIO(img_data))
                            tiles.append((tx, ty, img))
                    except Exception:
                        continue
            if not tiles:
                return
            tile_x_coords = [tx for tx, ty, _ in tiles]
            tile_y_coords = [ty for tx, ty, _ in tiles]
            min_tx, max_tx = min(tile_x_coords), max(tile_x_coords)
            min_ty, max_ty = min(tile_y_coords), max(tile_y_coords)
            composite_width = (max_tx - min_tx + 1) * tile_size
            composite_height = (max_ty - min_ty + 1) * tile_size
            composite = Image.new('RGB', (composite_width, composite_height))
            for tx, ty, img in tiles:
                x_off = (tx - min_tx) * tile_size
                y_off = (ty - min_ty) * tile_size
                composite.paste(img, (x_off, y_off))
            composite_array_3857 = np.array(composite)
            tile_x_min_merc = origin_x + min_tx * tile_size * res
            tile_x_max_merc = origin_x + (max_tx + 1) * tile_size * res
            tile_y_max_merc = origin_y - min_ty * tile_size * res
            tile_y_min_merc = origin_y - (max_ty + 1) * tile_size * res

            # Reproject stitched basemap from Web Mercator (EPSG:3857) to WGS84 (EPSG:4326)
            # so large-area views align correctly in lon/lat axes.
            src_transform = transform.from_bounds(
                tile_x_min_merc, tile_y_min_merc, tile_x_max_merc, tile_y_max_merc,
                composite_width, composite_height
            )
            # Destination grid follows the normalized request window used for tile fetch.
            dst_width = max(256, int(composite_width))
            dst_height = max(256, int(composite_height))
            dst_transform = transform.from_bounds(
                req_lon_min, req_lat_min, req_lon_max, req_lat_max,
                dst_width, dst_height
            )
            reprojected_bands = []
            num_bands = composite_array_3857.shape[2] if len(composite_array_3857.shape) == 3 else 1
            for band_idx in range(num_bands):
                source_band = composite_array_3857[:, :, band_idx] if len(composite_array_3857.shape) == 3 else composite_array_3857
                destination_band = np.zeros((dst_height, dst_width), dtype=source_band.dtype)
                reproject(
                    source=source_band,
                    destination=destination_band,
                    src_transform=src_transform,
                    src_crs='EPSG:3857',
                    dst_transform=dst_transform,
                    dst_crs='EPSG:4326',
                    resampling=Resampling.bilinear
                )
                reprojected_bands.append(destination_band)
            composite_array = np.stack(reprojected_bands, axis=-1) if len(reprojected_bands) > 1 else reprojected_bands[0]
            lon_min, lon_max = req_lon_min, req_lon_max
            lat_min, lat_max = req_lat_min, req_lat_max

            if hasattr(self, 'basemap_image_plot') and self.basemap_image_plot is not None:
                try:
                    self.basemap_image_plot.remove()
                except Exception:
                    pass
                self.basemap_image_plot = None
            base_extent = [lon_min, lon_max, lat_min, lat_max]

            # Draw enough wrapped copies to cover current visible longitude range.
            # This keeps imagery aligned when zoomed out or across the dateline.
            k_min = int(np.floor((xlim[0] - base_extent[1]) / 360.0))
            k_max = int(np.ceil((xlim[1] - base_extent[0]) / 360.0))
            extents_to_draw = []
            for k in range(k_min, k_max + 1):
                extents_to_draw.append([
                    base_extent[0] + 360.0 * k,
                    base_extent[1] + 360.0 * k,
                    base_extent[2],
                    base_extent[3]
                ])

            self.basemap_image_plots = []
            for extent_item in extents_to_draw:
                plot_item = self.ax.imshow(
                    composite_array,
                    extent=extent_item,
                    origin='upper',
                    zorder=-3,
                    interpolation='bilinear')
                self.basemap_image_plots.append(plot_item)
            self.basemap_image_plot = self.basemap_image_plots[0] if self.basemap_image_plots else None
            self.canvas.draw_idle()
        except Exception as e:
            print(f"Error loading basemap: {e}")
            traceback.print_exc()

    def _load_and_plot_noaa_charts(self, force_reload=False):
        """Load and display NOAA ENC Charts using ESRI REST MapServer export endpoint."""
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            if not xlim or not ylim or xlim[0] >= xlim[1] or ylim[0] >= ylim[1]:
                return
            extent_width = xlim[1] - xlim[0]
            extent_height = ylim[1] - ylim[0]
            current_extent = (xlim[0], xlim[1], ylim[0], ylim[1])
            if not force_reload:
                if hasattr(self, '_last_noaa_charts_extent') and self._last_noaa_charts_extent is not None:
                    last = self._last_noaa_charts_extent
                    lw, lh = last[1] - last[0], last[3] - last[2]
                    wc = abs(extent_width - lw) / max(extent_width, lw) if max(extent_width, lw) > 0 else 1.0
                    hc = abs(extent_height - lh) / max(extent_height, lh) if max(extent_height, lh) > 0 else 1.0
                    cx = (xlim[0] + xlim[1]) / 2
                    cy = (ylim[0] + ylim[1]) / 2
                    lcx, lcy = (last[0] + last[1]) / 2, (last[2] + last[3]) / 2
                    cxc = abs(cx - lcx) / extent_width if extent_width > 0 else 0
                    cyc = abs(cy - lcy) / extent_height if extent_height > 0 else 0
                    if wc < 0.15 and hc < 0.15 and cxc < 0.25 and cyc < 0.25:
                        return
            self._last_noaa_charts_extent = current_extent
            if not GEOSPATIAL_LIBS_AVAILABLE or pyproj is None:
                return
            wgs84_to_mercator = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
            width = xlim[1] - xlim[0]
            height = ylim[1] - ylim[0]
            base_padding = 0.15
            padded_xlim = (xlim[0] - width * base_padding, xlim[1] + width * base_padding)
            padded_ylim = (ylim[0] - height * base_padding, ylim[1] + height * base_padding)
            x_min_merc, y_min_merc = wgs84_to_mercator.transform(padded_xlim[0], padded_ylim[0])
            x_max_merc, y_max_merc = wgs84_to_mercator.transform(padded_xlim[1], padded_ylim[1])
            bbox_width_merc = x_max_merc - x_min_merc
            bbox_height_merc = y_max_merc - y_min_merc
            max_dimension_pixels = 2048
            if bbox_width_merc >= bbox_height_merc:
                cell_size = bbox_width_merc / max_dimension_pixels
                image_width = max_dimension_pixels
                image_height = int(bbox_height_merc / cell_size)
            else:
                cell_size = bbox_height_merc / max_dimension_pixels
                image_height = max_dimension_pixels
                image_width = int(bbox_width_merc / cell_size)
            if image_width > max_dimension_pixels:
                image_width = max_dimension_pixels
                cell_size = bbox_width_merc / image_width
                image_height = int(bbox_height_merc / cell_size)
            if image_height > max_dimension_pixels:
                image_height = max_dimension_pixels
                cell_size = bbox_height_merc / image_height
                image_width = int(bbox_width_merc / cell_size)
            image_width = max(image_width, 256)
            image_height = max(image_height, 256)
            bbox_3857 = f"{x_min_merc},{y_min_merc},{x_max_merc},{y_max_merc}"
            base_url = "https://gis.charttools.noaa.gov/arcgis/rest/services/MCS/NOAAChartDisplay/MapServer/exts/MaritimeChartService/MapServer/export"
            params = {
                'bbox': bbox_3857, 'bboxSR': '3857', 'imageSR': '3857',
                'size': f"{image_width},{image_height}", 'format': 'png',
                'transparent': 'true', 'f': 'image', 'layers': 'show:0,1,2,3,4,5,6,7'
            }
            url = f"{base_url}?{'&'.join([f'{k}={v}' for k, v in params.items()])}"
            req = Request(url, headers={'User-Agent': 'SAT_Planner'})
            with urlopen(req, timeout=10) as response:
                img_data = response.read()
            if len(img_data) < 8 or img_data[:8] != b'\x89PNG\r\n\x1a\n':
                return
            img = Image.open(BytesIO(img_data))
            img_array_3857 = np.array(img)
            if len(img_array_3857.shape) == 2:
                img_array_3857 = np.stack([img_array_3857] * 3, axis=-1)
            elif len(img_array_3857.shape) == 3 and img_array_3857.shape[2] == 4:
                if img_array_3857.dtype != np.uint8:
                    img_array_3857 = (img_array_3857 * 255).astype(np.uint8) if img_array_3857.max() <= 1.0 else img_array_3857.astype(np.uint8)
            elif len(img_array_3857.shape) == 3 and img_array_3857.shape[2] == 3:
                if img_array_3857.dtype != np.uint8:
                    img_array_3857 = (img_array_3857 * 255).astype(np.uint8) if img_array_3857.max() <= 1.0 else img_array_3857.astype(np.uint8)
            src_transform = transform.from_bounds(
                x_min_merc, y_min_merc, x_max_merc, y_max_merc,
                image_width, image_height)
            dst_transform = transform.from_bounds(
                padded_xlim[0], padded_ylim[0], padded_xlim[1], padded_ylim[1],
                image_width, image_height)
            num_bands = img_array_3857.shape[2] if len(img_array_3857.shape) == 3 else 1
            reprojected_bands = []
            for band_idx in range(num_bands):
                source_band = img_array_3857[:, :, band_idx] if len(img_array_3857.shape) == 3 else img_array_3857
                destination_band = np.zeros((image_height, image_width), dtype=source_band.dtype)
                reproject(
                    source=source_band, destination=destination_band,
                    src_transform=src_transform, src_crs='EPSG:3857',
                    dst_transform=dst_transform, dst_crs='EPSG:4326',
                    resampling=Resampling.bilinear)
                reprojected_bands.append(destination_band)
            img_array = np.stack(reprojected_bands, axis=-1) if len(reprojected_bands) > 1 else reprojected_bands[0]
            if img_array.dtype != np.uint8:
                img_array = img_array.astype(np.uint8)
            if hasattr(self, 'noaa_charts_image_plot') and self.noaa_charts_image_plot is not None:
                try:
                    self.noaa_charts_image_plot.remove()
                except Exception:
                    pass
            current_aspect = self.ax.get_aspect()
            current_xlim_before = self.ax.get_xlim()
            current_ylim_before = self.ax.get_ylim()
            alpha = self.noaa_charts_opacity / 100.0 if hasattr(self, 'noaa_charts_opacity') else 0.5
            self.noaa_charts_image_plot = self.ax.imshow(
                img_array,
                extent=[padded_xlim[0], padded_xlim[1], padded_ylim[0], padded_ylim[1]],
                origin='upper', zorder=10, alpha=alpha, interpolation='bilinear')
            current_xlim_after = self.ax.get_xlim()
            current_ylim_after = self.ax.get_ylim()
            current_aspect_after = self.ax.get_aspect()
            if (current_xlim_before != current_xlim_after or current_ylim_before != current_ylim_after or current_aspect != current_aspect_after):
                self.ax.set_xlim(current_xlim_before)
                self.ax.set_ylim(current_ylim_before)
                if current_aspect != 'auto':
                    self.ax.set_aspect(current_aspect, adjustable='datalim')
            self.canvas.draw_idle()
        except Exception as e:
            print(f"Error loading NOAA charts: {e}")
            traceback.print_exc()

    def _load_and_plot_eez(self, force_reload=False):
        """Load and display World EEZ polygons from ArcGIS REST as a transparent overlay image."""
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            if not xlim or not ylim or xlim[0] >= xlim[1] or ylim[0] >= ylim[1]:
                return
            extent_width = xlim[1] - xlim[0]
            extent_height = ylim[1] - ylim[0]
            current_extent = (xlim[0], xlim[1], ylim[0], ylim[1])
            if not force_reload:
                if hasattr(self, '_last_eez_extent') and self._last_eez_extent is not None:
                    last = self._last_eez_extent
                    lw, lh = last[1] - last[0], last[3] - last[2]
                    wc = abs(extent_width - lw) / max(extent_width, lw) if max(extent_width, lw) > 0 else 1.0
                    hc = abs(extent_height - lh) / max(extent_height, lh) if max(extent_height, lh) > 0 else 1.0
                    cx = (xlim[0] + xlim[1]) / 2
                    cy = (ylim[0] + ylim[1]) / 2
                    lcx, lcy = (last[0] + last[1]) / 2, (last[2] + last[3]) / 2
                    cxc = abs(cx - lcx) / extent_width if extent_width > 0 else 0
                    cyc = abs(cy - lcy) / extent_height if extent_height > 0 else 0
                    if wc < 0.15 and hc < 0.15 and cxc < 0.25 and cyc < 0.25:
                        return
            self._last_eez_extent = current_extent
            width = xlim[1] - xlim[0]
            height = ylim[1] - ylim[0]
            base_padding = 0.12
            padded_xlim = (xlim[0] - width * base_padding, xlim[1] + width * base_padding)
            padded_ylim = (ylim[0] - height * base_padding, ylim[1] + height * base_padding)

            # Normalize/clamp request bounds to valid WGS84 domain so the server bbox
            # and the plotted image extent stay consistent at large/global zoom scales.
            req_lon_min = max(-180.0, min(180.0, float(padded_xlim[0])))
            req_lon_max = max(-180.0, min(180.0, float(padded_xlim[1])))
            req_lat_min = max(-90.0, min(90.0, float(padded_ylim[0])))
            req_lat_max = max(-90.0, min(90.0, float(padded_ylim[1])))

            # Keep bbox ordered and non-degenerate for export requests.
            if req_lon_min > req_lon_max:
                req_lon_min, req_lon_max = req_lon_max, req_lon_min
            if req_lat_min > req_lat_max:
                req_lat_min, req_lat_max = req_lat_max, req_lat_min
            if abs(req_lon_max - req_lon_min) < 1e-9 or abs(req_lat_max - req_lat_min) < 1e-9:
                return

            request_extent = (req_lon_min, req_lon_max, req_lat_min, req_lat_max)
            width_px = 1400
            height_px = max(600, int(width_px * ((request_extent[3] - request_extent[2]) / max((request_extent[1] - request_extent[0]), 1e-9))))
            height_px = min(height_px, 2048)

            bbox_4326 = f"{request_extent[0]},{request_extent[2]},{request_extent[1]},{request_extent[3]}"
            base_url = "https://gis.ccom.unh.edu/server/rest/services/Global/World_EEZ_V12_20231025/MapServer/export"
            params = {
                'bbox': bbox_4326,
                'bboxSR': '4326',
                'imageSR': '4326',
                'size': f"{width_px},{height_px}",
                'format': 'png32',
                'transparent': 'true',
                'f': 'image',
                'layers': 'show:0'
            }
            url = f"{base_url}?{'&'.join([f'{k}={v}' for k, v in params.items()])}"
            req = Request(url, headers={'User-Agent': 'SAT_Planner'})
            with urlopen(req, timeout=10) as response:
                img_data = response.read()
            if len(img_data) < 8 or img_data[:8] != b'\x89PNG\r\n\x1a\n':
                return
            img = Image.open(BytesIO(img_data))
            img_array = np.array(img)
            if hasattr(self, 'eez_image_plots') and isinstance(self.eez_image_plots, list):
                for _plot in self.eez_image_plots:
                    try:
                        _plot.remove()
                    except Exception:
                        pass
                self.eez_image_plots = []
            if hasattr(self, 'eez_image_plot') and self.eez_image_plot is not None:
                try:
                    self.eez_image_plot.remove()
                except Exception:
                    pass
            current_aspect = self.ax.get_aspect()
            current_xlim_before = self.ax.get_xlim()
            current_ylim_before = self.ax.get_ylim()
            alpha = self.eez_opacity / 100.0 if hasattr(self, 'eez_opacity') else 0.5

            # Handle wrapped longitude views at large/global scale by drawing shifted copies.
            base_extent = [request_extent[0], request_extent[1], request_extent[2], request_extent[3]]
            extents_to_draw = [base_extent]
            if current_xlim_before[0] < -180.0:
                extents_to_draw.append([base_extent[0] - 360.0, base_extent[1] - 360.0, base_extent[2], base_extent[3]])
            if current_xlim_before[1] > 180.0:
                extents_to_draw.append([base_extent[0] + 360.0, base_extent[1] + 360.0, base_extent[2], base_extent[3]])

            self.eez_image_plots = []
            for extent_item in extents_to_draw:
                plot_item = self.ax.imshow(
                    img_array,
                    extent=extent_item,
                    origin='upper',
                    zorder=-0.8,
                    alpha=alpha,
                    interpolation='bilinear'
                )
                self.eez_image_plots.append(plot_item)
            self.eez_image_plot = self.eez_image_plots[0] if self.eez_image_plots else None
            if self.ax.get_xlim() != current_xlim_before or self.ax.get_ylim() != current_ylim_before:
                self.ax.set_xlim(current_xlim_before)
                self.ax.set_ylim(current_ylim_before)
                if current_aspect != 'auto':
                    self.ax.set_aspect(current_aspect, adjustable='datalim')
            self.canvas.draw_idle()
        except Exception as e:
            print(f"Error loading EEZ layer: {e}")

    def _get_eez_geoname_at_point(self, lat, lon):
        """Return EEZ GEONAME for a point in WGS84, or None."""
        try:
            if lat is None or lon is None:
                return None
            # Tiny cache to avoid repeated lookups on the same click point.
            cache_key = (round(float(lat), 5), round(float(lon), 5))
            if not hasattr(self, '_eez_point_query_cache'):
                self._eez_point_query_cache = {}
            if cache_key in self._eez_point_query_cache:
                return self._eez_point_query_cache[cache_key]

            base_url = "https://gis.ccom.unh.edu/server/rest/services/Global/World_EEZ_V12_20231025/MapServer/0/query"
            params = {
                "f": "json",
                "geometry": f"{float(lon)},{float(lat)}",
                "geometryType": "esriGeometryPoint",
                "inSR": "4326",
                "spatialRel": "esriSpatialRelIntersects",
                "outFields": "GEONAME",
                "returnGeometry": "false",
                "where": "1=1",
            }
            url = f"{base_url}?{urlencode(params)}"
            req = Request(url, headers={"User-Agent": "SAT_Planner"})
            with urlopen(req, timeout=6) as response:
                payload = response.read()
            data = json.loads(payload.decode("utf-8"))
            features = data.get("features") if isinstance(data, dict) else None
            geoname = None
            if isinstance(features, list) and features:
                attrs = features[0].get("attributes", {}) if isinstance(features[0], dict) else {}
                geoname = attrs.get("GEONAME") if isinstance(attrs, dict) else None
                if geoname is not None:
                    geoname = str(geoname).strip() or None

            # Keep cache small and bounded.
            if len(self._eez_point_query_cache) > 300:
                self._eez_point_query_cache.clear()
            self._eez_point_query_cache[cache_key] = geoname
            return geoname
        except Exception:
            return None

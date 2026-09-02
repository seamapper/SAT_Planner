"""
GMRT GridServer download and tiling helpers.

GMRT uses Lamont's custom REST APIs (not ArcGIS ImageServer):
  - Download: https://www.gmrt.org/services/GridServer
  - Preview:  https://www.gmrt.org/services/ImageServer
"""
import math
import os
import tempfile
import time
from io import BytesIO
from urllib.parse import urlencode

import numpy as np
import rasterio
import requests
from PyQt6.QtCore import QThread, pyqtSignal
from rasterio.crs import CRS
from rasterio.transform import from_bounds
from rasterio.warp import Resampling, reproject

GMRT_GRID_URL = "https://www.gmrt.org/services/GridServer"
GMRT_IMAGE_URL = "https://www.gmrt.org/services/ImageServer"
MAX_TILES_PER_DOWNLOAD = 253
TILE_REQUEST_TIMEOUT = 120
TILE_INTER_DELAY_SEC = 1.0
TILE_SIZE_BASE_MRES = 120.0
TILE_SIZE_BASE_DEGREES = 2.0
MIN_TILE_STRIP_DEGREES = 0.1
LAT_LIMIT = 85.0

_OVERLAP_BY_MRES = {
    960: 0.0192,
    480: 0.0096,
    240: 0.0048,
    120: 0.0024,
    60: 0.0012,
}


def tile_size_degrees_for_mres(mres):
    """Tile side length in degrees (60 m -> 1°, 120 m -> 2°, etc.)."""
    return TILE_SIZE_BASE_DEGREES * (float(mres) / TILE_SIZE_BASE_MRES)


def needs_tiling_for_spans(lon_span, lat_span, mres):
    """True when lon_span + lat_span exceeds two tile sides at this resolution."""
    tile_size = tile_size_degrees_for_mres(mres)
    return (lon_span + lat_span) > (2.0 * tile_size)


def overlap_degrees_for_mres(mres):
    """Overlap in degrees (~2 cells) for the given meter resolution."""
    mres_value = float(mres)
    if mres_value in _OVERLAP_BY_MRES:
        return _OVERLAP_BY_MRES[mres_value]
    return 2.0 * (mres_value / 111000.0)


def build_degree_strips(span_start, span_end, tile_size, overlap):
    """Return (start, end) pairs for tile strips along one axis."""
    if span_end <= span_start:
        return []
    strips = []
    current = span_start
    step = max(tile_size - overlap, 0.01)
    max_iters = max(1, int(math.ceil((span_end - span_start) / step)) + 2)
    iteration = 0
    while current < span_end and iteration < max_iters:
        tile_end = min(current + tile_size, span_end)
        if tile_end - current >= MIN_TILE_STRIP_DEGREES:
            strips.append((current, tile_end))
            current = tile_end - overlap
        else:
            break
        iteration += 1
    return strips


def generate_gmrt_tiles(west, east, south, north, mres):
    """Generate padded tile bounds for a GMRT AOI."""
    west = max(-180.0, min(180.0, west))
    east = max(-180.0, min(180.0, east))
    south = max(-LAT_LIMIT, min(LAT_LIMIT, south))
    north = max(-LAT_LIMIT, min(LAT_LIMIT, north))
    overlap = overlap_degrees_for_mres(mres)
    tile_size = tile_size_degrees_for_mres(mres)
    lon_tiles = build_degree_strips(west, east, tile_size, overlap)
    lat_tiles = build_degree_strips(south, north, tile_size, overlap)
    tiles = []
    for tw, te in lon_tiles:
        for ts, tn in lat_tiles:
            pw, pe, ps, pn = tw, te, ts, tn
            if pw > -180.0:
                pw = max(pw - overlap, -180.0)
            if pe < 180.0:
                pe = min(pe + overlap, 180.0)
            if ps > -LAT_LIMIT:
                ps = max(ps - overlap, -LAT_LIMIT)
            if pn < LAT_LIMIT:
                pn = min(pn + overlap, LAT_LIMIT)
            tiles.append((pw, pe, ps, pn))
    return tiles


def clamp_gmrt_bbox(west, south, east, north):
    """Clamp a GCS bbox to GMRT-supported limits."""
    west = max(-180.0, min(180.0, west))
    east = max(-180.0, min(180.0, east))
    south = max(-LAT_LIMIT, min(LAT_LIMIT, south))
    north = max(-LAT_LIMIT, min(LAT_LIMIT, north))
    return west, south, east, north


def estimate_gmrt_pixels(west, south, east, north, mres_meters):
    """Approximate output pixel dimensions for a GMRT AOI."""
    meters_per_deg = 111320.0
    width = int(((east - west) * meters_per_deg) / float(mres_meters))
    height = int(((north - south) * meters_per_deg) / float(mres_meters))
    return max(width, 1), max(height, 1)


def _download_grid_to_path(params, output_path, status_cb=None):
    """Download one GridServer GeoTIFF to output_path. Raises on failure."""
    request_params = {**params, "format": "geotiff"}
    url = f"{GMRT_GRID_URL}?{urlencode(request_params)}"
    if status_cb:
        status_cb(f"GMRT GridServer: {url}")
    with requests.get(GMRT_GRID_URL, params=request_params, stream=True, timeout=TILE_REQUEST_TIMEOUT) as response:
        if response.status_code != 200:
            detail = f"HTTP {response.status_code}"
            try:
                text = response.text.lower()
                if "invalid bounds" in text or "w/e/s/n" in text:
                    detail = "Invalid W/E/S/N bounds"
                elif "invalid resolution" in text:
                    detail = "Invalid resolution"
                elif "invalid layer" in text:
                    detail = "Invalid layer"
                elif response.status_code == 413:
                    detail = "Request too large for this resolution — enable tiling or coarsen cell size"
                elif response.status_code == 404:
                    detail = "No data returned for this area"
            except Exception:
                pass
            raise RuntimeError(f"GMRT GridServer error: {detail}")
        temp = tempfile.NamedTemporaryFile(suffix=".tif", delete=False)
        temp_path = temp.name
        temp.close()
        total = 0
        with open(temp_path, "wb") as out:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    out.write(chunk)
                    total += len(chunk)
        if total == 0:
            try:
                os.remove(temp_path)
            except OSError:
                pass
            raise RuntimeError("GMRT GridServer returned an empty GeoTIFF")
        os.replace(temp_path, output_path)
    return output_path


def _mosaic_geotiffs(tile_paths, output_path, bounds_4326, status_cb=None):
    """Mosaic GeoTIFF tiles with shallower-wins compositing into output_path."""
    if status_cb:
        status_cb(f"Mosaicking {len(tile_paths)} GMRT tiles...")
    datasets = []
    try:
        for path in tile_paths:
            if os.path.exists(path):
                datasets.append(rasterio.open(path))
        if not datasets:
            raise RuntimeError("No valid GMRT tiles to mosaic")

        cell_sizes = []
        for dataset in datasets:
            transform = dataset.transform
            bounds = dataset.bounds
            cell_size_x = abs(transform.a)
            cell_size_y = abs(transform.e)
            if cell_size_x <= 0 or cell_size_y <= 0:
                cell_size_x = (bounds.right - bounds.left) / max(dataset.width, 1)
                cell_size_y = abs(bounds.top - bounds.bottom) / max(dataset.height, 1)
            cell_sizes.append((cell_size_x, cell_size_y))

        min_cell_x = min(cs[0] for cs in cell_sizes)
        min_cell_y = min(cs[1] for cs in cell_sizes)
        min_x, min_y, max_x, max_y = bounds_4326
        width = max(1, int((max_x - min_x) / min_cell_x))
        height = max(1, int((max_y - min_y) / min_cell_y))
        output_transform = from_bounds(min_x, min_y, max_x, max_y, width, height)
        output_array = np.full((height, width), -99999, dtype=np.float32)
        dst_crs = datasets[0].crs or CRS.from_epsg(4326)

        for dataset in datasets:
            data = dataset.read(1).astype(np.float32)
            data = np.where(np.isnan(data) | np.isinf(data), -99999, data)
            tile_output = np.full((height, width), -99999, dtype=np.float32)
            reproject(
                source=data,
                destination=tile_output,
                src_transform=dataset.transform,
                src_crs=dataset.crs or dst_crs,
                dst_transform=output_transform,
                dst_crs=dst_crs,
                resampling=Resampling.nearest,
                src_nodata=-99999,
                dst_nodata=-99999,
            )
            valid_mask = (tile_output != -99999) & ~np.isnan(tile_output) & ~np.isinf(tile_output)
            if np.any(valid_mask):
                update_mask = valid_mask & ((output_array == -99999) | (tile_output > output_array))
                output_array[update_mask] = tile_output[update_mask]

        # Clamp absurd values to nodata
        output_array[(output_array > 9000) | (output_array < -12000)] = -99999

        profile = {
            "driver": "GTiff",
            "height": height,
            "width": width,
            "count": 1,
            "dtype": "float32",
            "crs": CRS.from_epsg(4326),
            "transform": output_transform,
            "compress": "lzw",
            "tiled": True,
            "nodata": -99999,
        }
        with rasterio.open(output_path, "w", **profile) as dst:
            dst.write(output_array, 1)
    finally:
        for dataset in datasets:
            try:
                dataset.close()
            except Exception:
                pass
    return output_path


class GMRTDownloader(QThread):
    """Download a GMRT GridServer GeoTIFF for a GCS bbox (with optional tiling)."""

    progress = pyqtSignal(int)
    status = pyqtSignal(str)
    finished = pyqtSignal(str)
    error = pyqtSignal(str)

    def __init__(self, bbox_4326, output_path, gmrt_layer="topo", mresolution=480,
                 use_tile_download=True):
        super().__init__()
        self.bbox = clamp_gmrt_bbox(*bbox_4326)
        self.output_path = output_path
        self.gmrt_layer = gmrt_layer
        self.mresolution = float(mresolution)
        self.use_tile_download = use_tile_download
        self.cancelled = False

    def cancel(self):
        self.cancelled = True

    def run(self):
        try:
            west, south, east, north = self.bbox
            if east <= west or north <= south:
                self.error.emit("Invalid GMRT bounds: East must be > West and North > South.")
                return

            lon_span = east - west
            lat_span = north - south
            do_tiles = self.use_tile_download and needs_tiling_for_spans(lon_span, lat_span, self.mresolution)

            if do_tiles:
                tiles = generate_gmrt_tiles(west, east, south, north, self.mresolution)
                if len(tiles) > MAX_TILES_PER_DOWNLOAD:
                    self.error.emit(
                        f"GMRT download would require {len(tiles)} tiles "
                        f"(limit {MAX_TILES_PER_DOWNLOAD}). Select a smaller area "
                        f"or a coarser cell size."
                    )
                    return
                self.status.emit(
                    f"GMRT tiled download: {len(tiles)} tiles at {self.mresolution:g} m "
                    f"(layer={self.gmrt_layer})"
                )
                temp_dir = tempfile.mkdtemp(prefix="gmrt_tiles_")
                tile_paths = []
                try:
                    for i, (tw, te, ts, tn) in enumerate(tiles):
                        if self.cancelled:
                            return
                        params = {
                            "west": tw,
                            "east": te,
                            "south": ts,
                            "north": tn,
                            "layer": self.gmrt_layer,
                            "mresolution": self.mresolution,
                        }
                        tile_path = os.path.join(temp_dir, f"tile_{i:04d}.tif")
                        self.status.emit(f"Downloading GMRT tile {i + 1}/{len(tiles)}...")
                        self.progress.emit(int(5 + (i / max(len(tiles), 1)) * 80))
                        _download_grid_to_path(params, tile_path, status_cb=self.status.emit)
                        tile_paths.append(tile_path)
                        if i < len(tiles) - 1:
                            time.sleep(TILE_INTER_DELAY_SEC)
                    if self.cancelled:
                        return
                    self.progress.emit(90)
                    _mosaic_geotiffs(
                        tile_paths,
                        self.output_path,
                        (west, south, east, north),
                        status_cb=self.status.emit,
                    )
                finally:
                    for path in tile_paths:
                        try:
                            os.remove(path)
                        except OSError:
                            pass
                    try:
                        os.rmdir(temp_dir)
                    except OSError:
                        pass
            else:
                params = {
                    "west": west,
                    "east": east,
                    "south": south,
                    "north": north,
                    "layer": self.gmrt_layer,
                    "mresolution": self.mresolution,
                }
                self.status.emit(
                    f"GMRT single download at {self.mresolution:g} m (layer={self.gmrt_layer})"
                )
                self.progress.emit(10)
                _download_grid_to_path(params, self.output_path, status_cb=self.status.emit)

            if self.cancelled:
                return
            self.progress.emit(100)
            self.status.emit(f"GMRT download complete: {self.output_path}")
            self.finished.emit(self.output_path)
        except Exception as exc:
            self.error.emit(str(exc))


class GMRTMapLoader(QThread):
    """Load a GMRT ImageServer JPEG preview for a GCS extent."""

    tileLoaded = pyqtSignal(object, float, float, float, float)  # QPixmap, west, south, east, north
    statusMessage = pyqtSignal(str)

    def __init__(self, bbox_4326, size, mask=False, purpose="GMRT preview"):
        super().__init__()
        self.bbox = clamp_gmrt_bbox(*bbox_4326)
        self.size = size
        self.mask = mask
        self.purpose = purpose

    def run(self):
        from PyQt6.QtGui import QPixmap, QImage

        try:
            west, south, east, north = self.bbox
            width, height = self.size
            params = {
                "minlongitude": west,
                "maxlongitude": east,
                "minlatitude": south,
                "maxlatitude": north,
                "width": max(int(width), 1),
                "mask": "1" if self.mask else "0",
            }
            # ImageServer uses width; approximate height via aspect if needed by client scaling
            url = f"{GMRT_IMAGE_URL}?{urlencode(params)}"
            self.statusMessage.emit(f"Map REST ({self.purpose}): {url}")
            response = requests.get(GMRT_IMAGE_URL, params=params, timeout=30)
            response.raise_for_status()
            image = QImage()
            if not image.loadFromData(response.content):
                self.tileLoaded.emit(QPixmap(), west, south, east, north)
                return
            pixmap = QPixmap.fromImage(image)
            # Keep native aspect (Web Mercator Y). MapWidget will KeepAspectRatio-scale.
            self.tileLoaded.emit(pixmap, west, south, east, north)
        except Exception as exc:
            print(f"GMRT map load error: {exc}")
            west, south, east, north = self.bbox
            from PyQt6.QtGui import QPixmap
            self.tileLoaded.emit(QPixmap(), west, south, east, north)

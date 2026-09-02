"""
Module for downloading bathymetry data from ArcGIS ImageServer and creating GeoTIFF files.

License: BSD 3-Clause License
Copyright (c) 2025, Center for Coastal and Ocean Mapping, University of New Hampshire
All rights reserved.

See LICENSE file for full license text.
"""
import requests
import numpy as np
import rasterio
from rasterio.transform import from_bounds
from rasterio.crs import CRS
from io import BytesIO
from PIL import Image
import pyproj
from PyQt6.QtCore import QThread, pyqtSignal
from urllib.parse import urlencode


def _format_rest_url(url, params):
    """Build a full REST request URL for activity log display."""
    return f"{url}?{urlencode(params)}"


class BathymetryDownloader(QThread):
    """Thread for downloading bathymetry data and creating GeoTIFF."""
    
    progress = pyqtSignal(int)  # Progress percentage
    status = pyqtSignal(str)  # Status message
    finished = pyqtSignal(str)  # Output file path
    error = pyqtSignal(str)  # Error message
    
    def __init__(self, base_url, bbox, output_path, output_crs="EPSG:3857", 
                 pixel_size=None, max_size=14000, use_tile_download=False,
                 bbox_in_4326=False, pixel_size_degrees=None, tid_url=None, download_mode="combined",
                 output_requests=None, ignore_source_nodata=False):
        super().__init__()
        self.base_url = base_url
        self.bbox = bbox  # (xmin, ymin, xmax, ymax) in EPSG:3857 or 4326 when bbox_in_4326
        self.output_path = output_path  # used when output_requests is None
        self.output_crs = output_crs
        self.pixel_size = pixel_size  # meters, or None; ignored when pixel_size_degrees is set
        self.bbox_in_4326 = bbox_in_4326
        self.pixel_size_degrees = pixel_size_degrees  # for services in 4326 (e.g. GEBCO)
        self.max_size = max_size
        self.use_tile_download = use_tile_download
        self.tile_overlap = 5
        self.tile_max_size = 2000
        self.cancelled = False
        self.tid_url = tid_url  # GEBCO TID ImageServer URL for bathymetry_only / land_only
        self.download_mode = download_mode  # used when output_requests is None
        self.ignore_source_nodata = ignore_source_nodata
        # output_requests: list of (mode, path) e.g. [("combined", path1), ("bathymetry_only", path2)]
        self.output_requests = output_requests if output_requests else None
        # Force int8 for GEBCO TID grids (signed 8-bit), int16 for other GEBCO bathymetry grids
        base_url_upper = base_url.upper()
        if "GEBCO" in base_url_upper and "TID" in base_url_upper:
            self._preserve_int8 = True
            self._preserve_int16 = False
        elif "GEBCO" in base_url_upper and ("2025" in base_url or "2026" in base_url):
            self._preserve_int8 = False
            self._preserve_int16 = True
        else:
            self._preserve_int8 = False
            self._preserve_int16 = False

    def _read_exported_tiff_array(self, content):
        """Read a downloaded export TIFF, preserving integer dtypes."""
        with rasterio.open(BytesIO(content)) as src:
            original_dtype = src.dtypes[0]
            img_array = src.read(1)
            transform = src.transform if src.transform else None
            crs = src.crs if src.crs else None

            if not self._preserve_int8 and not self._preserve_int16:
                if np.issubdtype(original_dtype, np.integer):
                    self._preserve_int8 = original_dtype == np.int8
                    self._preserve_int16 = original_dtype == np.int16

            if self.ignore_source_nodata:
                if np.issubdtype(img_array.dtype, np.integer):
                    return img_array, None, transform, crs
                return img_array.astype(np.float32), None, transform, crs

            source_nodata = src.nodata
            if source_nodata is not None:
                if np.issubdtype(img_array.dtype, np.integer):
                    img_array = img_array.astype(np.float32)
                    img_array = np.where(img_array == source_nodata, np.nan, img_array)
                else:
                    img_array = np.where(img_array == source_nodata, np.nan, img_array)
            elif not np.issubdtype(img_array.dtype, np.floating):
                img_array = img_array.astype(np.float32)

            return img_array, source_nodata, transform, crs
        
    def cancel(self):
        """Cancel the download."""
        self.cancelled = True

    def _export_image_params(self, xmin, ymin, xmax, ymax, width, height):
        """Build exportImage params for raw data download."""
        params = {
            "bbox": f"{xmin},{ymin},{xmax},{ymax}",
            "size": f"{width},{height}",
            "format": "tiff",
            "f": "image",
            "interpolation": "RSP_BilinearInterpolation",
        }
        if self.bbox_in_4326:
            params["bboxSR"] = "4326"
        if not self.ignore_source_nodata:
            params["noData"] = "true"
        return params
        
    def run(self):
        """Download data and create GeoTIFF."""
        try:
            xmin, ymin, xmax, ymax = self.bbox
            
            # Determine output size
            if self.bbox_in_4326 and self.pixel_size_degrees is not None:
                # Bbox is (lon_min, lat_min, lon_max, lat_max) in degrees
                width = int((xmax - xmin) / self.pixel_size_degrees)
                height = int((ymax - ymin) / abs(self.pixel_size_degrees))
            elif self.pixel_size:
                width = int((xmax - xmin) / self.pixel_size)
                height = int((ymax - ymin) / self.pixel_size)
            else:
                width = int((xmax - xmin) / 4.0)
                height = int((ymax - ymin) / 4.0)
                
            # Check if size exceeds maximum only if tiling is disabled
            # If tiling is enabled, we can handle any size
            if not self.use_tile_download and (width > self.max_size or height > self.max_size):
                error_msg = (
                    f"Selected area exceeds maximum download size ({self.max_size} x {self.max_size} pixels).\n"
                    f"Requested size: {width} x {height} pixels.\n\n"
                    f"Please either:\n"
                    f"1. Enable Tile Download option, or\n"
                    f"2. Select a smaller area, or\n"
                    f"3. Increase the cell size (currently {self.pixel_size if self.pixel_size else 4}m)"
                )
                self.error.emit(error_msg)
                return
            
            # Check if tiling is needed
            needs_tiling = self.use_tile_download and (width > self.tile_max_size or height > self.tile_max_size)
            
            if needs_tiling:
                # Use tiled download
                result = self._download_tiled(xmin, ymin, xmax, ymax, width, height)
                if result is None or result[0] is None:
                    return  # Error already emitted
                img_array, source_nodata, transform, downloaded_crs = result
            else:
                # Single download
                self.status.emit(f"Requesting data: {width}x{height} pixels...")
                self.progress.emit(10)
                
                if self.cancelled:
                    return
                    
                # Download raw data (no raster function for raw values)
                # Try multiple approaches to get raw data
                
                url = f"{self.base_url}/exportImage"
                
                # First, try to get raw TIFF data
                params = self._export_image_params(xmin, ymin, xmax, ymax, width, height)
                self.status.emit(f"Data REST (download): {_format_rest_url(url, params)}")
                
                # Don't specify rasterFunction to get raw values
                # The service should return raw F32 values
                
                self.progress.emit(30)
                self.status.emit("Downloading data...")
                
                if self.cancelled:
                    return
                    
                img_array = None
                source_nodata = None
                transform = None
                downloaded_crs = None
                
                try:
                    response = requests.get(url, params=params, timeout=300, stream=True)
                    
                    # Check for HTTP errors before processing
                    if response.status_code == 500:
                        # Server error - provide helpful message
                        error_msg = (
                            f"Server error (500) from REST endpoint.\n\n"
                            f"The server is unable to process this request. This may be due to:\n"
                            f"  - Requested area is too large\n"
                            f"  - Server is temporarily unavailable\n"
                            f"  - Request parameters are invalid\n\n"
                            f"Try:\n"
                            f"  - Selecting a smaller area\n"
                            f"  - Using a larger cell size (8m or 16m)\n"
                            f"  - Waiting a few moments and trying again\n\n"
                            f"Requested size: {width}x{height} pixels"
                        )
                        self.error.emit(error_msg)
                        return
                    
                    response.raise_for_status()
                    
                    # Check if we got a TIFF
                    content_type = response.headers.get('Content-Type', '')
                    if 'tiff' in content_type.lower() or response.content[:4] == b'II*\x00' or response.content[:4] == b'MM\x00*':
                        # We got a TIFF, try to read it with rasterio
                        try:
                            img_array, source_nodata, transform, downloaded_crs = self._read_exported_tiff_array(
                                response.content
                            )
                        except Exception as e:
                            self.status.emit(f"Could not read TIFF with rasterio: {e}. Trying PIL...")
                            # Fall back to PIL
                            img = Image.open(BytesIO(response.content))
                            img_array = np.array(img, dtype=np.float32)
                    else:
                        # Not a TIFF, try as image
                        img = Image.open(BytesIO(response.content))
                        if img.mode in ('RGB', 'RGBA'):
                            # If we got a colored image, we might need to use pixelBlock
                            # For now, convert to grayscale and warn
                            self.status.emit("Warning: Received RGB image. Using grayscale conversion.")
                            img_array = np.array(img.convert('L'), dtype=np.float32)
                        else:
                            img_array = np.array(img, dtype=np.float32)
                        
                except requests.exceptions.HTTPError as e:
                    # HTTP error (4xx, 5xx)
                    if hasattr(e.response, 'status_code'):
                        if e.response.status_code == 500:
                            error_msg = (
                                f"Server error (500) from REST endpoint.\n\n"
                                f"The server is unable to process this request. This may be due to:\n"
                                f"  - Requested area is too large\n"
                                f"  - Server is temporarily unavailable\n"
                                f"  - Request parameters are invalid\n\n"
                                f"Try:\n"
                                f"  - Selecting a smaller area\n"
                                f"  - Using a larger cell size (8m or 16m)\n"
                                f"  - Waiting a few moments and trying again\n\n"
                                f"Requested size: {width}x{height} pixels"
                            )
                        else:
                            error_msg = f"HTTP error {e.response.status_code} from REST endpoint: {str(e)}"
                    else:
                        error_msg = f"HTTP error from REST endpoint: {str(e)}"
                    self.error.emit(error_msg)
                    return
                except requests.exceptions.RequestException as e:
                    # Connection or network error
                    error_msg = f"Network error connecting to REST endpoint: {str(e)}"
                    self.error.emit(error_msg)
                    return
                except Exception as e:
                    self.status.emit(f"Error with TIFF format: {e}. Trying PNG format...")
                    # Fall back to PNG
                    try:
                        params["format"] = "png"
                        response = requests.get(url, params=params, timeout=300, stream=True)
                        response.raise_for_status()
                        
                        img = Image.open(BytesIO(response.content))
                        if img.mode in ('RGB', 'RGBA'):
                            self.status.emit("Warning: Received RGB PNG. Using grayscale conversion.")
                            img_array = np.array(img.convert('L'), dtype=np.float32)
                        else:
                            img_array = np.array(img, dtype=np.float32)
                    except Exception as e2:
                        error_msg = f"Failed to download data: {str(e2)}"
                        self.error.emit(error_msg)
                        return
                
            # For now, if we got a PNG, we don't have the actual bathymetry values
            # We need to use the pixelBlock operation or exportImage with proper format
            # Let's try a different approach - use exportImage with pixelType and no raster function
            
            self.progress.emit(50)
            self.status.emit("Processing data...")
            
            if self.cancelled:
                return
                
            # Try to get actual pixel values using pixelBlock or different export parameters
            # For now, we'll create the GeoTIFF with what we have
            # In a production version, you'd want to use the pixelBlock operation
            
            # Ensure img_array is defined and not None
            if img_array is None:
                error_msg = "Failed to load image data from server response"
                self.error.emit(error_msg)
                return
            
            # Multi-output path (GEBCO 2025: combined + bathymetry_only + land_only + direct_measurements_only)
            if self.output_requests:
                need_tid = any(m in ("bathymetry_only", "land_only", "direct_measurements_only", "direct_unknown_measurements_only") for m, _ in self.output_requests)
                tid_array = None
                if need_tid and self.tid_url:
                    if self.cancelled:
                        return
                    self.status.emit("Downloading TID grid for masking...")
                    tid_array = self._fetch_tid_grid(xmin, ymin, xmax, ymax, width, height)
                    if tid_array is None:
                        return
                # Build one array per output (copy + mask)
                arrays_to_write = []
                for mode, path in self.output_requests:
                    arr = img_array.copy()
                    if mode in ("bathymetry_only", "land_only", "direct_measurements_only", "direct_unknown_measurements_only") and tid_array is not None:
                        if not np.issubdtype(arr.dtype, np.floating):
                            arr = arr.astype(np.float32)
                        if mode == "bathymetry_only":
                            arr[tid_array == 0] = np.nan
                        elif mode == "land_only":
                            arr[tid_array != 0] = np.nan
                        elif mode == "direct_measurements_only":
                            arr[(tid_array < 10) | (tid_array > 20)] = np.nan
                        else:  # direct_unknown_measurements_only: TID in [10, 20] or 44 or 70
                            mask = ((tid_array >= 10) & (tid_array <= 20)) | (tid_array == 44) | (tid_array == 70)
                            arr[~mask] = np.nan
                    arrays_to_write.append((mode, path, arr))
                # Determine CRS and bbox (same for all outputs)
                transform = None
                if self.bbox_in_4326 and self.output_crs == "EPSG:4326":
                    crs = CRS.from_epsg(4326)
                    output_bbox = self.bbox
                else:
                    source_crs = CRS.from_epsg(3857)
                    if self.output_crs == "EPSG:3857":
                        crs = CRS.from_epsg(3857)
                        output_bbox = self.bbox
                    elif self.output_crs == "EPSG:4326":
                        from rasterio.warp import reproject, Resampling, calculate_default_transform
                        transformer = pyproj.Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)
                        lon_min, lat_min = transformer.transform(xmin, ymin)
                        lon_max, lat_max = transformer.transform(xmax, ymax)
                        output_bbox = (lon_min, lat_min, lon_max, lat_max)
                        crs = CRS.from_epsg(4326)
                        self.status.emit("Reprojecting to WGS84...")
                        transform, output_width, output_height = calculate_default_transform(
                            source_crs, crs, width, height, *self.bbox
                        )
                        source_transform = from_bounds(*self.bbox, width, height)
                        new_arrays = []
                        for mode, path, arr in arrays_to_write:
                            reprojected_array = np.zeros((output_height, output_width), dtype=arr.dtype)
                            reproject(
                                source=arr,
                                destination=reprojected_array,
                                src_transform=source_transform,
                                src_crs=source_crs,
                                dst_transform=transform,
                                dst_crs=crs,
                                resampling=Resampling.bilinear
                            )
                            new_arrays.append((mode, path, reprojected_array))
                        arrays_to_write = new_arrays
                        width = output_width
                        height = output_height
                    else:
                        crs = CRS.from_string(self.output_crs)
                        output_bbox = self.bbox
                if transform is None:
                    transform = from_bounds(*output_bbox, width, height)
                self.progress.emit(70)
                written_paths = []
                for mode, path, arr in arrays_to_write:
                    if self.cancelled:
                        return
                    self.status.emit(f"Writing {mode}...")
                    self._write_geotiff(arr, path, width, height, transform, crs, source_nodata)
                    written_paths.append(path)
                self.progress.emit(100)
                self.status.emit(f"Saved {len(written_paths)} file(s).")
                self.finished.emit("\n".join(written_paths))
                return
            # Single-output path: apply TID mask for bathymetry_only or land_only (GEBCO 2025)
            if self.download_mode in ("bathymetry_only", "land_only") and self.tid_url:
                if self.cancelled:
                    return
                self.status.emit("Downloading TID grid for masking...")
                tid_array = self._fetch_tid_grid(xmin, ymin, xmax, ymax, width, height)
                if tid_array is None:
                    return  # Error already emitted
                # Ensure float32 for NaN masking
                if not np.issubdtype(img_array.dtype, np.floating):
                    img_array = img_array.astype(np.float32)
                if self.download_mode == "bathymetry_only":
                    # Keep only cells where TID != 0 (water/bathymetry)
                    img_array[tid_array == 0] = np.nan
                else:  # land_only
                    # Keep only cells where TID == 0 (land)
                    img_array[tid_array != 0] = np.nan
            
            # Determine CRS and bbox
            transform = None
            if self.bbox_in_4326 and self.output_crs == "EPSG:4326":
                # Data and bbox already in WGS84; no reprojection
                crs = CRS.from_epsg(4326)
                output_bbox = self.bbox
            else:
                source_crs = CRS.from_epsg(3857)
                if self.output_crs == "EPSG:3857":
                    crs = CRS.from_epsg(3857)
                    output_bbox = self.bbox
                elif self.output_crs == "EPSG:4326":
                    transformer = pyproj.Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)
                    lon_min, lat_min = transformer.transform(xmin, ymin)
                    lon_max, lat_max = transformer.transform(xmax, ymax)
                    output_bbox = (lon_min, lat_min, lon_max, lat_max)
                    crs = CRS.from_epsg(4326)
                    self.status.emit("Reprojecting to WGS84...")
                    from rasterio.warp import reproject, Resampling, calculate_default_transform
                    transform, output_width, output_height = calculate_default_transform(
                        source_crs, crs, width, height, *self.bbox
                    )
                    reprojected_array = np.zeros((output_height, output_width), dtype=img_array.dtype)
                    source_transform = from_bounds(*self.bbox, width, height)
                    reproject(
                        source=img_array,
                        destination=reprojected_array,
                        src_transform=source_transform,
                        src_crs=source_crs,
                        dst_transform=transform,
                        dst_crs=crs,
                        resampling=Resampling.bilinear
                    )
                    img_array = reprojected_array
                    width = output_width
                    height = output_height
                else:
                    crs = CRS.from_string(self.output_crs)
                    output_bbox = self.bbox
                
            # Check if we should preserve int8 (for GEBCO 2025 TID) or int16 (for GEBCO 2025)
            # This must be checked before creating the GeoTIFF
            preserve_int8 = getattr(self, '_preserve_int8', False)
            preserve_int16 = getattr(self, '_preserve_int16', False)
                
            self.progress.emit(70)
            if preserve_int8:
                self.status.emit("Creating GeoTIFF (signed 8-bit integer)...")
            elif preserve_int16:
                self.status.emit("Creating GeoTIFF (signed 16-bit integer)...")
            else:
                self.status.emit("Creating GeoTIFF...")
            
            if self.cancelled:
                return
                
            # Create transform if not already set from reprojection
            if transform is None:
                transform = from_bounds(*output_bbox, width, height)
            
            self._write_geotiff(img_array, self.output_path, width, height, transform, crs, source_nodata)
                
            self.progress.emit(100)
            if preserve_int8:
                self.status.emit(f"GeoTIFF (signed 8-bit) saved to: {self.output_path}")
            elif preserve_int16:
                self.status.emit(f"GeoTIFF (signed 16-bit) saved to: {self.output_path}")
            else:
                self.status.emit(f"GeoTIFF saved to: {self.output_path}")
            self.finished.emit(self.output_path)
            
        except Exception as e:
            self.error.emit(f"Error: {str(e)}")
            import traceback
            traceback.print_exc()
    
    def _write_geotiff(self, img_array, path, width, height, transform, crs, source_nodata=None):
        """Write a single GeoTIFF from an array. Uses self._preserve_int8 / _preserve_int16 for dtype."""
        preserve_int8 = getattr(self, '_preserve_int8', False)
        preserve_int16 = getattr(self, '_preserve_int16', False)
        if len(img_array.shape) == 3:
            img_array = img_array[:, :, 0]
        elif len(img_array.shape) != 2:
            raise ValueError(f"Unexpected array shape: {img_array.shape}")

        has_masked_nans = np.issubdtype(img_array.dtype, np.floating) and np.any(np.isnan(img_array))

        if preserve_int8:
            nodata_value = -128
            output_dtype = np.int8
        elif preserve_int16:
            nodata_value = -32768
            output_dtype = np.int16
        else:
            nodata_value = -9999.0
            output_dtype = np.float32

        # GEBCO grids use 0 as a valid elevation/TID value; never adopt source nodata=0.
        if (
            source_nodata is not None
            and source_nodata != 0
            and not self.ignore_source_nodata
        ):
            if preserve_int8 or preserve_int16:
                nodata_value = int(source_nodata)
            else:
                nodata_value = float(source_nodata)

        img_array_for_write = img_array.copy()
        if has_masked_nans:
            img_array_for_write = np.nan_to_num(img_array_for_write, nan=nodata_value)

        if preserve_int8:
            img_array_for_write = np.round(img_array_for_write)
            img_array_for_write = np.clip(img_array_for_write, -127, 127)
            img_array_for_write = img_array_for_write.astype(np.int8)
            if has_masked_nans:
                img_array_for_write[np.isnan(img_array)] = nodata_value
        elif preserve_int16:
            img_array_for_write = np.round(img_array_for_write)
            img_array_for_write = np.clip(img_array_for_write, -32767, 32767)
            img_array_for_write = img_array_for_write.astype(np.int16)
            if has_masked_nans:
                img_array_for_write[np.isnan(img_array)] = nodata_value
        else:
            if img_array_for_write.dtype != np.float32:
                img_array_for_write = img_array_for_write.astype(np.float32)

        write_kwargs = {
            "driver": "GTiff",
            "height": height,
            "width": width,
            "count": 1,
            "dtype": output_dtype,
            "crs": crs,
            "transform": transform,
            "compress": "lzw",
        }
        if has_masked_nans:
            write_kwargs["nodata"] = nodata_value

        with rasterio.open(path, 'w', **write_kwargs) as dst:
            dst.write(img_array_for_write, 1)
    
    def _download_tiled(self, xmin, ymin, xmax, ymax, total_width, total_height):
        """Download data in tiles and reassemble.
        
        Returns:
            tuple: (img_array, source_nodata, transform, downloaded_crs)
        """
        # Calculate tile grid
        tiles_x = int(np.ceil(total_width / self.tile_max_size))
        tiles_y = int(np.ceil(total_height / self.tile_max_size))
        total_tiles = tiles_x * tiles_y
        
        self.status.emit(f"Data REST (download, {total_tiles} tiles): {self.base_url.rstrip('/')}/exportImage")
        self.status.emit(f"Downloading {total_tiles} tiles ({tiles_x}x{tiles_y})...")
        self.progress.emit(10)
        
        # Calculate pixel size in world coordinates
        pixel_size_x = (xmax - xmin) / total_width
        pixel_size_y = (ymax - ymin) / total_height
        
        # Initialize output array
        # Initialize array - check if we should preserve int16
        # Start with float32 for processing, we'll convert to int16 at the end if needed
        img_array = np.full((total_height, total_width), np.nan, dtype=np.float32)
        source_nodata = None
        transform = None
        downloaded_crs = None
        
        tile_num = 0
        
        for tile_y in range(tiles_y):
            for tile_x in range(tiles_x):
                if self.cancelled:
                    return None, None, None, None
                
                tile_num += 1
                progress = 10 + int(80 * tile_num / total_tiles)
                self.progress.emit(progress)
                self.status.emit(f"Downloading tile {tile_num}/{total_tiles}...")
                
                # Calculate tile bounds in pixels
                # Add overlap to ensure no gaps
                tile_start_x = tile_x * self.tile_max_size
                tile_start_y = tile_y * self.tile_max_size
                tile_end_x = min(tile_start_x + self.tile_max_size, total_width)
                tile_end_y = min(tile_start_y + self.tile_max_size, total_height)
                
                # Add overlap (extend bounds)
                overlap_start_x = max(0, tile_start_x - self.tile_overlap)
                overlap_start_y = max(0, tile_start_y - self.tile_overlap)
                overlap_end_x = min(total_width, tile_end_x + self.tile_overlap)
                overlap_end_y = min(total_height, tile_end_y + self.tile_overlap)
                
                # Calculate world bounds for this tile
                tile_xmin = xmin + overlap_start_x * pixel_size_x
                tile_ymin = ymin + (total_height - overlap_end_y) * pixel_size_y  # Y is inverted
                tile_xmax = xmin + overlap_end_x * pixel_size_x
                tile_ymax = ymin + (total_height - overlap_start_y) * pixel_size_y
                
                tile_width = overlap_end_x - overlap_start_x
                tile_height = overlap_end_y - overlap_start_y
                
                # Download this tile
                url = f"{self.base_url}/exportImage"
                params = self._export_image_params(tile_xmin, tile_ymin, tile_xmax, tile_ymax, tile_width, tile_height)
                if tile_num == 1:
                    self.status.emit(f"Data REST (download tile example): {_format_rest_url(url, params)}")
                
                try:
                    response = requests.get(url, params=params, timeout=300, stream=True)
                    
                    if response.status_code == 500:
                        error_msg = f"Server error (500) downloading tile {tile_num}/{total_tiles}"
                        self.error.emit(error_msg)
                        return None, None, None, None
                    
                    response.raise_for_status()
                    
                    # Read tile data
                    tile_array = None
                    content_type = response.headers.get('Content-Type', '')
                    
                    if 'tiff' in content_type.lower() or response.content[:4] == b'II*\x00' or response.content[:4] == b'MM\x00*':
                        try:
                            tile_array, tile_source_nodata, tile_transform, tile_crs = self._read_exported_tiff_array(
                                response.content
                            )
                            if source_nodata is None and tile_source_nodata is not None:
                                source_nodata = tile_source_nodata
                            if transform is None and tile_transform is not None:
                                transform = tile_transform
                            if downloaded_crs is None and tile_crs is not None:
                                downloaded_crs = tile_crs
                        except Exception:
                            img = Image.open(BytesIO(response.content))
                            tile_array = np.array(img, dtype=np.float32)
                    else:
                        img = Image.open(BytesIO(response.content))
                        if img.mode in ('RGB', 'RGBA'):
                            tile_array = np.array(img.convert('L'), dtype=np.float32)
                        else:
                            tile_array = np.array(img, dtype=np.float32)
                    
                    # Determine where to place this tile in the output array
                    # Account for overlap - use the non-overlap region for placement
                    output_start_x = tile_start_x
                    output_start_y = tile_start_y
                    output_end_x = tile_end_x
                    output_end_y = tile_end_y
                    
                    # Extract the non-overlap region from the tile
                    tile_offset_x = tile_start_x - overlap_start_x
                    tile_offset_y = tile_start_y - overlap_start_y
                    tile_region_x = tile_offset_x
                    tile_region_y = tile_offset_y
                    tile_region_width = output_end_x - output_start_x
                    tile_region_height = output_end_y - output_start_y
                    
                    # Place tile data into output array
                    tile_data = tile_array[tile_region_y:tile_region_y + tile_region_height,
                                          tile_region_x:tile_region_x + tile_region_width]
                    
                    # Handle overlap regions - use average or prefer non-NaN values
                    # img_array is float32 during processing, tile_data may be float32, int8, or int16
                    output_region = img_array[output_start_y:output_end_y, output_start_x:output_end_x]
                    
                    # Convert tile_data to float32 for overlap blending when needed
                    if np.issubdtype(tile_data.dtype, np.floating):
                        tile_data_float = tile_data.copy()
                    else:
                        tile_data_float = tile_data.astype(np.float32)
                    
                    # Combine: prefer non-NaN values, average if both have values
                    mask_both = ~np.isnan(output_region) & ~np.isnan(tile_data_float)
                    mask_new_only = np.isnan(output_region) & ~np.isnan(tile_data_float)
                    mask_old_only = ~np.isnan(output_region) & np.isnan(tile_data_float)
                    
                    # Average overlapping regions
                    combined = output_region.copy()
                    combined[mask_both] = (output_region[mask_both] + tile_data_float[mask_both]) / 2.0
                    combined[mask_new_only] = tile_data_float[mask_new_only]
                    # mask_old_only: keep existing values
                    
                    img_array[output_start_y:output_end_y, output_start_x:output_end_x] = combined
                    
                except Exception as e:
                    error_msg = f"Error downloading tile {tile_num}/{total_tiles}: {str(e)}"
                    self.error.emit(error_msg)
                    return None, None, None, None
        
        self.progress.emit(90)
        self.status.emit("Reassembling tiles...")
        
        # Convert to int8 or int16 if we detected integer data
        preserve_int8 = getattr(self, '_preserve_int8', False)
        preserve_int16 = getattr(self, '_preserve_int16', False)
        if preserve_int8:
            img_array = np.nan_to_num(img_array, nan=-128)
            img_array = np.round(img_array)
            img_array = np.clip(img_array, -127, 127).astype(np.int8)
        elif preserve_int16:
            img_array = np.nan_to_num(img_array, nan=-32768)
            img_array = np.round(img_array)
            img_array = np.clip(img_array, -32767, 32767).astype(np.int16)
        
        return img_array, source_nodata, transform, downloaded_crs
    
    def _fetch_tid_grid(self, xmin, ymin, xmax, ymax, width, height):
        """Fetch TID grid from TID ImageServer (same bbox and size). Returns 2D array (int8 or float32) or None on error."""
        url = f"{self.tid_url.rstrip('/')}/exportImage"
        params = self._export_image_params(xmin, ymin, xmax, ymax, width, height)
        self.status.emit(f"Data REST (TID grid): {_format_rest_url(url, params)}")
        try:
            if width <= self.tile_max_size and height <= self.tile_max_size:
                response = requests.get(url, params=params, timeout=300, stream=True)
                response.raise_for_status()
                content_type = response.headers.get('Content-Type', '')
                if 'tiff' in content_type.lower() or response.content[:4] in (b'II*\x00', b'MM\x00*'):
                    arr, _, _, _ = self._read_exported_tiff_array(response.content)
                    return arr
                else:
                    img = Image.open(BytesIO(response.content))
                    return np.array(img.convert('L') if img.mode in ('RGB', 'RGBA') else img, dtype=np.float32)
            else:
                # Tiled fetch for large areas
                tiles_x = int(np.ceil(width / self.tile_max_size))
                tiles_y = int(np.ceil(height / self.tile_max_size))
                pixel_size_x = (xmax - xmin) / width
                pixel_size_y = (ymax - ymin) / height
                out = np.full((height, width), -128, dtype=np.int8)  # TID nodata-like fill
                for tile_y in range(tiles_y):
                    for tile_x in range(tiles_x):
                        if self.cancelled:
                            return None
                        tx0 = tile_x * self.tile_max_size
                        ty0 = tile_y * self.tile_max_size
                        tx1 = min(tx0 + self.tile_max_size, width)
                        ty1 = min(ty0 + self.tile_max_size, height)
                        bx0 = xmin + tx0 * pixel_size_x
                        bx1 = xmin + tx1 * pixel_size_x
                        # Image row 0 = north (ymax), row increases southward
                        by0 = ymax - ty0 * pixel_size_y  # north
                        by1 = ymax - ty1 * pixel_size_y  # south
                        tw, th = tx1 - tx0, ty1 - ty0
                        resp = requests.get(url, params=self._export_image_params(bx0, by0, bx1, by1, tw, th), timeout=300, stream=True)
                        resp.raise_for_status()
                        with rasterio.open(BytesIO(resp.content)) as src:
                            tile = src.read(1)
                        out[ty0:ty1, tx0:tx1] = tile[:th, :tw]
                return out
        except Exception as e:
            self.error.emit(f"Failed to fetch TID grid: {str(e)}")
            return None


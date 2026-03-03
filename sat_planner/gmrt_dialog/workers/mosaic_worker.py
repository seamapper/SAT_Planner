# Copyright (c) 2026 Paul Johnson
# SPDX-License-Identifier: BSD-3-Clause
# Embedded in SAT Planner.

"""
Worker thread for mosaicking downloaded tiles.
"""

import os
import tempfile
from datetime import datetime

from PyQt6.QtCore import QThread, pyqtSignal
import numpy as np

try:
    import rasterio
except ImportError:
    rasterio = None
try:
    import netCDF4
except ImportError:
    netCDF4 = None


class MosaicWorker(QThread):
    """Worker thread for handling all mosaicking operations"""
    progress = pyqtSignal(str)
    finished = pyqtSignal(bool, str)

    def __init__(self, downloaded_tile_files, download_dir, layer_type,
                 west_spin, south_spin,
                 east_spin, north_spin, delete_tiles_checkbox, split_checkbox, format_type=None):
        super().__init__()
        self.downloaded_tile_files = downloaded_tile_files
        self.download_dir = download_dir
        self.layer_type = layer_type
        self.west_spin = west_spin
        self.south_spin = south_spin
        self.east_spin = east_spin
        self.north_spin = north_spin
        self.delete_tiles_checkbox = delete_tiles_checkbox
        self.split_checkbox = split_checkbox
        self.format_type = format_type
        self.mosaic_path = None

    def _open_raster_file(self, tile_file):
        try:
            dataset = rasterio.open(tile_file)
            return dataset
        except Exception:
            if tile_file.lower().endswith('.nc'):
                try:
                    dataset = rasterio.open(tile_file, driver='NetCDF')
                    return dataset
                except Exception:
                    gdal_patterns = [
                        f'NETCDF:"{tile_file}":z',
                        f'NETCDF:"{tile_file}":elevation',
                        f'NETCDF:"{tile_file}":topo',
                        f'NETCDF:"{tile_file}":bathy',
                        f'NETCDF:"{tile_file}"',
                    ]
                    for gdal_path in gdal_patterns:
                        try:
                            dataset = rasterio.open(gdal_path)
                            try:
                                src_bounds = dataset.bounds
                                src_crs = dataset.crs
                                needs_conversion = False
                                if src_crs is None:
                                    needs_conversion = True
                                elif not src_crs.is_geographic:
                                    needs_conversion = True
                                elif src_crs.to_epsg() != 4326:
                                    needs_conversion = True
                                elif (src_bounds[0] < -200 or src_bounds[2] > 200 or
                                      src_bounds[1] < -100 or src_bounds[3] > 100):
                                    needs_conversion = True
                                if needs_conversion:
                                    dataset.close()
                                    break
                                else:
                                    return dataset
                            except Exception:
                                dataset.close()
                                break
                        except Exception:
                            continue
                    if netCDF4 is not None:
                        try:
                            return self._open_netcdf_with_netcdf4(tile_file)
                        except Exception:
                            pass
                    return None
            return None
        return None

    def _open_netcdf_with_netcdf4(self, tile_file):
        nc = netCDF4.Dataset(tile_file, 'r')
        if 'z' in nc.variables and 'dimension' in nc.variables and 'x_range' in nc.variables and 'y_range' in nc.variables:
            try:
                dimension = nc.variables['dimension'][:]
                if len(dimension) >= 2:
                    ncols = int(dimension[0])
                    nrows = int(dimension[1])
                else:
                    nc.close()
                    raise Exception(f"Invalid dimension array: {dimension}")
                x_range = nc.variables['x_range'][:]
                y_range = nc.variables['y_range'][:]
                lon_min = float(x_range[0])
                lon_max = float(x_range[1])
                lat_min = float(y_range[0])
                lat_max = float(y_range[1])
                z_data = nc.variables['z'][:]
                if z_data.size == nrows * ncols:
                    data = z_data.reshape((nrows, ncols), order='C')
                else:
                    nc.close()
                    raise Exception(f"Data size mismatch: got {z_data.size}, expected {nrows * ncols}")
                from rasterio.transform import from_bounds
                from rasterio.crs import CRS
                transform = from_bounds(lon_min, lat_min, lon_max, lat_max, ncols, nrows)
                crs = CRS.from_epsg(4326)
                nodata_value = -99999 if 'z_range' in nc.variables else None
                nc.close()
            except Exception as e:
                nc.close()
                raise Exception(f"Error processing GMRT 1D format: {e}")
        else:
            data_var = None
            for var_name in nc.variables:
                var = nc.variables[var_name]
                if len(var.dimensions) >= 2:
                    data_var = var_name
                    break
            if data_var is None:
                nc.close()
                raise Exception("No 2D data variable found in NetCDF file")
            data = nc.variables[data_var][:]
            if len(data.shape) == 3:
                data = data[0, :, :]
            dims = nc.variables[data_var].dimensions
            lat_var = lon_var = None
            for dim in dims:
                if dim in nc.variables:
                    var = nc.variables[dim]
                    if hasattr(var, 'standard_name'):
                        if 'lat' in var.standard_name.lower():
                            lat_var = dim
                        elif 'lon' in var.standard_name.lower():
                            lon_var = dim
            if lat_var is None:
                for name in ['lat', 'latitude', 'y']:
                    if name in nc.variables:
                        lat_var = name
                        break
            if lon_var is None:
                for name in ['lon', 'longitude', 'x']:
                    if name in nc.variables:
                        lon_var = name
                        break
            if lat_var is None or lon_var is None:
                nc.close()
                raise Exception("Could not find latitude/longitude variables in NetCDF file")
            lats = nc.variables[lat_var][:]
            lons = nc.variables[lon_var][:]
            lat_min, lat_max = float(lats.min()), float(lats.max())
            lon_min, lon_max = float(lons.min()), float(lons.max())
            from rasterio.transform import from_bounds
            transform = from_bounds(lon_min, lat_min, lon_max, lat_max, data.shape[1], data.shape[0])
            crs = None
            if hasattr(nc, 'crs') or hasattr(nc, 'spatial_ref'):
                try:
                    import rasterio.crs
                    if hasattr(nc, 'crs_wkt'):
                        crs = rasterio.crs.CRS.from_wkt(nc.crs_wkt)
                    elif hasattr(nc, 'epsg'):
                        crs = rasterio.crs.CRS.from_epsg(int(nc.epsg))
                except Exception:
                    pass
            if crs is None:
                crs = rasterio.crs.CRS.from_epsg(4326)
            nodata_value = None
            if hasattr(nc.variables[data_var], '_FillValue'):
                nodata_value = nc.variables[data_var]._FillValue
            nc.close()

        temp_tiff = tempfile.NamedTemporaryFile(suffix='.tif', delete=False)
        temp_tiff.close()
        use_tiled = (data.shape[0] % 16 == 0) and (data.shape[1] % 16 == 0)
        profile = {
            'driver': 'GTiff',
            'height': data.shape[0],
            'width': data.shape[1],
            'count': 1,
            'dtype': data.dtype,
            'crs': crs,
            'transform': transform,
            'compress': 'lzw',
            'tiled': use_tiled
        }
        if nodata_value is not None:
            profile['nodata'] = nodata_value
        if use_tiled:
            profile['blockxsize'] = min(512, data.shape[1])
            profile['blockysize'] = min(512, data.shape[0])
            profile['blockxsize'] = (profile['blockxsize'] // 16) * 16
            profile['blockysize'] = (profile['blockysize'] // 16) * 16
        with rasterio.open(temp_tiff.name, 'w', **profile) as dst:
            dst.write(data, 1)
        dataset = rasterio.open(temp_tiff.name)
        dataset._temp_file = temp_tiff.name
        return dataset

    def run(self):
        try:
            self.progress.emit("Starting mosaicking process...")
            if not self.downloaded_tile_files:
                self.finished.emit(False, "No tiles to mosaic")
                return
            self.progress.emit(f"Mosaicking {len(self.downloaded_tile_files)} tiles...")
            missing_files = []
            for f in self.downloaded_tile_files:
                if not os.path.exists(f):
                    missing_files.append(f)
                else:
                    try:
                        with open(f, 'rb') as test_file:
                            test_file.read(1)
                    except Exception:
                        missing_files.append(f)
            if missing_files:
                self.finished.emit(False, f"Missing or unreadable tile files: {missing_files}")
                return
            if rasterio is None:
                self.finished.emit(False, "rasterio not available for mosaicking")
                return
            current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
            mosaic_filename = f"gmrt_{self.layer_type}_mosaic_{current_time}.tif"
            self.mosaic_path = os.path.join(self.download_dir, mosaic_filename)
            self.progress.emit(f"Output mosaic path: {self.mosaic_path}")
            datasets = []
            for i, tile_file in enumerate(self.downloaded_tile_files):
                if os.path.exists(tile_file):
                    try:
                        dataset = self._open_raster_file(tile_file)
                        if dataset is None:
                            continue
                        datasets.append(dataset)
                    except Exception:
                        continue
            if not datasets:
                self.finished.emit(False, "No valid tiles found for mosaicking")
                return
            self.progress.emit("Using rasterio for mosaicking...")
            fallback_datasets = []
            try:
                for tile_file in self.downloaded_tile_files:
                    if os.path.exists(tile_file):
                        try:
                            dataset = self._open_raster_file(tile_file)
                            if dataset is None:
                                continue
                            fallback_datasets.append(dataset)
                        except Exception:
                            continue
                if fallback_datasets:
                    original_bounds = (self.west_spin.value(), self.south_spin.value(),
                                     self.east_spin.value(), self.north_spin.value())
                    self._mosaic_with_rasterio_fallback(fallback_datasets, self.mosaic_path, original_bounds)
                    self.finished.emit(True, "Rasterio mosaicking completed successfully")
                else:
                    self.finished.emit(False, "No valid tiles found for fallback")
            except Exception as fallback_error:
                self.finished.emit(False, f"Rasterio fallback also failed: {str(fallback_error)}")
            finally:
                for dataset in fallback_datasets:
                    try:
                        temp_file = getattr(dataset, '_temp_file', None)
                        dataset.close()
                        if temp_file and os.path.exists(temp_file):
                            try:
                                os.remove(temp_file)
                            except Exception:
                                pass
                    except Exception:
                        pass
        except Exception as e:
            import traceback
            self.finished.emit(False, f"Mosaic worker error: {str(e)}\n{traceback.format_exc()}")

    def _mosaic_with_rasterio_fallback(self, datasets, mosaic_path, original_bounds):
        from rasterio.warp import reproject, Resampling
        if not datasets:
            raise Exception("No valid datasets provided for rasterio fallback")
        cell_sizes = []
        bounds_list = []
        for i, dataset in enumerate(datasets):
            transform = dataset.transform
            bounds = dataset.bounds
            cell_size_x = abs(transform[0])
            cell_size_y = abs(transform[4])
            if cell_size_x <= 0 or cell_size_y <= 0:
                if bounds[2] > bounds[0] and bounds[3] > bounds[1] and dataset.width > 0 and dataset.height > 0:
                    cell_size_x = (bounds[2] - bounds[0]) / dataset.width
                    cell_size_y = abs((bounds[3] - bounds[1]) / dataset.height)
            cell_sizes.append((cell_size_x, cell_size_y))
            bounds_list.append(bounds)
        min_cell_x = min(cs[0] for cs in cell_sizes)
        min_cell_y = min(cs[1] for cs in cell_sizes)
        first_tile_transform = datasets[0].transform
        from rasterio.transform import xy
        first_tile_first_pixel_lon, first_tile_first_pixel_lat = xy(first_tile_transform, 0, 0)
        min_x, min_y, max_x, max_y = original_bounds
        if max_x <= min_x or max_y <= min_y:
            raise Exception(f"Invalid bounds: max_x ({max_x}) <= min_x ({min_x}) or max_y ({max_y}) <= min_y ({min_y})")
        if min_cell_x <= 0 or min_cell_y <= 0:
            raise Exception(f"Invalid cell size: min_cell_x={min_cell_x}, min_cell_y={min_cell_y}")
        width = int((max_x - min_x) / min_cell_x)
        height = int((max_y - min_y) / min_cell_y)
        if width <= 0 or height <= 0:
            raise Exception(f"Invalid dimensions calculated: width={width}, height={height}")
        offset_from_first_pixel = min_x - first_tile_first_pixel_lon
        pixels_offset = round(offset_from_first_pixel / min_cell_x)
        aligned_west = first_tile_first_pixel_lon + pixels_offset * min_cell_x
        offset_from_first_pixel_lat = max_y - first_tile_first_pixel_lat
        pixels_offset_lat = round(offset_from_first_pixel_lat / min_cell_y)
        aligned_north = first_tile_first_pixel_lat + pixels_offset_lat * min_cell_y
        aligned_east = aligned_west + width * min_cell_x
        aligned_south = aligned_north - height * min_cell_y
        output_transform = rasterio.transform.from_bounds(aligned_west, aligned_south, aligned_east, aligned_north, width, height)
        output_array = np.full((height, width), -99999, dtype=np.float32)
        for dataset in datasets:
            data = dataset.read(1)
            data = np.where(np.isnan(data), -99999, data)
            data = np.where(np.isinf(data), -99999, data)
            src_transform = dataset.transform
            src_crs = dataset.crs
            tile_output = np.full((height, width), -99999, dtype=np.float32)
            reproject(
                source=data,
                destination=tile_output,
                src_transform=src_transform,
                src_crs=src_crs,
                dst_transform=output_transform,
                dst_crs=datasets[0].crs,
                resampling=Resampling.nearest,
                src_nodata=-99999,
                dst_nodata=-99999
            )
            valid_mask = (tile_output != -99999) & ~np.isnan(tile_output) & ~np.isinf(tile_output)
            if np.any(valid_mask):
                update_mask = valid_mask & ((output_array == -99999) | (tile_output > output_array))
                output_array[update_mask] = tile_output[update_mask]
        output_array = self.validate_bathymetry_data(output_array)
        profile = {
            'driver': 'GTiff',
            'height': output_array.shape[0],
            'width': output_array.shape[1],
            'count': 1,
            'dtype': output_array.dtype,
            'crs': datasets[0].crs,
            'transform': output_transform,
            'compress': 'lzw',
            'tiled': True,
            'nodata': -99999
        }
        with rasterio.open(mosaic_path, 'w', **profile) as dst:
            dst.write(output_array, 1)
        for dataset in datasets:
            try:
                temp_file = getattr(dataset, '_temp_file', None)
                dataset.close()
                if temp_file and os.path.exists(temp_file):
                    try:
                        os.remove(temp_file)
                    except Exception:
                        pass
            except Exception:
                pass

    def validate_bathymetry_data(self, data):
        max_elevation = 9000.0
        min_elevation = -12000.0
        validated_data = data.copy()
        nodata_value = -99999
        validated_data[validated_data > max_elevation] = nodata_value
        validated_data[validated_data < min_elevation] = nodata_value
        return validated_data

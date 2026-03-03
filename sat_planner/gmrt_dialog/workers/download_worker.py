# Copyright (c) 2026 Paul Johnson
# SPDX-License-Identifier: BSD-3-Clause
# Embedded in SAT Planner.

import os
import tempfile
import shutil

from PyQt6.QtCore import QThread, pyqtSignal
import requests

from ..config import GMRT_URL


class DownloadWorker(QThread):
    """Worker thread for downloading bathymetry data files from GMRT GridServer."""
    finished = pyqtSignal(bool, str)

    def __init__(self, params, filename, requested_format=None):
        super().__init__()
        self.params = params.copy()
        self.params['format'] = 'geotiff'
        self.filename = filename
        self.requested_format = requested_format

    def get_error_message(self, status_code):
        error_messages = {
            404: "No Data Returned - The requested area may be outside available data coverage",
            413: "Request Too Large - The requested area is too large for the specified resolution"
        }
        return error_messages.get(status_code, f"HTTP Error {status_code}")

    def run(self):
        try:
            param_str = ", ".join([f"{k}={v}" for k, v in self.params.items()])
            print(f"[DownloadWorker] Downloading bathymetry grid: {param_str}")
            with requests.get(GMRT_URL, params=self.params, stream=True, timeout=120) as r:
                if r.status_code == 200:
                    temp_geotiff = tempfile.NamedTemporaryFile(suffix='.tif', delete=False)
                    temp_geotiff_path = temp_geotiff.name
                    temp_geotiff.close()
                    total_bytes = 0
                    with open(temp_geotiff_path, 'wb') as f:
                        for chunk in r.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                                total_bytes += len(chunk)
                    print(f"[DownloadWorker] GeoTIFF download completed successfully ({total_bytes} bytes)")
                    if not os.path.exists(temp_geotiff_path):
                        self.finished.emit(False, f"Temporary GeoTIFF file not found: {temp_geotiff_path}")
                        return
                    file_size = os.path.getsize(temp_geotiff_path)
                    if file_size == 0:
                        try:
                            os.remove(temp_geotiff_path)
                        except Exception:
                            pass
                        self.finished.emit(False, "Downloaded GeoTIFF file is empty")
                        return
                    try:
                        if temp_geotiff_path != self.filename:
                            shutil.move(temp_geotiff_path, self.filename)
                    except Exception as e:
                        print(f"[DownloadWorker] Error moving GeoTIFF file: {e}")
                        try:
                            shutil.copy2(temp_geotiff_path, self.filename)
                            os.remove(temp_geotiff_path)
                        except Exception as e2:
                            self.finished.emit(False, f"Error saving GeoTIFF file: {e2}")
                            return
                    print(f"[DownloadWorker] Successfully downloaded GeoTIFF file")
                    self.finished.emit(True, self.filename)
                else:
                    error_msg = self.get_error_message(r.status_code)
                    if r.status_code == 404:
                        try:
                            content = r.text.lower()
                            if "invalid output format" in content:
                                error_msg = "Invalid Output Format Specified"
                            elif "invalid layer" in content:
                                error_msg = "Invalid Layer Specified"
                            elif "invalid bounds" in content or "w/e/s/n" in content:
                                error_msg = "Invalid W/E/S/N bounds specified"
                            elif "invalid resolution" in content:
                                error_msg = "Invalid Resolution"
                            else:
                                error_msg = "No Data Returned - The requested area may be outside available data coverage"
                        except Exception:
                            pass
                    detailed_error = f"Server Error {r.status_code}: {error_msg}"
                    print(f"[DownloadWorker] Grid download failed: {detailed_error}")
                    self.finished.emit(False, detailed_error)
        except requests.exceptions.Timeout:
            print(f"[DownloadWorker] Grid download timeout")
            self.finished.emit(False, "Request Timeout - The server took too long to respond")
        except requests.exceptions.ConnectionError:
            print(f"[DownloadWorker] Grid download connection error")
            self.finished.emit(False, "Connection Error - Unable to connect to GMRT server")
        except Exception as e:
            print(f"[DownloadWorker] Grid download error: {str(e)}")
            self.finished.emit(False, f"Download Error: {str(e)}")

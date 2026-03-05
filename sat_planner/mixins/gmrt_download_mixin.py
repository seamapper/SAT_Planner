"""
GMRT bathymetry grid download: prompt for save path, download GeoTIFF from GMRT GridServer, then load.
Reusable for Calibration, Reference, and Line Planning survey imports.
Expects host to provide _show_message(), and optionally last_geotiff_dir, _load_geotiff_from_path(), set_*_info_text().
"""
import os
import datetime

from PyQt6.QtWidgets import QFileDialog
from PyQt6.QtCore import QThread, pyqtSignal

try:
    import requests
except ImportError:
    requests = None

GMRT_GRIDSERVER_URL = "https://www.gmrt.org/services/GridServer"


class GMRTDownloadWorker(QThread):
    """Worker thread to download a GeoTIFF from GMRT GridServer and save to a file."""
    finished = pyqtSignal(bool, str)  # success, path_or_error_message

    def __init__(self, params, save_path):
        super().__init__()
        self.params = dict(params)
        self.params['format'] = 'geotiff'
        self.save_path = save_path

    def run(self):
        try:
            if requests is None:
                self.finished.emit(False, "requests library is required. Install with: pip install requests")
                return
            with requests.get(GMRT_GRIDSERVER_URL, params=self.params, stream=True, timeout=120) as r:
                if r.status_code != 200:
                    msg = r.text[:200] if r.text else f"HTTP {r.status_code}"
                    self.finished.emit(False, f"GMRT server error: {msg}")
                    return
                with open(self.save_path, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                if os.path.getsize(self.save_path) == 0:
                    try:
                        os.remove(self.save_path)
                    except OSError:
                        pass
                    self.finished.emit(False, "Downloaded file is empty.")
                    return
                self.finished.emit(True, self.save_path)
        except requests.exceptions.Timeout:
            self.finished.emit(False, "GMRT request timed out.")
        except requests.exceptions.ConnectionError:
            self.finished.emit(False, "Could not connect to GMRT server.")
        except Exception as e:
            self.finished.emit(False, str(e))


class GMRTDownloadMixin:
    """
    Mixin for downloading GMRT bathymetry GeoTIFFs and loading them.
    Use _download_gmrt_and_load(west, east, south, north, ...) after computing extent from survey lines.
    """

    def _download_gmrt_and_load(
        self,
        west,
        east,
        south,
        north,
        resolution=100,
        layer="topo",
        default_filename_prefix="GMRT_Bathy",
        log_func=None,
        default_directory=None,
    ):
        """
        Prompt for save path, download GMRT GeoTIFF for the given bounds, then load it.

        Args:
            west, east, south, north: Bounds in degrees (WGS84).
            resolution: Meter resolution (default 100).
            layer: GMRT layer (default "topo").
            default_filename_prefix: Used for default filename: {prefix}_{YYYYMMDD_HHMMSS}.tif.
            log_func: Optional callable(message, append=False) for activity log. If None, uses
                      set_cal_info_text, set_ref_info_text, or set_line_info_text if present.
            default_directory: Optional directory for the save dialog (e.g. directory of imported survey file).
                              If None, uses last_geotiff_dir or user home.
        """
        if requests is None:
            self._show_message(
                "warning",
                "GMRT Download",
                "The requests library is required to download GMRT grids. Install with: pip install requests",
            )
            return

        def _log(msg, append=True):
            if log_func is not None:
                log_func(msg, append=append)
                return
            if hasattr(self, "set_cal_info_text"):
                self.set_cal_info_text(msg, append=append)
            elif hasattr(self, "set_ref_info_text"):
                self.set_ref_info_text(msg, append=append)
            elif hasattr(self, "set_line_info_text"):
                self.set_line_info_text(msg, append=append)

        default_name = f"{default_filename_prefix}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.tif"
        if default_directory and os.path.isdir(default_directory):
            start_dir = default_directory
        else:
            start_dir = getattr(self, "last_geotiff_dir", None) or os.path.expanduser("~")
        suggested_path = os.path.join(start_dir, default_name)
        save_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save GMRT GeoTIFF",
            suggested_path,
            "GeoTIFF (*.tif *.tiff);;All files (*.*)",
        )
        if not save_path:
            return
        if not save_path.lower().endswith((".tif", ".tiff")):
            save_path = save_path + ".tif"

        _log("Downloading GMRT grid...", append=True)
        params = {
            "west": west,
            "east": east,
            "south": south,
            "north": north,
            "layer": layer,
            "mresolution": resolution,
        }
        worker = GMRTDownloadWorker(params, save_path)
        worker.finished.connect(
            lambda success, path_or_error: self._on_gmrt_download_finished(
                success, path_or_error, _log
            )
        )
        self._gmrt_download_worker = worker
        worker.start()

    def _on_gmrt_download_finished(self, success, path_or_error, log_func):
        if success:
            log_func("GMRT grid downloaded. Loading GeoTIFF...", append=True)
            if hasattr(self, "_load_geotiff_from_path"):
                self._load_geotiff_from_path(path_or_error)
            log_func("GMRT grid loaded.", append=True)
            # If this was a calibration-import GMRT download, set suggested export name from pitch line (offset + heading), same as when drawing
            if (hasattr(self, "param_notebook") and self.param_notebook.currentIndex() == 0
                    and hasattr(self, "cal_export_name_entry") and not self.cal_export_name_entry.text().strip()
                    and hasattr(self, "_update_cal_line_offset_from_pitch_line") and hasattr(self, "_update_cal_export_name_from_pitch_line")):
                self._update_cal_line_offset_from_pitch_line()
                self._update_cal_export_name_from_pitch_line()
        else:
            self._show_message("error", "GMRT Download", path_or_error)
            log_func(f"GMRT download failed: {path_or_error}", append=True)
        self._gmrt_download_worker = None

"""
GMRT bathymetry grid download: prompt for save path, download GeoTIFF from GMRT GridServer, then load.
Reusable for Calibration, Reference, and Line Planning survey imports.
Expects host to provide _show_message(), and optionally last_geotiff_dir, _load_geotiff_from_path(), set_*_info_text().
"""
import os
import datetime

from PyQt6.QtWidgets import QFileDialog
from PyQt6.QtCore import QThread, pyqtSignal

from ..gmrt_split import split_topo_bathy

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
        split_topo_depths=False,
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
            split_topo_depths: When True, after a successful download the GeoTIFF is split into
                              a topography file (values >= 0) and a bathymetry file (values < 0);
                              SAT Planner loads only the bathymetry file. When False, the combined
                              GeoTIFF is loaded as-is. If the bathymetry partition is empty the
                              downloaded files are deleted and the user is warned.
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

        # Use this directory as default for export dialogs (calibration, reference, line)
        gmrt_dir = os.path.dirname(save_path)
        if gmrt_dir and os.path.isdir(gmrt_dir):
            if hasattr(self, "last_export_dir"):
                self.last_export_dir = gmrt_dir
                if hasattr(self, "_save_last_export_dir"):
                    self._save_last_export_dir()
            if hasattr(self, "last_used_dir"):
                self.last_used_dir = gmrt_dir
                if hasattr(self, "_save_last_used_dir"):
                    self._save_last_used_dir()

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
                success, path_or_error, _log, split_topo_depths
            )
        )
        self._gmrt_download_worker = worker
        worker.start()

    def _on_gmrt_download_finished(self, success, path_or_error, log_func, split_topo_depths=False):
        if success:
            log_func("GMRT grid downloaded.", append=True)
            path_to_load = self._maybe_split_gmrt_grid(path_or_error, split_topo_depths, log_func)
            if path_to_load is None:
                # Split was requested but produced no bathymetry; files have been removed and
                # the user has been warned. Do not load anything.
                self._gmrt_download_worker = None
                return
            log_func("Loading GeoTIFF...", append=True)
            if hasattr(self, "_load_geotiff_from_path"):
                self._load_geotiff_from_path(path_to_load)
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

    def _maybe_split_gmrt_grid(self, downloaded_path, split_topo_depths, log_func):
        """Apply optional topo/bathy split to a freshly downloaded GMRT GeoTIFF.

        Returns the path the caller should load, or None if nothing should be loaded
        (e.g., split requested but no bathymetry in the download). When None is returned,
        all downloaded artifacts (combined, topo, bathy) have been deleted and the user
        has been informed.
        """
        if not split_topo_depths:
            return downloaded_path

        result = split_topo_bathy(
            downloaded_path,
            log=lambda msg: log_func(msg, append=True),
        )
        if result.error:
            self._show_message(
                "warning",
                "GMRT Split",
                f"{result.error}\n\nLoading the combined GeoTIFF instead.",
            )
            log_func(f"GMRT split failed: {result.error}", append=True)
            return downloaded_path

        if result.bathy_path is None:
            # No bathymetry in the downloaded extent. Per user request: warn and
            # delete everything that was downloaded so nothing is left behind.
            removed = []
            for path in (result.topo_path, downloaded_path):
                if path and os.path.exists(path):
                    try:
                        os.remove(path)
                        removed.append(os.path.basename(path))
                    except OSError as exc:
                        log_func(
                            f"Failed to delete {os.path.basename(path)} after empty-bathy download: {exc}",
                            append=True,
                        )
            if removed:
                log_func(
                    "Deleted downloaded GMRT grid(s) (no bathymetry): " + ", ".join(removed),
                    append=True,
                )
            self._show_message(
                "warning",
                "GMRT Download",
                "The downloaded GMRT grid contains no bathymetry (no values below 0 m). "
                "The downloaded file(s) have been deleted and nothing was loaded.",
            )
            return None

        # Split succeeded and we have a bathy file. Remove the original combined file.
        if os.path.exists(downloaded_path):
            try:
                os.remove(downloaded_path)
                log_func(
                    f"Deleted original combined GeoTIFF after split: {os.path.basename(downloaded_path)}",
                    append=True,
                )
            except OSError as exc:
                log_func(
                    f"Failed to delete combined GeoTIFF {os.path.basename(downloaded_path)}: {exc}",
                    append=True,
                )
        return result.bathy_path

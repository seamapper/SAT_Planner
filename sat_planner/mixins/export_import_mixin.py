"""
Export/import: save/load survey parameters (JSON), export survey files (CSV, Shapefile, GeoJSON, etc.).
_save_survey_parameters, _load_survey_parameters, _load_survey_parameters_dialog, _export_survey_files.
"""
import csv
import datetime
import json
import math
import os
import shutil

from PyQt6.QtWidgets import (
    QFileDialog,
    QDialog,
    QDialogButtonBox,
    QPushButton,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QSizePolicy,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from sat_planner.import_survey_dialog import ImportSurveyDialog

from sat_planner import GEOSPATIAL_LIBS_AVAILABLE, decimal_degrees_to_ddm
from sat_planner.constants import LineString, fiona, pyproj, rasterio, Window
from sat_planner import export_utils

try:
    from shapely.geometry import mapping
except ImportError:
    mapping = None

try:
    from rasterio.windows import from_bounds as window_from_bounds
except ImportError:
    window_from_bounds = None


class ExportImportMixin:
    """Mixin for save/load survey parameters and export survey files."""

    _SURVEY_IO_TAB_SPECS = (
        {
            "tab_label": "Calibration Survey",
            "import_label": "Import Calibration Survey",
            "export_label": "Export Calibration Survey",
            "import_handler": "_import_cal_survey_files",
            "export_handler": "_export_cal_survey_files",
            "export_name_attr": "cal_export_name_entry",
            "gmrt_prefix": "cal",
        },
        {
            "tab_label": "Accuracy Survey",
            "import_label": "Import Accuracy Survey",
            "export_label": "Export Accuracy Survey",
            "import_handler": "_import_survey_files",
            "export_handler": "_export_survey_files",
            "export_name_attr": "export_name_entry",
            "gmrt_prefix": "ref",
        },
        {
            "tab_label": "Line Plan",
            "import_label": "Import Line Plan",
            "export_label": "Export Line Plan",
            "import_handler": "_import_drawn_line",
            "export_handler": "_export_drawn_line",
            "export_name_attr": "line_export_name_entry",
            "gmrt_prefix": "line_plan",
        },
        {
            "tab_label": "Backscatter Line",
            "import_label": "Import Backscatter Line",
            "export_label": "Export Backscatter Line",
            "import_handler": "_import_backscatter_line",
            "export_handler": "_export_backscatter_line",
            "export_name_attr": "backscatter_export_name_entry",
            "gmrt_prefix": "backscatter",
            "requires_geospatial": True,
        },
        {
            "tab_label": "Performance Survey",
            "import_label": "Import Performance Survey",
            "export_label": "Export Performance Survey",
            "import_handler": "_import_performance_survey",
            "export_handler": "_export_performance_survey_files",
            "export_name_attr": "performance_export_name_entry",
            "gmrt_prefix": "performance",
        },
        {
            "tab_label": "ADCP Cal",
            "import_label": "Import ADCP Cal",
            "export_label": "Export ADCP Cal",
            "import_handler": "_import_adcp_cal",
            "export_handler": "_export_adcp_cal_files",
            "export_name_attr": "adcp_export_name_entry",
            "gmrt_prefix": "adcp",
        },
    )

    def _survey_io_spec(self, tab_index=None):
        if tab_index is None:
            tab_index = self.param_notebook.currentIndex()
        return self._SURVEY_IO_TAB_SPECS[tab_index]

    def _active_import_button(self):
        return getattr(self, "shared_import_btn", None)

    def _gmrt_widgets_for_prefix(self, prefix):
        return (
            getattr(self, f"{prefix}_download_gmrt_checkbox"),
            getattr(self, f"{prefix}_gmrt_buffer_spin"),
            getattr(self, f"{prefix}_gmrt_buffer_percent_spin"),
            getattr(self, f"{prefix}_gmrt_buffer_mode"),
            getattr(self, f"{prefix}_split_topo_depths_checkbox"),
            getattr(self, f"{prefix}_gmrt_cell_size_combo"),
        )

    def _create_hidden_gmrt_import_widgets(self, prefix, import_tooltip):
        """Create per-tab GMRT widgets (not shown on tab layouts; used by Import Survey dialog)."""
        from sat_planner.bathymetry_download.data_sources import (
            gmrt_cell_size_meters_options,
        )

        holder = QWidget(self)
        holder.hide()
        download_cb = QCheckBox(holder)
        download_cb.setChecked(False)
        download_cb.setToolTip(import_tooltip)
        buffer_spin = QDoubleSpinBox(holder)
        buffer_spin.setRange(0.01, 10.0)
        buffer_spin.setSingleStep(0.1)
        buffer_spin.setValue(0.5)
        buffer_spin.setDecimals(2)
        buffer_spin.setMinimumWidth(60)
        buffer_spin.setToolTip(
            "Fixed buffer in degrees around the survey center for GMRT download."
        )
        buffer_percent_spin = QDoubleSpinBox(holder)
        buffer_percent_spin.setRange(0.1, 200.0)
        buffer_percent_spin.setSingleStep(1.0)
        buffer_percent_spin.setValue(20.0)
        buffer_percent_spin.setDecimals(1)
        buffer_percent_spin.setMinimumWidth(60)
        buffer_percent_spin.setToolTip(
            "Expand the current map-frame view (the lon/lat area shown in the map "
            "window after import) by this percent of its width and height for the "
            "GMRT download. The map view itself is not changed."
        )
        buffer_mode = QComboBox(holder)
        buffer_mode.addItem("Degrees", "degrees")
        buffer_mode.addItem("Percent of survey", "percent")
        buffer_mode.setCurrentIndex(1)  # default: percent
        buffer_mode.setToolTip("Choose fixed degrees buffer or percent of map-frame extent.")
        cell_size_combo = QComboBox(holder)
        default_cell = 60
        for meters in gmrt_cell_size_meters_options():
            meters_i = int(meters)
            cell_size_combo.addItem(f"{meters_i} m", meters_i)
        idx = cell_size_combo.findData(default_cell)
        cell_size_combo.setCurrentIndex(idx if idx >= 0 else 0)
        cell_size_combo.setToolTip(
            "GMRT GridServer cell size (meters/pixel). Same presets as Download Bathymetry."
        )
        cell_size_combo.setMinimumWidth(80)
        split_cb = QCheckBox(holder)
        split_cb.setChecked(True)
        split_cb.setToolTip(
            "When on, the downloaded GMRT GeoTIFF is split into a topography file "
            "(values >= 0) and a bathymetry file (values < 0); SAT Planner loads only "
            "the bathymetry file. When off, a single combined topo+bathy GeoTIFF is loaded."
        )
        split_cb.setEnabled(download_cb.isChecked())
        download_cb.toggled.connect(split_cb.setEnabled)
        setattr(self, f"{prefix}_download_gmrt_checkbox", download_cb)
        setattr(self, f"{prefix}_gmrt_buffer_spin", buffer_spin)
        setattr(self, f"{prefix}_gmrt_buffer_percent_spin", buffer_percent_spin)
        setattr(self, f"{prefix}_gmrt_buffer_mode", buffer_mode)
        setattr(self, f"{prefix}_gmrt_cell_size_combo", cell_size_combo)
        setattr(self, f"{prefix}_split_topo_depths_checkbox", split_cb)

    def _gmrt_import_resolution_meters(self, prefix=None):
        """Return selected post-import GMRT cell size in meters."""
        defaults = (
            self._default_gmrt_import_options()
            if hasattr(self, "_default_gmrt_import_options")
            else {"cell_size_m": 60}
        )
        try:
            if prefix is None:
                prefix = self._survey_io_spec()["gmrt_prefix"]
            combo = getattr(self, f"{prefix}_gmrt_cell_size_combo", None)
            if combo is not None:
                data = combo.currentData()
                if data is not None:
                    return int(data)
                return int(float(combo.currentText().replace("m", "").strip()))
        except Exception:
            pass
        opts = getattr(self, "gmrt_import_options", None) or {}
        try:
            return int(float(opts.get("cell_size_m", defaults.get("cell_size_m", 60))))
        except (TypeError, ValueError):
            return int(defaults.get("cell_size_m", 60))

    def _apply_gmrt_import_options_to_widgets(self):
        """Push remembered GMRT import options onto all per-tab hidden widgets."""
        defaults = (
            self._default_gmrt_import_options()
            if hasattr(self, "_default_gmrt_import_options")
            else {
                "download": True,
                "buffer_mode": "percent",
                "buffer_deg": 0.5,
                "buffer_percent": 20.0,
                "cell_size_m": 60,
                "split_topo_depths": True,
            }
        )
        options = getattr(self, "gmrt_import_options", None) or defaults
        if hasattr(self, "_normalize_gmrt_import_options"):
            options = self._normalize_gmrt_import_options(options)
        download = bool(options.get("download", defaults["download"]))
        buffer_deg = float(options.get("buffer_deg", defaults["buffer_deg"]))
        buffer_percent = float(options.get("buffer_percent", defaults.get("buffer_percent", 20.0)))
        buffer_mode = options.get("buffer_mode", defaults.get("buffer_mode", "percent"))
        cell_size_m = int(options.get("cell_size_m", defaults.get("cell_size_m", 60)))
        split_topo = bool(options.get("split_topo_depths", defaults["split_topo_depths"]))
        for spec in self._SURVEY_IO_TAB_SPECS:
            download_cb, buffer_spin, buffer_percent_spin, mode_combo, split_cb, cell_combo = (
                self._gmrt_widgets_for_prefix(spec["gmrt_prefix"])
            )
            download_cb.blockSignals(True)
            download_cb.setChecked(download)
            download_cb.blockSignals(False)
            buffer_spin.setValue(buffer_deg)
            buffer_percent_spin.setValue(buffer_percent)
            idx = mode_combo.findData(buffer_mode)
            if idx < 0:
                idx = 0
            mode_combo.setCurrentIndex(idx)
            cell_idx = cell_combo.findData(cell_size_m)
            if cell_idx < 0:
                cell_idx = cell_combo.findData(int(defaults.get("cell_size_m", 60)))
            if cell_idx < 0:
                cell_idx = 0
            cell_combo.setCurrentIndex(cell_idx)
            split_cb.setChecked(split_topo)
            split_cb.setEnabled(download)

    def _capture_gmrt_import_options_from_widgets(self, prefix=None):
        """Update gmrt_import_options from the given (or current-tab) widgets."""
        if prefix is None:
            prefix = self._survey_io_spec()["gmrt_prefix"]
        download_cb, buffer_spin, buffer_percent_spin, mode_combo, split_cb, cell_combo = (
            self._gmrt_widgets_for_prefix(prefix)
        )
        defaults = (
            self._default_gmrt_import_options()
            if hasattr(self, "_default_gmrt_import_options")
            else {
                "download": True,
                "buffer_mode": "percent",
                "buffer_deg": 0.5,
                "buffer_percent": 20.0,
                "cell_size_m": 60,
                "split_topo_depths": True,
            }
        )
        try:
            buffer_deg = float(buffer_spin.value())
        except (TypeError, ValueError):
            buffer_deg = defaults["buffer_deg"]
        try:
            buffer_percent = float(buffer_percent_spin.value())
        except (TypeError, ValueError):
            buffer_percent = defaults.get("buffer_percent", 20.0)
        mode = mode_combo.currentData() if mode_combo is not None else "percent"
        if not mode:
            mode = "percent"
        try:
            cell_size_m = int(cell_combo.currentData())
        except (TypeError, ValueError):
            cell_size_m = int(defaults.get("cell_size_m", 60))
        options = {
            "download": bool(download_cb.isChecked()),
            "buffer_mode": "percent" if mode == "percent" else "degrees",
            "buffer_deg": max(0.01, min(10.0, buffer_deg)),
            "buffer_percent": max(0.1, min(200.0, buffer_percent)),
            "cell_size_m": cell_size_m,
            "split_topo_depths": bool(split_cb.isChecked()),
        }
        if hasattr(self, "_normalize_gmrt_import_options"):
            options = self._normalize_gmrt_import_options(options)
        self.gmrt_import_options = options

    def _gmrt_map_frame_extent(self):
        """Return (west, east, south, north) for the visible map frame, or None.

        Uses current axes limits and the same fill-the-frame math as
        ``_plot_survey_plan`` (geographic aspect vs figure aspect). If the
        axes are already filled, this is a no-op on the numbers; if they are
        still the tight survey/zoom box, this expands to the on-screen frame.
        Does not change the map view.
        """
        if not hasattr(self, "ax") or self.ax is None:
            return None
        try:
            if hasattr(self, "canvas") and self.canvas is not None:
                self.canvas.draw()
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            lon0, lon1 = sorted((float(xlim[0]), float(xlim[1])))
            lat0, lat1 = sorted((float(ylim[0]), float(ylim[1])))
            if lon1 <= lon0 or lat1 <= lat0:
                return None

            # Match _plot_survey_plan fill-frame adjustment.
            import numpy as np

            center_lat = (lat0 + lat1) / 2.0
            center_x = (lon0 + lon1) / 2.0
            center_y = (lat0 + lat1) / 2.0
            aspect_ratio = float(np.clip(1.0 / np.cos(np.radians(center_lat)), 0.1, 10.0))
            if hasattr(self, "figure") and self.figure is not None:
                fig_width, fig_height = self.figure.get_size_inches()
                plot_width = fig_width * (0.99 - 0.085)
                plot_height = fig_height * (0.95 - 0.08)
                figure_aspect = plot_width / plot_height if plot_height else 1.0
            else:
                figure_aspect = 1.0
            data_width = lon1 - lon0
            data_height = lat1 - lat0
            data_display_aspect = (data_width * aspect_ratio) / data_height
            if data_display_aspect > figure_aspect:
                new_height = (data_width * aspect_ratio) / figure_aspect
                lat0 = center_y - new_height / 2.0
                lat1 = center_y + new_height / 2.0
            elif data_display_aspect < figure_aspect:
                new_width = (data_height * figure_aspect) / aspect_ratio
                lon0 = center_x - new_width / 2.0
                lon1 = center_x + new_width / 2.0
            return (lon0, lon1, lat0, lat1)
        except Exception:
            return None

    def _gmrt_download_extent_from_points(self, points, prefix=None):
        """Return (west, east, south, north) for GMRT download from lat/lon points.

        Degrees mode: fixed box centered on the survey midpoint (existing behavior).
        Percent mode: start from the visible map-frame bounds (what fills the map
        window after import), then expand by buffer_percent of that frame's
        width and height. The map view is not changed.
        """
        if not points:
            raise ValueError("No points for GMRT extent")
        if prefix is None:
            prefix = self._survey_io_spec()["gmrt_prefix"]
        lats = [float(p[0]) for p in points]
        lons = [float(p[1]) for p in points]
        min_lat, max_lat = min(lats), max(lats)
        min_lon, max_lon = min(lons), max(lons)
        mid_lat = (min_lat + max_lat) / 2.0
        mid_lon = (min_lon + max_lon) / 2.0

        _, buffer_spin, buffer_percent_spin, mode_combo, _, _ = self._gmrt_widgets_for_prefix(prefix)
        mode = mode_combo.currentData() if mode_combo is not None else "percent"
        if not mode:
            mode = "percent"

        if mode == "percent":
            try:
                pct = float(buffer_percent_spin.value())
            except (TypeError, ValueError):
                pct = 20.0
            pct = max(0.1, min(200.0, pct))

            frame = self._gmrt_map_frame_extent()
            if frame is not None:
                view_min_lon, view_max_lon, view_min_lat, view_max_lat = frame
            else:
                # Fallback if axes are unavailable: survey bbox + Zoom-to-Plan padding.
                lat_span = max_lat - min_lat
                lon_span = max_lon - min_lon
                view_pad_lat = lat_span * 0.05 if lat_span != 0 else 0.01
                view_pad_lon = lon_span * 0.05 if lon_span != 0 else 0.01
                view_min_lat = min_lat - view_pad_lat
                view_max_lat = max_lat + view_pad_lat
                view_min_lon = min_lon - view_pad_lon
                view_max_lon = max_lon + view_pad_lon

            width = max(view_max_lon - view_min_lon, 1e-8)
            height = max(view_max_lat - view_min_lat, 1e-8)
            pad_x = width * (pct / 100.0) / 2.0
            pad_y = height * (pct / 100.0) / 2.0
            return (
                view_min_lon - pad_x,
                view_max_lon + pad_x,
                view_min_lat - pad_y,
                view_max_lat + pad_y,
            )

        try:
            buffer_deg = float(buffer_spin.value())
        except (TypeError, ValueError):
            buffer_deg = 0.5
        buffer_deg = max(0.01, min(10.0, buffer_deg))
        return (
            mid_lon - buffer_deg,
            mid_lon + buffer_deg,
            mid_lat - buffer_deg,
            mid_lat + buffer_deg,
        )

    def _gmrt_split_topo_depths_for_prefix(self, prefix=None):
        if prefix is None:
            prefix = self._survey_io_spec()["gmrt_prefix"]
        split_cb = getattr(self, f"{prefix}_split_topo_depths_checkbox", None)
        if split_cb is None:
            return True
        return bool(split_cb.isChecked())

    def _create_export_name_page(self, stack, default_text=""):
        page = QWidget()
        layout = QGridLayout(page)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(3)
        layout.setColumnStretch(0, 1)
        layout.setColumnStretch(1, 2)
        layout.addWidget(QLabel("Export Name:"), 0, 0)
        entry = QLineEdit()
        if default_text:
            entry.setText(default_text)
        layout.addWidget(entry, 0, 1)
        stack.addWidget(page)
        return entry

    def _create_shared_survey_io_widgets(self):
        """Create shared Import/Export widgets (call after all tab fields exist)."""
        gmrt_tooltip = (
            "When enabled, importing a survey will download a GMRT bathymetry GeoTIFF "
            "(buffer and 100 m resolution) and load it."
        )
        for spec in self._SURVEY_IO_TAB_SPECS:
            self._create_hidden_gmrt_import_widgets(spec["gmrt_prefix"], gmrt_tooltip)
        self._apply_gmrt_import_options_to_widgets()

        self._shared_io_button_row = QHBoxLayout()
        self.shared_import_btn = QPushButton("Import Survey")
        self.shared_import_btn.clicked.connect(self._on_shared_import_clicked)
        self._shared_io_button_row.addWidget(self.shared_import_btn, 1)
        self.shared_export_btn = QPushButton("Export Survey")
        self.shared_export_btn.clicked.connect(self._on_shared_export_clicked)
        self._shared_io_button_row.addWidget(self.shared_export_btn, 1)

        self.shared_export_name_stack = QStackedWidget()

        self.cal_export_name_entry = self._create_export_name_page(self.shared_export_name_stack)

        acc_export_default = "acc_depth0m_cross90deg"
        if hasattr(self, "heading_entry"):
            try:
                heading = float(self.heading_entry.text() or "0")
                cross = int(round((heading + 90) % 360))
                acc_export_default = f"acc_depth{0}m_cross{cross}deg".replace(
                    "acc_depth0m", f"acc_depth0m"
                )
                acc_export_default = f"acc_depth0m_cross{cross}deg"
            except Exception:
                pass
        self.export_name_entry = self._create_export_name_page(
            self.shared_export_name_stack, acc_export_default
        )

        line_default = f"Line_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
        self.line_export_name_entry = self._create_export_name_page(
            self.shared_export_name_stack, line_default
        )

        bs_default = f"BS_{datetime.datetime.now().strftime('%Y%m%d')}_0"
        self.backscatter_export_name_entry = self._create_export_name_page(
            self.shared_export_name_stack, bs_default
        )

        perf_default = (
            self._build_performance_export_basename()
            if hasattr(self, "_build_performance_export_basename")
            else "perf_swell0_depth0m"
        )
        self.performance_export_name_entry = self._create_export_name_page(
            self.shared_export_name_stack, perf_default
        )

        adcp_default = (
            self._build_adcp_export_basename()
            if hasattr(self, "_build_adcp_export_basename")
            else "adcp_cal"
        )
        self.adcp_export_name_entry = self._create_export_name_page(
            self.shared_export_name_stack, adcp_default
        )

    def _install_shared_survey_io(self, parent_layout):
        """Insert shared Import/Export controls above the tab widget."""
        if not hasattr(self, "shared_import_btn"):
            self._create_shared_survey_io_widgets()
        if not hasattr(self, "shared_survey_io_groupbox"):
            self.shared_survey_io_groupbox = QGroupBox("Import/Export")
            self.shared_survey_io_groupbox.setSizePolicy(
                QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
            )
            io_layout = QVBoxLayout(self.shared_survey_io_groupbox)
            io_layout.setContentsMargins(9, 9, 9, 9)
            io_layout.setSpacing(3)
            io_layout.addLayout(self._shared_io_button_row)
            io_layout.addWidget(self.shared_export_name_stack)
        parent_layout.insertWidget(0, self.shared_survey_io_groupbox)
        parent_layout.setStretch(0, 0)
        if parent_layout.count() > 1:
            parent_layout.setStretch(1, 1)
        self._update_shared_survey_io_ui()

    def _setup_shared_survey_io(self, parent_layout):
        """Backward-compatible alias."""
        self._install_shared_survey_io(parent_layout)

    def _shared_survey_io_import_enabled(self, tab_index):
        spec = self._survey_io_spec(tab_index)
        if spec.get("requires_geospatial") and not GEOSPATIAL_LIBS_AVAILABLE:
            return False
        return True

    def _shared_survey_io_export_enabled(self, tab_index):
        if tab_index == 5 and hasattr(self, "_adcp_plan_complete"):
            return self._adcp_plan_complete()
        if self._survey_io_spec(tab_index).get("requires_geospatial") and not GEOSPATIAL_LIBS_AVAILABLE:
            return False
        return True

    def _update_shared_survey_io_ui(self):
        if not hasattr(self, "param_notebook") or not hasattr(self, "shared_import_btn"):
            return
        tab_index = self.param_notebook.currentIndex()
        spec = self._survey_io_spec(tab_index)
        self.shared_export_name_stack.setCurrentIndex(tab_index)
        self.shared_import_btn.setText(spec["import_label"])
        self.shared_export_btn.setText(spec["export_label"])
        self.shared_import_btn.setEnabled(self._shared_survey_io_import_enabled(tab_index))
        self.shared_export_btn.setEnabled(self._shared_survey_io_export_enabled(tab_index))

    def _import_has_planning_geotiff(self, geotiff_path=None):
        """True when the imported survey references an existing planning GeoTIFF.

        A grid already loaded in the session (e.g. from a previous GMRT download)
        does not count — each import without its own bathymetry should still
        offer a GMRT download.
        """
        return bool(geotiff_path and os.path.isfile(geotiff_path))

    def _session_has_loaded_geotiff(self):
        """True when a planning bathymetry grid is currently loaded in the session."""
        return (
            getattr(self, "geotiff_data_array", None) is not None
            and getattr(self, "geotiff_extent", None) is not None
        ) or getattr(self, "geotiff_dataset_original", None) is not None

    def _maybe_prompt_gmrt_download_after_import(self, geotiff_path=None, download_callback=None):
        """After import, offer GMRT download when the survey has no planning bathymetry.

        If a session grid is already loaded, the dialog asks whether to keep it or
        download a new GMRT grid. Returns True when a GMRT download was started.
        """
        if not callable(download_callback):
            return False
        if self._import_has_planning_geotiff(geotiff_path):
            return False

        existing_grid = self._session_has_loaded_geotiff()
        spec = self._survey_io_spec()
        download_cb, buffer_spin, buffer_percent_spin, mode_combo, split_cb, cell_combo = (
            self._gmrt_widgets_for_prefix(spec["gmrt_prefix"])
        )
        dialog = ImportSurveyDialog(
            self,
            spec["tab_label"],
            download_cb,
            buffer_spin,
            buffer_percent_spin,
            mode_combo,
            split_cb,
            cell_combo,
            prompt_missing_geotiff=True,
            geotiff_path=geotiff_path,
            existing_grid_loaded=existing_grid,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return False
        self._capture_gmrt_import_options_from_widgets(spec["gmrt_prefix"])
        self._apply_gmrt_import_options_to_widgets()
        if hasattr(self, "_save_gmrt_import_options"):
            self._save_gmrt_import_options()
        if download_cb.isChecked():
            download_callback()
            return True
        return False

    def _on_shared_import_clicked(self):
        if hasattr(self, "_gmrt_is_downloading") and self._gmrt_is_downloading():
            self._gmrt_cancel_active_download()
            return
        spec = self._survey_io_spec()
        handler = getattr(self, spec["import_handler"])
        handler()

    def _on_shared_export_clicked(self):
        spec = self._survey_io_spec()
        handler = getattr(self, spec["export_handler"])
        handler()

    def _select_export_directory(self, start_dir=None, title="Select Export Directory"):
        """Directory picker with an Export Types button; returns path or None."""
        start = start_dir or getattr(self, "last_export_dir", None) or os.path.expanduser("~")
        if not os.path.isdir(start):
            start = os.path.expanduser("~")

        dialog = QFileDialog(self, title, start)
        dialog.setFileMode(QFileDialog.FileMode.Directory)
        dialog.setOption(QFileDialog.Option.ShowDirsOnly, True)
        dialog.setOption(QFileDialog.Option.DontUseNativeDialog, True)

        export_types_btn = QPushButton("Export Types")
        if hasattr(self, "_show_export_type_dialog"):
            export_types_btn.clicked.connect(self._show_export_type_dialog)
        else:
            export_types_btn.setEnabled(False)

        button_box = dialog.findChild(QDialogButtonBox)
        if button_box is not None:
            button_box.addButton(export_types_btn, QDialogButtonBox.ButtonRole.ActionRole)
        else:
            layout = dialog.layout()
            if layout is not None:
                layout.addWidget(export_types_btn)

        if dialog.exec() != QFileDialog.DialogCode.Accepted:
            return None
        selected = dialog.selectedFiles()
        if not selected:
            return None
        path = selected[0]
        return path if os.path.isdir(path) else None

    def _gpkg_path_for_shapefile(self, shapefile_path):
        """Return the .gpkg companion path for a given shapefile path."""
        if not shapefile_path:
            return None
        base, _ = os.path.splitext(shapefile_path)
        return base + ".gpkg"

    def _write_gpkg_if_enabled(self, shapefile_path, schema, features, *, crs="EPSG:4326", layer_name=None):
        """Write the same features that were exported to a shapefile to a GeoPackage.

        No-op (returns ``None``) when:
          * the ``gpkg`` export option is disabled,
          * ``fiona`` is not available, or
          * ``features`` is empty.

        Errors during write are swallowed (with a console message) so they
        never break the main shapefile export. Layer name defaults to the
        file basename, matching Fiona's default behavior.
        """
        if not self._export_type_enabled("gpkg"):
            return None
        if fiona is None:
            return None
        if not features:
            return None
        gpkg_path = self._gpkg_path_for_shapefile(shapefile_path)
        if not gpkg_path:
            return None
        try:
            export_utils.remove_export_file(gpkg_path)
            kwargs = {"driver": "GPKG", "crs": crs, "schema": schema}
            if layer_name:
                kwargs["layer"] = layer_name
            with fiona.open(gpkg_path, "w", **kwargs) as collection:
                collection.writerecords(features)
            return gpkg_path
        except Exception as e:
            print(f"Warning: failed to write GeoPackage {gpkg_path}: {e}")
            return None

    def _export_type_enabled(self, key):
        defaults = (
            self._default_export_type_options()
            if hasattr(self, "_default_export_type_options")
            else {}
        )
        options = getattr(self, "export_type_options", {}) or {}
        if key in options:
            return bool(options[key])
        legacy_png = {
            "map_png_high": "map_png",
            "map_png_low": "map_png",
            "profiles_png_high": "profiles_png",
            "profiles_png_low": "profiles_png",
        }
        if key in legacy_png and legacy_png[key] in options:
            return bool(options[legacy_png[key]])
        return bool(defaults.get(key, True))

    def _export_map_png_enabled(self):
        return self._export_type_enabled("map_png_high") or self._export_type_enabled("map_png_low")

    def _export_profiles_png_enabled(self):
        return self._export_type_enabled("profiles_png_high") or self._export_type_enabled("profiles_png_low")

    def _save_export_map_png(self, path, **kwargs):
        had_arrows = False
        if hasattr(self, "_set_travel_direction_arrows_visible"):
            had_arrows = self._set_travel_direction_arrows_visible(False)
        try:
            return export_utils.save_export_png(
                self.figure,
                path,
                save_high=self._export_type_enabled("map_png_high"),
                save_low=self._export_type_enabled("map_png_low"),
                **kwargs,
            )
        finally:
            if had_arrows and hasattr(self, "_set_travel_direction_arrows_visible"):
                self._set_travel_direction_arrows_visible(True)

    def _save_export_profile_png(self, path, **kwargs):
        return export_utils.save_export_png(
            self.profile_fig,
            path,
            save_high=self._export_type_enabled("profiles_png_high"),
            save_low=self._export_type_enabled("profiles_png_low"),
            **kwargs,
        )

    def _map_png_export_basenames(self, path):
        return export_utils.png_export_basenames(
            path,
            include_high=self._export_type_enabled("map_png_high"),
            include_low=self._export_type_enabled("map_png_low"),
        )

    def _profile_png_export_basenames(self, path):
        return export_utils.png_export_basenames(
            path,
            include_high=self._export_type_enabled("profiles_png_high"),
            include_low=self._export_type_enabled("profiles_png_low"),
        )

    def _geotiff_loaded_for_export(self):
        """True when a planning GeoTIFF is available to export."""
        if getattr(self, "geotiff_dataset_original", None) is not None:
            return True
        path = getattr(self, "current_geotiff_path", None)
        return bool(path and os.path.isfile(path))

    def _geotiff_cell_size_meters_int(self):
        """Larger of X/Y pixel size in meters (integer), or None if unavailable."""
        ds = getattr(self, "geotiff_dataset_original", None)
        if ds is None:
            path = getattr(self, "current_geotiff_path", None)
            if not (path and os.path.isfile(path) and rasterio is not None):
                return None
            try:
                with rasterio.open(path) as opened:
                    return self._geotiff_cell_size_meters_int_from_dataset(opened)
            except Exception:
                return None
        return self._geotiff_cell_size_meters_int_from_dataset(ds)

    def _geotiff_cell_size_meters_int_from_dataset(self, ds):
        try:
            transform = ds.transform
            res_x = abs(float(transform.a))
            res_y = abs(float(transform.e))
            crs = ds.crs
            geographic = False
            if crs is not None:
                try:
                    geographic = bool(crs.is_geographic)
                except Exception:
                    geographic = str(crs).upper() in ("EPSG:4326", "WGS84")
            if geographic:
                try:
                    center_lat = (float(ds.bounds.top) + float(ds.bounds.bottom)) / 2.0
                except Exception:
                    center_lat = 0.0
                m_per_deg_lat = 111320.0
                m_per_deg_lon = 111320.0 * math.cos(math.radians(center_lat))
                res_x *= m_per_deg_lon
                res_y *= m_per_deg_lat
            cell_m = max(res_x, res_y)
            if not math.isfinite(cell_m) or cell_m <= 0:
                return None
            return max(1, int(round(cell_m)))
        except Exception:
            return None

    def _path_for_params_sidecar(self, path, base_dir=None):
        """Return a portable path for ``*_params.json`` / GeoJSON sidecars.

        When ``path`` is in ``base_dir`` (same folder as the export package),
        store only the basename so the package can be moved. Otherwise store
        an absolute path. Returns ``None`` if ``path`` is empty.
        """
        if not path:
            return None
        try:
            abs_path = os.path.abspath(path)
        except Exception:
            return path
        if base_dir:
            try:
                abs_base = os.path.abspath(base_dir)
                if os.path.normcase(os.path.dirname(abs_path)) == os.path.normcase(abs_base):
                    return os.path.basename(abs_path)
            except Exception:
                pass
        return abs_path

    def _resolve_export_params_geotiff_path(self, exported_geotiff_path=None, export_dir=None):
        """geotiff_path for params: exported file, else loaded path, else None.

        When ``export_dir`` is set and the chosen file lives in that directory,
        the returned value is a relative basename for portability.
        """
        chosen = None
        if exported_geotiff_path and os.path.isfile(exported_geotiff_path):
            chosen = exported_geotiff_path
        else:
            path = getattr(self, "current_geotiff_path", None)
            if path and os.path.isfile(path):
                chosen = path
        if chosen is None:
            return None
        return self._path_for_params_sidecar(chosen, export_dir)

    def _resolve_import_sidecar_path(self, stored_path, base_dir=None):
        """Resolve a path stored in ``*_params.json`` / GeoJSON for loading.

        Relative paths (and bare filenames) resolve against ``base_dir`` (the
        directory of the imported survey or params file). Absolute paths that
        no longer exist also fall back to ``base_dir`` / basename so packages
        moved after an older absolute-path export still open.
        """
        if not stored_path or not isinstance(stored_path, str):
            return None
        stored_path = stored_path.strip()
        if not stored_path:
            return None

        candidates = []
        if os.path.isabs(stored_path):
            candidates.append(stored_path)
        if base_dir:
            candidates.append(os.path.normpath(os.path.join(base_dir, stored_path)))
            candidates.append(os.path.join(base_dir, os.path.basename(stored_path)))
        elif not os.path.isabs(stored_path):
            candidates.append(os.path.abspath(stored_path))

        seen = set()
        for candidate in candidates:
            if not candidate:
                continue
            try:
                key = os.path.normcase(os.path.abspath(candidate))
            except Exception:
                key = candidate
            if key in seen:
                continue
            seen.add(key)
            if os.path.isfile(candidate):
                try:
                    return os.path.abspath(candidate)
                except Exception:
                    return candidate
        if os.path.isabs(stored_path):
            return stored_path
        if base_dir:
            return os.path.normpath(os.path.join(base_dir, stored_path))
        return stored_path

    def _maybe_export_survey_geotiff(self, export_dir, export_name):
        """Export Full or View GeoTIFF when enabled. Returns output path or None."""
        export_full = self._export_type_enabled("geotiff_full")
        export_view = self._export_type_enabled("geotiff_view")
        if export_full and export_view:
            export_view = False
        if not export_full and not export_view:
            return None
        if not self._geotiff_loaded_for_export():
            return None
        if not GEOSPATIAL_LIBS_AVAILABLE or rasterio is None:
            return None

        cell_m = self._geotiff_cell_size_meters_int()
        if cell_m is None:
            cell_m = 0
        mode = "Full" if export_full else "View"
        source_tag = getattr(self, "current_geotiff_bathy_source_tag", None)
        source_suffix = f"_{source_tag}" if source_tag else ""
        out_name = f"{export_name}_{cell_m}m_{mode}{source_suffix}.tif"
        out_path = os.path.join(export_dir, out_name)
        try:
            export_utils.remove_export_file(out_path)
            if export_full:
                src_path = getattr(self, "current_geotiff_path", None)
                if src_path and os.path.isfile(src_path):
                    shutil.copy2(src_path, out_path)
                else:
                    ds = self.geotiff_dataset_original
                    profile = ds.profile.copy()
                    profile.update(driver="GTiff", count=1)
                    with rasterio.open(out_path, "w", **profile) as dst:
                        dst.write(ds.read(1), 1)
            else:
                if not self._export_geotiff_view_window(out_path):
                    return None
            return out_path if os.path.isfile(out_path) else None
        except Exception as e:
            print(f"Warning: GeoTIFF export failed: {e}")
            return None

    def _export_geotiff_view_window(self, out_path):
        """Write the map-view portion of the loaded GeoTIFF to out_path.

        Expands the current map view by 10% in width and height (5% each side),
        then clips to the original GeoTIFF bounds.
        """
        ds = getattr(self, "geotiff_dataset_original", None)
        if ds is None or not hasattr(self, "ax") or window_from_bounds is None:
            return False
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            lon_left, lon_right = float(xlim[0]), float(xlim[1])
            lat_bottom, lat_top = float(ylim[0]), float(ylim[1])
            if lon_right < lon_left:
                lon_left, lon_right = lon_right, lon_left
            if lat_top < lat_bottom:
                lat_bottom, lat_top = lat_top, lat_bottom

            # Expand view by 10% wider and taller (5% pad each side).
            pad_x = (lon_right - lon_left) * 0.05
            pad_y = (lat_top - lat_bottom) * 0.05
            lon_left -= pad_x
            lon_right += pad_x
            lat_bottom -= pad_y
            lat_top += pad_y

            if ds.crs is not None and not bool(getattr(ds.crs, "is_geographic", False)):
                if pyproj is None:
                    return False
                transformer = pyproj.Transformer.from_crs("EPSG:4326", ds.crs, always_xy=True)
                x0, y0 = transformer.transform(lon_left, lat_bottom)
                x1, y1 = transformer.transform(lon_right, lat_top)
                left, right = min(x0, x1), max(x0, x1)
                bottom, top = min(y0, y1), max(y0, y1)
            else:
                left, right = lon_left, lon_right
                bottom, top = lat_bottom, lat_top

            # Never exceed the original grid bounds.
            try:
                grid_left, grid_bottom, grid_right, grid_top = (
                    float(ds.bounds.left),
                    float(ds.bounds.bottom),
                    float(ds.bounds.right),
                    float(ds.bounds.top),
                )
                left = max(left, min(grid_left, grid_right))
                right = min(right, max(grid_left, grid_right))
                bottom = max(bottom, min(grid_bottom, grid_top))
                top = min(top, max(grid_bottom, grid_top))
            except Exception:
                pass
            if right <= left or top <= bottom:
                return False

            window = window_from_bounds(left, bottom, right, top, transform=ds.transform)
            if Window is not None:
                window = window.intersection(Window(0, 0, ds.width, ds.height))
            if window.width < 1 or window.height < 1:
                return False
            data = ds.read(1, window=window)
            out_transform = ds.window_transform(window)
            profile = ds.profile.copy()
            profile.update(
                driver="GTiff",
                height=data.shape[0],
                width=data.shape[1],
                transform=out_transform,
                count=1,
            )
            with rasterio.open(out_path, "w", **profile) as dst:
                dst.write(data, 1)
            return True
        except Exception as e:
            print(f"Warning: GeoTIFF view export failed: {e}")
            return False

    def _save_survey_parameters(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot save parameters.")
            return

        is_valid, values = self._validate_inputs()
        if not is_valid:
            return

        save_dir = QFileDialog.getExistingDirectory(self, "Select Directory to Save Survey Parameters", self.last_survey_params_dir)
        if not save_dir:
            return
        self.last_survey_params_dir = save_dir
        self._save_last_survey_params_dir()

        # Use the export name from the form, with .json extension
        export_name = values['export_name']
        if not export_name:
            # Fallback to default naming if export name is empty: depth at center
            # (m) + main line heading (deg). Depth comes from Pick Center or an
            # imported survey's params; falls back to 0 when not available.
            if hasattr(self, "_build_accuracy_export_basename"):
                export_name = self._build_accuracy_export_basename()
            else:
                picked_depth = getattr(self, "_depth_at_picked_point", None)
                try:
                    depth_int = int(round(abs(float(picked_depth)))) if picked_depth is not None else 0
                except (TypeError, ValueError):
                    depth_int = 0
                try:
                    cross_int = int(round((float(values["heading"]) + 90) % 360))
                except (TypeError, ValueError):
                    cross_int = 0
                export_name = f"acc_depth{depth_int}m_cross{cross_int}deg"
        filename = f"{export_name}_params.json"
        file_path = os.path.join(save_dir, filename)

        try:
            # Create parameters dictionary
            params = {
                'central_lat': float(self.central_lat_entry.text()),
                'central_lon': float(self.central_lon_entry.text()),
                'central_point_depth_m': float(getattr(self, '_depth_at_picked_point', 0.0)) if getattr(self, '_depth_at_picked_point', None) is not None else None,
                'line_length': float(self.line_length_entry.text()),
                'heading': float(self.heading_entry.text()),
                'dist_between_lines': float(self.dist_between_lines_entry.text()),
                'num_lines': int(self.num_lines_entry.text()),
                'bisect_lead': float(self.bisect_lead_entry.text()),
                'survey_speed': float(self.survey_speed_entry.text()),
                'export_name': self.export_name_entry.text().strip(),
                'geotiff_path': self._resolve_export_params_geotiff_path(export_dir=save_dir),
                # Backscatter survey uses BOTH:
                # - bathymetry grid (geotiff_path)
                # - optional separate backscatter raster (backscatter_geotiff_path)
                'backscatter_geotiff_path': self._path_for_params_sidecar(
                    (
                        self.backscatter_geotiff_path
                        if hasattr(self, 'backscatter_geotiff_path') and self.backscatter_geotiff_path
                        else None
                    ),
                    save_dir,
                ),
                'geotiff_nan_value': float(getattr(self, 'geotiff_nan_value', -11000.0)),
                'visualization_shapefile_paths': list(getattr(self, 'visualization_shapefile_paths', []) or []),
                'show_contours_var': bool(getattr(self, 'show_contours_var', False)),
                'contour_interval_m': (
                    float(self.contour_interval_entry.text())
                    if hasattr(self, 'contour_interval_entry') and self.contour_interval_entry.text()
                    else 200.0
                ),
                'offset_direction': self.offset_direction_var,
                'line_length_multiplier': self.line_length_multiplier,
                'dist_between_lines_multiplier': self.dist_between_lines_multiplier
            }
            self._add_geotiff_viz_params_to_params(params)

            # Save to JSON file
            with open(file_path, 'w') as f:
                json.dump(params, f, indent=4)

            self.set_ref_info_text(f"Survey parameters saved successfully to:\n{file_path}", append=False)

        except Exception as e:
            self._show_message("error","Save Error", f"Failed to save survey parameters: {e}")

    def _load_survey_parameters(self, file_path):
        """Load survey parameters from a JSON file."""
        try:
            with open(file_path, 'r') as f:
                params = json.load(f)

            # Update all input fields
            self.central_lat_entry.clear()
            self.central_lat_entry.setText(str(params['central_lat']))

            self.central_lon_entry.clear()
            self.central_lon_entry.setText(str(params['central_lon']))

            self.line_length_entry.setText(str(params['line_length']))

            self.heading_entry.setText(str(params['heading']))

            self.dist_between_lines_entry.setText(str(params['dist_between_lines']))

            self.num_lines_entry.setText(str(params['num_lines']))

            self.bisect_lead_entry.setText(str(params['bisect_lead']))

            self.survey_speed_entry.setText(str(params.get('survey_speed', '')))

            self.export_name_entry.clear()
            self.export_name_entry.setText(params['export_name'])

            self.offset_direction_var.set(params['offset_direction'])

            self.line_length_multiplier = float(params['line_length_multiplier'])
            self._update_multiplier_label_len(self.line_length_multiplier)

            self.dist_between_lines_multiplier = float(params['dist_between_lines_multiplier'])
            self._update_multiplier_label_dist(self.dist_between_lines_multiplier)

            if 'geotiff_nan_value' in params and hasattr(self, '_set_geotiff_nan_cutoff'):
                self._set_geotiff_nan_cutoff(params.get('geotiff_nan_value'), update_entry=True)
            if 'show_contours_var' in params:
                self.show_contours_var = bool(params.get('show_contours_var'))
                if hasattr(self, 'show_contours_checkbox'):
                    self.show_contours_checkbox.blockSignals(True)
                    self.show_contours_checkbox.setChecked(self.show_contours_var)
                    self.show_contours_checkbox.blockSignals(False)
            if params.get('contour_interval_m') is not None and hasattr(self, 'contour_interval_entry'):
                try:
                    self.contour_interval_entry.setText(f"{float(params.get('contour_interval_m')):g}")
                except Exception:
                    pass

            self._apply_geotiff_viz_params_from_params(params)

            # Restore bathymetry GeoTIFF and optional backscatter GeoTIFF (if present).
            # This is best-effort: if geotiff paths are missing or incompatible, we continue.
            params_dir = os.path.dirname(file_path)
            geotiff_path = self._resolve_import_sidecar_path(params.get('geotiff_path'), params_dir)
            backscatter_geotiff_path = self._resolve_import_sidecar_path(
                params.get('backscatter_geotiff_path'), params_dir
            )
            try:
                if geotiff_path and hasattr(self, '_load_geotiff_from_path') and os.path.exists(geotiff_path):
                    # Load bathymetry first (needed so backscatter can be reprojected/aligned).
                    self._load_geotiff_from_path(geotiff_path, use_background_loading=False)
                    # Re-apply stored NaN cutoff after load (GeoTIFF loader may detect its own nodata).
                    if 'geotiff_nan_value' in params and hasattr(self, '_set_geotiff_nan_cutoff'):
                        self._set_geotiff_nan_cutoff(params.get('geotiff_nan_value'), update_entry=True)
                    if hasattr(self, '_reload_geotiff_at_current_zoom'):
                        self._reload_geotiff_at_current_zoom()
            except Exception:
                pass

            try:
                if (
                    backscatter_geotiff_path
                    and hasattr(self, '_load_backscatter_geotiff_from_path')
                    and os.path.exists(backscatter_geotiff_path)
                ):
                    # Requires an already-loaded bathymetry grid (geotiff_data_array / extent).
                    self._load_backscatter_geotiff_from_path(backscatter_geotiff_path)
            except Exception:
                pass

            # Regenerate the plot with loaded parameters
            self._generate_and_plot()

            self.set_ref_info_text(f"Survey parameters loaded from: {file_path}", append=False)

        except Exception as e:
            self._show_message("error","Load Error", f"Failed to load survey parameters: {e}")

    def _load_survey_parameters_dialog(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot load parameters.")
            return

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Load Survey Parameters",
            self.last_survey_params_dir,
            "JSON files (*.json);;All files (*.*)"
        )
        if file_path:
            self.last_survey_params_dir = os.path.dirname(file_path)
            self._save_last_survey_params_dir()
            self._load_survey_parameters(file_path)

    def _export_survey_files(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot export survey files.")
            return

        if hasattr(self, "_commit_all_deferred_line_edits"):
            self._commit_all_deferred_line_edits()

        is_valid, values = self._validate_inputs()
        if not is_valid:
            return

        export_name = values['export_name']

        if not self.survey_lines_data and not self.cross_line_data:
            self._show_message("warning","No Data", "No survey lines to export. Generate them first.")
            return

        export_dir = self._select_export_directory(self.last_export_dir)
        if not export_dir:
            return
        self.last_export_dir = export_dir
        self._save_last_export_dir()

        try:
            profile_csv_path = None
            profile_png_paths = []
            exported_geotiff_path = self._maybe_export_survey_geotiff(export_dir, export_name)
            params_geotiff_path = self._resolve_export_params_geotiff_path(
                exported_geotiff_path, export_dir=export_dir
            )
            export_shapefile = self._export_type_enabled("esri_shapefile")
            export_gpkg = self._export_type_enabled("gpkg")
            export_sis = self._export_type_enabled("sis_asciiplan")
            export_gpx = self._export_type_enabled("gpx")
            export_text_csv = self._export_type_enabled("text_csv")
            export_text_txt = self._export_type_enabled("text_txt")
            export_hypack = self._export_type_enabled("hypack_lnw")
            export_map_png = self._export_map_png_enabled()
            export_profiles_png = self._export_profiles_png_enabled()
            if (export_shapefile or export_gpkg) and mapping is None:
                raise ImportError("shapely.geometry.mapping is required for shapefile/GeoPackage export")

            # --- Build common rows (line_num, line_name, point_label, lat, lon) for reference survey ---
            ref_rows = []
            for i, line in enumerate(self.survey_lines_data):
                if i % 2 == 0:
                    start, end = line[0], line[1]
                else:
                    start, end = line[1], line[0]
                line_name = f'ReferenceLine{i + 1}'
                ref_rows.append((i + 1, line_name, f'L{i+1}S', start[0], start[1]))
                ref_rows.append((i + 1, line_name, f'L{i+1}E', end[0], end[1]))
            if self.cross_line_data:
                ref_rows.append((0, 'Crossline', 'CLS', self.cross_line_data[0][0], self.cross_line_data[0][1]))
                ref_rows.append((0, 'Crossline', 'CLE', self.cross_line_data[1][0], self.cross_line_data[1][1]))

            csv_file_path = os.path.join(export_dir, f"{export_name}_DDD.csv")
            ddm_file_path = os.path.join(export_dir, f"{export_name}_DMM.csv")
            dms_file_path = os.path.join(export_dir, f"{export_name}_DMS.csv")
            ddm_txt_file_path = os.path.join(export_dir, f"{export_name}_DMM.txt")
            dms_txt_file_path = os.path.join(export_dir, f"{export_name}_DMS.txt")
            txt_file_path = os.path.join(export_dir, f"{export_name}_DDD.txt")
            if export_text_csv:
                export_utils.write_ddd_csv(csv_file_path, ref_rows, newline='')
                export_utils.write_dmm_csv(ddm_file_path, ref_rows)
                export_utils.write_dms_csv(dms_file_path, ref_rows)
            if export_text_txt:
                export_utils.write_dmm_txt(ddm_txt_file_path, ref_rows)
                export_utils.write_dms_txt(dms_txt_file_path, ref_rows)
                export_utils.write_ddd_txt(txt_file_path, ref_rows)

            # --- Export to ESRI Shapefile (.shp) and/or GeoPackage (.gpkg) ---
            shapefile_path = os.path.join(export_dir, f"{export_name}.shp")
            if export_shapefile or export_gpkg:
                schema = {
                    'geometry': 'LineString',
                    'properties': {'line_num': 'int', 'line_name': 'str'},
                }
                crs_epsg = 'EPSG:4326'  # WGS 84
                features = []
                # Add main survey lines (names match CSV / GPX: ReferenceLine1, …)
                for i, line_coords in enumerate(self.survey_lines_data):
                    shapely_line = LineString([(p[1], p[0]) for p in line_coords])
                    features.append({
                        'geometry': mapping(shapely_line),
                        'properties': {'line_num': i + 1, 'line_name': f'ReferenceLine{i + 1}'},
                    })
                # Add crossline
                if self.cross_line_data:
                    shapely_cross_line = LineString([(p[1], p[0]) for p in self.cross_line_data])
                    features.append({
                        'geometry': mapping(shapely_cross_line),
                        'properties': {'line_num': 0, 'line_name': 'Crossline'},
                    })
                if export_shapefile:
                    export_utils.remove_export_file(shapefile_path)
                    with fiona.open(shapefile_path, 'w', driver='ESRI Shapefile', crs=crs_epsg, schema=schema) as collection:
                        collection.writerecords(features)
                self._write_gpkg_if_enabled(shapefile_path, schema, features, crs=crs_epsg)

            # --- Export to GeoJSON ---
            geojson_file_path = os.path.join(export_dir, f"{export_name}.geojson")
            try:
                ref_export_speed = float(self.survey_speed_entry.text()) if hasattr(self, 'survey_speed_entry') and self.survey_speed_entry.text() else 8.0
            except (ValueError, TypeError):
                ref_export_speed = 8.0
            geojson_features = []
            # Main survey lines
            for i, line in enumerate(self.survey_lines_data):
                geojson_features.append({
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": [[line[0][1], line[0][0]], [line[1][1], line[1][0]]]
                    },
                    "properties": {
                        "line_num": i + 1,
                        "survey_speed": ref_export_speed,
                        "points": [
                            {"point_num": 1, "lat": line[0][0], "lon": line[0][1]},
                            {"point_num": 2, "lat": line[1][0], "lon": line[1][1]}
                        ]
                    }
                })
            # Crossline
            if self.cross_line_data:
                geojson_features.append({
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": [[self.cross_line_data[0][1], self.cross_line_data[0][0]], [self.cross_line_data[1][1], self.cross_line_data[1][0]]]
                    },
                    "properties": {
                        "line_num": 0,
                        "survey_speed": ref_export_speed,
                        "points": [
                            {"point_num": 1, "lat": self.cross_line_data[0][0], "lon": self.cross_line_data[0][1]},
                            {"point_num": 2, "lat": self.cross_line_data[1][0], "lon": self.cross_line_data[1][1]}
                        ]
                    }
                })
            geojson_collection = {
                "type": "FeatureCollection",
                "properties": {
                    "geotiff_path": params_geotiff_path,
                    "geotiff_nan_value": float(getattr(self, "geotiff_nan_value", -11000.0)),
                },
                "features": geojson_features
            }
            export_utils.remove_export_file(geojson_file_path)
            with open(geojson_file_path, 'w') as f:
                json.dump(geojson_collection, f, indent=2)

            # --- Export to Hypack LNW format (LIN/PTS/UTM), filename includes UTM zone ---
            lnw_file_path = None
            if export_hypack:
                lnw_lines = [(f"LINE{i+1:03d}", [line[0], line[1]]) for i, line in enumerate(self.survey_lines_data)]
                if self.cross_line_data:
                    lnw_lines.append(("CROSS", [self.cross_line_data[0], self.cross_line_data[1]]))
                if lnw_lines:
                    all_pts = [p for _name, pts in lnw_lines for p in pts]
                    zone, hem = export_utils.compute_utm_zone_from_points(all_pts)
                    utm_suffix = f"_UTM{zone}{'N' if hem == 'North' else 'S'}"
                    lnw_file_path = os.path.join(export_dir, f"{export_name}{utm_suffix}.lnw")
                    if not export_utils.write_lnw(lnw_file_path, lnw_lines):
                        lnw_file_path = None

            # --- Export to Kongsberg SIS ASCII Plan format ---
            sis_file_path = os.path.join(export_dir, f"{export_name}.asciiplan")
            ref_ascii_lines = []
            if self.cross_line_data:
                ref_ascii_lines.append(('Crossline', [self.cross_line_data[0], self.cross_line_data[1]]))
            for i, line in enumerate(self.survey_lines_data):
                ref_ascii_lines.append((f'Reference{i + 1}', [line[0], line[1]]))
            if export_sis:
                export_utils.write_asciiplan(sis_file_path, ref_ascii_lines)
            gpx_file_path = os.path.join(export_dir, f"{export_name}.gpx")
            ref_gpx_lines = []
            if self.cross_line_data:
                ref_gpx_lines.append(("Crossline", [self.cross_line_data[0], self.cross_line_data[1]]))
            for i, line in enumerate(self.survey_lines_data):
                if i % 2 == 0:
                    start, end = line[0], line[1]
                else:
                    start, end = line[1], line[0]
                ref_gpx_lines.append((f"ReferenceLine{i + 1}", [start, end]))
            gpx_written = False
            gpx_per_test_names = []
            if export_gpx:
                gpx_written = export_utils.write_gpx(gpx_file_path, ref_gpx_lines, creator="SAT Planner Accuracy")
                acc_gpx_tests = []
                if self.cross_line_data:
                    acc_gpx_tests.append(
                        ("Crossline", "Crossline", [self.cross_line_data[0], self.cross_line_data[1]])
                    )
                for i, line in enumerate(self.survey_lines_data):
                    if i % 2 == 0:
                        start, end = line[0], line[1]
                    else:
                        start, end = line[1], line[0]
                    suf = f"Line{i + 1:02d}"
                    acc_gpx_tests.append((suf, f"ReferenceLine{i + 1}", [start, end]))
                gpx_per_test_names = export_utils.write_gpx_per_test_files(
                    export_dir, export_name, acc_gpx_tests, creator="SAT Planner Accuracy"
                )

            # --- Export Accuracy Survey Information (same content as Accuracy Survey Info dialog) ---
            stats_file_path = os.path.join(export_dir, f"{export_name}_info.txt")
            info_text = None
            if hasattr(self, "_build_reference_planning_info_text"):
                info_text = self._build_reference_planning_info_text(
                    export_name=export_name,
                    export_date=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                )
            if not info_text:
                info_text = "ACCURACY SURVEY INFORMATION\nUnable to generate survey info text.\n"

            def _normalize_degree_symbols_for_export(text):
                normalized_lines = []
                for line in str(text).splitlines():
                    if line.endswith("°"):
                        line = line[:-1]
                    line = line.replace("°", " ")
                    normalized_lines.append(line)
                return "\n".join(normalized_lines) + "\n"

            export_utils.remove_export_file(stats_file_path)
            with open(stats_file_path, 'w', encoding='utf-8') as f:
                f.write(_normalize_degree_symbols_for_export(info_text))
            total_survey_time = self._calculate_total_survey_time()
            include_crossline = (
                total_survey_time.get('num_crossline_passes', 0) > 0
                and total_survey_time.get('crossline_total_distance_m', 0) > 0
            )

            # --- Export Survey Plan as PNG ---
            map_png_path = os.path.join(export_dir, f"{export_name}_map.png")
            if export_map_png:
                if hasattr(self, "_hide_map_hover_tooltip_for_export"):
                    self._hide_map_hover_tooltip_for_export()
                self._save_export_map_png(
                    map_png_path, dpi=300, bbox_inches='tight', facecolor='white'
                )
            
            # --- Export Profile Plots as PNG (crossline + each main line) ---
            if export_profiles_png and hasattr(self, 'profile_fig') and self.profile_fig is not None and hasattr(self, "_draw_segment_profile"):
                combo = getattr(self, "ref_profile_select_combo", None)
                previous_profile_selection = combo.currentText().strip() if combo is not None and combo.count() > 0 else "Crossline"

                # Crossline profile first (when crossline is enabled/present in the active plan)
                if include_crossline and self.cross_line_data and len(self.cross_line_data) == 2:
                    self._draw_segment_profile(self.cross_line_data, "Crossline Elevation Profile", "navy")
                    crossline_profile_png_path = os.path.join(export_dir, f"{export_name}_profile_crossline.png")
                    self._save_export_profile_png(
                        crossline_profile_png_path, dpi=300, bbox_inches='tight', facecolor='white'
                    )
                    profile_png_paths.append(crossline_profile_png_path)

                # Main line profiles
                for i, line_data in enumerate(self.survey_lines_data or []):
                    if not line_data or len(line_data) != 2:
                        continue
                    line_num = i + 1
                    self._draw_segment_profile(line_data, f"Main Line {line_num} Elevation Profile", "blue")
                    line_profile_png_path = os.path.join(
                        export_dir,
                        f"{export_name}_profile_main_line_{line_num:02d}.png"
                    )
                    self._save_export_profile_png(
                        line_profile_png_path, dpi=300, bbox_inches='tight', facecolor='white'
                    )
                    profile_png_paths.append(line_profile_png_path)

                # Restore previously selected profile for user continuity
                if combo is not None:
                    combo.blockSignals(True)
                    if combo.findText(previous_profile_selection) >= 0:
                        combo.setCurrentText(previous_profile_selection)
                    combo.blockSignals(False)
                if hasattr(self, "_draw_current_profile"):
                    self._draw_current_profile()

            # --- Export crossline elevation profile as CSV (Distance m, Elevation m, Slope deg) ---
            if export_text_csv and include_crossline and self.cross_line_data and len(self.cross_line_data) == 2:
                cl_a, cl_b = self.cross_line_data
                prof = self._profile_arrays_along_segment_endpoints(cl_a[0], cl_a[1], cl_b[0], cl_b[1])
                if prof[0] is not None:
                    profile_csv_path = os.path.join(export_dir, f"{export_name}_profile.csv")
                    export_utils.write_profile_csv(profile_csv_path, prof[0], prof[1], prof[2])

            # --- Export parameters metadata as JSON ---
            json_metadata_path = os.path.join(export_dir, f"{export_name}_params.json")
            try:
                # Get current parameter values (even if not all fields are filled, use what we have)
                params = {}
                try:
                    params['central_lat'] = float(self.central_lat_entry.text())
                except:
                    # Calculate from survey lines if not available
                    if self.survey_lines_data:
                        all_lats = [p[0] for line in self.survey_lines_data for p in line]
                        params['central_lat'] = (min(all_lats) + max(all_lats)) / 2.0
                    else:
                        params['central_lat'] = None
                
                try:
                    params['central_lon'] = float(self.central_lon_entry.text())
                except:
                    if self.survey_lines_data:
                        all_lons = [p[1] for line in self.survey_lines_data for p in line]
                        params['central_lon'] = (min(all_lons) + max(all_lons)) / 2.0
                    else:
                        params['central_lon'] = None
                
                try:
                    params['line_length'] = float(self.line_length_entry.text())
                except:
                    params['line_length'] = None
                
                try:
                    params['heading'] = float(self.heading_entry.text())
                except:
                    params['heading'] = None
                
                try:
                    params['dist_between_lines'] = float(self.dist_between_lines_entry.text())
                except:
                    params['dist_between_lines'] = None
                
                try:
                    params['num_lines'] = int(self.num_lines_entry.text())
                except:
                    params['num_lines'] = len(self.survey_lines_data) if self.survey_lines_data else None
                
                try:
                    params['bisect_lead'] = float(self.bisect_lead_entry.text())
                except:
                    params['bisect_lead'] = None
                
                try:
                    params['survey_speed'] = float(self.survey_speed_entry.text())
                except:
                    params['survey_speed'] = 8.0  # Default

                params['geotiff_path'] = params_geotiff_path
                params['geotiff_nan_value'] = float(getattr(self, 'geotiff_nan_value', -11000.0))
                params['show_contours_var'] = bool(getattr(self, 'show_contours_var', False))
                params['contour_interval_m'] = (
                    float(self.contour_interval_entry.text())
                    if hasattr(self, 'contour_interval_entry') and self.contour_interval_entry.text()
                    else 200.0
                )
                params['central_point_depth_m'] = (
                    float(self._depth_at_picked_point)
                    if hasattr(self, '_depth_at_picked_point') and self._depth_at_picked_point is not None
                    else None
                )
                params['visualization_shapefile_paths'] = list(getattr(self, 'visualization_shapefile_paths', []) or [])
                
                try:
                    params['crossline_passes'] = int(self.crossline_passes_entry.text())
                except:
                    params['crossline_passes'] = 2  # Default
                
                try:
                    params['export_name'] = self.export_name_entry.text().strip()
                except:
                    params['export_name'] = export_name
                
                try:
                    params['offset_direction'] = self.offset_direction_var
                except:
                    params['offset_direction'] = 'North'  # Default
                
                try:
                    params['line_length_multiplier'] = self.line_length_multiplier
                except:
                    params['line_length_multiplier'] = 8.0  # Default
                
                try:
                    params['dist_between_lines_multiplier'] = self.dist_between_lines_multiplier
                except:
                    params['dist_between_lines_multiplier'] = 1.0  # Default

                self._add_geotiff_viz_params_to_params(params)

                # Save metadata
                export_utils.remove_export_file(json_metadata_path)
                with open(json_metadata_path, 'w', encoding='utf-8') as f:
                    json.dump(params, f, indent=2)
            except Exception as e:
                # If metadata export fails, continue without it
                print(f"Warning: Could not export metadata: {e}")
            
            # Update success message to include stats file, text file, PNG files, and metadata
            success_files = [f"- {os.path.basename(geojson_file_path)}"]
            if export_text_csv:
                success_files.extend(
                    [
                        f"- {os.path.basename(csv_file_path)}",
                        f"- {os.path.basename(ddm_file_path)}",
                        f"- {os.path.basename(dms_file_path)}",
                    ]
                )
            if export_text_txt:
                success_files.extend(
                    [
                        f"- {os.path.basename(txt_file_path)}",
                        f"- {os.path.basename(ddm_txt_file_path)}",
                        f"- {os.path.basename(dms_txt_file_path)}",
                    ]
                )
            if export_shapefile:
                success_files.append(f"- {os.path.basename(shapefile_path)} (and associated files)")
            if lnw_file_path:
                success_files.append(f"- {os.path.basename(lnw_file_path)}")
            if gpx_written:
                success_files.append(f"- {os.path.basename(gpx_file_path)}")
            for bn in gpx_per_test_names:
                success_files.append(f"- {bn}")
            if export_sis:
                success_files.append(f"- {os.path.basename(sis_file_path)}")
            success_files.append(f"- {os.path.basename(stats_file_path)}")
            if export_map_png:
                for bn in self._map_png_export_basenames(map_png_path):
                    success_files.append(f"- {bn}")
            
            # Add metadata JSON if it was created
            try:
                if os.path.exists(json_metadata_path):
                    success_files.append(f"- {os.path.basename(json_metadata_path)}")
            except:
                pass
            
            # Add profile PNGs if they were created
            for profile_png_path in profile_png_paths:
                for bn in self._profile_png_export_basenames(profile_png_path):
                    success_files.append(f"- {bn}")
            if profile_csv_path and os.path.isfile(profile_csv_path):
                success_files.append(f"- {os.path.basename(profile_csv_path)}")
            if exported_geotiff_path and os.path.isfile(exported_geotiff_path):
                success_files.append(f"- {os.path.basename(exported_geotiff_path)}")

            # Per-file export status block for activity log (OK / FAILED).
            status_lines = []
            def _add_status(path):
                if not path:
                    return
                status_lines.append(
                    f"{'OK' if os.path.exists(path) else 'FAILED'}: {os.path.basename(path)}"
                )
            _add_status(geojson_file_path)
            if export_text_csv:
                _add_status(csv_file_path)
                _add_status(ddm_file_path)
                _add_status(dms_file_path)
            if export_text_txt:
                _add_status(txt_file_path)
                _add_status(ddm_txt_file_path)
                _add_status(dms_txt_file_path)
            if export_shapefile:
                _add_status(shapefile_path)
            _add_status(lnw_file_path)
            if export_sis:
                _add_status(sis_file_path)
            if gpx_written:
                _add_status(gpx_file_path)
            _add_status(stats_file_path)
            if export_map_png:
                _add_status(map_png_path)
            _add_status(json_metadata_path)
            for profile_png_path in profile_png_paths:
                _add_status(profile_png_path)
            _add_status(profile_csv_path)
            if exported_geotiff_path:
                _add_status(exported_geotiff_path)
            if status_lines:
                self.set_ref_info_text(
                    "Accuracy export results:\n" + "\n".join(status_lines),
                    append=True,
                )

            self.set_ref_info_text(
                f"Survey exported successfully to:\n" + "\n".join(success_files) + 
                f"\nin directory: {export_dir}", append=False)

        except Exception as e:
            self._show_message("error","Export Error", f"Failed to export survey files: {e}")

    def _export_performance_survey_files(self):
        if hasattr(self, "_commit_all_deferred_line_edits"):
            self._commit_all_deferred_line_edits()
        """Export performance swath lines and BIST segments (same products as Accuracy export)."""
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded. Cannot export.")
            return
        export_shapefile = self._export_type_enabled("esri_shapefile")
        export_gpkg = self._export_type_enabled("gpkg")
        export_sis = self._export_type_enabled("sis_asciiplan")
        export_gpx = self._export_type_enabled("gpx")
        export_text_csv = self._export_type_enabled("text_csv")
        export_text_txt = self._export_type_enabled("text_txt")
        export_hypack = self._export_type_enabled("hypack_lnw")
        export_map_png = self._export_map_png_enabled()
        export_profiles_png = self._export_profiles_png_enabled()
        if (export_shapefile or export_gpkg) and mapping is None:
            self._show_message("warning", "Export Error", "Shapely is required for shapefile/GeoPackage export.")
            return

        lines = getattr(self, "performance_test_lines_data", None) or []
        if len(lines) != 4:
            self._show_message("warning", "No Data", "Plot four performance test lines first (Plot Performance Lines).")
            return

        export_name = ""
        if hasattr(self, "performance_export_name_entry"):
            export_name = self.performance_export_name_entry.text().strip()
        if not export_name:
            export_name = (
                self._build_performance_export_basename()
                if hasattr(self, "_build_performance_export_basename")
                else "perf_swell0_depth0m"
            )
        bad = '<>:"/\\|?*'
        for c in bad:
            export_name = export_name.replace(c, "_")
        export_name = export_name.strip().strip(".")

        export_dir = self._select_export_directory(self.last_export_dir)
        if not export_dir:
            return
        self.last_export_dir = export_dir
        self._save_last_export_dir()

        bist_segs = getattr(self, "performance_bist_segments_data", None) or []
        has_bist = len(bist_segs) == 4

        try:
            exported_geotiff_path = self._maybe_export_survey_geotiff(export_dir, export_name)
            geotiff_path = self._resolve_export_params_geotiff_path(
                exported_geotiff_path, export_dir=export_dir
            )
            speed_kts = 8.0
            if hasattr(self, "performance_test_speed_entry"):
                try:
                    speed_kts = float(self.performance_test_speed_entry.text().strip())
                except (ValueError, TypeError):
                    speed_kts = 8.0

            perf_rows = []
            for i, line in enumerate(lines):
                n = i + 1
                start, end = line[0], line[1]
                lname = f"PerformanceLine{n}"
                perf_rows.append((n, lname, f"P{n}S", start[0], start[1]))
                perf_rows.append((n, lname, f"P{n}E", end[0], end[1]))
            if has_bist:
                for i, seg in enumerate(bist_segs):
                    bn = i + 1
                    line_num = 10 + bn
                    lname = f"BISTLine{bn}"
                    perf_rows.append((line_num, lname, f"B{bn}S", seg[0][0], seg[0][1]))
                    perf_rows.append((line_num, lname, f"B{bn}E", seg[1][0], seg[1][1]))

            csv_file_path = os.path.join(export_dir, f"{export_name}_DDD.csv")
            ddm_file_path = os.path.join(export_dir, f"{export_name}_DMM.csv")
            dms_file_path = os.path.join(export_dir, f"{export_name}_DMS.csv")
            ddm_txt_file_path = os.path.join(export_dir, f"{export_name}_DMM.txt")
            dms_txt_file_path = os.path.join(export_dir, f"{export_name}_DMS.txt")
            txt_file_path = os.path.join(export_dir, f"{export_name}_DDD.txt")
            if export_text_csv:
                export_utils.write_ddd_csv(csv_file_path, perf_rows, newline="")
                export_utils.write_dmm_csv(ddm_file_path, perf_rows)
                export_utils.write_dms_csv(dms_file_path, perf_rows)
            if export_text_txt:
                export_utils.write_dmm_txt(ddm_txt_file_path, perf_rows)
                export_utils.write_dms_txt(dms_txt_file_path, perf_rows)
                export_utils.write_ddd_txt(txt_file_path, perf_rows)

            shapefile_path = os.path.join(export_dir, f"{export_name}.shp")
            if export_shapefile or export_gpkg:
                schema = {"geometry": "LineString", "properties": {"line_num": "int", "line_name": "str"}}
                crs_epsg = "EPSG:4326"
                features = []
                for i, line_coords in enumerate(lines):
                    shapely_line = LineString([(p[1], p[0]) for p in line_coords])
                    n = i + 1
                    features.append({
                        "geometry": mapping(shapely_line),
                        "properties": {"line_num": n, "line_name": f"PerformanceLine{n}"},
                    })
                if has_bist:
                    for i, seg in enumerate(bist_segs):
                        shapely_line = LineString([(p[1], p[0]) for p in seg])
                        bn = i + 1
                        features.append({
                            "geometry": mapping(shapely_line),
                            "properties": {"line_num": 10 + bn, "line_name": f"BISTLine{bn}"},
                        })
                if export_shapefile:
                    export_utils.remove_export_file(shapefile_path)
                    with fiona.open(shapefile_path, "w", driver="ESRI Shapefile", crs=crs_epsg, schema=schema) as collection:
                        collection.writerecords(features)
                self._write_gpkg_if_enabled(shapefile_path, schema, features, crs=crs_epsg)

            geojson_features = []
            for i, line in enumerate(lines):
                geojson_features.append(
                    {
                        "type": "Feature",
                        "geometry": {
                            "type": "LineString",
                            "coordinates": [[line[0][1], line[0][0]], [line[1][1], line[1][0]]],
                        },
                        "properties": {
                            "line_num": i + 1,
                            "performance_segment": "swath",
                            "survey_speed": speed_kts,
                            "test_speed_kts": speed_kts,
                            "geotiff_path": geotiff_path,
                            "points": [
                                {"point_num": 1, "lat": line[0][0], "lon": line[0][1]},
                                {"point_num": 2, "lat": line[1][0], "lon": line[1][1]},
                            ],
                        },
                    }
                )
            if has_bist:
                for i, seg in enumerate(bist_segs):
                    geojson_features.append(
                        {
                            "type": "Feature",
                            "geometry": {
                                "type": "LineString",
                                "coordinates": [[seg[0][1], seg[0][0]], [seg[1][1], seg[1][0]]],
                            },
                            "properties": {
                                "line_num": 11 + i,
                                "performance_segment": "bist",
                                "survey_speed": speed_kts,
                                "test_speed_kts": speed_kts,
                                "geotiff_path": geotiff_path,
                                "points": [
                                    {"point_num": 1, "lat": seg[0][0], "lon": seg[0][1]},
                                    {"point_num": 2, "lat": seg[1][0], "lon": seg[1][1]},
                                ],
                            },
                        }
                    )
            geojson_collection = {
                "type": "FeatureCollection",
                "properties": {
                    "geotiff_path": geotiff_path,
                    "geotiff_nan_value": float(getattr(self, "geotiff_nan_value", -11000.0)),
                    "plan_type": "performance",
                },
                "features": geojson_features,
            }
            geojson_file_path = os.path.join(export_dir, f"{export_name}.geojson")
            export_utils.remove_export_file(geojson_file_path)
            with open(geojson_file_path, "w", encoding="utf-8") as f:
                json.dump(geojson_collection, f, indent=2)

            lnw_lines = [(f"PLN{i + 1:03d}", [line[0], line[1]]) for i, line in enumerate(lines)]
            if has_bist:
                for i, seg in enumerate(bist_segs):
                    lnw_lines.append((f"BIST{i + 1:03d}", [seg[0], seg[1]]))
            lnw_file_path = None
            if export_hypack and lnw_lines:
                all_pts = [p for _name, pts in lnw_lines for p in pts]
                zone, hem = export_utils.compute_utm_zone_from_points(all_pts)
                utm_suffix = f"_UTM{zone}{'N' if hem == 'North' else 'S'}"
                lnw_file_path = os.path.join(export_dir, f"{export_name}{utm_suffix}.lnw")
                if not export_utils.write_lnw(lnw_file_path, lnw_lines):
                    lnw_file_path = None

            sis_file_path = os.path.join(export_dir, f"{export_name}.asciiplan")
            perf_ascii_lines = [(f"Performance{i + 1}", [line[0], line[1]]) for i, line in enumerate(lines)]
            if has_bist:
                for i, seg in enumerate(bist_segs):
                    perf_ascii_lines.append((f"BIST{i + 1}", [seg[0], seg[1]]))
            if export_sis:
                export_utils.write_asciiplan(sis_file_path, perf_ascii_lines)

            gpx_file_path = os.path.join(export_dir, f"{export_name}.gpx")
            perf_gpx_lines = [(f"Performance{i + 1}", [line[0], line[1]]) for i, line in enumerate(lines)]
            if has_bist:
                for i, seg in enumerate(bist_segs):
                    perf_gpx_lines.append((f"BIST{i + 1}", [seg[0], seg[1]]))
            gpx_written = False
            gpx_per_test_names = []
            if export_gpx:
                gpx_written = export_utils.write_gpx(
                    gpx_file_path, perf_gpx_lines, creator="SAT Planner Performance"
                )
                perf_gpx_tests = [
                    (f"Performance{i + 1}", f"Performance{i + 1}", [line[0], line[1]])
                    for i, line in enumerate(lines)
                ]
                if has_bist:
                    for i, seg in enumerate(bist_segs):
                        perf_gpx_tests.append((f"BIST{i + 1}", f"BIST{i + 1}", [seg[0], seg[1]]))
                gpx_per_test_names = export_utils.write_gpx_per_test_files(
                    export_dir, export_name, perf_gpx_tests, creator="SAT Planner Performance"
                )

            stats_file_path = os.path.join(export_dir, f"{export_name}_info.txt")
            try:
                turn_min = 10.0
                if hasattr(self, "performance_turn_time_entry"):
                    raw_turn = self.performance_turn_time_entry.text().strip()
                    if raw_turn:
                        turn_min = float(raw_turn)
                if turn_min < 0:
                    turn_min = 0.0
            except Exception:
                turn_min = 10.0

            try:
                perf_info_text = self._build_performance_test_info_text(
                    lines, bist_segs if has_bist else [], float(speed_kts), float(turn_min)
                )
            except Exception as e:
                perf_info_text = (
                    "PERFORMANCE TEST SUMMARY\n"
                    + "=" * 24
                    + "\n"
                    + f"Could not build detailed performance info: {e}\n"
                )

            export_utils.remove_export_file(stats_file_path)
            with open(stats_file_path, "w", encoding="utf-8") as f:
                f.write("PERFORMANCE SURVEY INFORMATION\n")
                f.write("=" * 50 + "\n\n")
                f.write(f"Performance Survey: {export_name}\n")
                f.write(f"Export Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                f.write(perf_info_text)

            json_metadata_path = os.path.join(export_dir, f"{export_name}_performance_params.json")
            pmeta = {
                "export_name": export_name,
                "plan_type": "performance",
                "geotiff_path": geotiff_path,
                "geotiff_nan_value": float(getattr(self, "geotiff_nan_value", -11000.0)),
                "visualization_shapefile_paths": list(getattr(self, 'visualization_shapefile_paths', []) or []),
                "show_contours_var": bool(getattr(self, "show_contours_var", False)),
                "contour_interval_m": (
                    float(self.contour_interval_entry.text())
                    if hasattr(self, "contour_interval_entry") and self.contour_interval_entry.text()
                    else 200.0
                ),
                "test_speed_kts": speed_kts,
            }
            for key, attr in (
                ("central_lat", "performance_central_lat_entry"),
                ("central_lon", "performance_central_lon_entry"),
                ("test_depth_m", "performance_test_depth_entry"),
                ("swell_direction_deg", "performance_swell_direction_entry"),
                ("swath_angle_deg", "performance_swath_angle_entry"),
                ("bist_time_min", "performance_bist_time_entry"),
                ("turn_time_min", "performance_turn_time_entry"),
                ("num_pings", "performance_num_pings_entry"),
                ("sound_velocity", "performance_sound_velocity_entry"),
            ):
                w = getattr(self, attr, None)
                if w is not None and hasattr(w, "text"):
                    try:
                        pmeta[key] = w.text().strip()
                    except Exception:
                        pmeta[key] = None
                else:
                    pmeta[key] = None
            try:
                ll = getattr(self, "performance_line_length_m_entry", None)
                if ll is not None and hasattr(ll, "text"):
                    pmeta["line_length_m_label"] = ll.text().strip()
            except Exception:
                pass
            self._add_geotiff_viz_params_to_params(pmeta)
            export_utils.remove_export_file(json_metadata_path)
            with open(json_metadata_path, "w", encoding="utf-8") as f:
                json.dump(pmeta, f, indent=2)

            map_png_path = os.path.join(export_dir, f"{export_name}_map.png")
            if export_map_png:
                if hasattr(self, "_hide_map_hover_tooltip_for_export"):
                    self._hide_map_hover_tooltip_for_export()
                self._save_export_map_png(
                    map_png_path, dpi=300, bbox_inches="tight", facecolor="white"
                )

            profile_png_path = os.path.join(export_dir, f"{export_name}_profile.png")
            if export_profiles_png and hasattr(self, "profile_fig") and self.profile_fig is not None:
                self._save_export_profile_png(
                    profile_png_path, dpi=300, bbox_inches="tight", facecolor="white"
                )

            msg_lines = [f"- {os.path.basename(geojson_file_path)}"]
            if export_text_csv:
                msg_lines.extend(
                    [
                        f"- {os.path.basename(csv_file_path)}",
                        f"- {os.path.basename(ddm_file_path)}",
                        f"- {os.path.basename(dms_file_path)}",
                    ]
                )
            if export_text_txt:
                msg_lines.extend(
                    [
                        f"- {os.path.basename(txt_file_path)}",
                        f"- {os.path.basename(ddm_txt_file_path)}",
                        f"- {os.path.basename(dms_txt_file_path)}",
                    ]
                )
            if export_shapefile:
                msg_lines.append(f"- {os.path.basename(shapefile_path)} (and sidecars)")
            if gpx_written:
                msg_lines.append(f"- {os.path.basename(gpx_file_path)}")
            for bn in gpx_per_test_names:
                msg_lines.append(f"- {bn}")
            if lnw_file_path:
                msg_lines.append(f"- {os.path.basename(lnw_file_path)}")
            if export_sis:
                msg_lines.append(f"- {os.path.basename(sis_file_path)}")
            msg_lines.append(f"- {os.path.basename(stats_file_path)}")
            if export_map_png:
                for bn in self._map_png_export_basenames(map_png_path):
                    msg_lines.append(f"- {bn}")
            msg_lines.append(f"- {os.path.basename(json_metadata_path)}")
            if export_profiles_png and hasattr(self, "profile_fig") and self.profile_fig is not None:
                for bn in self._profile_png_export_basenames(profile_png_path):
                    msg_lines.append(f"- {bn}")
            if exported_geotiff_path and os.path.isfile(exported_geotiff_path):
                msg_lines.append(f"- {os.path.basename(exported_geotiff_path)}")
            status_lines = []
            def _add_status(path):
                if not path:
                    return
                status_lines.append(
                    f"{'OK' if os.path.exists(path) else 'FAILED'}: {os.path.basename(path)}"
                )
            _add_status(geojson_file_path)
            if export_text_csv:
                _add_status(csv_file_path)
                _add_status(ddm_file_path)
                _add_status(dms_file_path)
            if export_text_txt:
                _add_status(txt_file_path)
                _add_status(ddm_txt_file_path)
                _add_status(dms_txt_file_path)
            if export_shapefile:
                _add_status(shapefile_path)
            _add_status(lnw_file_path)
            if export_sis:
                _add_status(sis_file_path)
            if gpx_written:
                _add_status(gpx_file_path)
            _add_status(stats_file_path)
            if export_map_png:
                _add_status(map_png_path)
            _add_status(json_metadata_path)
            if export_profiles_png and hasattr(self, "profile_fig") and self.profile_fig is not None:
                _add_status(profile_png_path)
            if exported_geotiff_path:
                _add_status(exported_geotiff_path)
            if status_lines:
                if hasattr(self, "set_performance_activity_text"):
                    self.set_performance_activity_text(
                        "Performance export results:\n" + "\n".join(status_lines),
                        append=True,
                    )
                else:
                    self.set_ref_info_text(
                        "Performance export results:\n" + "\n".join(status_lines),
                        append=True,
                    )
            log_msg = "Performance survey exported successfully to:\n" + "\n".join(msg_lines) + f"\nin directory: {export_dir}"
            if hasattr(self, "set_performance_activity_text"):
                self.set_performance_activity_text(log_msg, append=False)
            else:
                self.set_ref_info_text(log_msg, append=False)

        except Exception as e:
            self._show_message("error", "Export Error", f"Failed to export performance survey: {e}")

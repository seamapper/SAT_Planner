"""Import Survey dialog with per-tab GMRT download options."""

from PyQt6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QRadioButton,
    QVBoxLayout,
)


class ImportSurveyDialog(QDialog):
    """Offer GMRT bathymetry download options after a survey import."""

    def __init__(
        self,
        parent,
        tab_label,
        download_checkbox,
        buffer_spin,
        buffer_percent_spin,
        buffer_mode_widget,
        split_checkbox,
        cell_size_combo=None,
        *,
        prompt_missing_geotiff=False,
        geotiff_path=None,
        existing_grid_loaded=False,
    ):
        super().__init__(parent)
        self.setWindowTitle(f"Import {tab_label}")
        self.setModal(True)
        self.resize(480, 300)

        self._download_checkbox = download_checkbox
        self._buffer_spin = buffer_spin
        self._buffer_percent_spin = buffer_percent_spin
        self._buffer_mode_widget = buffer_mode_widget
        self._split_checkbox = split_checkbox
        self._cell_size_combo = cell_size_combo
        self._existing_grid_loaded = bool(existing_grid_loaded and prompt_missing_geotiff)

        layout = QVBoxLayout(self)

        if prompt_missing_geotiff and self._existing_grid_loaded:
            if geotiff_path:
                intro_text = (
                    "The planning GeoTIFF referenced by this survey was not found:\n"
                    f"{geotiff_path}\n\n"
                    "A bathymetry grid is already loaded. Keep it, or download a new "
                    "GMRT grid for this survey area."
                )
            else:
                intro_text = (
                    "This survey does not include planning bathymetry (GeoTIFF).\n\n"
                    "A bathymetry grid is already loaded. Keep it, or download a new "
                    "GMRT grid for this survey area."
                )
            default_download = False
        elif prompt_missing_geotiff:
            if geotiff_path:
                intro_text = (
                    "The planning GeoTIFF referenced by this survey was not found:\n"
                    f"{geotiff_path}\n\n"
                    "Download GMRT bathymetry for the survey area, or skip to continue "
                    "without bathymetry."
                )
            else:
                intro_text = (
                    "This survey does not include planning bathymetry (GeoTIFF).\n\n"
                    "Download GMRT bathymetry for the survey area, or skip to continue "
                    "without bathymetry."
                )
            default_download = True
        else:
            intro_text = (
                "Choose GMRT options for this import, then click Continue to select a survey file."
            )
            default_download = download_checkbox.isChecked()

        intro = QLabel(intro_text)
        intro.setWordWrap(True)
        layout.addWidget(intro)

        gmrt_group = QGroupBox("GMRT Bathymetry Download")
        gmrt_layout = QVBoxLayout(gmrt_group)

        self.keep_existing_radio = None
        self.download_new_radio = None
        if self._existing_grid_loaded:
            choice_group = QButtonGroup(self)
            self.keep_existing_radio = QRadioButton("Keep existing grid")
            self.download_new_radio = QRadioButton("Download new GMRT grid")
            choice_group.addButton(self.keep_existing_radio)
            choice_group.addButton(self.download_new_radio)
            self.keep_existing_radio.setChecked(True)
            self.download_new_radio.setChecked(False)
            gmrt_layout.addWidget(self.keep_existing_radio)
            gmrt_layout.addWidget(self.download_new_radio)
            self.download_gmrt_checkbox = QCheckBox("Download GMRT bathymetry for survey area")
            self.download_gmrt_checkbox.setChecked(False)
            self.download_gmrt_checkbox.hide()
        else:
            self.download_gmrt_checkbox = QCheckBox("Download GMRT bathymetry for survey area")
            self.download_gmrt_checkbox.setChecked(default_download)
            self.download_gmrt_checkbox.setToolTip(download_checkbox.toolTip())
            gmrt_layout.addWidget(self.download_gmrt_checkbox)

        cell_row = QHBoxLayout()
        self.cell_size_label = QLabel("Cell Size (m):")
        cell_row.addWidget(self.cell_size_label)
        self.gmrt_cell_size_combo = QComboBox()
        if cell_size_combo is not None:
            for i in range(cell_size_combo.count()):
                self.gmrt_cell_size_combo.addItem(
                    cell_size_combo.itemText(i), cell_size_combo.itemData(i)
                )
            self.gmrt_cell_size_combo.setCurrentIndex(cell_size_combo.currentIndex())
            self.gmrt_cell_size_combo.setToolTip(cell_size_combo.toolTip())
        else:
            for meters in (60, 120, 240, 480, 960):
                self.gmrt_cell_size_combo.addItem(f"{meters} m", meters)
            self.gmrt_cell_size_combo.setCurrentIndex(0)  # 60 m
            self.gmrt_cell_size_combo.setToolTip(
                "GMRT GridServer cell size (meters/pixel). Same presets as Download Bathymetry."
            )
        self.gmrt_cell_size_combo.setMinimumWidth(100)
        cell_row.addWidget(self.gmrt_cell_size_combo)
        cell_row.addStretch()
        gmrt_layout.addLayout(cell_row)

        mode_row = QHBoxLayout()
        mode_row.addWidget(QLabel("Buffer:"))
        self.buffer_mode_degrees = QRadioButton("Degrees")
        self.buffer_mode_percent = QRadioButton("Percent of survey")
        self.buffer_mode_group = QButtonGroup(self)
        self.buffer_mode_group.addButton(self.buffer_mode_degrees)
        self.buffer_mode_group.addButton(self.buffer_mode_percent)
        mode_row.addWidget(self.buffer_mode_degrees)
        mode_row.addWidget(self.buffer_mode_percent)
        mode_row.addStretch()
        gmrt_layout.addLayout(mode_row)

        current_mode = "percent"
        if buffer_mode_widget is not None and hasattr(buffer_mode_widget, "currentData"):
            data = buffer_mode_widget.currentData()
            if data:
                current_mode = str(data)
        elif buffer_mode_widget is not None and hasattr(buffer_mode_widget, "currentText"):
            text = (buffer_mode_widget.currentText() or "").lower()
            if "degree" in text:
                current_mode = "degrees"
            elif "percent" in text:
                current_mode = "percent"
        if current_mode == "degrees":
            self.buffer_mode_degrees.setChecked(True)
        else:
            self.buffer_mode_percent.setChecked(True)

        deg_row = QHBoxLayout()
        self.buffer_deg_label = QLabel("Buffer (deg):")
        deg_row.addWidget(self.buffer_deg_label)
        self.gmrt_buffer_spin = QDoubleSpinBox()
        self.gmrt_buffer_spin.setRange(buffer_spin.minimum(), buffer_spin.maximum())
        self.gmrt_buffer_spin.setSingleStep(buffer_spin.singleStep())
        self.gmrt_buffer_spin.setDecimals(buffer_spin.decimals())
        self.gmrt_buffer_spin.setValue(buffer_spin.value())
        self.gmrt_buffer_spin.setToolTip(buffer_spin.toolTip())
        self.gmrt_buffer_spin.setMinimumWidth(80)
        deg_row.addWidget(self.gmrt_buffer_spin)
        deg_row.addStretch()
        gmrt_layout.addLayout(deg_row)

        pct_row = QHBoxLayout()
        self.buffer_percent_label = QLabel("Buffer (%):")
        pct_row.addWidget(self.buffer_percent_label)
        self.gmrt_buffer_percent_spin = QDoubleSpinBox()
        self.gmrt_buffer_percent_spin.setRange(
            buffer_percent_spin.minimum(), buffer_percent_spin.maximum()
        )
        self.gmrt_buffer_percent_spin.setSingleStep(buffer_percent_spin.singleStep())
        self.gmrt_buffer_percent_spin.setDecimals(buffer_percent_spin.decimals())
        self.gmrt_buffer_percent_spin.setValue(buffer_percent_spin.value())
        self.gmrt_buffer_percent_spin.setToolTip(buffer_percent_spin.toolTip())
        self.gmrt_buffer_percent_spin.setMinimumWidth(80)
        self.gmrt_buffer_percent_spin.setSuffix(" %")
        pct_row.addWidget(self.gmrt_buffer_percent_spin)
        pct_row.addStretch()
        gmrt_layout.addLayout(pct_row)

        self.split_topo_depths_checkbox = QCheckBox("Split Topo/Depths")
        self.split_topo_depths_checkbox.setChecked(split_checkbox.isChecked())
        self.split_topo_depths_checkbox.setToolTip(split_checkbox.toolTip())
        gmrt_layout.addWidget(self.split_topo_depths_checkbox)

        if self._existing_grid_loaded:
            self.keep_existing_radio.toggled.connect(self._update_enabled_state)
            self.download_new_radio.toggled.connect(self._update_enabled_state)
        else:
            self.download_gmrt_checkbox.toggled.connect(self._update_enabled_state)
        self.buffer_mode_degrees.toggled.connect(self._update_enabled_state)
        self.buffer_mode_percent.toggled.connect(self._update_enabled_state)
        self._update_enabled_state()

        layout.addWidget(gmrt_group)

        button_box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        button_box.button(QDialogButtonBox.StandardButton.Ok).setText(
            "Continue" if not prompt_missing_geotiff else "OK"
        )
        button_box.accepted.connect(self._on_accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def _buffer_mode(self):
        return "percent" if self.buffer_mode_percent.isChecked() else "degrees"

    def _download_selected(self):
        if self._existing_grid_loaded and self.download_new_radio is not None:
            return self.download_new_radio.isChecked()
        return self.download_gmrt_checkbox.isChecked()

    def _update_enabled_state(self):
        download_on = self._download_selected()
        mode = self._buffer_mode()
        self.gmrt_cell_size_combo.setEnabled(download_on)
        self.cell_size_label.setEnabled(download_on)
        self.buffer_mode_degrees.setEnabled(download_on)
        self.buffer_mode_percent.setEnabled(download_on)
        self.gmrt_buffer_spin.setEnabled(download_on and mode == "degrees")
        self.buffer_deg_label.setEnabled(download_on and mode == "degrees")
        self.gmrt_buffer_percent_spin.setEnabled(download_on and mode == "percent")
        self.buffer_percent_label.setEnabled(download_on and mode == "percent")
        self.split_topo_depths_checkbox.setEnabled(download_on)

    def _on_accept(self):
        self._download_checkbox.setChecked(self._download_selected())
        self._buffer_spin.setValue(self.gmrt_buffer_spin.value())
        self._buffer_percent_spin.setValue(self.gmrt_buffer_percent_spin.value())
        mode = self._buffer_mode()
        if self._buffer_mode_widget is not None and hasattr(self._buffer_mode_widget, "findData"):
            idx = self._buffer_mode_widget.findData(mode)
            if idx >= 0:
                self._buffer_mode_widget.setCurrentIndex(idx)
        if self._cell_size_combo is not None:
            idx = self._cell_size_combo.findData(self.gmrt_cell_size_combo.currentData())
            if idx >= 0:
                self._cell_size_combo.setCurrentIndex(idx)
            else:
                self._cell_size_combo.setCurrentIndex(self.gmrt_cell_size_combo.currentIndex())
        self._split_checkbox.setChecked(self.split_topo_depths_checkbox.isChecked())
        self.accept()

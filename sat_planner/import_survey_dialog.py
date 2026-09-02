"""Import Survey dialog with per-tab GMRT download options."""

from PyQt6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QGroupBox,
    QHBoxLayout,
    QLabel,
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
        split_checkbox,
        *,
        prompt_missing_geotiff=False,
        geotiff_path=None,
    ):
        super().__init__(parent)
        self.setWindowTitle(f"Import {tab_label}")
        self.setModal(True)
        self.resize(460, 200)

        self._download_checkbox = download_checkbox
        self._buffer_spin = buffer_spin
        self._split_checkbox = split_checkbox

        layout = QVBoxLayout(self)

        if prompt_missing_geotiff:
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

        self.download_gmrt_checkbox = QCheckBox("Download GMRT bathymetry for survey area")
        self.download_gmrt_checkbox.setChecked(default_download)
        self.download_gmrt_checkbox.setToolTip(download_checkbox.toolTip())
        gmrt_layout.addWidget(self.download_gmrt_checkbox)

        buffer_row = QHBoxLayout()
        buffer_row.addWidget(QLabel("Buffer (deg):"))
        self.gmrt_buffer_spin = QDoubleSpinBox()
        self.gmrt_buffer_spin.setRange(buffer_spin.minimum(), buffer_spin.maximum())
        self.gmrt_buffer_spin.setSingleStep(buffer_spin.singleStep())
        self.gmrt_buffer_spin.setDecimals(buffer_spin.decimals())
        self.gmrt_buffer_spin.setValue(buffer_spin.value())
        self.gmrt_buffer_spin.setToolTip(buffer_spin.toolTip())
        self.gmrt_buffer_spin.setMinimumWidth(80)
        buffer_row.addWidget(self.gmrt_buffer_spin)
        buffer_row.addStretch()
        gmrt_layout.addLayout(buffer_row)

        self.split_topo_depths_checkbox = QCheckBox("Split Topo/Depths")
        self.split_topo_depths_checkbox.setChecked(split_checkbox.isChecked())
        self.split_topo_depths_checkbox.setToolTip(split_checkbox.toolTip())
        self.split_topo_depths_checkbox.setEnabled(self.download_gmrt_checkbox.isChecked())
        self.download_gmrt_checkbox.toggled.connect(self.split_topo_depths_checkbox.setEnabled)
        self.download_gmrt_checkbox.toggled.connect(self.gmrt_buffer_spin.setEnabled)
        self.gmrt_buffer_spin.setEnabled(self.download_gmrt_checkbox.isChecked())
        gmrt_layout.addWidget(self.split_topo_depths_checkbox)

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

    def _on_accept(self):
        self._download_checkbox.setChecked(self.download_gmrt_checkbox.isChecked())
        self._buffer_spin.setValue(self.gmrt_buffer_spin.value())
        self._split_checkbox.setChecked(self.split_topo_depths_checkbox.isChecked())
        self.accept()

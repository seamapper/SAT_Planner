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
    """Configure GMRT import options, then continue to the file picker."""

    def __init__(self, parent, tab_label, download_checkbox, buffer_spin, split_checkbox):
        super().__init__(parent)
        self.setWindowTitle(f"Import {tab_label}")
        self.setModal(True)
        self.resize(420, 180)

        self._download_checkbox = download_checkbox
        self._buffer_spin = buffer_spin
        self._split_checkbox = split_checkbox

        layout = QVBoxLayout(self)

        intro = QLabel(
            "Choose GMRT options for this import, then click Continue to select a survey file."
        )
        intro.setWordWrap(True)
        layout.addWidget(intro)

        gmrt_group = QGroupBox("GMRT Bathymetry Download")
        gmrt_layout = QVBoxLayout(gmrt_group)

        self.download_gmrt_checkbox = QCheckBox("Download GMRT after import")
        self.download_gmrt_checkbox.setChecked(download_checkbox.isChecked())
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
        gmrt_layout.addWidget(self.split_topo_depths_checkbox)

        layout.addWidget(gmrt_group)

        button_box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        button_box.button(QDialogButtonBox.StandardButton.Ok).setText("Continue")
        button_box.accepted.connect(self._on_accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def _on_accept(self):
        self._download_checkbox.setChecked(self.download_gmrt_checkbox.isChecked())
        self._buffer_spin.setValue(self.gmrt_buffer_spin.value())
        self._split_checkbox.setChecked(self.split_topo_depths_checkbox.isChecked())
        self.accept()

"""Dialog for editing dynamic vertical exaggeration breakpoint tables."""

from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QDialogButtonBox,
    QMessageBox,
    QAbstractItemView,
    QHeaderView,
)
from PyQt6.QtCore import Qt


class DynVertExagDialog(QDialog):
    """Edit elevation-range breakpoints and vertical exaggeration values."""

    _COLUMNS = ("Elevation Range (m)", "Shaded Relief", "Shaded Relief Dyn")

    def __init__(self, parent, table_rows, default_rows):
        super().__init__(parent)
        self.setWindowTitle("Dynamic Vertical Exaggeration")
        self.setModal(True)
        self.resize(520, 320)
        self._default_rows = [dict(row) for row in default_rows]
        self._result_rows = None
        self._updating_table = False

        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)

        instructions = QLabel(
            "Each elevation range is the minimum value for that segment. "
            "The first range is always 0 m. Later ranges must increase. "
            "Edit elevation range and Shaded Relief; Shaded Relief Dyn updates automatically."
        )
        instructions.setWordWrap(True)
        layout.addWidget(instructions)

        self.table = QTableWidget(len(table_rows), len(self._COLUMNS), self)
        self.table.setHorizontalHeaderLabels(self._COLUMNS)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self.table.setEditTriggers(
            QAbstractItemView.EditTrigger.DoubleClicked
            | QAbstractItemView.EditTrigger.SelectedClicked
            | QAbstractItemView.EditTrigger.EditKeyPressed
        )
        self.table.itemChanged.connect(self._on_item_changed)
        layout.addWidget(self.table)

        self._populate_table(table_rows)

        button_row = QHBoxLayout()
        self.reset_btn = QPushButton("Reset to Defaults")
        self.reset_btn.clicked.connect(self._reset_to_defaults)
        button_row.addWidget(self.reset_btn)
        button_row.addStretch()
        layout.addLayout(button_row)

        self.button_box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.button_box.accepted.connect(self._on_accept)
        self.button_box.rejected.connect(self.reject)
        layout.addWidget(self.button_box)

    def _dyn_offset(self, row_idx):
        default = self._default_rows[row_idx]
        return default["shaded_relief_dyn"] - default["shaded_relief"]

    def _compute_dyn_value(self, row_idx, shaded_relief):
        return shaded_relief + self._dyn_offset(row_idx)

    def _read_only_item(self, text):
        item = QTableWidgetItem(text)
        item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
        return item

    def _populate_table(self, table_rows):
        self._updating_table = True
        try:
            for row_idx, row in enumerate(table_rows):
                elev_item = QTableWidgetItem(f"{row['elevation_range']:.4g}")
                elev_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                if row_idx == 0:
                    elev_item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
                self.table.setItem(row_idx, 0, elev_item)

                sr_item = QTableWidgetItem(f"{row['shaded_relief']:.4g}")
                sr_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                self.table.setItem(row_idx, 1, sr_item)

                dyn_value = self._compute_dyn_value(row_idx, row["shaded_relief"])
                self.table.setItem(row_idx, 2, self._read_only_item(f"{dyn_value:.4g}"))
        finally:
            self._updating_table = False

    def _reset_to_defaults(self):
        self._populate_table(self._default_rows)

    def _on_item_changed(self, item):
        if self._updating_table or item.column() != 1:
            return
        row_idx = item.row()
        try:
            shaded_relief = float(item.text().strip())
        except ValueError:
            return
        self._updating_table = True
        try:
            dyn_value = self._compute_dyn_value(row_idx, shaded_relief)
            self.table.setItem(row_idx, 2, self._read_only_item(f"{dyn_value:.4g}"))
        finally:
            self._updating_table = False

    def _parse_table(self):
        rows = []
        row_count = self.table.rowCount()
        for row_idx in range(row_count):
            elev_text = self.table.item(row_idx, 0).text().strip()
            sr_text = self.table.item(row_idx, 1).text().strip()
            try:
                elevation_range = float(elev_text)
                shaded_relief = float(sr_text)
            except ValueError as exc:
                raise ValueError(
                    f"Row {row_idx + 1}: enter valid numbers for elevation range and Shaded Relief."
                ) from exc
            shaded_relief_dyn = self._compute_dyn_value(row_idx, shaded_relief)
            rows.append({
                "elevation_range": elevation_range,
                "shaded_relief": shaded_relief,
                "shaded_relief_dyn": shaded_relief_dyn,
            })

        if rows[0]["elevation_range"] != 0.0:
            raise ValueError("The first elevation range must be 0 m.")

        for row_idx in range(1, len(rows)):
            if rows[row_idx]["elevation_range"] <= rows[row_idx - 1]["elevation_range"]:
                raise ValueError(
                    f"Row {row_idx + 1}: elevation range must be greater than "
                    f"{rows[row_idx - 1]['elevation_range']:.4g} m."
                )

        return rows

    def _on_accept(self):
        try:
            self._result_rows = self._parse_table()
        except ValueError as exc:
            QMessageBox.warning(self, "Invalid Values", str(exc))
            return
        self.accept()

    def get_table_rows(self):
        return self._result_rows

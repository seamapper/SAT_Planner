"""Dialog for assigning imported line segments to performance swath lines and BIST legs."""

from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QComboBox,
    QPushButton,
    QWidget,
    QMessageBox,
)

from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas


class PerformanceLineAssignmentDialog(QDialog):
    """Map each imported 2-point segment to Swath 1–4, BIST 1–4, or Ignore.

    Rules: Swath lines 1–4 must each be assigned exactly once.
    BIST: either assign all four BIST legs exactly once, or assign none (optional BIST).
    """

    _SWATH_LABELS = ("Swath 1", "Swath 2", "Swath 3", "Swath 4")
    _BIST_LABELS = ("BIST 1", "BIST 2", "BIST 3", "BIST 4")

    def __init__(self, parent, imported_lines):
        super().__init__(parent)
        self.setWindowTitle("Assign Performance Survey Lines")
        self.setMinimumSize(900, 700)
        self.setModal(True)
        self.imported_lines = imported_lines
        self.assignments = {}

        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(10)

        instructions = QLabel(
            "Assign each imported segment to a performance swath (1–4). "
            "Optionally assign four BIST segments (BIST 1–4), or leave BIST as Ignore. "
            "If any BIST is used, all four BIST slots must be filled."
        )
        instructions.setStyleSheet("font-weight: bold;")
        instructions.setWordWrap(True)
        main_layout.addWidget(instructions)

        self.figure = Figure(figsize=(8, 6))
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self.figure)
        main_layout.addWidget(self.canvas)

        self._plot_lines()

        assignment_widget = QWidget()
        assignment_layout = QVBoxLayout(assignment_widget)
        assignment_layout.setContentsMargins(5, 5, 5, 5)
        assign_title = QLabel("Line assignments:")
        assign_title.setStyleSheet("font-weight: bold;")
        assignment_layout.addWidget(assign_title)

        line_colors = [
            "red",
            "blue",
            "green",
            "purple",
            "orange",
            "cyan",
            "magenta",
            "goldenrod",
            "sienna",
            "olive",
        ]
        num_lines = len(imported_lines)
        self.comboboxes = []

        combo_options = ["", "Ignore"] + list(self._SWATH_LABELS) + list(self._BIST_LABELS)

        for i in range(num_lines):
            line_widget = QWidget()
            line_layout = QHBoxLayout(line_widget)
            line_layout.setContentsMargins(0, 0, 0, 0)
            line_layout.addWidget(QLabel(f"Segment {i + 1}:"))
            combo = QComboBox()
            for opt in combo_options:
                combo.addItem(opt)
            # Defaults: first four → swath 1–4; next four → BIST 1–4 if present
            if i < 4:
                combo.setCurrentIndex(2 + i)  # Swath 1 = index 2
            elif 4 <= i < 8:
                combo.setCurrentIndex(2 + 4 + (i - 4))  # BIST 1 starts after Swath 4
            else:
                combo.setCurrentIndex(1)  # Ignore
            self.comboboxes.append(combo)
            line_layout.addWidget(combo)
            color_label = QLabel("●")
            color_label.setStyleSheet(f"color: {line_colors[i % len(line_colors)]}; font-size: 16pt;")
            color_label.setMinimumWidth(24)
            line_layout.addWidget(color_label)
            assignment_layout.addWidget(line_widget)

        main_layout.addWidget(assignment_widget)

        button_layout = QHBoxLayout()
        button_layout.addStretch()
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        button_layout.addWidget(cancel_btn)
        ok_btn = QPushButton("OK")
        ok_btn.clicked.connect(self._on_ok_clicked)
        button_layout.addWidget(ok_btn)
        main_layout.addLayout(button_layout)

    def _plot_lines(self):
        self.ax.clear()
        line_colors = [
            "red",
            "blue",
            "green",
            "purple",
            "orange",
            "cyan",
            "magenta",
            "goldenrod",
            "sienna",
            "olive",
        ]
        all_lats, all_lons = [], []
        for i, line in enumerate(self.imported_lines):
            if len(line) != 2:
                continue
            (lat1, lon1), (lat2, lon2) = line[0], line[1]
            c = line_colors[i % len(line_colors)]
            self.ax.plot([lon1, lon2], [lat1, lat2], color=c, linewidth=2, label=f"Segment {i + 1}", zorder=10)
            mid_lat = (lat1 + lat2) / 2
            mid_lon = (lon1 + lon2) / 2
            self.ax.text(
                mid_lon,
                mid_lat,
                str(i + 1),
                fontsize=10,
                ha="center",
                va="center",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.7),
                zorder=11,
            )
            all_lats.extend([lat1, lat2])
            all_lons.extend([lon1, lon2])
        if all_lats and all_lons:
            lat_range = max(all_lats) - min(all_lats)
            lon_range = max(all_lons) - min(all_lons)
            pad_lat = max(lat_range * 0.1, 0.01)
            pad_lon = max(lon_range * 0.1, 0.01)
            self.ax.set_xlim(min(all_lons) - pad_lon, max(all_lons) + pad_lon)
            self.ax.set_ylim(min(all_lats) - pad_lat, max(all_lats) + pad_lat)
        self.ax.set_xlabel("Longitude")
        self.ax.set_ylabel("Latitude")
        self.ax.grid(True, alpha=0.3)
        self.ax.legend(loc="upper right")
        self.figure.tight_layout()
        self.canvas.draw()

    def _role_from_text(self, text):
        if not text or text == "Ignore":
            return None
        if text in self._SWATH_LABELS:
            return "swath", self._SWATH_LABELS.index(text) + 1
        if text in self._BIST_LABELS:
            return "bist", self._BIST_LABELS.index(text) + 1
        return None

    def _on_ok_clicked(self):
        swath_slots = {1: None, 2: None, 3: None, 4: None}
        bist_slots = {1: None, 2: None, 3: None, 4: None}

        for i, combo in enumerate(self.comboboxes):
            text = combo.currentText()
            parsed = self._role_from_text(text)
            if parsed is None:
                continue
            kind, num = parsed
            if kind == "swath":
                if swath_slots[num] is not None:
                    QMessageBox.warning(self, "Assignment Error", f"Swath {num} is assigned more than once.")
                    return
                swath_slots[num] = i
            else:
                if bist_slots[num] is not None:
                    QMessageBox.warning(self, "Assignment Error", f"BIST {num} is assigned more than once.")
                    return
                bist_slots[num] = i

        for n in range(1, 5):
            if swath_slots[n] is None:
                QMessageBox.warning(self, "Assignment Error", f"Swath line {n} must be assigned exactly once.")
                return

        bist_assigned = sum(1 for v in bist_slots.values() if v is not None)
        if bist_assigned not in (0, 4):
            QMessageBox.warning(
                self,
                "Assignment Error",
                "Either assign all four BIST segments or leave all BIST as Ignore.",
            )
            return

        self.assignments.clear()
        for n in range(1, 5):
            self.assignments[swath_slots[n]] = f"swath_{n}"
        if bist_assigned == 4:
            for n in range(1, 5):
                self.assignments[bist_slots[n]] = f"bist_{n}"

        self.accept()

    def get_assignments(self):
        return self.assignments.copy()

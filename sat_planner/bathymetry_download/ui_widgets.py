"""Reusable PyQt widgets."""

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QMouseEvent
from PyQt6.QtWidgets import QLabel


class ClickableLabel(QLabel):
    """A QLabel that emits a clicked signal when clicked."""

    clicked = pyqtSignal()

    def mousePressEvent(self, event: QMouseEvent):
        if event.button() == Qt.MouseButton.LeftButton:
            self.clicked.emit()
        super().mousePressEvent(event)

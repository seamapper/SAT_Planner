# Copyright (c) 2026 Paul Johnson
# SPDX-License-Identifier: BSD-3-Clause
# Embedded in SAT Planner.

from PyQt6.QtCore import QThread, pyqtSignal
from PyQt6.QtGui import QPixmap, QImage
import requests

from ..config import GMRT_IMAGE_URL


class MapWorker(QThread):
    """Worker thread for loading map preview images from GMRT ImageServer."""
    map_loaded = pyqtSignal(QPixmap)
    map_error = pyqtSignal(str)

    def __init__(self, west, east, south, north, width=800, mask=True):
        super().__init__()
        self.west = west
        self.east = east
        self.south = south
        self.north = north
        self.width = width
        self.mask = mask

    def run(self):
        try:
            params = {
                "minlongitude": self.west,
                "maxlongitude": self.east,
                "minlatitude": self.south,
                "maxlatitude": self.north,
                "width": self.width,
                "mask": "1" if self.mask else "0"
            }
            response = requests.get(GMRT_IMAGE_URL, params=params, timeout=30)
            if response.status_code == 200:
                image = QImage()
                image.loadFromData(response.content)
                pixmap = QPixmap.fromImage(image)
                self.map_loaded.emit(pixmap)
            else:
                self.map_error.emit(f"Failed to load map: HTTP {response.status_code}")
        except Exception as e:
            self.map_error.emit(f"Map loading error: {str(e)}")

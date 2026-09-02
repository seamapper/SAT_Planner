"""Async loader for ArcGIS ImageServer metadata."""

import requests
from PyQt6.QtCore import QThread, pyqtSignal


class ServiceInfoLoader(QThread):
    """Thread for loading service information asynchronously."""

    loaded = pyqtSignal(dict)
    error = pyqtSignal(str)
    statusMessage = pyqtSignal(str)

    def __init__(self, base_url):
        super().__init__()
        self.base_url = base_url

    def run(self):
        try:
            url = f"{self.base_url}?f=json"
            self.statusMessage.emit(f"Service REST (metadata): {url}")
            response = requests.get(url, timeout=15)
            response.raise_for_status()
            data = response.json()

            extent = data.get("extent", {})
            extent_dict = {
                "xmin": extent.get("xmin", -8254538.5),
                "ymin": extent.get("ymin", 4898563.25),
                "xmax": extent.get("xmax", -7411670.5),
                "ymax": extent.get("ymax", 5636075.25),
            }

            raster_functions = ["None"]
            for rf_info in data.get("rasterFunctionInfos", []):
                name = rf_info.get("name", "")
                if name and name != "None":
                    raster_functions.append(name)

            result = {
                "extent": extent_dict,
                "raster_functions": raster_functions,
                "pixel_size_x": data.get("pixelSizeX", None),
                "pixel_size_y": data.get("pixelSizeY", None),
            }
            self.loaded.emit(result)

        except requests.exceptions.Timeout:
            self.error.emit("Connection timeout. Using default extent.")
        except requests.exceptions.RequestException as exc:
            self.error.emit(
                f"Network error connecting to REST endpoint: {exc}. Using default extent."
            )
        except Exception as exc:
            self.error.emit(f"Error loading service info: {exc}. Using default extent.")

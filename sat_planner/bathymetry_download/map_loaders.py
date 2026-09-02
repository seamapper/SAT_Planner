"""Background threads for loading map preview images from REST services."""

import json
from io import BytesIO
from urllib.parse import urlencode

import numpy as np
import requests
from PIL import Image
from PyQt6.QtCore import QThread, pyqtSignal
from PyQt6.QtGui import QImage, QPixmap


def format_rest_url(url, params):
    """Build a full REST request URL for activity log display."""
    return f"{url}?{urlencode(params)}"


class BasemapLoader(QThread):
    """Load World Imagery via MapServer export for an exact map extent."""

    tileLoaded = pyqtSignal(QPixmap)
    statusMessage = pyqtSignal(str)

    def __init__(self, bbox, size, bbox_sr=None):
        super().__init__()
        self.bbox = bbox
        self.size = size
        self.bbox_sr = bbox_sr
        self.basemap_url = (
            "https://services.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/export"
        )

    def run(self):
        try:
            xmin, ymin, xmax, ymax = self.bbox
            width, height = self.size
            params = {
                "bbox": f"{xmin},{ymin},{xmax},{ymax}",
                "size": f"{width},{height}",
                "format": "png",
                "f": "image",
                "transparent": "false",
            }
            if self.bbox_sr:
                params["bboxSR"] = self.bbox_sr
                params["imageSR"] = self.bbox_sr
            full_url = format_rest_url(self.basemap_url, params)
            self.statusMessage.emit(f"Map REST (World Imagery): {full_url}")
            response = requests.get(self.basemap_url, params=params, timeout=30)
            response.raise_for_status()

            img = Image.open(BytesIO(response.content))
            img_bytes = BytesIO()
            img.save(img_bytes, format="PNG")
            img_bytes.seek(0)
            pixmap = QPixmap()
            pixmap.loadFromData(img_bytes.getvalue(), "PNG")
            self.tileLoaded.emit(pixmap)
        except Exception as exc:
            print(f"Error loading basemap: {exc}")
            self.tileLoaded.emit(QPixmap())


class MapTileLoader(QThread):
    """Thread for loading map tiles asynchronously from ArcGIS ImageServer."""

    tileLoaded = pyqtSignal(QPixmap, float, float, float, float)
    statusMessage = pyqtSignal(str)

    def __init__(
        self,
        base_url,
        bbox,
        size,
        raster_function="Haxby Percent Clip DRA",
        bbox_sr=None,
        purpose="map display",
    ):
        super().__init__()
        self.base_url = base_url
        self.bbox = bbox
        self.size = size
        self.raster_function = raster_function
        self.bbox_sr = bbox_sr
        self.purpose = purpose

    def run(self):
        try:
            xmin, ymin, xmax, ymax = self.bbox
            width, height = self.size

            url = f"{self.base_url}/exportImage"
            params = {
                "bbox": f"{xmin},{ymin},{xmax},{ymax}",
                "size": f"{width},{height}",
                "format": "png",
                "f": "image",
            }
            if self.bbox_sr:
                params["bboxSR"] = self.bbox_sr
                params["imageSR"] = self.bbox_sr

            if self.raster_function and self.raster_function != "None":
                rendering_rule = {"rasterFunction": self.raster_function}
                params["renderingRule"] = json.dumps(rendering_rule)

            full_url = format_rest_url(url, params)
            self.statusMessage.emit(f"Map REST ({self.purpose}): {full_url}")

            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()

            img = Image.open(BytesIO(response.content))
            img_bytes = BytesIO()
            img.save(img_bytes, format="PNG")
            img_bytes.seek(0)

            pixmap = QPixmap()
            pixmap.loadFromData(img_bytes.getvalue(), "PNG")

            if pixmap.isNull():
                img_rgb = img.convert("RGB")
                img_array = np.array(img_rgb, dtype=np.uint8)
                height, width, _channel = img_array.shape
                if not img_array.flags["C_CONTIGUOUS"]:
                    img_array = np.ascontiguousarray(img_array)
                bytes_per_line = 3 * width
                q_image = QImage(
                    img_array.data, width, height, bytes_per_line, QImage.Format.Format_RGB888
                )
                pixmap = QPixmap.fromImage(q_image.copy())

            self.tileLoaded.emit(pixmap, xmin, ymin, xmax, ymax)

        except Exception as exc:
            print(f"Error loading tile: {exc}")
            import traceback

            traceback.print_exc()
            self.tileLoaded.emit(QPixmap(), *self.bbox)


class MapServerLoader(QThread):
    """Load a single image from an ArcGIS MapServer (e.g. GEBCO Haxby). Bbox in GCS (4326)."""

    tileLoaded = pyqtSignal(QPixmap, float, float, float, float)
    statusMessage = pyqtSignal(str)

    def __init__(self, map_server_url, bbox_4326, size, transparent=False, purpose="map display"):
        super().__init__()
        self.map_server_url = map_server_url.rstrip("/")
        self.bbox_4326 = bbox_4326
        self.size = size
        self.transparent = transparent
        self.purpose = purpose
        self.max_retries = 3
        self.retry_delay_seconds = 1.0

    def run(self):
        try:
            west, south, east, north = self.bbox_4326
            width, height = self.size
            max_side = 4096
            if width > max_side or height > max_side:
                scale = min(max_side / width, max_side / height)
                width = int(width * scale)
                height = int(height * scale)
            url = f"{self.map_server_url}/export"
            params = {
                "bbox": f"{west},{south},{east},{north}",
                "bboxSR": "4326",
                "size": f"{width},{height}",
                "format": "png",
                "f": "image",
                "transparent": "true" if self.transparent else "false",
            }
            self.statusMessage.emit(f"Map REST ({self.purpose}): {format_rest_url(url, params)}")
            response = None
            for attempt in range(1, self.max_retries + 1):
                try:
                    response = requests.get(url, params=params, timeout=60)
                    response.raise_for_status()
                    break
                except requests.exceptions.HTTPError as http_err:
                    status_code = getattr(http_err.response, "status_code", None)
                    is_retryable = status_code in (500, 502, 503, 504)
                    if is_retryable and attempt < self.max_retries:
                        self.statusMessage.emit(
                            f'<span style="color: orange;">Warning: Map server returned {status_code}. '
                            f"Retrying ({attempt + 1}/{self.max_retries})...</span>"
                        )
                        QThread.msleep(int(self.retry_delay_seconds * 1000 * attempt))
                        continue
                    raise
                except (requests.exceptions.Timeout, requests.exceptions.ConnectionError):
                    if attempt < self.max_retries:
                        self.statusMessage.emit(
                            f'<span style="color: orange;">Warning: Map server request timed out/connection failed. '
                            f"Retrying ({attempt + 1}/{self.max_retries})...</span>"
                        )
                        QThread.msleep(int(self.retry_delay_seconds * 1000 * attempt))
                        continue
                    raise

            if response is None:
                raise RuntimeError("Map server request failed without a response.")

            img = Image.open(BytesIO(response.content))
            img_bytes = BytesIO()
            img.save(img_bytes, format="PNG")
            img_bytes.seek(0)
            pixmap = QPixmap()
            pixmap.loadFromData(img_bytes.getvalue(), "PNG")
            if pixmap.isNull():
                img_rgb = img.convert("RGB")
                img_array = np.array(img_rgb, dtype=np.uint8)
                h, w = img_array.shape[:2]
                bytes_per_line = 3 * w
                q_image = QImage(img_array.data, w, h, bytes_per_line, QImage.Format.Format_RGB888)
                pixmap = QPixmap.fromImage(q_image.copy())
            self.tileLoaded.emit(pixmap, west, south, east, north)
        except Exception as exc:
            print(f"MapServerLoader error: {exc}")
            self.statusMessage.emit(
                f'<span style="color: orange;">Warning: Map load problem ({exc}). Please try again.</span>'
            )
            self.tileLoaded.emit(QPixmap(), *self.bbox_4326)

"""
Local MapLibre map embedded in the main window.

MapLibre GL JS is served from this package over 127.0.0.1 so the page loads
offline and its WebGL worker is allowed to start. Survey lines are copied
onto this map from the matplotlib axes after each redraw. World imagery is a
remote tile layer. NOAA ENC charts are the same Maritime Chart Service export
the desktop app already uses, served here as {z}/{x}/{y} tiles.
"""
import base64
import json
import math
import mimetypes
import os
import re
import sys
import threading
from collections import OrderedDict
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from urllib.request import Request, urlopen

from PyQt6.QtCore import QEventLoop, QFile, QObject, QTimer, QUrl, pyqtSignal, pyqtSlot
from PyQt6.QtWebChannel import QWebChannel
from PyQt6.QtWebEngineCore import QWebEnginePage, QWebEngineSettings
from PyQt6.QtWebEngineWidgets import QWebEngineView

from sat_planner.constants import PLANNING_PLACEHOLDER_TEXT
from sat_planner.utils_ui import get_about_program_info, media_path

mimetypes.add_type("application/javascript", ".js")
mimetypes.add_type("text/css", ".css")

_server_lock = threading.Lock()
_server = None

# Same ArcGIS World Imagery template basemap_mixin downloads for the hidden canvas.
IMAGERY_TILE_URL = "https://services.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}"
_NOAA_EXPORT_URL = (
    "https://gis.charttools.noaa.gov/arcgis/rest/services/MCS/NOAAChartDisplay/"
    "MapServer/exts/MaritimeChartService/MapServer/export"
)
_NOAA_TILE_PATH = re.compile(r"^/noaa/(\d+)/(\d+)/(\d+)\.png$")
_NOAA_CACHE = OrderedDict()
_NOAA_CACHE_LOCK = threading.Lock()
_NOAA_CACHE_LIMIT = 128
_MERCATOR_LIMIT = 20037508.342789244
# 1×1 transparent PNG, returned when a chart tile cannot be fetched.
_TRANSPARENT_PNG = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAAC0lEQVQI12NgAAIAAAUAAeImBZsAAAAASUVORK5CYII="
)


def _asset_dir():
    """Directory that holds index.html and the MapLibre files."""
    if getattr(sys, "frozen", False):
        base = getattr(sys, "_MEIPASS", os.path.dirname(sys.executable))
        return os.path.join(base, "sat_planner", "maplibre")
    return os.path.dirname(os.path.abspath(__file__))


def _qwebchannel_source():
    """Qt's qwebchannel.js, shipped as a Qt resource."""
    qfile = QFile(":/qtwebchannel/qwebchannel.js")
    if not qfile.open(QFile.OpenModeFlag.ReadOnly):
        return b""
    return bytes(qfile.readAll())


def _tile_bbox_3857(z, x, y):
    """Web Mercator bounds of one XYZ tile, as an ArcGIS export bbox string."""
    n = 2 ** z
    west = x / n * 360.0 - 180.0
    east = (x + 1) / n * 360.0 - 180.0
    north = math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * y / n))))
    south = math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * (y + 1) / n))))

    def lon_to_m(lon):
        return lon * _MERCATOR_LIMIT / 180.0

    def lat_to_m(lat):
        lat = max(-85.05112878, min(85.05112878, lat))
        rad = math.radians(lat)
        return math.log(math.tan(math.pi / 4.0 + rad / 2.0)) * _MERCATOR_LIMIT / math.pi

    return f"{lon_to_m(west)},{lat_to_m(south)},{lon_to_m(east)},{lat_to_m(north)}"


def _noaa_tile_png(z, x, y):
    """Fetch one NOAA ENC tile. Failures become a transparent tile."""
    key = (z, x, y)
    with _NOAA_CACHE_LOCK:
        cached = _NOAA_CACHE.get(key)
        if cached is not None:
            _NOAA_CACHE.move_to_end(key)
            return cached
    data = _TRANSPARENT_PNG
    try:
        bbox = _tile_bbox_3857(z, x, y)
        query = (
            f"bbox={bbox}&bboxSR=3857&imageSR=3857&size=256,256&format=png"
            "&transparent=true&f=image&layers=show:0,1,2,3,4,5,6,7"
        )
        req = Request(f"{_NOAA_EXPORT_URL}?{query}", headers={"User-Agent": "SAT_Planner"})
        with urlopen(req, timeout=12) as response:
            payload = response.read()
        if payload.startswith(b"\x89PNG"):
            data = payload
    except Exception:
        data = _TRANSPARENT_PNG
    with _NOAA_CACHE_LOCK:
        _NOAA_CACHE[key] = data
        _NOAA_CACHE.move_to_end(key)
        while len(_NOAA_CACHE) > _NOAA_CACHE_LIMIT:
            _NOAA_CACHE.popitem(last=False)
    return data


def _handler_class(directory):
    class _QuietHandler(SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=directory, **kwargs)

        def do_GET(self):
            path = self.path.split("?", 1)[0]
            if path == "/qwebchannel.js":
                data = _qwebchannel_source()
                self._send_bytes(data, "application/javascript", data != b"")
                return
            match = _NOAA_TILE_PATH.match(path)
            if match:
                z, x, y = (int(part) for part in match.groups())
                limit = 2 ** z
                if z > 18 or x >= limit or y >= limit:
                    self._send_bytes(_TRANSPARENT_PNG, "image/png", True)
                    return
                self._send_bytes(_noaa_tile_png(z, x, y), "image/png", True, cache=True)
                return
            super().do_GET()

        def _send_bytes(self, data, content_type, ok, cache=False):
            self.send_response(200 if ok else 404)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(data)))
            if cache:
                self.send_header("Cache-Control", "public, max-age=86400")
            self.end_headers()
            self.wfile.write(data)

        def log_message(self, format, *args):
            return

    return _QuietHandler


class _AssetServer(threading.Thread):
    """Serve the map page on localhost. One server is shared by the process."""

    def __init__(self, directory):
        super().__init__(name="sat-maplibre", daemon=True)
        self.httpd = ThreadingHTTPServer(("127.0.0.1", 0), _handler_class(directory))
        self.port = self.httpd.server_address[1]

    def run(self):
        self.httpd.serve_forever()


def _startup_logo_data_url():
    """CCOM mark used on the empty-map screen, as a data URL."""
    path = media_path("CCOM.png")
    if not os.path.isfile(path):
        return ""
    try:
        with open(path, "rb") as handle:
            encoded = base64.b64encode(handle.read()).decode("ascii")
    except OSError:
        return ""
    return "data:image/png;base64," + encoded


def _startup_overlay_info():
    """About text and logo shown until a plan or bathymetry grid is on the map."""
    try:
        import sat_planner.app_core as app_core

        about = get_about_program_info(source_path=app_core.__file__)
    except Exception:
        about = get_about_program_info()
    return {
        "title": about.get("title", ""),
        "frozen": about.get("compile_date", ""),
        "logo": _startup_logo_data_url(),
        "author": about.get("author", ""),
        "email": about.get("email", ""),
        "institution": about.get("institution", ""),
        "grant": about.get("grant", ""),
        "license": about.get("license", ""),
        "message": PLANNING_PLACEHOLDER_TEXT,
    }


def map_page_url():
    """Start the asset server if needed and return the map page URL."""
    global _server
    with _server_lock:
        if _server is None:
            directory = _asset_dir()
            if not os.path.isfile(os.path.join(directory, "index.html")):
                raise FileNotFoundError(f"Map page is missing: {directory}")
            _server = _AssetServer(directory)
            _server.start()
    return f"http://127.0.0.1:{_server.port}/index.html"


class _MapPage(QWebEnginePage):
    def javaScriptConsoleMessage(self, level, message, line_number, source_id):
        if level == QWebEnginePage.JavaScriptConsoleMessageLevel.InfoMessageLevel:
            return
        print(f"maplibre: {message} ({source_id}:{line_number})")


class _MapBridge(QObject):
    """Object the map page calls when the user clicks, moves, or finishes a pan."""

    def __init__(self, view):
        super().__init__(view)
        self._view = view

    @pyqtSlot(float, float, float)
    def mapClick(self, lon, lat, button):
        self._view.mapClicked.emit(float(lon), float(lat), int(button))

    @pyqtSlot(float, float, float)
    def mapPress(self, lon, lat, button):
        self._view.mapPressed.emit(float(lon), float(lat), int(button))

    @pyqtSlot(float, float, float)
    def mapRelease(self, lon, lat, button):
        self._view.mapReleased.emit(float(lon), float(lat), int(button))

    @pyqtSlot(float, float)
    def mapMove(self, lon, lat):
        self._view.mapMoved.emit(float(lon), float(lat))

    @pyqtSlot(float, float, float, float)
    def viewChanged(self, west, south, east, north):
        self._view.viewChanged.emit(float(west), float(south), float(east), float(north))


class MapLibreView(QWebEngineView):
    """MapLibre map with a Qt channel for pointer events and drawing commands."""

    mapClicked = pyqtSignal(float, float, int)
    mapMoved = pyqtSignal(float, float)
    mapPressed = pyqtSignal(float, float, int)
    mapReleased = pyqtSignal(float, float, int)
    viewChanged = pyqtSignal(float, float, float, float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._startup_applied = False
        self._pending_plan = None
        self._pending_bathymetry = None
        self._pending_cursor = None
        self._pending_profile_cursor = None
        self._pending_hover = None
        self._pending_colorbar = None
        self._pending_legend = None
        self._pending_web_layers = None
        self._pending_fit = None
        page = _MapPage(self)
        self._bridge = _MapBridge(self)
        self._channel = QWebChannel(page)
        self._channel.registerObject("bridge", self._bridge)
        page.setWebChannel(self._channel)
        self.setPage(page)
        settings = self.settings()
        settings.setAttribute(QWebEngineSettings.WebAttribute.JavascriptEnabled, True)
        settings.setAttribute(QWebEngineSettings.WebAttribute.LocalContentCanAccessRemoteUrls, True)
        settings.setAttribute(QWebEngineSettings.WebAttribute.WebGLEnabled, True)
        self.loadFinished.connect(self._on_load_finished)
        self.load(QUrl(map_page_url()))

    def _run(self, script):
        if self._startup_applied:
            self.page().runJavaScript(script)

    def _on_load_finished(self, ok):
        if not ok or self._startup_applied:
            return
        self._startup_applied = True
        self._run(
            "if (window.setStartup) window.setStartup({info});".format(
                info=json.dumps(_startup_overlay_info()),
            )
        )
        if self._pending_web_layers is not None:
            self.set_web_layers(*self._pending_web_layers)
        if self._pending_bathymetry is not None:
            self.show_bathymetry(*self._pending_bathymetry)
        if self._pending_plan is not None:
            collection, bounds = self._pending_plan
            self.show_plan(collection, bounds)
        if self._pending_fit is not None:
            self.fit_bounds(*self._pending_fit)
        if self._pending_cursor is not None:
            self.set_cursor(self._pending_cursor)
        if self._pending_profile_cursor is not None:
            self.set_profile_cursor(*self._pending_profile_cursor)
        if self._pending_hover is not None:
            self.set_hover(self._pending_hover)
        if self._pending_colorbar is not None:
            self.set_colorbar(*self._pending_colorbar)
        if self._pending_legend is not None:
            self.set_legend(*self._pending_legend)

    def show_plan(self, collection, bounds=None):
        """Draw a GeoJSON feature collection and optionally fit the view to bounds."""
        self._pending_plan = (collection, bounds)
        self._run(
            "if (window.setPlan) window.setPlan({collection}, {bounds});".format(
                collection=json.dumps(collection),
                bounds=json.dumps(bounds),
            )
        )

    def set_web_layers(self, imagery_on, noaa_on, noaa_opacity):
        """Show or hide World Imagery tiles and the NOAA ENC tile overlay."""
        opacity = max(0.0, min(1.0, float(noaa_opacity)))
        self._pending_web_layers = (bool(imagery_on), bool(noaa_on), opacity)
        self._run(
            "if (window.setWebLayers) window.setWebLayers({imagery}, {noaa}, {opacity}, {tiles});".format(
                imagery="true" if imagery_on else "false",
                noaa="true" if noaa_on else "false",
                opacity=json.dumps(opacity),
                tiles=json.dumps([IMAGERY_TILE_URL]),
            )
        )

    def fit_bounds(self, west, south, east, north):
        """Move the map to a longitude/latitude window without changing the drawing."""
        self._pending_fit = (west, south, east, north)
        self._run(
            "if (window.fitView) window.fitView({west}, {south}, {east}, {north});".format(
                west=json.dumps(west),
                south=json.dumps(south),
                east=json.dumps(east),
                north=json.dumps(north),
            )
        )

    def show_bathymetry(self, data_url, west, south, east, north):
        """Replace the bathymetry image. An empty url removes it."""
        self._pending_bathymetry = (data_url, west, south, east, north)
        self._run(
            "if (window.setBathymetry) window.setBathymetry({url}, {west}, {south}, {east}, {north});".format(
                url=json.dumps(data_url or ""),
                west=json.dumps(west),
                south=json.dumps(south),
                east=json.dumps(east),
                north=json.dumps(north),
            )
        )

    def set_cursor(self, css_cursor):
        """Set the map canvas cursor. An empty string restores the default."""
        self._pending_cursor = css_cursor or ""
        self._run("if (window.setCursor) window.setCursor({cursor});".format(cursor=json.dumps(self._pending_cursor)))

    def set_profile_cursor(self, lon, lat):
        """Show or clear the marker for the profile pointer's distance along the line."""
        if lon is None or lat is None:
            self._pending_profile_cursor = None
            self._run("if (window.setProfileCursor) window.setProfileCursor(NaN, NaN);")
            return
        lon = float(lon)
        lat = float(lat)
        self._pending_profile_cursor = (lon, lat)
        self._run(
            "if (window.setProfileCursor) window.setProfileCursor({lon}, {lat});".format(
                lon=json.dumps(lon),
                lat=json.dumps(lat),
            )
        )

    def set_hover(self, text):
        """Show or clear the hover readout on the map."""
        self._pending_hover = text or ""
        self._run("if (window.setHover) window.setHover({text});".format(text=json.dumps(self._pending_hover)))

    def set_legend(self, data_url, visible):
        """Show or hide the survey-line legend on the map."""
        self._pending_legend = (data_url or "", bool(visible))
        self._run(
            "if (window.setLegend) window.setLegend({url}, {visible});".format(
                url=json.dumps(data_url or ""),
                visible="true" if visible else "false",
            )
        )

    def set_colorbar(self, data_url, visible):
        """Show or hide the depth/slope colorbar image on the map."""
        self._pending_colorbar = (data_url or "", bool(visible))
        self._run(
            "if (window.setColorbar) window.setColorbar({url}, {visible});".format(
                url=json.dumps(data_url or ""),
                visible="true" if visible else "false",
            )
        )

    def capture_png(self, timeout_ms=4000):
        """Return PNG bytes of the current on-screen map, or None if the page is not ready."""
        if not self._startup_applied or getattr(self, "_capturing_png", False):
            return None
        self._capturing_png = True
        box = {"payload": None}
        loop = QEventLoop()
        timer = QTimer()
        timer.setSingleShot(True)
        timer.timeout.connect(loop.quit)

        def finished(payload):
            box["payload"] = payload
            if loop.isRunning():
                loop.quit()

        script = (
            "(function(){ try {"
            " if (window.captureMapPngNow) return window.captureMapPngNow();"
            " return '';"
            " } catch (e) { return ''; } })()"
        )
        try:
            self.page().runJavaScript(script, finished)
            if box["payload"] is None:
                timer.start(int(timeout_ms))
                loop.exec()
        finally:
            timer.stop()
            self._capturing_png = False
        payload = box["payload"]
        if isinstance(payload, str) and "," in payload and payload.startswith("data:image/png"):
            try:
                self._last_capture_includes_overlays = False
                return base64.b64decode(payload.split(",", 1)[1])
            except Exception:
                pass
        self._last_capture_includes_overlays = True
        return self._grab_widget_png()

    def _grab_widget_png(self):
        """Screenshot the map widget when the page canvas cannot be read."""
        from PyQt6.QtCore import QBuffer, QIODevice

        pixmap = self.grab()
        if pixmap.isNull() or pixmap.width() < 2 or pixmap.height() < 2:
            return None
        image = pixmap.toImage()
        if image.isNull():
            return None
        corners = (
            image.pixelColor(0, 0),
            image.pixelColor(image.width() - 1, 0),
            image.pixelColor(0, image.height() - 1),
            image.pixelColor(image.width() // 2, image.height() // 2),
        )
        if len({(color.red(), color.green(), color.blue()) for color in corners}) < 2:
            return None
        buffer = QBuffer()
        buffer.open(QIODevice.OpenModeFlag.WriteOnly)
        if not image.save(buffer, "PNG"):
            return None
        return bytes(buffer.data())

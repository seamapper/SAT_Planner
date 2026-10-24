"""MapLibre preview of imported survey lines for the assignment dialogs."""

import math

from PyQt6.QtWidgets import QHBoxLayout, QPushButton, QSizePolicy, QVBoxLayout, QWidget

from sat_planner.maplibre.map_view import MapLibreView


def _mercator_y(lat):
    """Web Mercator y for a latitude, in the same units as longitude degrees."""
    limited = max(-85.0, min(85.0, float(lat)))
    return math.degrees(math.log(math.tan(math.pi / 4.0 + math.radians(limited) / 2.0)))


def _label_rotation_deg(lon1, lat1, lon2, lat2):
    """Screen angle of the segment, kept within -90 to 90 so the text stays upright."""
    dx = float(lon2) - float(lon1)
    dy = -(_mercator_y(lat2) - _mercator_y(lat1))
    if abs(dx) < 1e-12 and abs(dy) < 1e-12:
        return 0.0
    angle = math.degrees(math.atan2(dy, dx))
    while angle > 90.0:
        angle -= 180.0
    while angle < -90.0:
        angle += 180.0
    return round(angle, 2)


def _line_collection(lines, colors, label_for_index):
    """GeoJSON for the imported segments, plus a view framed on those lines."""
    features = []
    lats = []
    lons = []
    for index, line in enumerate(lines or []):
        if not line or len(line) != 2:
            continue
        (lat1, lon1), (lat2, lon2) = line
        try:
            lat1, lon1, lat2, lon2 = float(lat1), float(lon1), float(lat2), float(lon2)
        except (TypeError, ValueError):
            continue
        color = colors[index % len(colors)]
        features.append({
            "type": "Feature",
            "properties": {"color": color, "width": 3, "dash": False},
            "geometry": {"type": "LineString", "coordinates": [[lon1, lat1], [lon2, lat2]]},
        })
        features.append({
            "type": "Feature",
            "properties": {
                "label": label_for_index(index),
                "marker": False,
                "labelRotation": _label_rotation_deg(lon1, lat1, lon2, lat2),
            },
            "geometry": {"type": "Point", "coordinates": [(lon1 + lon2) / 2.0, (lat1 + lat2) / 2.0]},
        })
        lats.extend((lat1, lat2))
        lons.extend((lon1, lon2))
    bounds = None
    if lats and lons:
        lat_span = max(lats) - min(lats)
        lon_span = max(lons) - min(lons)
        pad_lat = max(lat_span * 0.12, 0.01)
        pad_lon = max(lon_span * 0.12, 0.01)
        bounds = [min(lons) - pad_lon, min(lats) - pad_lat, max(lons) + pad_lon, max(lats) + pad_lat]
    return {"type": "FeatureCollection", "features": features}, bounds


class AssignmentMapPreview(QWidget):
    """Imported lines on a MapLibre map, with a button for the imagery basemap."""

    def __init__(self, lines, colors, label_for_index, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)

        bar = QHBoxLayout()
        self.imagery_btn = QPushButton("Imagery")
        self.imagery_btn.setCheckable(True)
        self.imagery_btn.setChecked(False)
        self.imagery_btn.setToolTip("Show or hide the satellite imagery basemap")
        self.imagery_btn.toggled.connect(self._on_imagery_toggled)
        bar.addWidget(self.imagery_btn)
        bar.addStretch()
        layout.addLayout(bar)

        self.map_view = MapLibreView(self)
        self.map_view.setMinimumHeight(320)
        self.map_view.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(self.map_view, 1)

        collection, bounds = _line_collection(lines, colors, label_for_index)
        self.map_view.show_plan(collection, bounds)
        self.map_view.set_web_layers(False, False, 0)

    def _on_imagery_toggled(self, checked):
        self.imagery_btn.setText("Imagery On" if checked else "Imagery")
        self.map_view.set_web_layers(bool(checked), False, 0)

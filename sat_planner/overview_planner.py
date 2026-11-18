"""
Overview of several finished test plans: order, transit, and schedule.

Times come from each plan's exported *_info.txt. Lines and the arrival and
departure points come from the GeoJSON saved beside that file.
"""
import json
import os
import re
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timedelta


PLAN_COLORS = (
    "#1f4e79",
    "#c0392b",
    "#1e8449",
    "#b9770e",
    "#6c3483",
    "#0e6655",
    "#1a5276",
    "#922b21",
)

# First matching line wins. Calibration publishes both "Total Survey Time"
# and "Total Time"; the second includes turns and the steam inside the site.
# Older line plans published the whole run as "Estimated Time".
# New line plans use "Total Time", the same label as the other surveys.
_DURATION_PATTERNS = (
    re.compile(r"Grand total \(speed \+ turns\):\s*([0-9]+(?:\.[0-9]+)?)\s*min", re.IGNORECASE),
    re.compile(r"Grand total:\s*([0-9]+(?:\.[0-9]+)?)\s*min", re.IGNORECASE),
    re.compile(r"(?m)^[ \t]*Total Time:\s*([0-9]+(?:\.[0-9]+)?)\s*min", re.IGNORECASE),
    re.compile(r"(?m)^[ \t]*Total Survey Time:\s*([0-9]+(?:\.[0-9]+)?)\s*min", re.IGNORECASE),
    re.compile(r"(?m)^[ \t]*Estimated Time:\s*([0-9]+(?:\.[0-9]+)?)\s*min", re.IGNORECASE),
)


@dataclass
class Stop:
    """One row in the order of operations."""

    kind: str
    name: str
    uid: str = field(default_factory=lambda: uuid.uuid4().hex[:12])
    duration_min: float = 0.0
    extra_min: float = 0.0
    # Gap before this stop, after the previous activity. The first stop ignores it.
    lead_min: float = 0.0
    # Added to the plan time and extra time when the end clock is edited.
    adjust_min: float = 0.0
    info_path: str = ""
    geojson_path: str = ""
    test_type: str = ""
    lines: list = field(default_factory=list)
    flipped: bool = False
    lat: float | None = None
    lon: float | None = None
    note: str = ""
    arrive_lat: float | None = None
    arrive_lon: float | None = None
    depart_lat: float | None = None
    depart_lon: float | None = None

    def endpoints(self):
        """(arrival lat/lon, departure lat/lon). A point uses the same coordinate for both.

        A survey uses the line ends from its GeoJSON. Overview-only overrides replace
        either end without changing that file.
        """
        if self.kind == "point":
            if self.lat is None or self.lon is None:
                return None, None
            point = (float(self.lat), float(self.lon))
            return point, point
        arrive, depart = self._file_endpoints()
        if self.arrive_lat is not None and self.arrive_lon is not None:
            arrive = (float(self.arrive_lat), float(self.arrive_lon))
        if self.depart_lat is not None and self.depart_lon is not None:
            depart = (float(self.depart_lat), float(self.depart_lon))
        return arrive, depart

    def _file_endpoints(self):
        if not self.lines:
            return None, None
        first = self.lines[0]
        last = self.lines[-1]
        if len(first) < 1 or len(last) < 1:
            return None, None
        arrive = (float(first[0][1]), float(first[0][0]))
        depart = (float(last[-1][1]), float(last[-1][0]))
        if self.flipped:
            return depart, arrive
        return arrive, depart

    def occupied_minutes(self):
        """Plan time, extra time, and any change made by editing the end clock."""
        total = float(self.duration_min) + float(self.extra_min) + float(self.adjust_min)
        return max(0.0, total)


@dataclass
class ScheduleRow:
    kind: str
    name: str
    duration_min: float
    distance_m: float = 0.0
    speed_kts: float | None = None
    lead_min: float = 0.0
    start: datetime | None = None
    end: datetime | None = None


def duration_minutes_from_info(text):
    """Minutes for the whole test, including that test's own transits and turns."""
    if not text:
        return None
    for pattern in _DURATION_PATTERNS:
        match = pattern.search(text)
        if match:
            return float(match.group(1))
    return None


def test_type_from_info(text):
    head = (text or "")[:800].upper()
    if "ADCP" in head:
        return "ADCP"
    if "PERFORMANCE" in head:
        return "Performance"
    if "BACKSCATTER" in head:
        return "Backscatter"
    if "CALIBRATION" in head:
        return "Calibration"
    if "ACCURACY" in head or "REFERENCE SURVEY" in head:
        return "Accuracy"
    if "LINE" in head:
        return "Line"
    return "Test"


def geojson_path_for_info(info_path):
    """SiteA_info.txt is exported next to SiteA.geojson."""
    folder, name = os.path.split(os.path.abspath(info_path))
    stem = name
    if stem.lower().endswith("_info.txt"):
        stem = stem[: -len("_info.txt")]
    elif stem.lower().endswith(".txt"):
        stem = stem[:-4]
    return os.path.join(folder, stem + ".geojson")


def lines_from_geojson(path):
    """LineStrings as lists of [lon, lat]."""
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    lines = []
    for feature in data.get("features") or []:
        geometry = feature.get("geometry") or {}
        kind = geometry.get("type")
        coords = geometry.get("coordinates") or []
        if kind == "LineString" and len(coords) >= 2:
            lines.append(coords)
        elif kind == "MultiLineString":
            for part in coords:
                if len(part) >= 2:
                    lines.append(part)
    return lines


def load_plan(info_path):
    """Read one exported test. Returns (Stop, warning). warning is empty when the files are complete."""
    info_path = os.path.abspath(info_path)
    name = os.path.splitext(os.path.basename(info_path))[0]
    if name.lower().endswith("_info"):
        name = name[: -len("_info")]
    stop = Stop(kind="plan", name=name, info_path=info_path)
    if not os.path.isfile(info_path):
        stop.note = "Info file not found"
        return stop, f"{name}: info file not found"
    try:
        with open(info_path, encoding="utf-8") as handle:
            text = handle.read()
    except OSError as exc:
        stop.note = "Could not read info file"
        return stop, f"{name}: {exc}"
    stop.test_type = test_type_from_info(text)
    minutes = duration_minutes_from_info(text)
    warnings = []
    if minutes is None:
        stop.note = "No total time in the info file"
        warnings.append(f"{name}: no total time in the info file")
    else:
        stop.duration_min = minutes
    geojson_path = geojson_path_for_info(info_path)
    stop.geojson_path = geojson_path
    if not os.path.isfile(geojson_path):
        stop.note = (stop.note + "; " if stop.note else "") + "GeoJSON not found"
        warnings.append(f"{name}: no GeoJSON next to the info file")
    else:
        try:
            stop.lines = lines_from_geojson(geojson_path)
        except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
            stop.note = (stop.note + "; " if stop.note else "") + "Could not read GeoJSON"
            warnings.append(f"{name}: {exc}")
        else:
            if not stop.lines:
                stop.note = (stop.note + "; " if stop.note else "") + "GeoJSON has no lines"
                warnings.append(f"{name}: GeoJSON has no lines")
    return stop, "\n".join(warnings)


def line_ends(stop):
    """(lat, lon) at the first and last vertex of each survey line."""
    points = []
    for line in stop.lines or []:
        if not line:
            continue
        for coord in (line[0], line[-1]):
            if len(coord) >= 2:
                points.append((float(coord[1]), float(coord[0])))
    return points


def nearest_line_end(stop, lat, lon):
    """Closest line end to a click. Returns ((lat, lon) or None, meters or None)."""
    best = None
    best_m = None
    for plat, plon in line_ends(stop):
        meters = geodesic_meters(lat, lon, plat, plon)
        if best_m is None or meters < best_m:
            best = (plat, plon)
            best_m = meters
    return best, best_m


def end_snap_limit_m(stop):
    """How far a click may be and still count as choosing a line end."""
    points = line_ends(stop)
    farthest = 0.0
    for index, origin in enumerate(points):
        for other in points[index + 1:]:
            farthest = max(farthest, geodesic_meters(origin[0], origin[1], other[0], other[1]))
    return max(2000.0, farthest * 0.6)


def swap_endpoints(stop):
    """Exchange the overview arrival and departure. The survey file is untouched."""
    if stop.kind != "plan":
        return False
    arrive, depart = stop.endpoints()
    if arrive is None or depart is None:
        return False
    stop.flipped = not stop.flipped
    stop.arrive_lat, stop.depart_lat = stop.depart_lat, stop.arrive_lat
    stop.arrive_lon, stop.depart_lon = stop.depart_lon, stop.arrive_lon
    return True


def clear_end_overrides(stop):
    """Go back to the start and end stored in the survey file."""
    stop.arrive_lat = None
    stop.arrive_lon = None
    stop.depart_lat = None
    stop.depart_lon = None
    stop.flipped = False


def geodesic_meters(lat1, lon1, lat2, lon2):
    try:
        import pyproj

        _az12, _az21, distance = pyproj.Geod(ellps="WGS84").inv(lon1, lat1, lon2, lat2)
        return float(distance)
    except Exception:
        from math import asin, cos, radians, sin, sqrt

        radius = 6371000.0
        dlat = radians(lat2 - lat1)
        dlon = radians(lon2 - lon1)
        lat1 = radians(lat1)
        lat2 = radians(lat2)
        a = sin(dlat / 2.0) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2.0) ** 2
        return 2.0 * radius * asin(sqrt(a))


def transit_between(depart, arrive, speed_kts):
    """(distance meters, minutes) at the overview transit speed."""
    if depart is None or arrive is None:
        return 0.0, 0.0
    distance = geodesic_meters(depart[0], depart[1], arrive[0], arrive[1])
    speed = float(speed_kts)
    if speed <= 0 or distance <= 0:
        return distance, 0.0
    meters_per_minute = speed * 1852.0 / 60.0
    return distance, distance / meters_per_minute


def activity_rows(stops):
    """Table rows: each stop, then the transit to the next stop."""
    rows = []
    for index in range(len(stops)):
        rows.append(("stop", index))
        if index + 1 < len(stops):
            rows.append(("transit", index))
    return rows


def overview_export_lines(stops):
    """Survey lines and transits as (kind, name, [(lat, lon), ...])."""
    lines = []
    for index, stop in enumerate(stops):
        if stop.kind == "plan":
            usable = [line for line in stop.lines if len(line) >= 2]
            for number, line in enumerate(usable, start=1):
                points = [(float(coord[1]), float(coord[0])) for coord in line if len(coord) >= 2]
                if len(points) < 2:
                    continue
                label = stop.name if len(usable) == 1 else f"{stop.name} {number}"
                lines.append(("Plan", f"{index + 1}. {label}", points))
        if index + 1 >= len(stops):
            continue
        depart = stop.endpoints()[1]
        arrive = stops[index + 1].endpoints()[0]
        if depart is None or arrive is None:
            continue
        lines.append((
            "Transit",
            f"Transit {index + 1} to {index + 2}",
            [(float(depart[0]), float(depart[1])), (float(arrive[0]), float(arrive[1]))],
        ))
    return lines


def overview_export_points(stops):
    """Arrival of each stop as (name, lat, lon)."""
    points = []
    for index, stop in enumerate(stops):
        arrive, _depart = stop.endpoints()
        if arrive is None:
            continue
        points.append((f"{index + 1}. {stop.name}", float(arrive[0]), float(arrive[1])))
    return points


def overview_kml_text(stops, document_name):
    """KML document for the survey lines, transits, and stop positions."""
    lines = overview_export_lines(stops)
    points = overview_export_points(stops)
    styles = [
        '    <Style id="transit">',
        "      <LineStyle>",
        "        <color>ff4a4a4a</color>",
        "        <width>2</width>",
        "      </LineStyle>",
        "    </Style>",
    ]
    for index, color in enumerate(PLAN_COLORS):
        styles.extend([
            f'    <Style id="plan{index}">',
            "      <LineStyle>",
            f"        <color>{_kml_color(color)}</color>",
            "        <width>3</width>",
            "      </LineStyle>",
            "    </Style>",
            f'    <Style id="stop{index}">',
            "      <IconStyle>",
            f"        <color>{_kml_color(color)}</color>",
            "        <scale>0.8</scale>",
            "      </IconStyle>",
            "    </Style>",
        ])
    folders = []
    plan_marks = []
    transit_marks = []
    for kind, label, path in lines:
        if kind == "Transit":
            transit_marks.append(_kml_line(label, path, "#transit"))
        else:
            plan_marks.append(_kml_line(label, path, f"#plan{_kml_style_index(label)}"))
    if plan_marks:
        folders.append(_kml_folder("Survey lines", plan_marks))
    if transit_marks:
        folders.append(_kml_folder("Transits", transit_marks))
    if points:
        marks = [
            _kml_point(label, lat, lon, f"#stop{_kml_style_index(label)}")
            for label, lat, lon in points
        ]
        folders.append(_kml_folder("Stops", marks))
    body = "\n".join(styles + folders)
    title = _kml_escape(document_name or "Overview")
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<kml xmlns="http://www.opengis.net/kml/2.2">\n'
        "  <Document>\n"
        f"    <name>{title}</name>\n"
        f"{body}\n"
        "  </Document>\n"
        "</kml>\n"
    )


def _kml_style_index(label):
    head = (label or "").split(".", 1)[0].strip()
    try:
        index = int(head) - 1
    except ValueError:
        index = 0
    if index < 0:
        index = 0
    return index % len(PLAN_COLORS)


def _kml_color(hex_color):
    text = (hex_color or "#000000").lstrip("#")
    if len(text) != 6:
        text = "000000"
    return "ff" + text[4:6] + text[2:4] + text[0:2]


def _kml_escape(text):
    return (
        (text or "")
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )


def _kml_coordinates(points):
    return " ".join(f"{lon:.6f},{lat:.6f}" for lat, lon in points)


def _kml_line(name, points, style_url):
    return (
        "      <Placemark>\n"
        f"        <name>{_kml_escape(name)}</name>\n"
        f"        <styleUrl>{style_url}</styleUrl>\n"
        "        <LineString>\n"
        "          <tessellate>1</tessellate>\n"
        f"          <coordinates>{_kml_coordinates(points)}</coordinates>\n"
        "        </LineString>\n"
        "      </Placemark>"
    )


def _kml_point(name, lat, lon, style_url):
    return (
        "      <Placemark>\n"
        f"        <name>{_kml_escape(name)}</name>\n"
        f"        <styleUrl>{style_url}</styleUrl>\n"
        "        <Point>\n"
        f"          <coordinates>{lon:.6f},{lat:.6f}</coordinates>\n"
        "        </Point>\n"
        "      </Placemark>"
    )


def _kml_folder(name, placemarks):
    body = "\n".join(placemarks)
    return (
        "    <Folder>\n"
        f"      <name>{_kml_escape(name)}</name>\n"
        f"{body}\n"
        "    </Folder>"
    )


def transit_leg(left, right, default_speed, overrides=None):
    """Distance and time from the end of one stop to the start of the next.

    A speed stored for this pair of stops replaces the default. The survey files
    are not involved.
    """
    overrides = overrides or {}
    key = (left.uid, right.uid)
    custom_speed = overrides.get(key)
    custom = custom_speed is not None
    speed = float(custom_speed) if custom else float(default_speed)
    _arrive, depart = left.endpoints()
    arrive, _depart = right.endpoints()
    distance, minutes = transit_between(depart, arrive, speed)
    return {
        "distance_m": distance,
        "minutes": minutes,
        "speed_kts": speed,
        "custom": custom,
        "key": key,
    }


def build_schedule(
    stops,
    speed_kts,
    start=None,
    transit_speeds=None,
    transit_leads=None,
    transit_durations=None,
):
    """One row per stop and one row per transit.

    Clock times follow `start` when it is set. A lead is a gap before that
    activity. A stored transit duration replaces the distance-and-speed time.
    """
    leads = transit_leads or {}
    durations = transit_durations or {}
    rows = []
    cursor = start
    for index, stop in enumerate(stops):
        lead = float(stop.lead_min or 0.0) if index else 0.0
        if cursor is not None and lead:
            cursor = cursor + timedelta(minutes=lead)
        minutes = stop.occupied_minutes()
        row = ScheduleRow(
            kind=stop.kind,
            name=stop.name,
            duration_min=minutes,
            lead_min=lead,
            start=cursor,
        )
        if cursor is not None:
            row.end = cursor + timedelta(minutes=minutes)
            cursor = row.end
        rows.append(row)
        if index + 1 >= len(stops):
            continue
        leg = transit_leg(stop, stops[index + 1], speed_kts, transit_speeds)
        key = leg["key"]
        lead = float(leads.get(key, 0.0) or 0.0)
        if cursor is not None and lead:
            cursor = cursor + timedelta(minutes=lead)
        minutes = leg["minutes"]
        if key in durations:
            minutes = max(0.0, float(durations[key]))
        transit = ScheduleRow(
            kind="transit",
            name=f"Transit to {stops[index + 1].name}",
            duration_min=minutes,
            distance_m=leg["distance_m"],
            speed_kts=leg["speed_kts"],
            lead_min=lead,
            start=cursor,
        )
        if cursor is not None:
            transit.end = cursor + timedelta(minutes=minutes)
            cursor = transit.end
        rows.append(transit)
    return rows


def format_schedule(rows, start=None):
    if not rows:
        return "Load test plans, then set the order of operations."
    lines = []
    total = 0.0
    current_day = None
    for row in rows:
        total += row.lead_min + row.duration_min
        if start is not None and row.start is not None and row.end is not None:
            day = row.start.date()
            if day != current_day:
                if lines:
                    lines.append("")
                lines.append(row.start.strftime("%A %d %b %Y"))
                current_day = day
            if row.end.date() != row.start.date():
                when = f"{row.start:%Y-%m-%d %H:%M}–{row.end:%Y-%m-%d %H:%M}  "
            else:
                when = f"{row.start:%H:%M}–{row.end:%H:%M}  "
        else:
            when = ""
        if row.kind == "transit":
            speed = f", {row.speed_kts:.1f} kn" if row.speed_kts else ""
            detail = f"{row.distance_m / 1852.0:.1f} nm{speed}, {_format_minutes(row.duration_min)}"
        else:
            detail = _format_minutes(row.duration_min)
        lines.append(f"{when}{row.name}  {detail}")
    lines.append("")
    lines.append(f"Total {_format_minutes(total)}")
    return "\n".join(lines)


def _format_minutes(minutes):
    minutes = float(minutes)
    if minutes >= 60.0:
        return f"{minutes / 60.0:.2f} hr"
    return f"{minutes:.1f} min"


def legend_items(stops):
    """(label, color, dashed) for the map legend, in stop order."""
    items = []
    for index, stop in enumerate(stops):
        color = PLAN_COLORS[index % len(PLAN_COLORS)]
        name = (stop.name or "").strip() or "Stop"
        items.append((f"{index + 1}. {name}", color, False))
    for index in range(len(stops) - 1):
        depart = stops[index].endpoints()[1]
        arrive = stops[index + 1].endpoints()[0]
        if depart is None or arrive is None:
            continue
        items.append(("Transit", "#4a4a4a", True))
        break
    return items


def activity_time_features(stops, speed_kts, transit_speeds=None, transit_durations=None):
    """A centered time label for each survey, point, and transit.

    A stop's time is its plan time plus extra time, including a clock adjustment.
    A transit uses its duration, or the duration stored for that leg.
    """
    durations = transit_durations or {}
    features = []
    for index, stop in enumerate(stops):
        center = _stop_label_center(stop)
        if center is not None:
            lon, lat, place = center
            features.append(_time_label_feature(lon, lat, stop.occupied_minutes(), place))
        if index + 1 >= len(stops):
            continue
        nxt = stops[index + 1]
        _arrive, depart = stop.endpoints()
        arrive, _depart = nxt.endpoints()
        if depart is None or arrive is None:
            continue
        leg = transit_leg(stop, nxt, speed_kts, transit_speeds)
        minutes = leg["minutes"]
        if leg["key"] in durations:
            minutes = max(0.0, float(durations[leg["key"]]))
        lon = (float(depart[1]) + float(arrive[1])) / 2.0
        lat = (float(depart[0]) + float(arrive[0])) / 2.0
        features.append(_time_label_feature(lon, lat, minutes, "center"))
    return features


def _stop_label_center(stop):
    """(lon, lat, place) at the middle of a stop, or None when it has no position."""
    if stop.kind == "point":
        if stop.lat is None or stop.lon is None:
            return None
        return float(stop.lon), float(stop.lat), "below"
    coords = []
    for line in stop.lines or []:
        for coord in line:
            if len(coord) >= 2:
                try:
                    lon = float(coord[0])
                    lat = float(coord[1])
                except (TypeError, ValueError):
                    continue
                coords.append((lon, lat))
    if coords:
        lons = [item[0] for item in coords]
        lats = [item[1] for item in coords]
        return (min(lons) + max(lons)) / 2.0, (min(lats) + max(lats)) / 2.0, "center"
    arrive, _depart = stop.endpoints()
    if arrive is None:
        return None
    return float(arrive[1]), float(arrive[0]), "below"


def _time_label_feature(lon, lat, minutes, place):
    return {
        "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [lon, lat]},
        "properties": {
            "label": _format_minutes(minutes),
            "labelPlace": place,
        },
    }


def map_features(stops, selected_index=None):
    """GeoJSON for the imagery map: plan lines, stop numbers, and dashed transits."""
    features = []
    for index, stop in enumerate(stops):
        color = PLAN_COLORS[index % len(PLAN_COLORS)]
        width = 5 if index == selected_index else 2.5
        for line in stop.lines:
            features.append({
                "type": "Feature",
                "geometry": {"type": "LineString", "coordinates": line},
                "properties": {"color": color, "width": width, "dash": False, "activity": f"stop:{index}"},
            })
        arrive, depart = stop.endpoints()
        if arrive is not None:
            features.append({
                "type": "Feature",
                "geometry": {"type": "Point", "coordinates": [arrive[1], arrive[0]]},
                "properties": {"color": color, "marker": True, "label": str(index + 1), "activity": f"stop:{index}"},
            })
        if index == selected_index and stop.kind == "plan" and depart is not None:
            features.append({
                "type": "Feature",
                "geometry": {"type": "Point", "coordinates": [depart[1], depart[0]]},
                "properties": {"color": color, "marker": True, "pick": True, "label": "End", "activity": f"stop:{index}"},
            })
        if index + 1 >= len(stops):
            continue
        next_arrive, _next_depart = stops[index + 1].endpoints()
        if depart is None or next_arrive is None:
            continue
        features.append({
            "type": "Feature",
            "geometry": {
                "type": "LineString",
                "coordinates": [
                    [depart[1], depart[0]],
                    [next_arrive[1], next_arrive[0]],
                ],
            },
            "properties": {"color": "#4a4a4a", "width": 2, "dash": True, "activity": f"transit:{index}"},
        })
    return {"type": "FeatureCollection", "features": features}


def feature_bounds(collection):
    """[west, south, east, north] with a small pad, or None."""
    lons = []
    lats = []
    for feature in collection.get("features") or []:
        geometry = feature.get("geometry") or {}
        coords = geometry.get("coordinates") or []
        if geometry.get("type") == "Point" and len(coords) >= 2:
            lons.append(float(coords[0]))
            lats.append(float(coords[1]))
        elif geometry.get("type") == "LineString":
            for coord in coords:
                if len(coord) >= 2:
                    lons.append(float(coord[0]))
                    lats.append(float(coord[1]))
    if not lons:
        return None
    west, east = min(lons), max(lons)
    south, north = min(lats), max(lats)
    pad_x = max((east - west) * 0.08, 0.02)
    pad_y = max((north - south) * 0.08, 0.02)
    return [west - pad_x, south - pad_y, east + pad_x, north + pad_y]


_CLOCK_LIMIT_MIN = 100000.0


def format_clock(moment):
    """Date and time for one activity boundary."""
    if moment is None:
        return ""
    return moment.strftime("%Y-%m-%d %H:%M")


def parse_clock(text, default=None):
    """Read a date and time. A bare time keeps the date from `default`."""
    text = (text or "").strip()
    if not text:
        return None
    for fmt in (
        "%Y-%m-%d %H:%M",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%dT%H:%M",
        "%Y-%m-%dT%H:%M:%S",
        "%m/%d/%Y %H:%M",
        "%m/%d/%Y %I:%M %p",
    ):
        try:
            return datetime.strptime(text, fmt).replace(second=0, microsecond=0)
        except ValueError:
            continue
    for fmt in ("%H:%M", "%H:%M:%S", "%I:%M %p"):
        try:
            clock = datetime.strptime(text, fmt)
        except ValueError:
            continue
        base = default or datetime.now()
        return base.replace(hour=clock.hour, minute=clock.minute, second=0, microsecond=0)
    return None


def apply_activity_clock(
    stops,
    kind,
    index,
    which,
    new_dt,
    start,
    speed_kts,
    transit_speeds,
    transit_leads,
    transit_durations,
):
    """Set one activity's start or end.

    With a start already set, a new start keeps that activity's length and moves
    it and every later activity. A new end changes that activity's length, and
    later activities follow. With no start set, the typed time places the whole
    sequence on the clock and leaves each length as it is.

    Returns (overview start, changed).
    """
    new_dt = _naive_clock(new_dt)
    start = _naive_clock(start)
    if new_dt is None or which not in ("start", "end"):
        return start, False
    if kind == "stop" and not (0 <= index < len(stops)):
        return start, False
    if kind == "transit" and not (0 <= index < len(stops) - 1):
        return start, False
    leads = transit_leads if transit_leads is not None else {}
    durations = transit_durations if transit_durations is not None else {}
    speeds = transit_speeds if transit_speeds is not None else {}
    anchor = start if start is not None else datetime(2000, 1, 1)
    rows = build_schedule(stops, speed_kts, anchor, speeds, leads, durations)
    row = _schedule_row(stops, rows, kind, index)
    if row is None or row.start is None or row.end is None:
        return start, False
    target = row.start if which == "start" else row.end
    if abs((new_dt - target).total_seconds()) < 30:
        return start, False
    if start is None:
        placed = new_dt - (target - anchor)
        return placed.replace(second=0, microsecond=0), True
    if which == "start":
        if kind == "stop" and index == 0:
            return new_dt, True
        changed = _add_lead(stops, kind, index, (new_dt - row.start).total_seconds() / 60.0, leads)
        return start, changed
    if kind == "stop":
        changed = _set_occupied(stops[index], (new_dt - row.start).total_seconds() / 60.0)
        return start, changed
    changed = _set_transit_duration(
        stops, index, (new_dt - row.start).total_seconds() / 60.0, speed_kts, speeds, durations,
    )
    return start, changed


def _schedule_row(stops, rows, kind, index):
    cursor = 0
    for row_kind, stop_index in activity_rows(stops):
        if row_kind == kind and stop_index == index:
            if cursor < len(rows):
                return rows[cursor]
            return None
        cursor += 1
    return None


def _add_lead(stops, kind, index, delta_min, transit_leads):
    if kind == "stop":
        updated = _clamp_clock_offset(float(stops[index].lead_min) + delta_min)
        if abs(updated - float(stops[index].lead_min)) < 0.05:
            return False
        stops[index].lead_min = 0.0 if abs(updated) < 0.05 else updated
        return True
    key = (stops[index].uid, stops[index + 1].uid)
    old = float(transit_leads.get(key, 0.0))
    updated = _clamp_clock_offset(old + delta_min)
    if abs(updated - old) < 0.05:
        return False
    if abs(updated) < 0.05:
        transit_leads.pop(key, None)
    else:
        transit_leads[key] = updated
    return True


def _set_occupied(stop, minutes):
    old = stop.occupied_minutes()
    cap = max(0.0, float(stop.duration_min)) + _CLOCK_LIMIT_MIN
    minutes = min(max(float(minutes), 0.0), cap)
    if abs(minutes - old) < 0.05:
        return False
    stop.adjust_min = minutes - float(stop.duration_min) - float(stop.extra_min)
    if abs(stop.adjust_min) < 0.05:
        stop.adjust_min = 0.0
    return True


def _set_transit_duration(stops, index, minutes, speed_kts, transit_speeds, transit_durations):
    key = (stops[index].uid, stops[index + 1].uid)
    leg_minutes = transit_leg(stops[index], stops[index + 1], speed_kts, transit_speeds)["minutes"]
    current = min(max(float(minutes), 0.0), _CLOCK_LIMIT_MIN)
    previous = transit_durations.get(key)
    if abs(current - leg_minutes) < 0.05:
        if previous is None:
            return False
        transit_durations.pop(key, None)
        return True
    if previous is not None and abs(float(previous) - current) < 0.05:
        return False
    transit_durations[key] = current
    return True


def _naive_clock(moment):
    if moment is None:
        return None
    if moment.tzinfo is not None:
        moment = moment.replace(tzinfo=None)
    return moment.replace(second=0, microsecond=0)


def _clamp_clock_offset(minutes):
    return min(max(float(minutes), -_CLOCK_LIMIT_MIN), _CLOCK_LIMIT_MIN)


def stops_to_document(stops, speed_kts, start, transit_speeds=None, transit_leads=None, transit_durations=None):
    ids = {stop.uid for stop in stops}
    speeds = []
    for (src, dst), value in (transit_speeds or {}).items():
        if src in ids and dst in ids:
            speeds.append({"from": src, "to": dst, "speed_kts": float(value)})
    payload = {
        "version": 1,
        "transit_speed_kts": float(speed_kts),
        "start": start.isoformat(timespec="minutes") if start is not None else "",
        "transit_speeds": speeds,
        "transit_leads": _stored_pairs(ids, transit_leads, "minutes"),
        "transit_durations": _stored_pairs(ids, transit_durations, "minutes"),
        "stops": [],
    }
    for stop in stops:
        payload["stops"].append({
            "uid": stop.uid,
            "kind": stop.kind,
            "name": stop.name,
            "extra_min": float(stop.extra_min),
            "lead_min": float(stop.lead_min),
            "adjust_min": float(stop.adjust_min),
            "info_path": stop.info_path,
            "flipped": bool(stop.flipped),
            "lat": stop.lat,
            "lon": stop.lon,
            "arrive_lat": stop.arrive_lat,
            "arrive_lon": stop.arrive_lon,
            "depart_lat": stop.depart_lat,
            "depart_lon": stop.depart_lon,
            "duration_min": float(stop.duration_min) if stop.kind == "point" else None,
        })
    return payload


def document_to_stops(payload):
    """Reload plans from an overview file.

    Returns (stops, default speed, start, warnings, transit speeds, transit leads,
    transit durations).
    """
    warnings = []
    stops = []
    for item in payload.get("stops") or []:
        kind = item.get("kind") or "plan"
        if kind == "plan":
            info_path = item.get("info_path") or ""
            stop, warning = load_plan(info_path)
            if warning:
                warnings.append(warning)
            stop.extra_min = _float_or(item.get("extra_min"), 0.0)
            stop.lead_min = _float_or(item.get("lead_min"), 0.0)
            stop.adjust_min = _float_or(item.get("adjust_min"), 0.0)
            stop.flipped = bool(item.get("flipped"))
            stop.arrive_lat, stop.arrive_lon = _coord_pair(item, "arrive_lat", "arrive_lon")
            stop.depart_lat, stop.depart_lon = _coord_pair(item, "depart_lat", "depart_lon")
            _apply_saved_uid(stop, item)
            stops.append(stop)
            continue
        stop = Stop(
            kind="point",
            name=item.get("name") or "Point",
            duration_min=_float_or(item.get("duration_min"), 0.0),
            extra_min=_float_or(item.get("extra_min"), 0.0),
            lead_min=_float_or(item.get("lead_min"), 0.0),
            adjust_min=_float_or(item.get("adjust_min"), 0.0),
            lat=_optional_float(item.get("lat")),
            lon=_optional_float(item.get("lon")),
        )
        _apply_saved_uid(stop, item)
        stops.append(stop)
    speed = _float_or(payload.get("transit_speed_kts"), 10.0)
    start = _parse_start(payload.get("start") or "")
    speeds = {}
    for item in payload.get("transit_speeds") or []:
        src = str(item.get("from") or "")
        dst = str(item.get("to") or "")
        if src and dst:
            speeds[(src, dst)] = _float_or(item.get("speed_kts"), speed)
    leads = _pair_map(payload.get("transit_leads"), "minutes")
    durations = _pair_map(payload.get("transit_durations"), "minutes")
    return stops, speed, start, warnings, speeds, leads, durations


def _stored_pairs(ids, mapping, field):
    stored = []
    for (src, dst), value in (mapping or {}).items():
        if src in ids and dst in ids:
            stored.append({"from": src, "to": dst, field: float(value)})
    return stored


def _pair_map(items, field):
    result = {}
    for item in items or []:
        src = str(item.get("from") or "")
        dst = str(item.get("to") or "")
        if src and dst:
            result[(src, dst)] = _float_or(item.get(field), 0.0)
    return result


def _apply_saved_uid(stop, item):
    saved = str(item.get("uid") or "").strip()
    if saved:
        stop.uid = saved


def _float_or(value, default):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _coord_pair(item, lat_key, lon_key):
    lat = _optional_float(item.get(lat_key))
    lon = _optional_float(item.get(lon_key))
    if lat is None or lon is None:
        return None, None
    return lat, lon


def _optional_float(value):
    try:
        if value is None or value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _parse_start(text):
    text = (text or "").strip()
    if not text:
        return None
    for fmt in ("%Y-%m-%dT%H:%M", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M"):
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            continue
    return None

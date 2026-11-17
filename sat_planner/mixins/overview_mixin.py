"""
Overview tab: several finished test plans, their order, transits, and a schedule.
"""
import csv
import json
import os
from datetime import datetime

from PyQt6.QtCore import QDate, QDateTime, QPoint, QTime, Qt
from PyQt6.QtWidgets import (
    QCheckBox,
    QDateTimeEdit,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QAbstractItemView,
    QHeaderView,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QStyledItemDelegate,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from sat_planner.overview_planner import (
    Stop,
    activity_rows,
    apply_activity_clock,
    build_schedule,
    document_to_stops,
    clear_end_overrides,
    end_snap_limit_m,
    feature_bounds,
    format_clock,
    format_schedule,
    legend_items,
    load_plan,
    map_features,
    nearest_line_end,
    overview_kml_text,
    overview_export_lines,
    overview_export_points,
    parse_clock,
    stops_to_document,
    swap_endpoints,
    transit_leg,
)

_OVERVIEW_COL_ORDER = 0
_OVERVIEW_COL_TYPE = 1
_OVERVIEW_COL_NAME = 2
_OVERVIEW_COL_TIME = 3
_OVERVIEW_COL_EXTRA = 4
_OVERVIEW_COL_SPEED = 5
_OVERVIEW_COL_TOTAL = 6
_OVERVIEW_COL_START = 7
_OVERVIEW_COL_END = 8
_OVERVIEW_COL_NOTE = 9

_OVERVIEW_EXPORT_CHOICES = (
    ("esri_shapefile", "ESRI Shapefile (.shp)"),
    ("gpkg", "GeoPackage (.gpkg)"),
    ("kml", "KML (.kml)"),
    ("table_csv", "Table (.csv)"),
    ("info_txt", "Info (.txt)"),
    ("map_png_high", "Map PNG — high resolution"),
    ("map_png_low", "Map PNG — low resolution (email)"),
)


def _as_qdatetime(moment):
    return QDateTime(
        QDate(moment.year, moment.month, moment.day),
        QTime(moment.hour, moment.minute),
    )


class _OverviewClockDelegate(QStyledItemDelegate):
    """Date and time editor for a Start or End cell."""

    def __init__(self, parent, fallback):
        super().__init__(parent)
        self._fallback = fallback

    def createEditor(self, parent, option, index):
        editor = QDateTimeEdit(parent)
        editor.setDisplayFormat("yyyy-MM-dd HH:mm")
        editor.setCalendarPopup(True)
        editor.setKeyboardTracking(False)
        return editor

    def setEditorData(self, editor, index):
        parsed = parse_clock(index.data(Qt.ItemDataRole.DisplayRole) or "")
        if parsed is None:
            parsed = self._fallback()
        editor.setDateTime(_as_qdatetime(parsed))

    def setModelData(self, editor, model, index):
        model.setData(
            index,
            editor.dateTime().toString("yyyy-MM-dd HH:mm"),
            Qt.ItemDataRole.EditRole,
        )


class OverviewMixin:
    """Load exported plans and show the operation order on the imagery map."""

    def _build_overview_tab(self):
        frame = getattr(self, "overview_frame", None)
        if frame is None:
            return
        self.overview_stops = []
        self.overview_transit_speeds = {}
        self.overview_transit_leads = {}
        self.overview_transit_durations = {}
        self._overview_row_refs = []
        self._overview_map_active = False
        self._overview_filling = False
        self._overview_end_pick = None

        layout = QVBoxLayout(frame)
        layout.setSpacing(4)

        note = QLabel(
            "Load finished test plans. Their times come from the info files. "
            "The map uses imagery and shows the order, including the steam between stops. "
            "Start and End are the date and time of each activity."
        )
        note.setWordWrap(True)
        layout.addWidget(note)

        file_row = QHBoxLayout()
        load_btn = QPushButton("Load Test Plans")
        load_btn.clicked.connect(self._overview_load_plans)
        add_btn = QPushButton("Add Point")
        add_btn.setToolTip("Type a coordinate, or click the map while the point dialog is open.")
        add_btn.clicked.connect(self._overview_add_point)
        open_btn = QPushButton("Open Overview")
        open_btn.clicked.connect(self._overview_open)
        save_btn = QPushButton("Save Overview")
        save_btn.setToolTip("Choose a folder and a name, then pick which files to write.")
        save_btn.clicked.connect(self._overview_save)
        for button in (load_btn, add_btn, open_btn, save_btn):
            button.setMinimumWidth(0)
            file_row.addWidget(button)
        layout.addLayout(file_row)

        speed_row = QHBoxLayout()
        speed_row.addWidget(QLabel("Default transit speed"))
        self.overview_speed_spin = QDoubleSpinBox()
        self.overview_speed_spin.setRange(0.1, 40.0)
        self.overview_speed_spin.setDecimals(1)
        self.overview_speed_spin.setValue(10.0)
        self.overview_speed_spin.setSuffix(" kn")
        self.overview_speed_spin.setToolTip("Used for a transit row unless that row has its own speed.")
        self.overview_speed_spin.valueChanged.connect(self._overview_inputs_changed)
        speed_row.addWidget(self.overview_speed_spin)
        layout.addLayout(speed_row)

        start_row = QHBoxLayout()
        self.overview_start_check = QCheckBox("Start")
        self.overview_start_check.toggled.connect(self._overview_inputs_changed)
        self.overview_start_edit = QDateTimeEdit(QDateTime.currentDateTime())
        self.overview_start_edit.setDisplayFormat("yyyy-MM-dd HH:mm")
        self.overview_start_edit.setCalendarPopup(True)
        self.overview_start_edit.dateTimeChanged.connect(self._overview_inputs_changed)
        self.overview_start_edit.setEnabled(False)
        start_row.addWidget(self.overview_start_check)
        start_row.addWidget(self.overview_start_edit, 1)
        layout.addLayout(start_row)

        self.overview_list = QTableWidget(0, 10)
        self.overview_list.setHorizontalHeaderLabels(
            ["#", "Type", "Name", "Time", "Extra", "Speed", "Total", "Start", "End", "Note"]
        )
        self.overview_list.verticalHeader().setVisible(False)
        self.overview_list.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.overview_list.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.overview_list.setEditTriggers(
            QAbstractItemView.EditTrigger.DoubleClicked
            | QAbstractItemView.EditTrigger.SelectedClicked
            | QAbstractItemView.EditTrigger.EditKeyPressed
        )
        self.overview_list.setAlternatingRowColors(True)
        self.overview_list.setMinimumHeight(160)
        self.overview_list.setWordWrap(False)
        self.overview_list.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.overview_list.setHorizontalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        header = self.overview_list.horizontalHeader()
        header.setStretchLastSection(False)
        header.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        for column, width in (
            (_OVERVIEW_COL_ORDER, 36),
            (_OVERVIEW_COL_TYPE, 96),
            (_OVERVIEW_COL_NAME, 150),
            (_OVERVIEW_COL_TIME, 64),
            (_OVERVIEW_COL_EXTRA, 64),
            (_OVERVIEW_COL_SPEED, 64),
            (_OVERVIEW_COL_TOTAL, 64),
            (_OVERVIEW_COL_START, 156),
            (_OVERVIEW_COL_END, 156),
            (_OVERVIEW_COL_NOTE, 180),
        ):
            self.overview_list.setColumnWidth(column, width)
        clock_delegate = _OverviewClockDelegate(self.overview_list, self._overview_clock_fallback)
        self.overview_list.setItemDelegateForColumn(_OVERVIEW_COL_START, clock_delegate)
        self.overview_list.setItemDelegateForColumn(_OVERVIEW_COL_END, clock_delegate)
        self.overview_list.horizontalHeaderItem(_OVERVIEW_COL_TIME).setToolTip(
            "Minutes from the info file, or the time entered for a point."
        )
        self.overview_list.horizontalHeaderItem(_OVERVIEW_COL_EXTRA).setToolTip(
            "Extra minutes at this stop. Edit this cell. The survey file is not changed."
        )
        self.overview_list.horizontalHeaderItem(_OVERVIEW_COL_SPEED).setToolTip(
            "Speed for this transit. Clear the cell to use the default transit speed. "
            "Editing this speed replaces an end time set on this transit."
        )
        self.overview_list.horizontalHeaderItem(_OVERVIEW_COL_TOTAL).setToolTip(
            "Plan time plus extra time, including a change made by editing the end time."
        )
        self.overview_list.horizontalHeaderItem(_OVERVIEW_COL_START).setToolTip(
            "Date and time this activity starts. Changing it moves this activity and "
            "every activity after it. The length of this activity stays the same."
        )
        self.overview_list.horizontalHeaderItem(_OVERVIEW_COL_END).setToolTip(
            "Date and time this activity ends. Changing it makes this activity longer "
            "or shorter, and every activity after it moves to follow."
        )
        self.overview_list.itemChanged.connect(self._overview_extra_edited)
        self.overview_list.currentCellChanged.connect(self._overview_selection_changed)
        layout.addWidget(self.overview_list)

        order_row = QHBoxLayout()
        up_btn = QPushButton("Up")
        up_btn.clicked.connect(lambda: self._overview_move(-1))
        down_btn = QPushButton("Down")
        down_btn.clicked.connect(lambda: self._overview_move(1))
        flip_btn = QPushButton("Swap ends")
        flip_btn.setToolTip("Swap this survey's overview start and end. The survey file stays unchanged.")
        flip_btn.clicked.connect(self._overview_flip)
        edit_btn = QPushButton("Edit Point")
        edit_btn.setToolTip("Change this point by typing, or click the map while the dialog is open.")
        edit_btn.clicked.connect(self._overview_edit_point)
        remove_btn = QPushButton("Remove")
        remove_btn.clicked.connect(self._overview_remove)
        for button in (up_btn, down_btn, flip_btn, edit_btn, remove_btn):
            button.setMinimumWidth(0)
            order_row.addWidget(button)
        layout.addLayout(order_row)

        ends_row = QHBoxLayout()
        start_btn = QPushButton("Set start")
        start_btn.setToolTip("Click a line end to choose where this survey starts in the overview.")
        start_btn.clicked.connect(lambda: self._overview_begin_end_pick("arrive"))
        end_btn = QPushButton("Set end")
        end_btn.setToolTip("Click a line end to choose where this survey ends in the overview.")
        end_btn.clicked.connect(lambda: self._overview_begin_end_pick("depart"))
        reset_btn = QPushButton("Use plan ends")
        reset_btn.setToolTip("Use the start and end from the survey file.")
        reset_btn.clicked.connect(self._overview_reset_ends)
        for button in (start_btn, end_btn, reset_btn):
            button.setMinimumWidth(0)
            ends_row.addWidget(button)
        layout.addLayout(ends_row)
        self.overview_end_hint = QLabel("")
        self.overview_end_hint.setWordWrap(True)
        self.overview_end_hint.hide()
        layout.addWidget(self.overview_end_hint)

        self.overview_schedule = QPlainTextEdit()
        self.overview_schedule.setReadOnly(True)
        self.overview_schedule.setMinimumHeight(180)
        self.overview_schedule.setPlaceholderText("The schedule appears here.")
        layout.addWidget(self.overview_schedule, 1)
        self._fill_overview_schedule()

    def _overview_directory(self):
        return getattr(self, "last_overview_dir", None) or getattr(self, "last_used_dir", "") or ""

    def _overview_remember_directory(self, path):
        folder = os.path.dirname(path)
        if not folder or not os.path.isdir(folder):
            return
        self.last_overview_dir = folder
        if hasattr(self, "_save_last_overview_dir"):
            self._save_last_overview_dir()

    def _overview_message(self, text):
        if hasattr(self, "_show_message"):
            self._show_message("info", "Overview", text)
            return
        from PyQt6.QtWidgets import QMessageBox

        QMessageBox.information(self, "Overview", text)

    def _overview_load_plans(self):
        paths, _selected = QFileDialog.getOpenFileNames(
            self,
            "Load Test Plans",
            self._overview_directory(),
            "Test info files (*_info.txt);;Text files (*.txt);;All files (*.*)",
        )
        if not paths:
            return
        self.last_overview_dir = os.path.dirname(paths[0])
        warnings = []
        for path in paths:
            stop, warning = load_plan(path)
            self.overview_stops.append(stop)
            if warning:
                warnings.append(warning)
        self._fill_overview_list(select=len(self.overview_stops) - 1)
        self._refresh_overview(fit=True)
        if warnings:
            self._overview_message("\n".join(warnings))

    def _overview_add_point(self):
        def commit(values):
            self.overview_stops.append(Stop(kind="point", **values))
            self._fill_overview_list(select=len(self.overview_stops) - 1)
            self._refresh_overview(fit=True)

        self._open_overview_point_dialog(None, commit)

    def _overview_current_ref(self):
        row = self.overview_list.currentRow() if hasattr(self, "overview_list") else -1
        refs = getattr(self, "_overview_row_refs", [])
        if 0 <= row < len(refs):
            return refs[row]
        return ("stop", -1)

    def _overview_selected_stop_index(self):
        kind, index = self._overview_current_ref()
        if kind == "stop" and 0 <= index < len(self.overview_stops):
            return index
        return -1

    def _overview_table_row_for(self, kind, index):
        for row, ref in enumerate(getattr(self, "_overview_row_refs", [])):
            if ref == (kind, index):
                return row
        return -1

    def _overview_edit_point(self):
        index = self._overview_selected_stop_index()
        if index < 0 or index >= len(self.overview_stops):
            self._overview_message("Select a point to edit.")
            return
        stop = self.overview_stops[index]
        if stop.kind != "point":
            self._overview_message("Select a point to edit. Test plans keep the time from their info file.")
            return

        def commit(values, index=index):
            if index < 0 or index >= len(self.overview_stops):
                return
            edited = self.overview_stops[index]
            if edited.kind != "point":
                return
            edited.name = values["name"]
            edited.lat = values["lat"]
            edited.lon = values["lon"]
            edited.duration_min = values["duration_min"]
            self._fill_overview_list(select=index)
            self._refresh_overview(fit=False)

        self._open_overview_point_dialog(stop, commit)

    def _open_overview_point_dialog(self, stop, on_accept):
        """Keep the point dialog open so a map click can fill the coordinate."""
        previous = getattr(self, "_overview_point_dialog_widget", None)
        if previous is not None:
            previous.reject()
        self._overview_cancel_end_pick(redraw=False)
        dialog = QDialog(self)
        dialog.setWindowTitle("Overview point")
        dialog.setModal(False)
        dialog.setWindowModality(Qt.WindowModality.NonModal)
        dialog.setWindowFlags(
            Qt.WindowType.Tool | Qt.WindowType.WindowTitleHint | Qt.WindowType.WindowCloseButtonHint
        )
        dialog.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        form = QFormLayout(dialog)
        hint = QLabel("Click the map to fill the coordinate. Typing still works.")
        hint.setWordWrap(True)
        form.addRow(hint)
        name_edit = QLineEdit(stop.name if stop is not None else "Port")
        lat_edit = QLineEdit("" if stop is None or stop.lat is None else f"{stop.lat:.6f}")
        lon_edit = QLineEdit("" if stop is None or stop.lon is None else f"{stop.lon:.6f}")
        time_spin = QDoubleSpinBox()
        time_spin.setRange(0.0, 100000.0)
        time_spin.setDecimals(1)
        time_spin.setSuffix(" min")
        time_spin.setValue(0.0 if stop is None else float(stop.duration_min))
        form.addRow("Name", name_edit)
        form.addRow("Latitude", lat_edit)
        form.addRow("Longitude", lon_edit)
        form.addRow("Time at point", time_spin)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        form.addRow(buttons)

        def accept_point():
            values = self._overview_point_values(name_edit, lat_edit, lon_edit, time_spin)
            if values is None:
                return
            self._overview_point_pick = None
            self._overview_pick_preview = None
            on_accept(values)
            dialog.done(QDialog.DialogCode.Accepted)

        buttons.accepted.connect(accept_point)
        buttons.rejected.connect(dialog.reject)
        dialog.finished.connect(self._overview_point_dialog_closed)
        self._overview_point_dialog_widget = dialog
        self._overview_point_pick = (lat_edit, lon_edit)
        if stop is not None and stop.lat is not None and stop.lon is not None:
            self._overview_pick_preview = (float(stop.lat), float(stop.lon))
            self._show_overview_on_map(fit=False, legend=False)
        else:
            self._overview_pick_preview = None
        dialog.adjustSize()
        dialog.move(self.mapToGlobal(QPoint(16, 72)))
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()
        view = getattr(self, "map_view", None)
        if view is not None and hasattr(view, "set_cursor"):
            view.set_cursor("crosshair")

    def _overview_point_values(self, name_edit, lat_edit, lon_edit, time_spin):
        name = name_edit.text().strip() or "Point"
        try:
            lat = float(lat_edit.text().strip())
            lon = float(lon_edit.text().strip())
        except ValueError:
            self._overview_message("Enter a numeric latitude and longitude.")
            return None
        if not -90.0 <= lat <= 90.0 or not -180.0 <= lon <= 180.0:
            self._overview_message("Latitude must be between -90 and 90, and longitude between -180 and 180.")
            return None
        return {"name": name, "lat": lat, "lon": lon, "duration_min": float(time_spin.value())}

    def _overview_apply_map_click(self, lon, lat):
        """Fill the open point dialog, or set a survey end, from a map click."""
        if getattr(self, "_overview_point_pick", None):
            self._overview_apply_point_click(lon, lat)
            return
        if getattr(self, "_overview_end_pick", None):
            self._overview_apply_end_click(lon, lat)

    def _overview_apply_point_click(self, lon, lat):
        """Fill the open point dialog from a map click and mark that spot."""
        fields = getattr(self, "_overview_point_pick", None)
        if not fields:
            return
        lat_edit, lon_edit = fields
        try:
            lat_edit.setText(f"{float(lat):.6f}")
            lon_edit.setText(f"{float(lon):.6f}")
        except RuntimeError:
            self._overview_point_pick = None
            return
        self._overview_pick_preview = (float(lat), float(lon))
        if self._overview_tab_is_current():
            self._show_overview_on_map(fit=False, legend=False)

    def _overview_point_dialog_closed(self, _result=None):
        self._overview_point_pick = None
        self._overview_pick_preview = None
        self._overview_point_dialog_widget = None
        view = getattr(self, "map_view", None)
        if view is not None and hasattr(view, "set_cursor"):
            view.set_cursor("")
        if self._overview_tab_is_current() or getattr(self, "_overview_map_active", False):
            self._show_overview_on_map(fit=False)

    def _overview_config_path(self):
        return getattr(self, "CONFIG_FILENAME", None)

    def _overview_load_config(self):
        """Return (settings, readable). readable is false when the file is there but cannot be read."""
        path = self._overview_config_path()
        if not path or not os.path.isfile(path):
            return {}, True
        for encoding in ("utf-8-sig", None):
            try:
                with open(path, encoding=encoding) as handle:
                    data = json.load(handle)
            except (OSError, json.JSONDecodeError, UnicodeError):
                continue
            return (data if isinstance(data, dict) else {}), True
        return {}, False

    def _overview_read_config(self):
        data, _readable = self._overview_load_config()
        return data

    def _overview_update_config(self, updates):
        path = self._overview_config_path()
        if not path:
            return
        config, readable = self._overview_load_config()
        if not readable:
            return
        config.update(updates)
        try:
            with open(path, "w", encoding="utf-8") as handle:
                json.dump(config, handle)
        except OSError:
            pass

    def _default_overview_export_options(self):
        return {key: True for key, _label in _OVERVIEW_EXPORT_CHOICES}

    def _overview_export_options(self):
        options = self._default_overview_export_options()
        saved = getattr(self, "overview_export_options", None) or {}
        for key in options:
            if key in saved:
                options[key] = bool(saved[key])
        return options

    def _overview_clean_name(self, text):
        cleaned = "".join(ch if ch not in '<>:"/\\|?*' and ord(ch) >= 32 else "_" for ch in (text or "").strip())
        return cleaned.strip(" .")

    def _overview_save(self):
        directory = QFileDialog.getExistingDirectory(
            self,
            "Select Overview Folder",
            self._overview_directory(),
        )
        if not directory:
            return
        chosen = self._overview_save_dialog()
        if chosen is None:
            return
        name, options = chosen
        self.overview_export_options = dict(options)
        if hasattr(self, "_save_overview_export_options"):
            self._save_overview_export_options()
        written, errors = self._overview_write_export(directory, name, options)
        self.last_overview_dir = directory
        self._overview_update_config({
            "last_overview_dir": directory,
            "last_overview_name": name,
            "overview_export_options": options,
        })
        lines = [f"Saved the overview in {directory}:"]
        lines.extend(f"- {os.path.basename(path)}" for path in written)
        if errors:
            lines.append("")
            lines.append("Some files could not be written:")
            lines.extend(f"- {os.path.basename(path)}: {message}" for path, message in errors)
        self._overview_message("\n".join(lines))

    def _overview_save_dialog(self):
        """Name and export choices. Returns (name, options) or None."""
        dialog = QDialog(self)
        dialog.setWindowTitle("Save Overview")
        dialog.setMinimumWidth(420)
        layout = QVBoxLayout(dialog)
        form = QFormLayout()
        name_edit = QLineEdit(str(self._overview_read_config().get("last_overview_name") or ""))
        name_edit.setPlaceholderText("Overview name")
        form.addRow("Name", name_edit)
        layout.addLayout(form)
        saved = self._overview_export_options()
        boxes = {}
        for key, label in _OVERVIEW_EXPORT_CHOICES:
            box = QCheckBox(label)
            box.setChecked(bool(saved.get(key)))
            layout.addWidget(box)
            boxes[key] = box
        note = QLabel(
            "The overview is always saved as a .json file using this name. "
            "Table is the activity table. Info is the schedule at the bottom of this tab."
        )
        note.setWordWrap(True)
        layout.addWidget(note)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )

        def accept_named():
            if not self._overview_clean_name(name_edit.text()):
                self._overview_message("Enter a name for the overview.")
                return
            dialog.accept()

        buttons.accepted.connect(accept_named)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return None
        name = self._overview_clean_name(name_edit.text())
        options = {key: box.isChecked() for key, box in boxes.items()}
        return name, options

    def _overview_write_export(self, directory, name, options):
        written = []
        errors = []
        json_path = os.path.join(directory, f"{name}.json")
        payload = stops_to_document(
            self.overview_stops,
            self.overview_speed_spin.value(),
            self._overview_start(),
            self.overview_transit_speeds,
            self.overview_transit_leads,
            self.overview_transit_durations,
        )
        try:
            with open(json_path, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2)
            written.append(json_path)
        except OSError as exc:
            errors.append((json_path, str(exc)))
            return written, errors
        if options.get("table_csv"):
            self._overview_write_table_csv(directory, name, written, errors)
        if options.get("info_txt"):
            self._overview_write_info_txt(directory, name, written, errors)
        if options.get("esri_shapefile") or options.get("gpkg"):
            self._overview_write_vectors(directory, name, options, written, errors)
        if options.get("kml"):
            self._overview_write_kml(directory, name, written, errors)
        if options.get("map_png_high") or options.get("map_png_low"):
            self._overview_write_map_png(directory, name, options, written, errors)
        return written, errors

    def _overview_write_table_csv(self, directory, name, written, errors):
        path = os.path.join(directory, f"{name}.csv")
        headers = []
        table = self.overview_list
        for column in range(table.columnCount()):
            item = table.horizontalHeaderItem(column)
            headers.append(item.text() if item is not None else "")
        try:
            with open(path, "w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(headers)
                for row in range(table.rowCount()):
                    values = []
                    for column in range(table.columnCount()):
                        item = table.item(row, column)
                        values.append(item.text() if item is not None else "")
                    writer.writerow(values)
            written.append(path)
        except OSError as exc:
            errors.append((path, str(exc)))

    def _overview_write_info_txt(self, directory, name, written, errors):
        path = os.path.join(directory, f"{name}_info.txt")
        self._fill_overview_schedule()
        text = self.overview_schedule.toPlainText() if hasattr(self, "overview_schedule") else ""
        try:
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(text)
                if text and not text.endswith("\n"):
                    handle.write("\n")
            written.append(path)
        except OSError as exc:
            errors.append((path, str(exc)))

    def _overview_vector_records(self):
        line_schema = {
            "geometry": "LineString",
            "properties": {"seq": "int", "kind": "str:20", "name": "str:80"},
        }
        point_schema = {
            "geometry": "Point",
            "properties": {"seq": "int", "name": "str:80"},
        }
        line_features = []
        for seq, (kind, label, points) in enumerate(overview_export_lines(self.overview_stops), start=1):
            line_features.append({
                "geometry": {
                    "type": "LineString",
                    "coordinates": [[lon, lat] for lat, lon in points],
                },
                "properties": {"seq": seq, "kind": kind, "name": label[:80]},
            })
        point_features = []
        for seq, (label, lat, lon) in enumerate(overview_export_points(self.overview_stops), start=1):
            point_features.append({
                "geometry": {"type": "Point", "coordinates": [lon, lat]},
                "properties": {"seq": seq, "name": label[:80]},
            })
        return line_schema, line_features, point_schema, point_features

    def _overview_write_vectors(self, directory, name, options, written, errors):
        line_schema, line_features, point_schema, point_features = self._overview_vector_records()
        if not line_features and not point_features:
            errors.append((os.path.join(directory, f"{name}.shp"), "No positions to export"))
            return
        try:
            import fiona
        except ImportError:
            errors.append((os.path.join(directory, f"{name}.shp"), "Fiona is not available"))
            return
        from sat_planner.export_utils import remove_export_file

        jobs = []
        if line_features:
            jobs.append((f"{name}.shp", f"{name}.gpkg", "lines", line_schema, line_features))
        if point_features:
            jobs.append((f"{name}_points.shp", f"{name}_points.gpkg", "points", point_schema, point_features))
        for shp_name, gpkg_name, layer, schema, features in jobs:
            if options.get("esri_shapefile"):
                path = os.path.join(directory, shp_name)
                try:
                    remove_export_file(path)
                    with fiona.open(path, "w", driver="ESRI Shapefile", crs="EPSG:4326", schema=schema) as collection:
                        collection.writerecords(features)
                    written.append(path)
                except Exception as exc:
                    errors.append((path, str(exc)))
            if options.get("gpkg"):
                path = os.path.join(directory, gpkg_name)
                try:
                    remove_export_file(path)
                    with fiona.open(
                        path, "w", driver="GPKG", crs="EPSG:4326", schema=schema, layer=layer,
                    ) as collection:
                        collection.writerecords(features)
                    written.append(path)
                except Exception as exc:
                    errors.append((path, str(exc)))

    def _overview_write_kml(self, directory, name, written, errors):
        path = os.path.join(directory, f"{name}.kml")
        text = overview_kml_text(self.overview_stops, name)
        if "<Placemark>" not in text:
            errors.append((path, "No positions to export"))
            return
        try:
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(text)
            written.append(path)
        except OSError as exc:
            errors.append((path, str(exc)))

    def _overview_write_map_png(self, directory, name, options, written, errors):
        from sat_planner.export_utils import low_resolution_png_path, remove_export_file, save_existing_png_with_low_res

        view = getattr(self, "map_view", None)
        path = os.path.join(directory, f"{name}_map.png")
        if view is None or not hasattr(view, "capture_png"):
            errors.append((path, "The map is not available"))
            return
        data = view.capture_png()
        if not data:
            errors.append((path, "The map picture was not ready"))
            return
        save_high = bool(options.get("map_png_high"))
        save_low = bool(options.get("map_png_low"))
        try:
            remove_export_file(path)
            with open(path, "wb") as handle:
                handle.write(data)
            if save_low:
                low_path = save_existing_png_with_low_res(path)
                if low_path:
                    written.append(low_path)
                else:
                    errors.append((low_resolution_png_path(path), "Could not write the smaller map picture"))
            if save_high:
                written.append(path)
            else:
                remove_export_file(path)
        except OSError as exc:
            errors.append((path, str(exc)))

    def _overview_open(self):
        path, _selected = QFileDialog.getOpenFileName(
            self,
            "Open Overview",
            self._overview_directory(),
            "Overview (*.json);;All files (*.*)",
        )
        if not path:
            return
        self._overview_remember_directory(path)
        try:
            with open(path, encoding="utf-8") as handle:
                payload = json.load(handle)
        except (OSError, json.JSONDecodeError) as exc:
            self._overview_message(f"Could not open the overview.\n{exc}")
            return
        stops, speed, start, warnings, speeds, leads, durations = document_to_stops(payload)
        self.overview_stops = stops
        self.overview_transit_speeds = speeds
        self.overview_transit_leads = leads
        self.overview_transit_durations = durations
        self.overview_speed_spin.blockSignals(True)
        self.overview_speed_spin.setValue(speed if speed > 0 else 10.0)
        self.overview_speed_spin.blockSignals(False)
        self.overview_start_check.blockSignals(True)
        self.overview_start_edit.blockSignals(True)
        if start is None:
            self.overview_start_check.setChecked(False)
        else:
            self.overview_start_check.setChecked(True)
            self.overview_start_edit.setDateTime(QDateTime(
                QDate(start.year, start.month, start.day),
                QTime(start.hour, start.minute),
            ))
        self.overview_start_check.blockSignals(False)
        self.overview_start_edit.blockSignals(False)
        self._fill_overview_list(select=0 if stops else -1)
        self._refresh_overview(fit=True)
        if warnings:
            self._overview_message("\n".join(warnings))

    def _overview_move(self, delta):
        index = self._overview_selected_stop_index()
        target = index + delta
        if index < 0 or target < 0 or target >= len(self.overview_stops):
            if index < 0:
                self._overview_message("Select a survey or point to move.")
            return
        stops = self.overview_stops
        stops[index], stops[target] = stops[target], stops[index]
        self._fill_overview_list(select=target)
        self._refresh_overview(fit=False)

    def _overview_flip(self):
        index = self._overview_selected_stop_index()
        if index < 0 or index >= len(self.overview_stops):
            return
        stop = self.overview_stops[index]
        if not swap_endpoints(stop):
            return
        self._fill_overview_list(select=index)
        self._refresh_overview(fit=False)

    def _overview_begin_end_pick(self, mode):
        index = self._overview_selected_stop_index()
        if index < 0 or index >= len(self.overview_stops):
            self._overview_message("Select a survey first.")
            return
        stop = self.overview_stops[index]
        if stop.kind != "plan":
            self._overview_message("Select a survey. Points already use one coordinate.")
            return
        if not stop.lines:
            self._overview_message("This survey has no lines to pick an end from.")
            return
        if getattr(self, "_overview_end_pick", None) == (mode, index):
            self._overview_cancel_end_pick()
            return
        dialog = getattr(self, "_overview_point_dialog_widget", None)
        if dialog is not None:
            dialog.reject()
        self._overview_end_pick = (mode, index)
        if mode == "arrive":
            hint = "Click a line end on this survey to set the start. The survey file stays unchanged."
        else:
            hint = "Click a line end on this survey to set the end. The survey file stays unchanged."
        self.overview_end_hint.setText(hint)
        self.overview_end_hint.show()
        view = getattr(self, "map_view", None)
        if view is not None and hasattr(view, "set_cursor"):
            view.set_cursor("crosshair")

    def _overview_apply_end_click(self, lon, lat):
        pick = getattr(self, "_overview_end_pick", None)
        if not pick:
            return
        mode, index = pick
        if index < 0 or index >= len(self.overview_stops):
            self._overview_cancel_end_pick()
            return
        stop = self.overview_stops[index]
        point, meters = nearest_line_end(stop, float(lat), float(lon))
        if point is None or meters is None:
            self._overview_message("This survey has no line ends to use.")
            return
        if meters > end_snap_limit_m(stop):
            self._overview_message("Click nearer a line end of the selected survey.")
            return
        if mode == "arrive":
            stop.arrive_lat, stop.arrive_lon = point
        else:
            stop.depart_lat, stop.depart_lon = point
        self._overview_cancel_end_pick(redraw=False)
        self._fill_overview_list(select=index)
        self._refresh_overview(fit=False)

    def _overview_reset_ends(self):
        index = self._overview_selected_stop_index()
        if index < 0 or index >= len(self.overview_stops):
            return
        stop = self.overview_stops[index]
        if stop.kind != "plan":
            return
        clear_end_overrides(stop)
        self._overview_cancel_end_pick(redraw=False)
        self._fill_overview_list(select=index)
        self._refresh_overview(fit=False)

    def _overview_cancel_end_pick(self, redraw=True):
        was = getattr(self, "_overview_end_pick", None) is not None
        self._overview_end_pick = None
        hint = getattr(self, "overview_end_hint", None)
        if hint is not None:
            hint.hide()
        if not was:
            return
        view = getattr(self, "map_view", None)
        if view is not None and hasattr(view, "set_cursor"):
            view.set_cursor("")
        if redraw and (self._overview_tab_is_current() or getattr(self, "_overview_map_active", False)):
            self._show_overview_on_map(fit=False)

    def _overview_remove(self):
        index = self._overview_selected_stop_index()
        if index < 0 or index >= len(self.overview_stops):
            self._overview_message("Select a survey or point to remove.")
            return
        removed = self.overview_stops[index].uid
        for name in ("overview_transit_speeds", "overview_transit_leads", "overview_transit_durations"):
            store = getattr(self, name, {})
            setattr(self, name, {
                key: value for key, value in store.items() if removed not in key
            })
        del self.overview_stops[index]
        next_index = min(index, len(self.overview_stops) - 1)
        self._fill_overview_list(select=next_index)
        self._refresh_overview(fit=True)

    def _overview_inputs_changed(self, *_args):
        if self._overview_filling:
            return
        self.overview_start_edit.setEnabled(self.overview_start_check.isChecked())
        self._fill_overview_list()
        self._refresh_overview(fit=False)

    def _overview_selection_changed(self, current_row=None, _column=None, previous_row=None, _previous_column=None):
        if self._overview_filling:
            return
        if current_row is not None and previous_row is not None and current_row == previous_row:
            return
        index = self._overview_selected_stop_index()
        pick = getattr(self, "_overview_end_pick", None)
        if pick is not None and pick[1] != index:
            self._overview_cancel_end_pick(redraw=False)
        if getattr(self, "_overview_map_active", False):
            self._show_overview_on_map(fit=False)

    def _overview_start(self):
        if not self.overview_start_check.isChecked():
            return None
        value = self.overview_start_edit.dateTime().toPyDateTime()
        if value.tzinfo is not None:
            value = value.replace(tzinfo=None)
        return value.replace(second=0, microsecond=0)

    def _overview_note_text(self, stop):
        parts = []
        if stop.flipped:
            parts.append("ends swapped")
        if stop.arrive_lat is not None:
            parts.append("custom start")
        if stop.depart_lat is not None:
            parts.append("custom end")
        if abs(float(stop.adjust_min)) >= 0.05:
            parts.append("time adjusted")
        if stop.note:
            parts.append(stop.note)
        return "; ".join(parts)

    def _overview_text_item(self, text, editable=False, align_right=False):
        item = QTableWidgetItem(text)
        flags = Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled
        if editable:
            flags |= Qt.ItemFlag.ItemIsEditable
        item.setFlags(flags)
        if align_right:
            item.setTextAlignment(int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        return item

    def _overview_default_speed(self):
        spin = getattr(self, "overview_speed_spin", None)
        return spin.value() if spin is not None else 10.0

    def _overview_extra_edited(self, item):
        if self._overview_filling:
            return
        refs = getattr(self, "_overview_row_refs", [])
        row = item.row()
        if row < 0 or row >= len(refs):
            return
        kind, index = refs[row]
        column = item.column()
        text = item.text()
        if column == _OVERVIEW_COL_EXTRA and kind == "stop":
            self._overview_apply_extra(index, text)
        elif column == _OVERVIEW_COL_NAME and kind == "stop":
            self._overview_apply_point_name(index, text)
        elif column == _OVERVIEW_COL_SPEED and kind == "transit":
            self._overview_apply_transit_speed(index, text)
        elif column in (_OVERVIEW_COL_START, _OVERVIEW_COL_END):
            which = "start" if column == _OVERVIEW_COL_START else "end"
            self._overview_apply_clock(kind, index, which, text)

    def _overview_apply_extra(self, index, text):
        if index < 0 or index >= len(self.overview_stops):
            return
        stop = self.overview_stops[index]
        try:
            value = float(str(text).strip() or "0")
        except ValueError:
            value = float(stop.extra_min)
        stop.extra_min = min(max(value, 0.0), 100000.0)
        self._fill_overview_list(select=index)
        self._refresh_overview(fit=False)

    def _overview_apply_point_name(self, index, text):
        if index < 0 or index >= len(self.overview_stops):
            return
        stop = self.overview_stops[index]
        if stop.kind != "point":
            return
        stop.name = str(text).strip() or "Point"
        self._fill_overview_list(select=index, column=_OVERVIEW_COL_NAME)
        self._refresh_overview(fit=False)

    def _overview_apply_transit_speed(self, index, text):
        if index < 0 or index + 1 >= len(self.overview_stops):
            return
        left = self.overview_stops[index]
        right = self.overview_stops[index + 1]
        key = (left.uid, right.uid)
        default = self._overview_default_speed()
        raw = str(text).strip()
        if not raw:
            self.overview_transit_speeds.pop(key, None)
            self.overview_transit_durations.pop(key, None)
        else:
            try:
                speed = float(raw)
            except ValueError:
                speed = None
            if speed is not None:
                self.overview_transit_durations.pop(key, None)
                speed = min(max(speed, 0.1), 40.0)
                if abs(speed - default) < 0.05:
                    self.overview_transit_speeds.pop(key, None)
                else:
                    self.overview_transit_speeds[key] = speed
        self._fill_overview_list()
        self._refresh_overview(fit=False)

    def _overview_clock_fallback(self):
        current = self._overview_start() if hasattr(self, "overview_start_check") else None
        if current is not None:
            return current
        return datetime.now().replace(second=0, microsecond=0)

    def _overview_build_schedule(self):
        start = self._overview_start() if hasattr(self, "overview_start_check") else None
        return build_schedule(
            getattr(self, "overview_stops", []),
            self._overview_default_speed(),
            start,
            getattr(self, "overview_transit_speeds", None),
            getattr(self, "overview_transit_leads", None),
            getattr(self, "overview_transit_durations", None),
        )

    def _overview_apply_clock(self, kind, index, which, text):
        column = _OVERVIEW_COL_START if which == "start" else _OVERVIEW_COL_END
        row = self._overview_boundary(kind, index)
        current = None
        if row is not None:
            current = row.start if which == "start" else row.end
        parsed = parse_clock(text, current)
        if parsed is None:
            self._fill_overview_list(column=column)
            return
        new_start, changed = apply_activity_clock(
            self.overview_stops,
            kind,
            index,
            which,
            parsed,
            self._overview_start(),
            self._overview_default_speed(),
            self.overview_transit_speeds,
            self.overview_transit_leads,
            self.overview_transit_durations,
        )
        if changed and new_start is not None:
            self._overview_set_start_widget(new_start)
        self._fill_overview_list(column=column)
        self._refresh_overview(fit=False)

    def _overview_boundary(self, kind, index):
        rows = self._overview_build_schedule()
        for row, ref in enumerate(activity_rows(self.overview_stops)):
            if ref == (kind, index) and row < len(rows):
                return rows[row]
        return None

    def _overview_set_start_widget(self, moment):
        self.overview_start_check.blockSignals(True)
        self.overview_start_edit.blockSignals(True)
        self.overview_start_check.setChecked(True)
        self.overview_start_edit.setDateTime(_as_qdatetime(moment))
        self.overview_start_edit.setEnabled(True)
        self.overview_start_check.blockSignals(False)
        self.overview_start_edit.blockSignals(False)

    def _overview_stop_cells(self, index, stop, clock):
        start_text, end_text = self._overview_clock_texts(clock)
        return (
            (_OVERVIEW_COL_ORDER, str(index + 1), False, False),
            (_OVERVIEW_COL_TYPE, stop.test_type or "Point", False, False),
            (_OVERVIEW_COL_NAME, stop.name, stop.kind == "point", False),
            (_OVERVIEW_COL_TIME, f"{stop.duration_min:.1f}", False, True),
            (_OVERVIEW_COL_EXTRA, f"{stop.extra_min:.1f}", True, True),
            (_OVERVIEW_COL_SPEED, "", False, True),
            (_OVERVIEW_COL_TOTAL, f"{stop.occupied_minutes():.1f}", False, True),
            (_OVERVIEW_COL_START, start_text, True, False),
            (_OVERVIEW_COL_END, end_text, True, False),
            (_OVERVIEW_COL_NOTE, self._overview_note_text(stop), False, False),
        )

    def _overview_transit_cells(self, index, clock):
        left = self.overview_stops[index]
        right = self.overview_stops[index + 1]
        leg = transit_leg(left, right, self._overview_default_speed(), self.overview_transit_speeds)
        minutes = clock.duration_min if clock is not None else leg["minutes"]
        if left.endpoints()[1] is None or right.endpoints()[0] is None:
            note = "no positions"
        else:
            note = f"{leg['distance_m'] / 1852.0:.1f} nm"
            note += "; custom speed" if leg["custom"] else "; default speed"
        if (left.uid, right.uid) in self.overview_transit_durations:
            note += "; time adjusted"
        start_text, end_text = self._overview_clock_texts(clock)
        return (
            (_OVERVIEW_COL_ORDER, "", False, False),
            (_OVERVIEW_COL_TYPE, "Transit", False, False),
            (_OVERVIEW_COL_NAME, f"Transit to {right.name}", False, False),
            (_OVERVIEW_COL_TIME, f"{minutes:.1f}", False, True),
            (_OVERVIEW_COL_EXTRA, "", False, True),
            (_OVERVIEW_COL_SPEED, f"{leg['speed_kts']:.1f}", True, True),
            (_OVERVIEW_COL_TOTAL, f"{minutes:.1f}", False, True),
            (_OVERVIEW_COL_START, start_text, True, False),
            (_OVERVIEW_COL_END, end_text, True, False),
            (_OVERVIEW_COL_NOTE, note, False, False),
        )

    def _overview_clock_texts(self, clock):
        if clock is None or clock.start is None or clock.end is None:
            return "", ""
        return format_clock(clock.start), format_clock(clock.end)

    def _fill_overview_list(self, select=None, column=None):
        if select is None:
            kind, index = self._overview_current_ref()
        else:
            kind, index = "stop", select
        rows = activity_rows(self.overview_stops)
        schedule = self._overview_build_schedule()
        chosen = column
        self._overview_filling = True
        try:
            self._overview_row_refs = rows
            self.overview_list.setRowCount(len(rows))
            for row, (row_kind, stop_index) in enumerate(rows):
                clock = schedule[row] if row < len(schedule) else None
                if row_kind == "stop":
                    cells = self._overview_stop_cells(stop_index, self.overview_stops[stop_index], clock)
                else:
                    cells = self._overview_transit_cells(stop_index, clock)
                for cell_column, text, editable, align_right in cells:
                    self.overview_list.setItem(
                        row,
                        cell_column,
                        self._overview_text_item(text, editable=editable, align_right=align_right),
                    )
            table_row = self._overview_table_row_for(kind, index)
            if table_row < 0 and rows:
                table_row = 0
            if table_row >= 0:
                if chosen is None:
                    chosen = _OVERVIEW_COL_SPEED if kind == "transit" else _OVERVIEW_COL_EXTRA
                self.overview_list.setCurrentCell(table_row, chosen)
        finally:
            self._overview_filling = False

    def _fill_overview_schedule(self):
        start = self._overview_start() if hasattr(self, "overview_start_check") else None
        rows = self._overview_build_schedule()
        self.overview_schedule.setPlainText(format_schedule(rows, start))

    def _overview_tab_is_current(self):
        notebook = getattr(self, "param_notebook", None)
        index = getattr(self, "overview_tab_index", None)
        return notebook is not None and index is not None and notebook.currentIndex() == index

    def _refresh_overview(self, fit=False):
        self._fill_overview_schedule()
        if self._overview_tab_is_current() or getattr(self, "_overview_map_active", False):
            self._show_overview_on_map(fit=fit)

    def _on_overview_tab_changed(self, current_tab):
        index = getattr(self, "overview_tab_index", None)
        if index is None:
            return
        active = current_tab == index
        was = getattr(self, "_overview_map_active", False)
        if active:
            self._overview_set_bottom_strip_visible(False)
            self._show_overview_on_map(fit=not was)
            return
        self._overview_set_bottom_strip_visible(True)
        self._overview_cancel_end_pick(redraw=False)
        if was:
            self._overview_map_active = False
            if hasattr(self, "_sync_maplibre_plan"):
                self._sync_maplibre_plan(fit=True, terrain=True)

    def _overview_set_bottom_strip_visible(self, visible):
        """Hide the profile band so the map fills the right side while Overview is open."""
        strip = getattr(self, "bottom_strip_widget", None)
        if strip is None:
            return
        strip.setVisible(bool(visible))

    def _zoom_overview_home(self):
        """Fit the map to every stop and transit in the overview."""
        view = getattr(self, "map_view", None)
        if view is None or not hasattr(view, "fit_bounds"):
            return
        collection = map_features(getattr(self, "overview_stops", []) or [])
        bounds = feature_bounds(collection)
        if bounds is None:
            return
        view.fit_bounds(bounds[0], bounds[1], bounds[2], bounds[3])

    def _show_overview_on_map(self, fit=False, legend=True):
        """Draw the sequence on the imagery map. The survey grid stays put underneath."""
        self._overview_set_bottom_strip_visible(False)
        view = getattr(self, "map_view", None)
        if view is None or not hasattr(view, "show_plan"):
            self._overview_map_active = True
            return
        self._overview_map_active = True
        selected = self._overview_selected_stop_index() if hasattr(self, "overview_list") else -1
        if selected < 0:
            kind, index = self._overview_current_ref()
            selected = index if kind == "transit" and index >= 0 else None
        collection = map_features(self.overview_stops, selected)
        bounds = feature_bounds(collection) if fit else None
        if hasattr(view, "set_web_layers"):
            view.set_web_layers(True, False, 0.0)
        if hasattr(view, "show_bathymetry"):
            view.show_bathymetry("", None, None, None, None)
        if hasattr(view, "set_colorbar"):
            view.set_colorbar("", False)
        if legend:
            self._sync_overview_legend()
        preview = getattr(self, "_overview_pick_preview", None)
        if preview is not None:
            collection["features"].append({
                "type": "Feature",
                "geometry": {"type": "Point", "coordinates": [preview[1], preview[0]]},
                "properties": {"color": "#111111", "marker": True, "pick": True, "label": "Pick"},
            })
        # Draw after the imagery change. An earlier plan would be replaced by
        # the empty picture that clearing the depth grid can leave behind.
        view.show_plan(collection, bounds)
        if fit and bounds is not None and hasattr(view, "fit_bounds"):
            view.fit_bounds(bounds[0], bounds[1], bounds[2], bounds[3])

    def _sync_overview_legend(self):
        """Show each stop, in its map color, when Map Options has Legend on."""
        import base64

        from matplotlib.lines import Line2D

        view = getattr(self, "map_view", None)
        if view is None or not hasattr(view, "set_legend"):
            return
        items = legend_items(getattr(self, "overview_stops", []) or [])
        if not items or not hasattr(self, "_render_legend_png"):
            view.set_legend("", False)
            return
        handles = []
        labels = []
        for label, color, dashed in items:
            handles.append(Line2D(
                [0], [0],
                color=color,
                linewidth=2.2,
                linestyle=(0, (4, 2)) if dashed else "solid",
            ))
            labels.append(label)
        png = self._render_legend_png(handles, labels)
        self._map_legend_png_bytes = png
        if not png:
            view.set_legend("", False)
            return
        url = "data:image/png;base64," + base64.b64encode(png).decode("ascii")
        view.set_legend(url, bool(getattr(self, "show_map_legend_var", False)))

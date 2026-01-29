# Refactoring Plan: SAT_Planner_PyQt.py into Smaller Modules

## Would splitting help?

**Yes.** A single 10,000+ line file with one giant class is hard to:

- **Navigate** – finding the right method or section takes a long time
- **Edit safely** – changes in one area can accidentally affect another
- **Test** – individual features are difficult to unit test
- **Review** – diffs and code reviews are overwhelming
- **Share work** – multiple people can’t easily work on different features in parallel

Breaking the app into logical modules makes each of these easier.

---

## Current structure (summary)

- **One class:** `SurveyPlanApp(QMainWindow)` (~10,100 lines)
- **117+ methods** covering:
  - GeoTIFF load/display, dynamic resolution, contours
  - Reference planning (survey lines, crossline, export/import)
  - Calibration (pitch/roll lines, heading, export/import)
  - Line planning (draw line, profile, statistics)
  - Map interaction (click, scroll, pan, zoom)
  - Profile plots (crossline, pitch, line planning)
  - UI setup (`_create_widgets`), dialogs, activity log
  - Utilities (messages, coordinate conversion, config paths)

---

## Recommended approach: package + mixins

Keep `SurveyPlanApp` as the main window, but **split its behavior into mixins** and **shared code into small modules**. The app class then inherits from the mixins and uses the shared modules. No need to pass `self` everywhere or rewrite the app as many separate classes.

### 1. Package layout

```
MBES_Test_Planner/
├── SAT_Planner_PyQt.py      # optional: thin launcher that imports from sat_planner
├── sat_planner/
│   ├── __init__.py          # __version__, GEOSPATIAL_LIBS_AVAILABLE, SurveyPlanApp
│   ├── main.py              # if __name__ == "__main__": QApplication, window.show()
│   ├── constants.py         # CONFIG_FILENAME, version, geospatial import flag
│   ├── utils_ui.py          # _show_message, _ask_yes_no, _ask_ok_cancel (need QApplication)
│   ├── utils_geo.py         # decimal_degrees_to_ddm, depth/slope helpers (pure)
│   ├── app_core.py          # SurveyPlanApp __init__, _setup_layout, _create_widgets (or import)
│   ├── mixins/
│   │   ├── __init__.py
│   │   ├── geotiff_mixin.py      # Load/remove GeoTIFF, resolution, display mode, contours
│   │   ├── plotting_mixin.py     # _plot_survey_plan, _clear_plot, colorbars
│   │   ├── reference_mixin.py    # Reference tab, survey lines, export/import reference
│   │   ├── calibration_mixin.py  # Pitch/roll lines, calibration export/import, cal stats
│   │   ├── line_planning_mixin.py# Line planning tab, draw/edit, profile, line stats
│   │   ├── profiles_mixin.py      # Crossline/pitch/line planning profile drawing
│   │   ├── map_interaction_mixin.py # Click, scroll, pan, zoom, pick center/pitch/roll
│   │   ├── export_import_mixin.py   # Save/load params, export/import survey files
│   │   └── config_mixin.py        # Last-used dirs, config load/save
│   └── widgets/
│       ├── __init__.py
│       └── tab_reference.py   # optional: build reference tab given parent app
│       └── tab_calibration.py
│       └── tab_line_planning.py
```

- **constants.py** – single source of truth for version and config path.
- **utils_ui.py** – small helpers that need Qt (message boxes, etc.); can stay as methods on the app or become functions that take `parent` if you prefer.
- **utils_geo.py** – pure coordinate/depth/slope helpers; easy to unit test.
- **mixins/** – each file holds one logical group of methods. `SurveyPlanApp` inherits from all of them (e.g. `class SurveyPlanApp(GeoTIFFMixin, PlottingMixin, ...):`). No change to how methods call each other or use `self`.
- **widgets/** – optional; you can later move tab-building code into functions or small classes that take `self` (the app) and attach widgets to it.

### 2. Main app class (after refactor)

```python
# sat_planner/app_core.py
from PyQt6.QtWidgets import QMainWindow, ...
from .constants import CONFIG_FILENAME, __version__
from .mixins.geotiff_mixin import GeoTIFFMixin
from .mixins.plotting_mixin import PlottingMixin
from .mixins.reference_mixin import ReferenceMixin
# ... other mixins

class SurveyPlanApp(GeoTIFFMixin, PlottingMixin, ReferenceMixin, CalibrationMixin,
                    LinePlanningMixin, ProfilesMixin, MapInteractionMixin,
                    ExportImportMixin, ConfigMixin, QMainWindow):
    def __init__(self):
        super().__init__()
        # Only init state and layout here; mixins provide the methods.
        ...
```

### 3. What to move where (method → module)

| Module | Methods (examples) |
|--------|---------------------|
| **constants** | CONFIG_FILENAME, __version__, GEOSPATIAL_LIBS_AVAILABLE |
| **utils_ui** | _show_message, _ask_yes_no, _ask_ok_cancel |
| **utils_geo** | _decimal_degrees_to_ddm, _get_depth_at_point, _calculate_slope_at_point |
| **geotiff_mixin** | _load_geotiff, _remove_geotiff, _load_geotiff_at_resolution, _on_geotiff_display_mode_changed, _on_contour_*, _toggle_dynamic_resolution |
| **plotting_mixin** | _plot_survey_plan, _clear_plot, _remove_colorbar, _generate_and_plot |
| **reference_mixin** | _export_survey_data, _import_survey_files (ref), _validate_inputs (ref), set_ref_info_text, _reset_reference_tab, _show_reference_planning_info |
| **calibration_mixin** | _toggle_pick_pitch_line_mode, _toggle_edit_pitch_line_mode, _toggle_pick_roll_line_mode, pitch/roll handle events, _export_cal_survey_files, _import_cal_survey_files, _show_calibration_statistics |
| **line_planning_mixin** | _toggle_line_planning_mode, _clear_line_planning, line planning click/handle events, _draw_line_planning_profile, _show_line_information |
| **profiles_mixin** | _get_profile_data_from_geotiff, _draw_crossline_profile, _draw_pitch_line_profile, _draw_current_profile |
| **map_interaction_mixin** | _on_plot_click, _on_scroll, _on_middle_*, _zoom_to_*, _toggle_pick_center_mode |
| **export_import_mixin** | _save_survey_parameters, _load_survey_parameters, _export_survey_files, etc. |
| **config_mixin** | _load_last_used_dir, _save_last_used_dir, _load_last_geotiff_dir, ... |

You can refine this table as you split (e.g. move a method to a different mixin if it fits better).

### 4. Migration strategy (low risk)

1. **Create the package** – Add `sat_planner/`, `constants.py`, `utils_geo.py`, and empty mixin files.
2. **Extract constants and pure utils** – Move version, config path, and geospatial flag to `constants.py`. Move coordinate/depth/slope helpers to `utils_geo.py` and have the main file import them. Run the app and tests to confirm nothing breaks.
3. **Extract one mixin** – Pick the smallest coherent block (e.g. GeoTIFF or Config “last used dirs”). Move those methods into a mixin, add the mixin as a base class of `SurveyPlanApp`, and remove those methods from the monolith. Run again.
4. **Repeat** – Extract one mixin at a time, run after each step. Fix any missing `self` or cross-mixin calls.
5. **Optional** – Move `_create_widgets` into smaller functions in `widgets/` (e.g. one function per tab) that take `self` and create the widgets. The main app still calls them from `_create_widgets`.

### 5. Alternative: composition instead of mixins

Instead of mixins, you could have:

- `SurveyPlanApp` hold instances like `self.geotiff_manager`, `self.reference_manager`, etc.
- Each manager in its own file and responsible for one feature.

That’s cleaner OOP but requires more refactoring: the current code often does `self._load_geotiff()` and uses `self.figure`, `self.ax`, etc. With composition you’d either pass the app (or figure/axes) into the manager or expose callbacks. Mixins get you to “smaller files” with minimal change to how the rest of the code runs.

---

## Summary

- **Yes, breaking the program into smaller pieces will make editing (and navigating, testing, and reviewing) easier.**
- **Recommended first step:** Introduce the `sat_planner` package, move constants and pure geo utils, then extract one mixin (e.g. config or GeoTIFF) and wire it up. After that, repeat for the other logical blocks using the table above.

### 6. Quick start (package created)

A starter package is in place:

- **`sat_planner/constants.py`** – `__version__`, `CONFIG_FILENAME`, `GEOSPATIAL_LIBS_AVAILABLE`
- **`sat_planner/utils_geo.py`** – `decimal_degrees_to_ddm(decimal_deg, is_latitude=True)` (pure function)
- **`sat_planner/__init__.py`** – re-exports the above
- **`sat_planner/mixins/`** – placeholder for future mixins

To use it from `SAT_Planner_PyQt.py` without moving the whole app:

1. At the top, add:  
   `from sat_planner import __version__, CONFIG_FILENAME, GEOSPATIAL_LIBS_AVAILABLE, decimal_degrees_to_ddm`
2. Remove the local `__version__`, `CONFIG_FILENAME`, and the geospatial try/except block (they’re now in `sat_planner.constants`).
3. Replace every `self._decimal_degrees_to_ddm(...)` with `decimal_degrees_to_ddm(...)` (and keep the same two arguments).

The app can still be run as `python SAT_Planner_PyQt.py`. After that, extract one mixin at a time following the table in section 3.

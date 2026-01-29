# Refactoring Plan: SAT_Planner_PyQt.py into Smaller Modules

## Progress summary

| Phase | Status |
|-------|--------|
| Package + constants + utils | **Done** – `sat_planner/`, `constants.py`, `utils_geo.py`, `utils_ui.py` |
| All mixins extracted | **Done** – Basemap, GeoTIFF, Plotting, Reference, Calibration, Line Planning, Profiles, Map Interaction, Export/Import, Config |
| GeoTIFF/plotting helpers moved | **Done** – resolution, zoom/pan, _reset_to_consistent_view, _calculate_consistent_plot_limits, _calculate_slope_at_point, _get_depth_at_point in mixins |
| Basemap/NOAA web services | **Done** – `BasemapMixin` (_toggle_imagery_basemap, _load_and_plot_basemap, _load_and_plot_noaa_charts, etc.) |
| Survey plan DDM axis labels | **Done** – degrees–decimal minutes in `plotting_mixin` via `utils_geo.decimal_degrees_to_ddm` |
| Move app class to package | **Remaining** – put `SurveyPlanApp` in `sat_planner/app_core.py`, launcher in `SAT_Planner_PyQt.py` |
| Move remaining methods from main | **Remaining** – ~22 methods still in main file (see table below) |
| Optional: `main.py`, `widgets/` | **Optional** – `python -m sat_planner`; split `_create_widgets` by tab |

---

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

- **Entry point:** `SAT_Planner_PyQt.py` – creates `SurveyPlanApp`, runs app (~1,800 lines; down from ~10,100).
- **Package:** `sat_planner/` – constants, utils_geo, utils_ui, and 10 mixins (Basemap, GeoTIFF, Plotting, Reference, Calibration, Line Planning, Profiles, Map Interaction, Export/Import, Config).
- **Still in main file:** ~22 methods – `__init__`, `_setup_layout`, `_create_widgets`, message wrappers, parameter/tab handlers, calibration helpers, colormap/temp-line events, `_on_mouse_motion`, `_show_statistics_dialog`, etc.

---

## Recommended approach: package + mixins

Keep `SurveyPlanApp` as the main window, but **split its behavior into mixins** and **shared code into small modules**. The app class then inherits from the mixins and uses the shared modules. No need to pass `self` everywhere or rewrite the app as many separate classes.

### 1. Package layout

**In place:**  
`SAT_Planner_PyQt.py` (launcher + `SurveyPlanApp`), `sat_planner/` with:

```
sat_planner/
├── __init__.py              # __version__, CONFIG_FILENAME, GEOSPATIAL_LIBS_AVAILABLE, decimal_degrees_to_ddm, show_message, ask_yes_no, ask_ok_cancel
├── constants.py             # CONFIG_FILENAME, __version__, geospatial imports (rasterio, pyproj, etc.)
├── utils_ui.py               # show_message, ask_yes_no, ask_ok_cancel (take parent)
├── utils_geo.py             # decimal_degrees_to_ddm, decimal_degrees_to_dms_string (pure)
├── mixins/
│   ├── __init__.py
│   ├── basemap_mixin.py     # Imagery basemap, NOAA ENC Charts (toggles, load, opacity)
│   ├── geotiff_mixin.py     # Load/remove GeoTIFF, resolution, display mode, contours, slope/depth at point, reset view
│   ├── plotting_mixin.py    # _plot_survey_plan, _clear_plot, colorbars, _calculate_consistent_plot_limits, DDM axis labels
│   ├── reference_mixin.py   # Reference tab, survey lines, export/import reference
│   ├── calibration_mixin.py # Pitch/roll lines, calibration export/import, cal stats
│   ├── line_planning_mixin.py # Line planning tab, draw/edit, profile, line stats
│   ├── profiles_mixin.py    # Crossline/pitch/line planning profile drawing
│   ├── map_interaction_mixin.py # Click, scroll, pan, zoom, pick center/pitch/roll, basemap/NOAA reload timers
│   ├── export_import_mixin.py   # Save/load params, export/import survey files
│   └── config_mixin.py      # Last-used dirs, config load/save
```

**Remaining / optional:**

- **app_core.py** (remaining) – Move `SurveyPlanApp` class here (with `__init__`, `_setup_layout`, `_create_widgets`, `_on_tab_changed`, `_quit_app`, `_show_about_dialog`). Launcher stays in `SAT_Planner_PyQt.py` and imports `SurveyPlanApp` from `sat_planner`.
- **main.py** (optional) – `if __name__ == "__main__":` so you can run `python -m sat_planner`.
- **widgets/** (optional) – `tab_reference.py`, `tab_calibration.py`, `tab_line_planning.py`; each builds one tab given `self` (app). `_create_widgets` in app_core just calls them.

- **constants.py** – single source of truth for version and config path.
- **utils_ui.py** – Qt helpers (message boxes, etc.); used as functions with `parent`.
- **utils_geo.py** – pure coordinate helpers (DDM, DMS); easy to unit test.
- **mixins/** – `SurveyPlanApp` inherits from all 10 mixins; methods use `self` as before.

### 2. Main app class (current and target)

**Current:** `SurveyPlanApp` is defined in `SAT_Planner_PyQt.py` and inherits from all 10 mixins:

```python
# SAT_Planner_PyQt.py (current)
class SurveyPlanApp(BasemapMixin, GeoTIFFMixin, PlottingMixin, ReferenceMixin, CalibrationMixin,
                    LinePlanningMixin, ProfilesMixin, MapInteractionMixin,
                    ExportImportMixin, ConfigMixin, QMainWindow):
    ...
```

**Target (remaining):** Move the class into `sat_planner/app_core.py` so `SAT_Planner_PyQt.py` becomes a thin launcher that imports `SurveyPlanApp` and runs the app. Same inheritance list; `__init__`, `_setup_layout`, `_create_widgets`, and any remaining “core” methods live in app_core.

### 3. What to move where (method → module)

**Done:** constants, utils_ui, utils_geo; all 10 mixins exist and contain the methods below (plus BasemapMixin for basemap/NOAA). GeoTIFF and Plotting mixins also have: _load_geotiff_at_resolution, _reload_geotiff_at_current_zoom, _clear_rapid_zoom_mode, _clear_panning_mode, _toggle_dynamic_resolution, _reset_to_consistent_view, _calculate_consistent_plot_limits, _calculate_slope_at_point, _get_depth_at_point. Plotting has DDM axis labels via utils_geo.decimal_degrees_to_ddm.

**Remaining (still in SAT_Planner_PyQt.py):**

| Still in main file | Target |
|--------------------|--------|
| _show_message, _ask_yes_no, _ask_ok_cancel | Remove wrappers; call utils_ui show_message, ask_yes_no, ask_ok_cancel directly where used. |
| _show_about_dialog | Keep in app (app_core). |
| _update_multiplier_label_len, _update_multiplier_label_dist | reference_mixin |
| _setup_layout | app_core |
| _on_parameter_changed, _auto_regenerate_survey_plan | reference_mixin |
| _quit_app | app_core |
| _on_line_length_or_speed_change, _update_export_name | reference_mixin |
| _update_cal_line_offset_from_pitch_line, _add_heading_lines_from_pitch_line, _update_cal_export_name_from_pitch_line | calibration_mixin |
| _on_draw_event_update_colormap, _on_temp_line_motion | geotiff_mixin or plotting_mixin |
| _on_mouse_motion | map_interaction_mixin |
| _on_tab_changed | app_core |
| _create_widgets | app_core (or later split into widgets/) |
| _show_statistics_dialog | utils_ui (function with parent) or keep in app |

**Original mapping (for reference):**

| Module | Methods (examples) |
|--------|---------------------|
| **constants** | CONFIG_FILENAME, __version__, GEOSPATIAL_LIBS_AVAILABLE ✓ |
| **utils_ui** | show_message, ask_yes_no, ask_ok_cancel ✓ |
| **utils_geo** | decimal_degrees_to_ddm, (get_depth/slope in geotiff_mixin) ✓ |
| **basemap_mixin** | _toggle_imagery_basemap, _toggle_noaa_charts, _load_and_plot_basemap, _load_and_plot_noaa_charts, _update_noaa_charts_opacity ✓ |
| **geotiff_mixin** | _load_geotiff, _remove_geotiff, _load_geotiff_at_resolution, _toggle_dynamic_resolution, _reset_to_consistent_view, _calculate_slope_at_point, _get_depth_at_point, … ✓ |
| **plotting_mixin** | _plot_survey_plan, _clear_plot, _calculate_consistent_plot_limits, DDM axis labels, … ✓ |
| **reference_mixin** | _export_survey_data, _validate_inputs (ref), set_ref_info_text, … (+ remaining: _on_parameter_changed, _auto_regenerate_survey_plan, multiplier labels, _update_export_name, _on_line_length_or_speed_change) |
| **calibration_mixin** | pitch/roll/heading, export/import cal (+ remaining: _update_cal_line_offset_from_pitch_line, _add_heading_lines_from_pitch_line, _update_cal_export_name_from_pitch_line) |
| **line_planning_mixin** | _toggle_line_planning_mode, _draw_line_planning_profile, … ✓ |
| **profiles_mixin** | _draw_crossline_profile, _draw_pitch_line_profile, _draw_current_profile ✓ |
| **map_interaction_mixin** | _on_plot_click, _on_scroll, pan, zoom, pick center (+ remaining: _on_mouse_motion) |
| **export_import_mixin** | _save_survey_parameters, _load_survey_parameters, _export_survey_files, … ✓ |
| **config_mixin** | _load_last_used_dir, _save_last_used_dir, … ✓ |

### 4. Migration strategy (low risk)

**Done:**

1. ~~Create the package~~ – `sat_planner/`, `constants.py`, `utils_geo.py`, `utils_ui.py`, and all 10 mixin files.
2. ~~Extract constants and pure utils~~ – Version, config path, geospatial imports in `constants.py`; DDM/DMS in `utils_geo.py`; message/confirm in `utils_ui.py`.
3. ~~Extract mixins~~ – All 10 mixins extracted (Basemap, GeoTIFF, Plotting, Reference, Calibration, Line Planning, Profiles, Map Interaction, Export/Import, Config). GeoTIFF/plotting helpers and basemap/NOAA moved into mixins. Survey plan axes use DDM labels.

**Next (recommended order):**

4. **Move remaining methods from main file** – Move the ~22 methods still in `SAT_Planner_PyQt.py` into the target mixins or app_core (see table in §3). Remove message wrappers and call `utils_ui` directly. Run after each move.
5. **Move app class to package** – Create `sat_planner/app_core.py` with `SurveyPlanApp` (and `__init__`, `_setup_layout`, `_create_widgets`, `_on_tab_changed`, `_quit_app`, `_show_about_dialog`). Have `SAT_Planner_PyQt.py` import `SurveyPlanApp` from `sat_planner` and run the app. Main file becomes a thin launcher.
6. **Optional** – Add `sat_planner/main.py` for `python -m sat_planner`. Optionally split `_create_widgets` into `widgets/tab_*.py` (one function per tab).

### 5. Alternative: composition instead of mixins

Instead of mixins, you could have:

- `SurveyPlanApp` hold instances like `self.geotiff_manager`, `self.reference_manager`, etc.
- Each manager in its own file and responsible for one feature.

That’s cleaner OOP but requires more refactoring: the current code often does `self._load_geotiff()` and uses `self.figure`, `self.ax`, etc. With composition you’d either pass the app (or figure/axes) into the manager or expose callbacks. Mixins get you to “smaller files” with minimal change to how the rest of the code runs.

---

## Summary

- **Yes, breaking the program into smaller pieces will make editing (and navigating, testing, and reviewing) easier.**
- **Done:** Package, constants, utils, all 10 mixins (including Basemap), geotiff/plotting helpers, DDM axis labels. Main file reduced from ~10,100 to ~1,800 lines.
- **Next:** Move the ~22 remaining methods from `SAT_Planner_PyQt.py` into the target mixins (see §3). Then move `SurveyPlanApp` into `sat_planner/app_core.py` so the main file is a thin launcher. Optional: `main.py`, `widgets/`.

# Refactoring Plan: SAT_Planner_PyQt.py into Smaller Modules

*Last updated: 2026-01-29*

## Progress summary

| Phase | Status |
|-------|--------|
| Package + constants + utils | **Done** – `sat_planner/`, `constants.py`, `utils_geo.py`, `utils_ui.py` |
| All mixins extracted | **Done** – Basemap, GeoTIFF, Plotting, Reference, Calibration, Line Planning, Profiles, Map Interaction, Export/Import, Config |
| GeoTIFF/plotting helpers moved | **Done** – resolution, zoom/pan, _reset_to_consistent_view, _calculate_consistent_plot_limits, _calculate_slope_at_point, _get_depth_at_point in mixins |
| Basemap/NOAA web services | **Done** – `BasemapMixin` (_toggle_imagery_basemap, _load_and_plot_basemap, _load_and_plot_noaa_charts, etc.) |
| Survey plan DDM axis labels | **Done** – degrees–decimal minutes in `plotting_mixin` via `utils_geo.decimal_degrees_to_ddm` |
| Move remaining methods from main | **Done** – multiplier/param/export → reference_mixin; cal pitch/heading/export → calibration_mixin; colormap/temp_line → geotiff_mixin; _on_mouse_motion → map_interaction_mixin; show_statistics_dialog → utils_ui |
| Move app class to package | **Done** – `SurveyPlanApp` in `sat_planner/app_core.py`; `SAT_Planner_PyQt.py` is thin launcher |
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

- **Entry point:** `SAT_Planner_PyQt.py` – thin launcher: imports `SurveyPlanApp` from `sat_planner`, runs app (~25 lines).
- **Package:** `sat_planner/` – `app_core.py` (SurveyPlanApp + 9 core methods), `constants.py`, `utils_geo.py`, `utils_ui.py`, and 10 mixins in `mixins/`.
- **App class:** `sat_planner/app_core.py` – `SurveyPlanApp` with `__init__`, `_show_message`, `_ask_yes_no`, `_ask_ok_cancel`, `_show_about_dialog`, `_setup_layout`, `_quit_app`, `_on_tab_changed`, `_create_widgets`. `_create_widgets` contains nested callbacks that are candidates for the optional `widgets/` split.

---

## Recommended approach: package + mixins

Keep `SurveyPlanApp` as the main window, but **split its behavior into mixins** and **shared code into small modules**. The app class then inherits from the mixins and uses the shared modules. No need to pass `self` everywhere or rewrite the app as many separate classes.

### 1. Package layout

**In place:**  
`SAT_Planner_PyQt.py` (thin launcher only), `sat_planner/` with:

```
sat_planner/
├── __init__.py              # __version__, CONFIG_FILENAME, GEOSPATIAL_LIBS_AVAILABLE, decimal_degrees_to_ddm, show_message, ask_yes_no, ask_ok_cancel, SurveyPlanApp
├── app_core.py              # SurveyPlanApp (__init__, _setup_layout, _create_widgets, _on_tab_changed, _quit_app, _show_about_dialog, message wrappers)
├── constants.py             # CONFIG_FILENAME, __version__, geospatial imports (rasterio, pyproj, etc.)
├── utils_ui.py              # show_message, ask_yes_no, ask_ok_cancel, show_statistics_dialog (take parent)
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

**Optional:**

- **main.py** (optional) – `if __name__ == "__main__":` so you can run `python -m sat_planner`.
- **widgets/** (optional) – `tab_reference.py`, `tab_calibration.py`, `tab_line_planning.py`; each builds one tab given `self` (app). `_create_widgets` in app_core just calls them.

- **constants.py** – single source of truth for version and config path.
- **utils_ui.py** – Qt helpers (message boxes, etc.); used as functions with `parent`.
- **utils_geo.py** – pure coordinate helpers (DDM, DMS); easy to unit test.
- **mixins/** – `SurveyPlanApp` inherits from all 10 mixins; methods use `self` as before.

### 2. Main app class (current and target)

**Current:** `SurveyPlanApp` is defined in `sat_planner/app_core.py` and inherits from all 10 mixins. `SAT_Planner_PyQt.py` is a thin launcher that imports `SurveyPlanApp` from `sat_planner` and runs the app.

```python
# SAT_Planner_PyQt.py (current)
class SurveyPlanApp(BasemapMixin, GeoTIFFMixin, PlottingMixin, ReferenceMixin, CalibrationMixin,
                    LinePlanningMixin, ProfilesMixin, MapInteractionMixin,
                    ExportImportMixin, ConfigMixin, QMainWindow):
    ...
```

**Target (remaining):** Move the class into `sat_planner/app_core.py` so `SAT_Planner_PyQt.py` becomes a thin launcher that imports `SurveyPlanApp` and runs the app. Same inheritance list; `__init__`, `_setup_layout`, `_create_widgets`, and any remaining “core” methods live in app_core.

### 3. What to move where (method → module)

**Done:** constants, utils_ui (including show_statistics_dialog), utils_geo; all 10 mixins exist. Reference_mixin has _update_multiplier_label_len/dist, _on_parameter_changed, _auto_regenerate_survey_plan, _on_line_length_or_speed_change, _update_export_name. Calibration_mixin has _update_cal_line_offset_from_pitch_line, _add_heading_lines_from_pitch_line, _update_cal_export_name_from_pitch_line. Geotiff_mixin has _on_draw_event_update_colormap, _on_temp_line_motion. Map_interaction_mixin has _on_mouse_motion. Mixins call show_statistics_dialog(self, title, text) from utils_ui.

**In app_core (done):** All core methods live in `sat_planner/app_core.py`: `__init__`, `_setup_layout`, `_create_widgets`, `_show_message`, `_ask_yes_no`, `_ask_ok_cancel`, `_show_about_dialog`, `_quit_app`, `_on_tab_changed`.

**Original mapping (for reference):**

| Module | Methods (examples) |
|--------|---------------------|
| **constants** | CONFIG_FILENAME, __version__, GEOSPATIAL_LIBS_AVAILABLE ✓ |
| **utils_ui** | show_message, ask_yes_no, ask_ok_cancel, show_statistics_dialog ✓ |
| **utils_geo** | decimal_degrees_to_ddm, (get_depth/slope in geotiff_mixin) ✓ |
| **basemap_mixin** | _toggle_imagery_basemap, _toggle_noaa_charts, _load_and_plot_basemap, _load_and_plot_noaa_charts, _update_noaa_charts_opacity ✓ |
| **geotiff_mixin** | _load_geotiff, _remove_geotiff, _load_geotiff_at_resolution, _toggle_dynamic_resolution, _reset_to_consistent_view, _calculate_slope_at_point, _get_depth_at_point, _on_draw_event_update_colormap, _on_temp_line_motion ✓ |
| **plotting_mixin** | _plot_survey_plan, _clear_plot, _calculate_consistent_plot_limits, DDM axis labels ✓ |
| **reference_mixin** | _export_survey_data, _validate_inputs (ref), set_ref_info_text, _update_multiplier_label_len/dist, _on_parameter_changed, _auto_regenerate_survey_plan, _on_line_length_or_speed_change, _update_export_name ✓ |
| **calibration_mixin** | pitch/roll/heading, export/import cal, _update_cal_line_offset_from_pitch_line, _add_heading_lines_from_pitch_line, _update_cal_export_name_from_pitch_line ✓ |
| **line_planning_mixin** | _toggle_line_planning_mode, _draw_line_planning_profile, … ✓ |
| **profiles_mixin** | _draw_crossline_profile, _draw_pitch_line_profile, _draw_current_profile ✓ |
| **map_interaction_mixin** | _on_plot_click, _on_scroll, pan, zoom, pick center, _on_mouse_motion ✓ |
| **export_import_mixin** | _save_survey_parameters, _load_survey_parameters, _export_survey_files, … ✓ |
| **config_mixin** | _load_last_used_dir, _save_last_used_dir, … ✓ |

### 4. Migration strategy (low risk)

**Done:**

1. ~~Create the package~~ – `sat_planner/`, `constants.py`, `utils_geo.py`, `utils_ui.py`, and all 10 mixin files.
2. ~~Extract constants and pure utils~~ – Version, config path, geospatial imports in `constants.py`; DDM/DMS in `utils_geo.py`; message/confirm in `utils_ui.py`.
3. ~~Extract mixins~~ – All 10 mixins extracted (Basemap, GeoTIFF, Plotting, Reference, Calibration, Line Planning, Profiles, Map Interaction, Export/Import, Config). GeoTIFF/plotting helpers and basemap/NOAA moved into mixins. Survey plan axes use DDM labels.
4. ~~Move remaining methods from main file~~ – Multiplier/param/export → reference_mixin; cal pitch/heading/export → calibration_mixin; colormap/temp_line → geotiff_mixin; _on_mouse_motion → map_interaction_mixin; show_statistics_dialog → utils_ui. Main file now ~1,290 lines with 9 core methods only.

**Next (recommended order):**

5. ~~**Move app class to package**~~ – **Done.** `sat_planner/app_core.py` contains `SurveyPlanApp`; `SAT_Planner_PyQt.py` imports from `sat_planner` and runs the app (thin launcher).
6. **Optional** – Add `sat_planner/main.py` for `python -m sat_planner`. Optionally split `_create_widgets` into `widgets/tab_*.py` (one builder per tab); move nested callbacks such as `_on_pitch_line_motion`, `_on_pick_center_motion`, `_on_geotiff_hover_motion` into those modules or keep as closures.

### 5. Alternative: composition instead of mixins

Instead of mixins, you could have:

- `SurveyPlanApp` hold instances like `self.geotiff_manager`, `self.reference_manager`, etc.
- Each manager in its own file and responsible for one feature.

That’s cleaner OOP but requires more refactoring: the current code often does `self._load_geotiff()` and uses `self.figure`, `self.ax`, etc. With composition you’d either pass the app (or figure/axes) into the manager or expose callbacks. Mixins get you to “smaller files” with minimal change to how the rest of the code runs.

---

## Summary

- **Yes, breaking the program into smaller pieces will make editing (and navigating, testing, and reviewing) easier.**
- **Done:** Package, constants, utils, all 10 mixins, geotiff/plotting helpers, DDM axis labels, and app class in `sat_planner/app_core.py`. Main file `SAT_Planner_PyQt.py` is a thin launcher (~25 lines).
- **Next (optional):** `sat_planner/main.py` for `python -m sat_planner`; `sat_planner/widgets/` to split `_create_widgets` and nested callbacks by tab.

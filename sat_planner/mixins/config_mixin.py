"""
Config: load/save last-used directories to/from JSON config file.
_load_last_used_dir, _save_last_used_dir, _load_last_geotiff_dir, _save_last_geotiff_dir,
_load_last_backscatter_dir, _save_last_backscatter_dir,
_load_last_backscatter_import_dir, _save_last_backscatter_import_dir,
_load_last_backscatter_export_dir, _save_last_backscatter_export_dir,
_load_last_survey_params_dir, _save_last_survey_params_dir, _load_last_export_dir, _save_last_export_dir,
_load_last_ref_import_dir, _save_last_ref_import_dir, _load_last_line_import_dir, _save_last_line_import_dir,
_load_last_perf_import_dir, _save_last_perf_import_dir, _load_last_adcp_import_dir, _save_last_adcp_import_dir,
_load_last_shapefile_dir, _save_last_shapefile_dir,
_load_export_type_options, _save_export_type_options,
_load_vert_exag_table, _save_vert_exag_table, _default_vert_exag_table,
_vert_exag_table_for_params, _add_vert_exag_table_to_params, _apply_vert_exag_table_from_params,
_load_shaded_relief_cmap, _save_shaded_relief_cmap, _normalize_shaded_relief_cmap,
_shaded_relief_cmap_label, _load_slope_overlay_bands, _save_slope_overlay_bands,
_add_geotiff_viz_params_to_params, _apply_geotiff_viz_params_from_params.
"""
import json
import os

from sat_planner.constants import (
    DEFAULT_SHADED_RELIEF_CMAP,
    DEFAULT_SLOPE_OVERLAY_BANDS,
    SHADED_RELIEF_CMAP_OPTIONS,
)


class ConfigMixin:
    """Mixin for loading/saving last-used directories from config file."""

    def _default_export_type_options(self):
        return {
            "esri_shapefile": True,
            "gpkg": False,
            "sis_asciiplan": True,
            "gpx": True,
            "text_csv": True,
            "text_txt": True,
            "hypack_lnw": True,
            "map_png_high": True,
            "map_png_low": True,
            "profiles_png_high": True,
            "profiles_png_low": True,
        }

    def _load_last_used_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_used_dir' in config and os.path.isdir(config['last_used_dir']):
                    self.last_used_dir = config['last_used_dir']
        except Exception:
            pass

    def _save_last_used_dir(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
            config['last_used_dir'] = self.last_used_dir
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _load_last_geotiff_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_geotiff_dir' in config and os.path.isdir(config['last_geotiff_dir']):
                    self.last_geotiff_dir = config['last_geotiff_dir']
        except Exception:
            pass

    def _save_last_geotiff_dir(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
            config['last_geotiff_dir'] = self.last_geotiff_dir
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _load_last_backscatter_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_backscatter_dir' in config and os.path.isdir(config['last_backscatter_dir']):
                    self.last_backscatter_dir = config['last_backscatter_dir']
        except Exception:
            pass

    def _save_last_backscatter_dir(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
            config['last_backscatter_dir'] = self.last_backscatter_dir
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _load_last_backscatter_import_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_backscatter_import_dir' in config and os.path.isdir(config['last_backscatter_import_dir']):
                    self.last_backscatter_import_dir = config['last_backscatter_import_dir']
        except Exception:
            pass

    def _save_last_backscatter_import_dir(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
            config['last_backscatter_import_dir'] = self.last_backscatter_import_dir
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _load_last_backscatter_export_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_backscatter_export_dir' in config and os.path.isdir(config['last_backscatter_export_dir']):
                    self.last_backscatter_export_dir = config['last_backscatter_export_dir']
        except Exception:
            pass

    def _save_last_backscatter_export_dir(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
            config['last_backscatter_export_dir'] = self.last_backscatter_export_dir
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _load_last_survey_params_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_survey_params_dir' in config and os.path.isdir(config['last_survey_params_dir']):
                    self.last_survey_params_dir = config['last_survey_params_dir']
        except Exception:
            pass

    def _save_last_survey_params_dir(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
            config['last_survey_params_dir'] = self.last_survey_params_dir
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _load_last_export_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_export_dir' in config and os.path.isdir(config['last_export_dir']):
                    self.last_export_dir = config['last_export_dir']
        except Exception:
            pass

    def _save_last_export_dir(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
            config['last_export_dir'] = self.last_export_dir
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _load_last_ref_import_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_ref_import_dir' in config and os.path.isdir(config['last_ref_import_dir']):
                    self.last_ref_import_dir = config['last_ref_import_dir']
        except Exception:
            pass

    def _save_last_ref_import_dir(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
            config['last_ref_import_dir'] = self.last_ref_import_dir
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _load_last_line_import_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_line_import_dir' in config and os.path.isdir(config['last_line_import_dir']):
                    self.last_line_import_dir = config['last_line_import_dir']
        except Exception:
            pass

    def _save_last_line_import_dir(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
            config['last_line_import_dir'] = self.last_line_import_dir
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _load_last_perf_import_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_perf_import_dir' in config and os.path.isdir(config['last_perf_import_dir']):
                    self.last_perf_import_dir = config['last_perf_import_dir']
        except Exception:
            pass

    def _save_last_perf_import_dir(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
            config['last_perf_import_dir'] = self.last_perf_import_dir
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _load_last_adcp_import_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_adcp_import_dir' in config and os.path.isdir(config['last_adcp_import_dir']):
                    self.last_adcp_import_dir = config['last_adcp_import_dir']
        except Exception:
            pass

    def _save_last_adcp_import_dir(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
            config['last_adcp_import_dir'] = self.last_adcp_import_dir
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _load_last_shapefile_dir(self):
        try:
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
                if 'last_shapefile_dir' in config and os.path.isdir(config['last_shapefile_dir']):
                    self.last_shapefile_dir = config['last_shapefile_dir']
        except Exception:
            pass

    def _save_last_shapefile_dir(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, 'r') as f:
                    config = json.load(f)
            config['last_shapefile_dir'] = self.last_shapefile_dir
            with open(self.CONFIG_FILENAME, 'w') as f:
                json.dump(config, f)
        except Exception:
            pass

    def _load_export_type_options(self):
        try:
            defaults = self._default_export_type_options()
            self.export_type_options = dict(defaults)
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, "r") as f:
                    config = json.load(f)
                saved = config.get("export_type_options", {})
                if isinstance(saved, dict):
                    migrated = dict(saved)
                    if "map_png" in migrated:
                        if "map_png_high" not in migrated:
                            migrated["map_png_high"] = bool(migrated["map_png"])
                        if "map_png_low" not in migrated:
                            migrated["map_png_low"] = bool(migrated["map_png"])
                    if "profiles_png" in migrated:
                        if "profiles_png_high" not in migrated:
                            migrated["profiles_png_high"] = bool(migrated["profiles_png"])
                        if "profiles_png_low" not in migrated:
                            migrated["profiles_png_low"] = bool(migrated["profiles_png"])
                    for key in defaults:
                        if key in migrated:
                            self.export_type_options[key] = bool(migrated[key])
        except Exception:
            self.export_type_options = dict(self._default_export_type_options())

    def _save_export_type_options(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, "r") as f:
                    config = json.load(f)
            defaults = self._default_export_type_options()
            current = getattr(self, "export_type_options", defaults) or defaults
            config["export_type_options"] = {
                key: bool(current.get(key, defaults[key])) for key in defaults
            }
            with open(self.CONFIG_FILENAME, "w") as f:
                json.dump(config, f)
        except Exception:
            pass

    def _default_vert_exag_table(self):
        return [
            {"elevation_range": 0.0, "shaded_relief": 1.5, "shaded_relief_dyn": 3.0},
            {"elevation_range": 20.0, "shaded_relief": 0.8, "shaded_relief_dyn": 1.5},
            {"elevation_range": 50.0, "shaded_relief": 0.4, "shaded_relief_dyn": 0.8},
            {"elevation_range": 200.0, "shaded_relief": 0.1, "shaded_relief_dyn": 0.4},
            {"elevation_range": 1000.0, "shaded_relief": 0.05, "shaded_relief_dyn": 0.05},
        ]

    def _normalize_vert_exag_table(self, table):
        defaults = self._default_vert_exag_table()
        if not isinstance(table, list) or len(table) != len(defaults):
            return [dict(row) for row in defaults]

        normalized = []
        for default_row, row in zip(defaults, table):
            if not isinstance(row, dict):
                normalized.append(dict(default_row))
                continue
            try:
                elevation_range = float(row.get("elevation_range", default_row["elevation_range"]))
                shaded_relief = float(row.get("shaded_relief", default_row["shaded_relief"]))
                shaded_relief_dyn = float(row.get("shaded_relief_dyn", default_row["shaded_relief_dyn"]))
            except (TypeError, ValueError):
                normalized.append(dict(default_row))
                continue
            normalized.append({
                "elevation_range": elevation_range,
                "shaded_relief": shaded_relief,
                "shaded_relief_dyn": shaded_relief_dyn,
            })

        normalized[0]["elevation_range"] = 0.0
        return normalized

    def _vert_exag_table_for_params(self):
        table = getattr(self, "vert_exag_table", None) or self._default_vert_exag_table()
        return [dict(row) for row in self._normalize_vert_exag_table(table)]

    def _add_vert_exag_table_to_params(self, params):
        if isinstance(params, dict):
            params["vert_exag_table"] = self._vert_exag_table_for_params()
        return params

    def _apply_vert_exag_table_from_params(self, params):
        if not isinstance(params, dict):
            return False
        saved = params.get("vert_exag_table")
        if saved is None:
            return False
        self.vert_exag_table = self._normalize_vert_exag_table(saved)
        return True

    def _valid_shaded_relief_cmap_names(self):
        return [cmap_name for _label, cmap_name in SHADED_RELIEF_CMAP_OPTIONS]

    def _normalize_shaded_relief_cmap(self, cmap_name):
        if cmap_name == "Spectral":
            cmap_name = "Spectral_r"
        elif cmap_name == "RdYlBu":
            cmap_name = "RdYlBu_r"
        elif cmap_name == "hsv_shift":
            cmap_name = "hsv_r"
        elif cmap_name == "ice_shift":
            cmap_name = "ice"
        valid = self._valid_shaded_relief_cmap_names()
        if cmap_name in valid:
            return cmap_name
        return DEFAULT_SHADED_RELIEF_CMAP

    def _shaded_relief_cmap_label(self, cmap_name=None):
        cmap_name = self._normalize_shaded_relief_cmap(
            cmap_name if cmap_name is not None else getattr(self, "shaded_relief_cmap", DEFAULT_SHADED_RELIEF_CMAP)
        )
        for label, name in SHADED_RELIEF_CMAP_OPTIONS:
            if name == cmap_name:
                return label
        return cmap_name

    def _load_shaded_relief_cmap(self):
        try:
            self.shaded_relief_cmap = DEFAULT_SHADED_RELIEF_CMAP
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, "r") as f:
                    config = json.load(f)
                saved = config.get("shaded_relief_cmap")
                if saved is not None:
                    self.shaded_relief_cmap = self._normalize_shaded_relief_cmap(saved)
        except Exception:
            self.shaded_relief_cmap = DEFAULT_SHADED_RELIEF_CMAP

    def _save_shaded_relief_cmap(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, "r") as f:
                    config = json.load(f)
            config["shaded_relief_cmap"] = self._normalize_shaded_relief_cmap(
                getattr(self, "shaded_relief_cmap", DEFAULT_SHADED_RELIEF_CMAP)
            )
            with open(self.CONFIG_FILENAME, "w") as f:
                json.dump(config, f, indent=2)
        except Exception:
            pass

    def _add_shaded_relief_cmap_to_params(self, params):
        if isinstance(params, dict):
            params["shaded_relief_cmap"] = self._normalize_shaded_relief_cmap(
                getattr(self, "shaded_relief_cmap", DEFAULT_SHADED_RELIEF_CMAP)
            )
        return params

    def _apply_shaded_relief_cmap_from_params(self, params):
        if not isinstance(params, dict):
            return False
        saved = params.get("shaded_relief_cmap")
        if saved is None:
            return False
        self.shaded_relief_cmap = self._normalize_shaded_relief_cmap(saved)
        return True

    def _default_slope_overlay_bands(self):
        return [dict(band) for band in DEFAULT_SLOPE_OVERLAY_BANDS]

    def _normalize_slope_overlay_bands(self, bands):
        defaults = self._default_slope_overlay_bands()
        if not isinstance(bands, (list, tuple)) or not bands:
            return defaults
        normalized = []
        for i, default in enumerate(defaults):
            src = bands[i] if i < len(bands) and isinstance(bands[i], dict) else {}
            min_val = src.get("min", default.get("min"))
            max_val = src.get("max", default.get("max"))
            try:
                min_val = None if min_val is None or str(min_val).strip() in ("", "-") else float(min_val)
            except (TypeError, ValueError):
                min_val = default.get("min")
            try:
                max_val = None if max_val is None or str(max_val).strip() in ("", "-") else float(max_val)
            except (TypeError, ValueError):
                max_val = default.get("max")
            color_fallback = default.get("color_hex", "#00ff00")
            color_hex = src.get("color_hex", color_fallback)
            if hasattr(self, "_normalize_slope_overlay_color_hex"):
                color_hex = self._normalize_slope_overlay_color_hex(color_hex, fallback=color_fallback)
            else:
                text = str(color_hex or color_fallback).strip()
                if not text.startswith("#"):
                    text = f"#{text}"
                color_hex = text.lower() if len(text) == 7 else color_fallback
            normalized.append({"min": min_val, "max": max_val, "color_hex": color_hex})
        return normalized

    def _load_slope_overlay_bands(self):
        try:
            self.slope_overlay_bands = self._default_slope_overlay_bands()
            self.slope_overlay_opacity = 40
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, "r") as f:
                    config = json.load(f)
                saved_bands = config.get("slope_overlay_bands")
                if saved_bands is not None:
                    self.slope_overlay_bands = self._normalize_slope_overlay_bands(saved_bands)
                else:
                    # Migrate legacy single-color key into range 1
                    legacy_color = config.get("slope_overlay_color_hex")
                    if legacy_color is not None:
                        self.slope_overlay_bands[0]["color_hex"] = (
                            self._normalize_slope_overlay_color_hex(legacy_color)
                            if hasattr(self, "_normalize_slope_overlay_color_hex")
                            else str(legacy_color)
                        )
                saved_opacity = config.get("slope_overlay_opacity")
                if saved_opacity is not None:
                    try:
                        self.slope_overlay_opacity = max(0, min(100, int(saved_opacity)))
                    except (TypeError, ValueError):
                        self.slope_overlay_opacity = 40
            if hasattr(self, "_sync_slope_overlay_legacy_vars"):
                self._sync_slope_overlay_legacy_vars()
            if hasattr(self, "_sync_slope_overlay_opacity_widget"):
                self._sync_slope_overlay_opacity_widget()
        except Exception:
            self.slope_overlay_bands = self._default_slope_overlay_bands()
            self.slope_overlay_opacity = 40
            if hasattr(self, "_sync_slope_overlay_legacy_vars"):
                self._sync_slope_overlay_legacy_vars()

    def _save_slope_overlay_bands(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, "r") as f:
                    config = json.load(f)
            bands = self._normalize_slope_overlay_bands(
                getattr(self, "slope_overlay_bands", None) or self._default_slope_overlay_bands()
            )
            self.slope_overlay_bands = bands
            config["slope_overlay_bands"] = bands
            # Keep legacy key in sync with range 1 color for older readers
            config["slope_overlay_color_hex"] = bands[0]["color_hex"]
            try:
                opacity = max(0, min(100, int(getattr(self, "slope_overlay_opacity", 40))))
            except (TypeError, ValueError):
                opacity = 40
            self.slope_overlay_opacity = opacity
            config["slope_overlay_opacity"] = opacity
            with open(self.CONFIG_FILENAME, "w") as f:
                json.dump(config, f, indent=2)
        except Exception:
            pass

    def _add_slope_overlay_bands_to_params(self, params):
        if isinstance(params, dict):
            bands = self._normalize_slope_overlay_bands(
                getattr(self, "slope_overlay_bands", None) or self._default_slope_overlay_bands()
            )
            params["slope_overlay_bands"] = bands
            params["slope_overlay_color_hex"] = bands[0]["color_hex"]
            try:
                params["slope_overlay_opacity"] = max(
                    0, min(100, int(getattr(self, "slope_overlay_opacity", 40)))
                )
            except (TypeError, ValueError):
                params["slope_overlay_opacity"] = 40
        return params

    def _apply_slope_overlay_bands_from_params(self, params):
        if not isinstance(params, dict):
            return False
        applied = False
        saved_bands = params.get("slope_overlay_bands")
        if saved_bands is not None:
            self.slope_overlay_bands = self._normalize_slope_overlay_bands(saved_bands)
            if hasattr(self, "_sync_slope_overlay_legacy_vars"):
                self._sync_slope_overlay_legacy_vars()
            applied = True
        else:
            # Legacy single-color fallback into range 1 only
            saved_color = params.get("slope_overlay_color_hex")
            if saved_color is not None:
                if not getattr(self, "slope_overlay_bands", None):
                    self.slope_overlay_bands = self._default_slope_overlay_bands()
                self.slope_overlay_bands[0]["color_hex"] = (
                    self._normalize_slope_overlay_color_hex(saved_color)
                    if hasattr(self, "_normalize_slope_overlay_color_hex")
                    else str(saved_color)
                )
                if hasattr(self, "_sync_slope_overlay_legacy_vars"):
                    self._sync_slope_overlay_legacy_vars()
                applied = True
        if "slope_overlay_opacity" in params:
            try:
                self.slope_overlay_opacity = max(0, min(100, int(params.get("slope_overlay_opacity"))))
            except (TypeError, ValueError):
                pass
            else:
                applied = True
        return applied

    def _add_geotiff_viz_params_to_params(self, params):
        self._add_vert_exag_table_to_params(params)
        self._add_shaded_relief_cmap_to_params(params)
        self._add_slope_overlay_bands_to_params(params)
        return params

    def _apply_geotiff_viz_params_from_params(self, params):
        applied_ve = self._apply_vert_exag_table_from_params(params)
        applied_cmap = self._apply_shaded_relief_cmap_from_params(params)
        applied_slope = self._apply_slope_overlay_bands_from_params(params)
        if applied_cmap and hasattr(self, "_update_shaded_relief_cmap_button"):
            self._update_shaded_relief_cmap_button()
        if applied_slope and hasattr(self, "_sync_slope_overlay_band_widgets"):
            self._sync_slope_overlay_band_widgets()
        if applied_slope and hasattr(self, "_sync_slope_overlay_opacity_widget"):
            self._sync_slope_overlay_opacity_widget()
        return applied_ve or applied_cmap or applied_slope

    def _load_vert_exag_table(self):
        try:
            defaults = self._default_vert_exag_table()
            self.vert_exag_table = [dict(row) for row in defaults]
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, "r") as f:
                    config = json.load(f)
                saved = config.get("vert_exag_table")
                if saved is not None:
                    self.vert_exag_table = self._normalize_vert_exag_table(saved)
        except Exception:
            self.vert_exag_table = [dict(row) for row in self._default_vert_exag_table()]

    def _save_vert_exag_table(self):
        try:
            config = {}
            if os.path.exists(self.CONFIG_FILENAME):
                with open(self.CONFIG_FILENAME, "r") as f:
                    config = json.load(f)
            table = getattr(self, "vert_exag_table", None) or self._default_vert_exag_table()
            config["vert_exag_table"] = self._normalize_vert_exag_table(table)
            with open(self.CONFIG_FILENAME, "w") as f:
                json.dump(config, f, indent=2)
        except Exception:
            pass

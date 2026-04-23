"""
Config: load/save last-used directories to/from JSON config file.
_load_last_used_dir, _save_last_used_dir, _load_last_geotiff_dir, _save_last_geotiff_dir,
_load_last_survey_params_dir, _save_last_survey_params_dir, _load_last_export_dir, _save_last_export_dir,
_load_last_ref_import_dir, _save_last_ref_import_dir, _load_last_line_import_dir, _save_last_line_import_dir,
_load_last_perf_import_dir, _save_last_perf_import_dir, _load_last_shapefile_dir, _save_last_shapefile_dir,
_load_export_type_options, _save_export_type_options.
"""
import json
import os


class ConfigMixin:
    """Mixin for loading/saving last-used directories from config file."""

    def _default_export_type_options(self):
        return {
            "esri_shapefile": True,
            "sis_asciiplan": True,
            "gpx": True,
            "text_csv": True,
            "text_txt": True,
            "hypack_lnw": True,
            "map_png": True,
            "profiles_png": True,
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
                    for key in defaults:
                        if key in saved:
                            self.export_type_options[key] = bool(saved[key])
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

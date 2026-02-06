"""
Config: load/save last-used directories to/from JSON config file.
_load_last_used_dir, _save_last_used_dir, _load_last_geotiff_dir, _save_last_geotiff_dir,
_load_last_survey_params_dir, _save_last_survey_params_dir, _load_last_export_dir, _save_last_export_dir,
_load_last_ref_import_dir, _save_last_ref_import_dir, _load_last_line_import_dir, _save_last_line_import_dir.
"""
import json
import os


class ConfigMixin:
    """Mixin for loading/saving last-used directories from config file."""

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
            config = {'last_used_dir': self.last_used_dir}
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

# -*- mode: python ; coding: utf-8 -*-
# PyInstaller spec for SAT Planner. Based on working old.spec; exe name/version from sat_planner.constants.

import sys
import os

block_cipher = None

# Ensure project root is on path so we can import bundle helpers
_spec_dir = os.path.dirname(os.path.abspath(SPEC))
sys.path.insert(0, _spec_dir)

# Read version without importing sat_planner.constants (that imports geospatial libs).
import re
_constants_path = os.path.join(_spec_dir, 'sat_planner', 'constants.py')
with open(_constants_path, encoding='utf-8') as _f:
    _version_match = re.search(r'__version__\s*=\s*["\']([^"\']+)', _f.read())
if not _version_match:
    raise RuntimeError(f'Could not read __version__ from {_constants_path}')
__version__ = _version_match.group(1)
exe_name = 'SAT_Planner_v' + __version__  # e.g. SAT_Planner_v2026.01

# Icon path relative to project root
icon_path = os.path.join(_spec_dir, 'media', 'CCOM.ico')

# Data files: package data plus conda GDAL/PROJ share trees (not inside site-packages).
from PyInstaller.utils.hooks import collect_data_files
from pyi_geospatial_bundle import (
    conda_library_dirs,
    collect_gdal_dlls,
    collect_gdal_proj_data,
)

datas = []
for pkg in ('pyproj', 'fiona', 'shapely', 'rasterio', 'matplotlib', 'colormaps'):
    try:
        datas += collect_data_files(pkg)
    except Exception:
        pass
datas += [(os.path.join(_spec_dir, 'media'), 'media')]

binaries = []
_conda_bin, _conda_share = conda_library_dirs()
if _conda_share:
    _gdal_proj_data = collect_gdal_proj_data(_conda_share)
    datas += _gdal_proj_data
    print(f'Bundling GDAL/PROJ data: {len(_gdal_proj_data)} director{"y" if len(_gdal_proj_data) == 1 else "ies"}')
else:
    print('WARNING: conda Library/share not found; GDAL/PROJ data will not be bundled.')
if _conda_bin:
    binaries += collect_gdal_dlls(_conda_bin)
    print(f'Bundling GDAL/PROJ DLLs: {len(binaries)} file(s)')
else:
    print('WARNING: conda Library/bin not found; GDAL/PROJ DLLs will not be bundled.')

a = Analysis(
    ['SAT_Planner_PyQt.py'],
    pathex=[_spec_dir],
    binaries=binaries,
    datas=datas,
    hiddenimports=[
        # Refactored package
        'sat_planner',
        'sat_planner.app_core',
        'sat_planner.constants',
        'sat_planner.utils_geo',
        'sat_planner.utils_ui',
        'sat_planner.mixins.basemap_mixin',
        'sat_planner.mixins.geotiff_mixin',
        'sat_planner.mixins.plotting_mixin',
        'sat_planner.mixins.reference_mixin',
        'sat_planner.mixins.survey_parsers_mixin',
        'sat_planner.mixins.gmrt_download_mixin',
        'sat_planner.mixins.calibration_mixin',
        'sat_planner.mixins.line_planning_mixin',
        'sat_planner.mixins.performance_mixin',
        'sat_planner.mixins.adcp_mixin',
        'sat_planner.performance_import_dialog',
        'sat_planner.mixins.profiles_mixin',
        'sat_planner.mixins.map_interaction_mixin',
        'sat_planner.mixins.export_import_mixin',
        'sat_planner.mixins.config_mixin',
        # GMRT Download dialog (still used by import-path tooling / legacy)
        'sat_planner.gmrt_dialog',
        'sat_planner.gmrt_dialog.main_window',
        'sat_planner.gmrt_dialog.map_widget',
        'sat_planner.gmrt_dialog.config',
        'sat_planner.gmrt_dialog.workers',
        'sat_planner.gmrt_dialog.workers.download_worker',
        'sat_planner.gmrt_dialog.workers.map_worker',
        'sat_planner.gmrt_dialog.workers.mosaic_worker',
        # Unified bathymetry download dialog (vendored from Bathymetry_Downloader)
        'sat_planner.bathymetry_download',
        'sat_planner.bathymetry_download.data_sources',
        'sat_planner.bathymetry_download.main_window',
        'sat_planner.bathymetry_download.ui_layout',
        'sat_planner.bathymetry_download.map_widget',
        'sat_planner.bathymetry_download.map_loaders',
        'sat_planner.bathymetry_download.service_loader',
        'sat_planner.bathymetry_download.ui_widgets',
        'sat_planner.bathymetry_download.download_module',
        'sat_planner.bathymetry_download.gmrt_module',
        'sat_planner.bathymetry_download.geo_utils',
        'sat_planner.gmrt_split',
        'sat_planner.import_survey_dialog',
        'sat_planner.map_options_dialog',
        'sat_planner.dyn_vert_exag_dialog',
        # Geospatial: from working old.spec (do not add osgeo)
        'fiona', 'shapely', 'pyproj', 'rasterio',
        'rasterio.sample', 'rasterio.io', 'rasterio.warp', 'rasterio.transform',
        'rasterio.crs', 'rasterio.features', 'rasterio.mask', 'rasterio.plot',
        'rasterio.windows', 'rasterio.errors', 'rasterio.dtypes', 'rasterio.profiles',
        'rasterio.env', 'rasterio.vrt', 'rasterio._base', 'rasterio._io',
        'rasterio._warp', 'rasterio._transform', 'rasterio._crs', 'rasterio._features',
        'rasterio._mask', 'rasterio._plot', 'rasterio._windows', 'rasterio._errors',
        'rasterio._dtypes', 'rasterio._profiles', 'rasterio._env', 'rasterio._vrt',
        'pyproj.datadir', 'pyproj.crs', 'pyproj.transformer', 'pyproj.enums',
        'pyproj.exceptions',
        'fiona.crs', 'fiona.errors', 'fiona.schema',
        'shapely.geometry', 'shapely.ops', 'shapely.validation', 'shapely.speedups',
        'shapely.algorithms',
        'matplotlib.backends.backend_qtagg', 'matplotlib.backends.backend_qt5agg',
        'matplotlib.figure', 'matplotlib.colors', 'matplotlib.pyplot',
        'PyQt6', 'PyQt6.QtCore', 'PyQt6.QtGui', 'PyQt6.QtWidgets',
        'requests',  # optional: for GMRT grid download on calibration import
        'colormaps', 'colormaps.cmaps', 'colormaps._registry', 'colormaps.colormap',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[os.path.join(_spec_dir, 'pyi_rth_geospatial_dll.py')],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

# One-file exe (same as before refactoring). Geospatial libs are imported in SAT_Planner_PyQt.py so PyInstaller bundles them.
exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name=exe_name,
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=icon_path,
)

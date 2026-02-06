# -*- mode: python ; coding: utf-8 -*-
# PyInstaller spec for SAT Planner. Based on working old.spec; exe name/version from sat_planner.constants.

import sys
import os

block_cipher = None

# Ensure project root is on path so we can import version
_spec_dir = os.path.dirname(os.path.abspath(SPEC))
sys.path.insert(0, _spec_dir)

from sat_planner.constants import __version__
exe_name = 'SAT_Planner_v' + __version__  # e.g. SAT_Planner_v2026.01

# Icon path relative to project root
icon_path = os.path.join(_spec_dir, 'media', 'CCOM.ico')

# Data files: same pattern as working old.spec (pyproj/fiona/shapely/rasterio need their data)
from PyInstaller.utils.hooks import collect_data_files

datas = []
for pkg in ('pyproj', 'fiona', 'shapely', 'rasterio', 'matplotlib'):
    try:
        datas += collect_data_files(pkg)
    except Exception:
        pass
datas += [(os.path.join(_spec_dir, 'media'), 'media')]

a = Analysis(
    ['SAT_Planner_PyQt.py'],
    pathex=[_spec_dir],
    binaries=[],
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
        'sat_planner.mixins.calibration_mixin',
        'sat_planner.mixins.line_planning_mixin',
        'sat_planner.mixins.profiles_mixin',
        'sat_planner.mixins.map_interaction_mixin',
        'sat_planner.mixins.export_import_mixin',
        'sat_planner.mixins.config_mixin',
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
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=icon_path,
)

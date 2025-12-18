# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_data_files

block_cipher = None

# Collect data files for geospatial libraries
datas = []
datas += collect_data_files('pyproj')
datas += collect_data_files('fiona')
datas += collect_data_files('shapely')
datas += collect_data_files('rasterio')
datas += collect_data_files('matplotlib')
datas += [('media/CCOM.ico', 'media')]
# Include CCOM.png if it exists for the About dialog
import os
if os.path.exists('media/CCOM.png'):
    datas += [('media/CCOM.png', 'media')]

a = Analysis(
    ['SAT_Planner_PyQt.py'],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=[
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
        'fiona._shim', 'fiona.crs', 'fiona.errors', 'fiona.schema',
        'shapely.geometry', 'shapely.ops', 'shapely.validation', 'shapely.speedups', 
        'shapely.algorithms',
        'matplotlib.backends.backend_qtagg', 'matplotlib.backends.backend_qt5agg',
        'matplotlib.figure', 'matplotlib.colors', 'matplotlib.pyplot',
        'PyQt6', 'PyQt6.QtCore', 'PyQt6.QtGui', 'PyQt6.QtWidgets'
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='Sat_Planner_v2025.10',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,  # No console window
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon='media/CCOM.ico',
)






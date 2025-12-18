# -*- mode: python ; coding: utf-8 -*-
# macOS spec file for SAT Planner
# Build with: pyinstaller Sat_Planner_macOS.spec

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
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

# For macOS, use BUNDLE instead of EXE to create an .app bundle
# Note: Icon file should be in .icns format for macOS
# Convert .ico to .icns using: sips -s format icns media/CCOM.ico --out media/CCOM.icns
# Or use online converters or iconutil command
app = BUNDLE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='SAT_Planner',
    icon='media/CCOM.icns',  # Use .icns format for macOS (create from .ico if needed)
    bundle_identifier='edu.unh.ccom.satplanner',
    version=None,
    info_plist={
        'NSPrincipalClass': 'NSApplication',
        'NSHighResolutionCapable': 'True',
        'CFBundleShortVersionString': '2025.11',
        'CFBundleVersion': '2025.11',
        'NSHumanReadableCopyright': 'Copyright (c) 2025, University of New Hampshire Center for Coastal and Ocean Mapping / Joint Hydrographic Center (UNH/CCOM-JHC)',
    },
)


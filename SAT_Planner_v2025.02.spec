# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_data_files

datas = []
datas += collect_data_files('pyproj')
datas += collect_data_files('fiona')
datas += collect_data_files('shapely')
datas += collect_data_files('rasterio')
datas += collect_data_files('matplotlib')


a = Analysis(
    ['C:\\Users\\pjohnson\\PycharmProjects\\MultibeamToolsAI\\MultibeamToolsAI\\multibeam_tools\\apps\\sat_qat_planner.py'],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=['fiona', 'shapely', 'pyproj', 'rasterio', 'rasterio.sample', 'rasterio.io', 'rasterio.warp', 'rasterio.transform', 'rasterio.crs', 'rasterio.features', 'rasterio.mask', 'rasterio.plot', 'rasterio.windows', 'rasterio.errors', 'rasterio.dtypes', 'rasterio.profiles', 'rasterio.env', 'rasterio.vrt', 'rasterio._base', 'rasterio._io', 'rasterio._warp', 'rasterio._transform', 'rasterio._crs', 'rasterio._features', 'rasterio._mask', 'rasterio._plot', 'rasterio._windows', 'rasterio._errors', 'rasterio._dtypes', 'rasterio._profiles', 'rasterio._env', 'rasterio._vrt', 'pyproj.datadir', 'pyproj.crs', 'pyproj.transformer', 'pyproj.enums', 'pyproj.exceptions', 'fiona._shim', 'fiona.crs', 'fiona.errors', 'fiona.schema', 'shapely.geometry', 'shapely.ops', 'shapely.validation', 'shapely.speedups', 'shapely.algorithms', 'matplotlib.backends.backend_tkagg', 'matplotlib.backends.backend_tk', 'matplotlib.figure', 'matplotlib.colors', 'matplotlib.pyplot', 'tkinter', 'tkinter.filedialog', 'tkinter.messagebox', 'tkinter.ttk'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='SAT_QAT_Planner_v2025.02',
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
    icon=['C:\\Users\\pjohnson\\PycharmProjects\\MultibeamToolsAI\\MultibeamToolsAI\\media\\CCOM.ico'],
)

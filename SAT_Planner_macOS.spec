# -*- mode: python ; coding: utf-8 -*-
# PyInstaller spec for the macOS SAT Planner app.
# Run build_app.sh. This spec produces dist/SAT_Planner_v<version>.app.

import os
import re
import sys

block_cipher = None

_spec_dir = os.path.dirname(os.path.abspath(SPEC))
sys.path.insert(0, _spec_dir)

if not sys.platform.startswith("darwin"):
    raise SystemExit("SAT_Planner_macOS.spec must be built on macOS. Use build_exe.bat on Windows.")

_constants_path = os.path.join(_spec_dir, "sat_planner", "constants.py")
with open(_constants_path, encoding="utf-8") as _f:
    _version_match = re.search(r'__version__\s*=\s*["\']([^"\']+)', _f.read())
if not _version_match:
    raise RuntimeError(f"Could not read __version__ from {_constants_path}")
__version__ = _version_match.group(1)
exe_name = "SAT_Planner_v" + __version__

icon_path = os.path.join(_spec_dir, "media", "CCOM.icns")
if not os.path.isfile(icon_path):
    raise SystemExit(
        "media/CCOM.icns is missing. Run ./build_app.sh, which builds it from media/CCOM.png."
    )
entitlements = os.path.join(_spec_dir, "macos", "entitlements.plist")

from PyInstaller.utils.hooks import collect_data_files, collect_submodules
from pyi_geospatial_bundle import (
    collect_macos_dylibs,
    collect_macos_gdal_data,
    conda_macos_prefix,
)

datas = []
for pkg in ("pyproj", "fiona", "shapely", "rasterio", "matplotlib", "colormaps"):
    try:
        datas += collect_data_files(pkg)
    except Exception:
        pass
datas += [(os.path.join(_spec_dir, "media"), "media")]
_maplibre_dir = os.path.join(_spec_dir, "sat_planner", "maplibre")
for _map_asset in ("index.html", "maplibre-gl.js", "maplibre-gl.css"):
    datas.append((os.path.join(_maplibre_dir, _map_asset), os.path.join("sat_planner", "maplibre")))

binaries = []
_prefix = conda_macos_prefix()
if _prefix:
    _gdal_proj_data = collect_macos_gdal_data(_prefix)
    datas += _gdal_proj_data
    binaries += collect_macos_dylibs(_prefix)
    print(f"Bundling GDAL/PROJ data: {len(_gdal_proj_data)} director{'y' if len(_gdal_proj_data) == 1 else 'ies'}")
    print(f"Bundling conda dylibs: {len(binaries)} file(s)")
else:
    print("WARNING: conda share/gdal or share/proj not found; GeoTIFF support may fail in the app.")

hiddenimports = collect_submodules("sat_planner")
# Cython extensions (rasterio._base and friends) import pure-Python modules
# such as rasterio.serde. PyInstaller cannot see those imports.
for _pkg in ("rasterio", "fiona", "shapely", "pyproj"):
    try:
        hiddenimports += collect_submodules(_pkg)
    except Exception as _exc:
        print(f"WARNING: could not collect {_pkg} submodules: {_exc}")
hiddenimports += [
    "PyQt6.QtWebEngineWidgets",
    "PyQt6.QtWebEngineCore",
    "PyQt6.QtWebChannel",
    "matplotlib.backends.backend_qtagg",
    "requests",
    "colormaps",
]

a = Analysis(
    ["SAT_Planner_PyQt.py"],
    pathex=[_spec_dir],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[os.path.join(_spec_dir, "pyi_rth_geospatial_dll.py")],
    excludes=[],
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name=exe_name,
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity="-",
    entitlements_file=entitlements,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name=exe_name,
)

app = BUNDLE(
    coll,
    name=exe_name + ".app",
    icon=icon_path,
    bundle_identifier="edu.unh.ccom.satplanner",
    version=__version__,
    info_plist={
        "NSHighResolutionCapable": True,
        "CFBundleName": "SAT Planner",
        "CFBundleDisplayName": "SAT Planner",
        "CFBundleShortVersionString": __version__,
        "CFBundleVersion": __version__,
    },
)

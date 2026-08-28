"""
Shared constants for SAT Planner.
Single source of truth for version, config path, and geospatial availability.
"""
import os
import sys
import traceback

__version__ = "2026.38"

CONFIG_FILENAME = os.path.join(os.path.expanduser("~"), ".cal_ref_planner_config.json")

# Register third-party colormaps with matplotlib (via colormaps package).
_OPTIONAL_COLORMAP_NAMES = ("ice", "arctic", "sapphire", "torch")
_COLORMAPS_EXTRA = ()
try:
    import colormaps as cmaps

    _available = []
    for _name in _OPTIONAL_COLORMAP_NAMES:
        try:
            getattr(cmaps, _name)
            _available.append((_name, _name))
        except (AttributeError, OSError, ValueError):
            pass
    _COLORMAPS_EXTRA = tuple(_available)
except ImportError:
    pass

# Shaded Relief elevation overlay colormaps: (button label, matplotlib cmap name)
SHADED_RELIEF_CMAP_OPTIONS = (
    ("rainbow", "rainbow"),
    ("viridis", "viridis"),
    ("cividis", "cividis"),
    ("turbo", "turbo"),
    ("CnBu", "BuGn_r"),
    ("Greys", "Greys_r"),
    *_COLORMAPS_EXTRA,
    ("RdYlBu", "RdYlBu_r"),
    ("Spectral", "Spectral_r"),
    ("hsv", "hsv_r"),
    ("jet", "jet"),
    ("winter", "winter"),
)
DEFAULT_SHADED_RELIEF_CMAP = "rainbow"

# GeoTIFF slope overlay bands: min/max in degrees (None = undefined), color as #rrggbb
DEFAULT_SLOPE_OVERLAY_BANDS = (
    {"min": 10.0, "max": 20.0, "color_hex": "#00ff00"},
    {"min": None, "max": None, "color_hex": "#1e90ff"},  # dodgerblue
    {"min": None, "max": None, "color_hex": "#ff4500"},  # orangered
)

# Conditional imports for geospatial libraries (re-exported for use by SAT_Planner_PyQt)
GEOSPATIAL_LIBS_AVAILABLE = True
rasterio = None
transform = None
RasterioIOError = None
Window = None
window_transform = None
window_bounds = None
rowcol = None
reproject = None
Resampling = None
pyproj = None
CRSError = None
LineString = None
fiona = None
try:
    import rasterio
    from rasterio import transform
    from rasterio.errors import RasterioIOError
    from rasterio.windows import Window
    from rasterio.windows import transform as window_transform
    from rasterio.windows import bounds as window_bounds
    from rasterio.transform import rowcol
    from rasterio.warp import reproject, Resampling
    import pyproj
    from pyproj.exceptions import CRSError
    from shapely.geometry import LineString
    import fiona
except (ImportError, OSError) as e:
    GEOSPATIAL_LIBS_AVAILABLE = False
    _geo_msg = (
        f"Warning: A geospatial library not found: {e}. "
        "GeoTIFF/Shapefile features will be disabled."
    )
    print(_geo_msg)
    if getattr(sys, "frozen", False):
        try:
            _log_path = os.path.join(os.path.expanduser("~"), "sat_planner_geospatial_error.log")
            with open(_log_path, "w", encoding="utf-8") as _log:
                _log.write(_geo_msg + "\n\n")
                traceback.print_exc(file=_log)
        except OSError:
            pass

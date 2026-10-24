"""
Shared constants for SAT Planner.
Single source of truth for version, config path, and geospatial availability.
"""
import os
import sys
import traceback

__version__ = "2026.46"

CONFIG_FILENAME = os.path.join(os.path.expanduser("~"), ".cal_ref_planner_config.json")

PLANNING_PLACEHOLDER_TEXT = "Load a Test Plan or Bathymetry\nto begin planning"

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

# GeoTIFF depth overlay bands: min/max elevation/depth (None = undefined), color as #rrggbb
DEFAULT_DEPTH_OVERLAY_BANDS = (
    {"min": None, "max": None, "color_hex": "#fa8072"},  # salmon
    {"min": None, "max": None, "color_hex": "#3cb371"},  # mediumseagreen
    {"min": None, "max": None, "color_hex": "#6a5acd"},  # slateblue
)

# Kongsberg depth modes. Depths are positive down, in meters.
# min_m is inclusive. A shared boundary belongs to the deeper mode, so max_m is
# exclusive. max_m None means the mode continues to full ocean depth.
# Color is chosen by mode name, so the same name matches across systems.
DEPTH_MODE_NONE_LABEL = "No system selected"

DEPTH_MODE_COLORS = {
    "Very Shallow": "#EEBABA",
    "Shallow": "#FF4B4B",
    "Medium": "#FFB84F",
    "Deep": "#FFF529",
    "Deeper": "#33CC33",
    "Very Deep": "#1DD4FF",
    "Extra Deep": "#FF40FF",
    "Extreme Deep": "#9966FF",
}


def _kongsberg_modes(specs):
    return tuple(
        {
            "name": name,
            "min_m": min_m,
            "max_m": max_m,
            "color_hex": DEPTH_MODE_COLORS[name],
        }
        for name, min_m, max_m in specs
    )


KONGSBERG_DEPTH_MODES = {
    "EM712": _kongsberg_modes((
        ("Very Shallow", 2, 100),
        ("Shallow", 100, 200),
        ("Medium", 200, 300),
        ("Deep", 300, 500),
        ("Very Deep", 500, 1000),
        ("Extra Deep", 1000, None),
    )),
    "EM304": _kongsberg_modes((
        ("Shallow", 10, 250),
        ("Medium", 250, 750),
        ("Deep", 750, 1000),
        ("Deeper", 1000, 3300),
        ("Very Deep", 3300, 5000),
        ("Extra Deep", 5000, 7000),
        ("Extreme Deep", 7000, None),
    )),
    "EM124": _kongsberg_modes((
        ("Shallow", 10, 450),
        ("Medium", 450, 1000),
        ("Deep", 1000, 2600),
        ("Deeper", 2600, 9000),
        ("Very Deep", 9000, None),
    )),
}
KONGSBERG_DEPTH_MODE_SYSTEMS = ("EM712", "EM304", "EM124")
KONGSBERG_DEPTH_MODE_DISPLAY_LABELS = {
    "EM712": "EM712 Depth Modes",
    "EM304": "EM304 Depth Modes",
    "EM124": "EM124 Depth Modes",
}


def kongsberg_depth_mode_indices(depth_m, bands):
    """Mode index for each positive-down depth. -1 means no mode (too shallow or NaN).

    Bands must be shallow-to-deep. A depth equal to a boundary takes the deeper mode.
    """
    import numpy as np

    depth = np.asarray(depth_m, dtype=float)
    indices = np.full(depth.shape, -1, dtype=np.int16)
    finite = np.isfinite(depth)
    for i, band in enumerate(bands):
        indices = np.where(finite & (depth >= band["min_m"]), i, indices)
    return indices

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

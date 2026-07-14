"""
Shared constants for SAT Planner.
Single source of truth for version, config path, and geospatial availability.
"""
import os

__version__ = "2026.36"

CONFIG_FILENAME = os.path.join(os.path.expanduser("~"), ".cal_ref_planner_config.json")

# Shaded Relief elevation overlay colormaps: (button label, matplotlib cmap name)
SHADED_RELIEF_CMAP_OPTIONS = (
    ("rainbow", "rainbow"),
    ("viridis", "viridis"),
    ("cividis", "cividis"),
    ("turbo", "turbo"),
    ("CnBu", "BuGn_r"),
    ("Greys", "Greys_r"),
    ("RdYlBu", "RdYlBu"),
    ("Spectral", "Spectral_r"),
    ("hsv", "hsv_r"),
    ("jet", "jet"),
    ("winter", "winter"),
)
DEFAULT_SHADED_RELIEF_CMAP = "rainbow"

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
except ImportError as e:
    GEOSPATIAL_LIBS_AVAILABLE = False
    print(f"Warning: A geospatial library not found: {e}. GeoTIFF/Shapefile features will be disabled.")

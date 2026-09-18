"""Data source registry for SAT Planner bathymetry downloads (vendored).

GEBCO 2026 TID remains in DATA_SOURCES for internal GEBCO masking but is
excluded from UI_DATA_SOURCES so it is never offered as a user-facing source.
"""

from .gmrt_module import LAT_LIMIT as GMRT_LAT_LIMIT

WORLD_EXTENT_4326 = (-180.0, -90.0, 180.0, 90.0)
DEFAULT_DATA_SOURCE = "GMRT Topo-Bathy"
HIDDEN_UI_SOURCES = frozenset({"GEBCO 2026 TID"})


def build_data_sources():
    """Return the configured data source registry."""
    return {
        "GEBCO 2026": {
            "url": "https://gis.ccom.unh.edu/server/rest/services/GEBCO2026/gebco_2026_IS/ImageServer",
            "display_url": "https://gis.ccom.unh.edu/server/rest/services/GEBCO/GEBCO_2026_Depths_Haxby_GCS/MapServer",
            "land_display_url": "https://gis.ccom.unh.edu/server/rest/services/GEBCO/GEBCO_2026_Land_Grey_GCS/MapServer",
            "bathymetry_raster_function": "None",
            "hillshade_raster_function": "None",
            "default_extent": WORLD_EXTENT_4326,
            "service_crs": "EPSG:4326",
            "native_resolution_only": True,
            "native_pixel_size_degrees": 0.004166666666666667,
            "show_output_data_types": False,
            "ignore_source_nodata": True,
            "download_filename_prefix": "GEBCO_2026",
            "attribution": (
                "GEBCO Bathymetric Compilation Group 2026. The GEBCO_2026 Grid - a continuous "
                "terrain model for oceans and land at 15 arc-second intervals."
            ),
            "attribution_url": (
                "https://www.bodc.ac.uk/data/published_data_library/catalogue/"
                "10.5285/4f68d5c7-45eb-f999-e063-7086abc036fa"
            ),
        },
        "GEBCO 2026 TID": {
            "url": "https://gis.ccom.unh.edu/server/rest/services/GEBCO2026/gebco_2026_tid_IS/ImageServer",
            "display_url": "https://gis.ccom.unh.edu/server/rest/services/GEBCO2026/GEBCO_2026_TID_GCS/MapServer",
            "land_display_url": "https://gis.ccom.unh.edu/server/rest/services/GEBCO/GEBCO_2026_Land_Grey_GCS/MapServer",
            "bathymetry_raster_function": "None",
            "hillshade_raster_function": "None",
            "default_extent": WORLD_EXTENT_4326,
            "service_crs": "EPSG:4326",
            "native_resolution_only": True,
            "native_pixel_size_degrees": 0.004166666666666667,
            "show_output_data_types": False,
            "ignore_source_nodata": True,
            "download_filename_prefix": "GEBCO_2026_TID",
            "attribution": (
                "GEBCO Bathymetric Compilation Group 2026. The GEBCO_2026 Grid - a continuous "
                "terrain model for oceans and land at 15 arc-second intervals."
            ),
            "attribution_url": (
                "https://www.bodc.ac.uk/data/published_data_library/catalogue/"
                "10.5285/4f68d5c7-45eb-f999-e063-7086abc036fa"
            ),
        },
        "GMRT Topo-Bathy": {
            "api": "gmrt",
            "url": "https://www.gmrt.org/services/GridServer",
            "preview_url": "https://www.gmrt.org/services/ImageServer",
            "gmrt_layer": "topo",
            "bathymetry_raster_function": "None",
            "hillshade_raster_function": "None",
            "default_extent": (-180.0, -GMRT_LAT_LIMIT, 180.0, GMRT_LAT_LIMIT),
            "service_crs": "EPSG:4326",
            "native_resolution_only": False,
            "show_output_data_types": False,
            "ignore_source_nodata": False,
            "cell_size_meters_options": [60, 120, 240, 480, 960],
            "default_cell_size_meters": 120,
            "download_filename_prefix": "gmrt_topo",
            "attribution": (
                "Ryan, W.B.F. et al. (2009), Global Multi-Resolution Topography (GMRT) synthesis, "
                "Geochem. Geophys. Geosyst., doi:10.1029/2008GC002332"
            ),
            "attribution_url": "https://doi.org/10.1029/2008GC002332",
            "lat_limit": GMRT_LAT_LIMIT,
        },
        "GMRT Topo-Bathy (Observed Only)": {
            "api": "gmrt",
            "url": "https://www.gmrt.org/services/GridServer",
            "preview_url": "https://www.gmrt.org/services/ImageServer",
            "gmrt_layer": "topo-mask",
            "bathymetry_raster_function": "None",
            "hillshade_raster_function": "None",
            "default_extent": (-180.0, -GMRT_LAT_LIMIT, 180.0, GMRT_LAT_LIMIT),
            "service_crs": "EPSG:4326",
            "native_resolution_only": False,
            "show_output_data_types": False,
            "ignore_source_nodata": False,
            "cell_size_meters_options": [60, 120, 240, 480, 960],
            "default_cell_size_meters": 120,
            "download_filename_prefix": "gmrt_topo-mask",
            "attribution": (
                "Ryan, W.B.F. et al. (2009), Global Multi-Resolution Topography (GMRT) synthesis, "
                "Geochem. Geophys. Geosyst., doi:10.1029/2008GC002332"
            ),
            "attribution_url": "https://doi.org/10.1029/2008GC002332",
            "lat_limit": GMRT_LAT_LIMIT,
        },
        "NCEI Multibeam Mosaic Raw": {
            "url": "https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_mosaics/multibeam_mosaic_raw/ImageServer",
            "land_display_url": "https://gis.ccom.unh.edu/server/rest/services/GEBCO/GEBCO_2026_Land_Grey_GCS/MapServer",
            "bathymetry_raster_function": "ColorHillshadeHaxby_8000-0",
            "hillshade_raster_function": "None",
            "default_extent": WORLD_EXTENT_4326,
            "service_crs": "EPSG:4326",
            "native_resolution_only": True,
            "native_pixel_size_degrees": 8.333333333333334e-4,
            "show_output_data_types": False,
            "ignore_source_nodata": False,
            "configurable_cell_size_degrees": True,
            "show_bathymetry_only_notice": True,
            "download_filename_prefix": "multibeam_mosaic_raw",
            "attribution": "NOAA National Centers for Environmental Information (NCEI) Multibeam Mosaic",
            "attribution_url": (
                "https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_mosaics/"
                "multibeam_mosaic_raw/ImageServer"
            ),
        },
        "NCEI Multibeam Mosaic Proc": {
            "url": "https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_mosaics/multibeam_mosaic_processed/ImageServer",
            "land_display_url": "https://gis.ccom.unh.edu/server/rest/services/GEBCO/GEBCO_2026_Land_Grey_GCS/MapServer",
            "bathymetry_raster_function": "ColorHillshadeHaxby_8000-0",
            "hillshade_raster_function": "None",
            "default_extent": WORLD_EXTENT_4326,
            "service_crs": "EPSG:4326",
            "native_resolution_only": True,
            "native_pixel_size_degrees": 8.333333333333334e-4,
            "show_output_data_types": False,
            "ignore_source_nodata": False,
            "configurable_cell_size_degrees": True,
            "show_bathymetry_only_notice": True,
            "download_filename_prefix": "multibeam_mosaic_processed",
            "attribution": "NOAA National Centers for Environmental Information (NCEI) Multibeam Mosaic (Processed)",
            "attribution_url": (
                "https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_mosaics/"
                "multibeam_mosaic_processed/ImageServer"
            ),
        },
        "WGOM-LI-SNE Hi Resolution": {
            "url": (
                "https://gis.ccom.unh.edu/server/rest/services/WGOM_LI_SNE/"
                "WGOM_LI_SNE_BTY_4m_20231005_WMAS_2_IS/ImageServer"
            ),
            "bathymetry_raster_function": "StdDev - BlueGreen",
            "hillshade_raster_function": "Multidirectional Hillshade 3x",
            "default_extent": (-8254538.5, 4898559.25, -7411670.5, 5636075.25),
            "service_crs": "EPSG:3857",
            "native_resolution_only": False,
            "show_output_data_types": False,
            "ignore_source_nodata": False,
            "show_basemap": True,
            "show_hillshade": True,
            "default_cell_size_meters": 4.0,
            "download_filename_prefix": "WGOM_LI_SNE_4m",
            "attribution": (
                "Center for Coastal and Ocean Mapping (CCOM), University of New Hampshire — "
                "WGOM-LI-SNE Bathymetry (4 m)"
            ),
            "attribution_url": (
                "https://gis.ccom.unh.edu/server/rest/services/WGOM_LI_SNE/"
                "WGOM_LI_SNE_BTY_4m_20231005_WMAS_2_IS/ImageServer"
            ),
        },
        "WGOM-LI-SNE Regional": {
            "url": (
                "https://gis.ccom.unh.edu/server/rest/services/WGOM_LI_SNE/"
                "WGOM_LI_SNE_BTY_20231004_16m_2_WMAS_IS/ImageServer"
            ),
            "bathymetry_raster_function": "StdDev - BlueGreen",
            "hillshade_raster_function": "Multidirectional Hillshade 3x",
            "default_extent": (-8313630.50001078, 4898555.25001255, -7411662.50001078, 5636075.25001255),
            "service_crs": "EPSG:3857",
            "native_resolution_only": False,
            "show_output_data_types": False,
            "ignore_source_nodata": False,
            "show_basemap": True,
            "show_hillshade": True,
            "default_cell_size_meters": 16.0,
            "download_filename_prefix": "WGOM_LI_SNE_16m",
            "attribution": (
                "Center for Coastal and Ocean Mapping (CCOM), University of New Hampshire — "
                "WGOM-LI-SNE Bathymetry (16 m)"
            ),
            "attribution_url": (
                "https://gis.ccom.unh.edu/server/rest/services/WGOM_LI_SNE/"
                "WGOM_LI_SNE_BTY_20231004_16m_2_WMAS_IS/ImageServer"
            ),
        },
    }


DATA_SOURCES = build_data_sources()
UI_DATA_SOURCES = {
    name: cfg for name, cfg in DATA_SOURCES.items() if name not in HIDDEN_UI_SOURCES
}
UI_DATA_SOURCE_ORDER = (
    "GMRT Topo-Bathy",
    "GMRT Topo-Bathy (Observed Only)",
    "GEBCO 2026",
    "NCEI Multibeam Mosaic Raw",
    "NCEI Multibeam Mosaic Proc",
    "WGOM-LI-SNE Hi Resolution",
    "WGOM-LI-SNE Regional",
)

# Short tags appended to exported Full/View GeoTIFF filenames.
EXPORT_SOURCE_TAGS = {
    "GMRT Topo-Bathy": "GMRT",
    "GMRT Topo-Bathy (Observed Only)": "GMRT",
    "GEBCO 2026": "GEBCO",
    "NCEI Multibeam Mosaic Raw": "NCEI",
    "NCEI Multibeam Mosaic Proc": "NCEI",
    "WGOM-LI-SNE Hi Resolution": "CCOM",
    "WGOM-LI-SNE Regional": "CCOM",
}


def export_source_tag_for_data_source(name):
    """Return filename source tag (e.g. GMRT) for a Download Data source, or None."""
    if not name:
        return None
    return EXPORT_SOURCE_TAGS.get(name)


def get_source(data_sources, name):
    """Return config dict for a named data source."""
    return data_sources.get(name, {})


def is_gmrt_source(data_sources, name):
    """Return True for Lamont GMRT GridServer sources."""
    return get_source(data_sources, name).get("api") == "gmrt"


def uses_web_mercator_map(data_sources, name):
    """Return True when map extents/selection are in EPSG:3857 (WGOM)."""
    return get_source(data_sources, name).get("service_crs") == "EPSG:3857"


def uses_meter_cell_size(data_sources, name):
    """Return True when the UI shows meter cell-size options (WGOM or GMRT)."""
    ds = get_source(data_sources, name)
    return ds.get("service_crs") == "EPSG:3857" or ds.get("api") == "gmrt"


def download_filename_prefix(data_sources, name):
    """Return the filename prefix for downloaded GeoTIFFs."""
    ds = get_source(data_sources, name)
    prefix = ds.get("download_filename_prefix")
    if prefix:
        return prefix
    return name.replace(" ", "_")


def format_cell_size_degrees(value):
    """Format a cell size in degrees for display in the UI."""
    return f"{value:.12g}"


def shows_output_data_types(data_sources, name):
    """Return True when the source supports multiple output grid types."""
    return get_source(data_sources, name).get("show_output_data_types", False)


def bbox_sr_for_data_source(data_sources, name):
    """Return bboxSR for ImageServer export, or None for native service CRS."""
    ds = get_source(data_sources, name)
    if ds.get("api") == "gmrt":
        return None
    if ds.get("service_crs") == "EPSG:3857":
        return None
    if ds.get("display_url"):
        return None
    return "4326"


def output_mode_filename_suffix(mode):
    """Map an output mode to its filename suffix."""
    if mode == "bathymetry_only":
        return "bathymetry"
    if mode == "land_only":
        return "land"
    if mode == "direct_measurements_only":
        return "direct"
    if mode == "direct_unknown_measurements_only":
        return "direct_unknown"
    return mode

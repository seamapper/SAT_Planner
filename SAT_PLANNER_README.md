# SAT/QAT Planner

![Example Plot](media/SAT_Planner.jpg)

A comprehensive Shipboard Acceptance Testing (SAT) and Quality Assurance Testing (QAT) planning tool with GeoTIFF support, built with PyQt6. The GUI uses a consistent dark theme.

## Overview

The SAT/QAT Planner is a desktop application designed for planning and visualizing multibeam testing operations. It supports three main planning modes:
- **Calibration Survey Planning**: Plan pitch, roll, and heading calibration lines
- **Accuracy Survey Planning**: Generate parallel survey lines with customizable parameters
- **Line Planning**: Interactive line drawing with real-time elevation profiles

## Features

### Core Functionality
- **Multi-tab Interface**: Separate tabs for Calibration, Accuracy, and Line planning (left panel)
- **Dark Theme**: Qt GUI uses a dark theme; map (matplotlib) keeps default styling
- **GeoTIFF Support**: Load and visualize elevation data from GeoTIFF files
- **GMRT Download**: Optional download of GMRT bathymetry GeoTIFF when importing surveys (Calibration, Accuracy, Line tabs; configurable buffer)
- **Download Data**: Source selector in GeoTIFF Controls (default `Select Source`; current source `GMRT`) opens the "Download GMRT Grid" dialog immediately on selection and keeps the selected source after successful download
- **Dynamic Resolution**: Automatically adjust GeoTIFF resolution based on zoom level
- **Interactive Plotting**: Pan (middle mouse), zoom (scroll), no toolbar
- **Real-time Statistics**: Survey distances, times, and comprehensive statistics
- **Elevation Profiles**: View elevation and slope profiles for drawn lines
- **Activity Log**: Collapsible side panel below the map
- **Export**: DDD/DMM/DMS CSV and TXT, SIS asciiplan, Hypack LNW, Shapefile; statistics and *_info.txt with waypoint sections (DMM then DDD)
- **Survey Import**: DDD, DMS, DMM, LNW, CSV, GeoJSON
- **EEZ Overlay**: EEZ layer with opacity control (default 80%) and hover `GEONAME` tooltip lookup

### Calibration Survey Planning
- Draw pitch and roll calibration lines; generate heading lines from pitch line
- **Reverse Line Direction**: Flip start/end of Pitch, Roll, Heading1, Heading2 (checkboxes)
- Turn time parameter; optional GMRT download on import (checkbox + buffer)
- Import from DDD/DMS/DMM/LNW, CSV, GeoJSON; **suggested line assignment** (metadata labels or geometry: Pitch = middle parallel, Roll = fourth line, Heading1/2 = outer two)
- **Calibration Survey Info** dialog and *_info.txt: **Calibration Waypoints (DMM)** and **Calibration Waypoints (DDD)** (Pitch/Roll/Heading1/Heading2 start and end)
- Validation when heading offset > 2× shallowest depth
- Export via shared `export_utils` (DDD/DMM/DMS CSV/TXT, asciiplan, LNW)

### Accuracy Survey Planning
- Parallel lines with line length, spacing, heading, speed, turn time
- Optional GMRT download on import; import from DDD/DMS/DMM/LNW, CSV, GeoJSON; **suggested crossline and reference line order**
- **Accuracy Survey Info** dialog and *_info.txt: **Accuracy Waypoints (DMM)** and **Accuracy Waypoints (DDD)** (L1S/L1E, …, CLS/CLE)
- Survey time breakdown (main lines, crossline, transit, turn)
- Export via shared `export_utils`

### Line Planning
- Interactive line drawing; **Reverse Line Direction** (flip line start/end); real-time elevation profiles; edit by dragging waypoints
- Optional GMRT download on import; import/export (DDD, DMS, DMM, LNW, CSV, GeoJSON)
- **Survey Info** dialog and *_info.txt: **Line Plan Waypoints (DMM)** and **Line Plan Waypoints (DDD)** (WP1, WP2, …)
- Survey statistics for drawn lines; export via shared `export_utils`

### GeoTIFF Visualization
- Elevation/slope display; hillshade; dynamic resolution; multiple CRS
- Axis labels in degrees–decimal minutes (DDM)
- Contour interval and slope min/max updates are debounced while typing
- Survey legend is drawn above map overlays/layers
- EEZ overlay reloads on pan/zoom and supports paused-hover name lookup

### GeoJSON Metadata
- Calibration, Accuracy, and Line GeoJSON exports include `survey_speed`
- Calibration, Accuracy, and Line GeoJSON exports include saved `geotiff_path`
- On import, SAT Planner attempts to reload the saved GeoTIFF path; if unavailable, import continues with a warning

### Download GMRT Grid Dialog
- **Download Data** -> **GMRT** opens the "Download GMRT Grid" window (separate from main app)
- Set area (North/South/East/West or draw on map)
- Cell resolution: 100 m, 200 m, 400 m, or Custom (e.g. 50 m); default 100 m
- **Split Grid Into Bathymetry and Topography** is enabled by default; SAT Planner then loads the bathy grid
- Warning (orange) when estimated pixels > 16,000,000
- Download progress displays in-dialog (tile mode: `x of y`; single-grid mode: indeterminate "Downloading...")
- **Close GMRT Downloader** button at bottom
- Config: `~/.gmrtgrab_sat_planner_config.json`

## Requirements

- **Python**: 3.7+
- **Core**: PyQt6, matplotlib, numpy
- **Geospatial**: rasterio, pyproj, shapely, fiona
- **Optional**: Pillow (Imagery Basemap, NOAA ENC Charts), requests (GMRT download)

## Installation

### Pre-built Executable
- Download from [Releases](https://github.com/seamapper/SAT_Planner/releases): e.g. `SAT_Planner_v2026.14.exe` (Windows) or newer; version in filename from `sat_planner/constants.py`.

### From Source
1. Clone: `git clone https://github.com/seamapper/SAT_Planner.git && cd SAT_Planner`
2. Install: `pip install PyQt6 matplotlib numpy rasterio pyproj shapely fiona Pillow requests`
3. Run: `python SAT_Planner_PyQt.py`

## Building from Source

- **Windows**: `pip install pyinstaller` then run `build_exe.bat` or `pyinstaller SAT_Planner.spec`; output in `dist/` as `SAT_Planner_v<version>.exe` (version from `sat_planner/constants.py`).
- **macOS**: Create or use a PyInstaller spec that includes `sat_planner` and `SAT_Planner_PyQt.py`; icon `media/CCOM.icns` if available.

## Usage

1. **Load GeoTIFF** (optional): "Load GeoTIFF" or **Download Data** -> **GMRT** (opens GMRT dialog; if split used, app loads bathy grid).
2. **Planning mode**: Calibration, Accuracy, or Line tab.
3. **Configure** parameters, generate or draw, view statistics, export.

### Map Controls (no toolbar)
- **Left click**: Add waypoint (drawing); pick points in Calibration/Accuracy/Line.
- **Right click**: Finish drawing line.
- **Middle mouse**: Pan.
- **Scroll**: Zoom.

## Configuration

- `~/.cal_ref_planner_config.json` – user preferences, last directories, survey parameters.
- GMRT dialog: `~/.gmrtgrab_sat_planner_config.json`.

## Project Structure

- **`SAT_Planner_PyQt.py`** – Entry point.
- **`sat_planner/`** – Core package:
  - **`constants.py`** – Version, config path.
  - **`export_utils.py`** – Shared export: DDD/DMM/DMS CSV and TXT, SIS asciiplan, Hypack LNW; UTM zone from points (used by Calibration, Accuracy, Line planning).
  - **`utils_geo.py`**, **`utils_ui.py`** – Helpers.
  - **`gmrt_dialog/`** – Embedded GMRT Download dialog (GeoTIFF-only, optional split).
  - **`mixins/`** – Basemap, GeoTIFF, Plotting, SurveyParsers, GMRTDownload, Reference/Accuracy, Calibration, LinePlanning, Profiles, MapInteraction, ExportImport, Config.

## Version History

- **v2026.14**: Added GeoTIFF Controls **Download Data** source dropdown (`Select Source`, `GMRT`) that opens the source flow immediately and keeps selection after successful download. Added GMRT dialog download progress bar (`x of y` for tiled downloads, indeterminate for single download). Added EEZ pan/zoom-driven EEZ refresh, paused-hover EEZ `GEONAME` tooltip lookups (including without a loaded GeoTIFF), and default EEZ opacity 80%. Improved large-area/world-scale alignment for EEZ and Imagery Basemap overlays. Updated GMRT dialog so **Split Grid Into Bathymetry and Topography** is enabled by default. Fixed line-plan zoom GeoTIFF coverage refresh, line-plan profile refresh after import, and calibration Pitch Line Info refresh after survey import.
- **v2026.11**: UI and workflow updates. Accuracy tab naming in UI/docs (formerly Reference). Import/Export button labels updated (Import/Export Accuracy Survey, Import/Export Calibration Survey, Import/Export Line Survey) and reordered on tabs. "Download GMRT" defaults to unchecked and placement updated in import/export groups. Activity Log changed to collapsible side panel. GeoJSON export/import now includes `survey_speed` and saved `geotiff_path` for Calibration/Accuracy/Line; missing GeoTIFF path does not block import (warning + continue). Contour interval and slope min/max redraws are debounced while typing. Survey legend now draws above all overlays.
- **v2026.10**: Shared `export_utils` (DDD/DMM/DMS CSV/TXT, asciiplan, LNW; UTM from points). Calibration: Reverse Line Direction, import suggestion (metadata/geometry), Calibration Survey Info + Calibration Waypoints (DMM/DDD) in dialog and *_info.txt. Reference: import suggestion, Reference Survey Info + Reference Waypoints (DMM/DDD). Line planning: Reverse Line Direction, Survey Info + Line Plan Waypoints (DMM/DDD). Exe build uses version from constants; output `SAT_Planner_v2026.10.exe`.
- **v2026.09+**: GMRT dialog ("Download GMRT GeoTIFF"), "Download GMRT Grid" window, GeoTIFF-only, cell resolution 100/200/400/Custom (50 m), large-area warning >16M pixels, split → load bathy, "Close GMRT Downloader" button, Activity Log 320 px.
- **v2026.08**: Dark theme, Activity Log on right, GMRT on import (Calibration/Reference/Line), Line import DDD/DMS/DMM/LNW/CSV/GeoJSON, no nav toolbar (scroll zoom, middle-mouse pan).
- **v2026.04**: Hover coordinates in DDM; default window height 1110 px.
- **v2026.02**: Turn time; statistics breakdowns; heading offset validation; export with full statistics.
- **v2026.01**: Package refactor with mixins; DDM axis labels.
- **v2025.11**: Dynamic resolution fixes, About dialog.
- **v2025.10**: Imagery Basemap, NOAA ENC overlays, navigation toolbar (later removed in v2026.08).
- **v2025.09–v2025.01**: Line planning, PyQt6, metadata, contours, initial release.

## Troubleshooting

- **GeoTIFF**: Ensure valid file/CRS; use Dynamic Resolution for large files.
- **Missing modules**: `pip install rasterio pyproj shapely fiona` (and optionally Pillow, requests).
- **Limited mode**: App runs without geospatial libs but GeoTIFF disabled.

## License

BSD 3-Clause. See [LICENSE](LICENSE). Copyright (c) 2025, University of New Hampshire Center for Coastal and Ocean Mapping / Joint Hydrographic Center (UNH/CCOM-JHC).

## Contact

- pjohnson@ccom.unh.edu  
- UNH/CCOM-JHC

## Acknowledgments

Developed at UNH/CCOM-JHC under grant NA25NOSX400C0001-T1-01 from NOAA.

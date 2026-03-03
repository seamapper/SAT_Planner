# SAT/QAT Planner

![Example Plot](media/SAT_Planner.jpg)

A comprehensive Shipboard Acceptance Testing (SAT) and Quality Assurance Testing (QAT) planning tool with GeoTIFF support, built with PyQt6. The GUI uses a consistent dark theme.

## Overview

The SAT/QAT Planner is a desktop application designed for planning and visualizing multibeam testing operations. It supports three main planning modes:
- **Calibration Survey Planning**: Plan pitch, roll, and heading calibration lines
- **Reference Survey Planning**: Generate parallel survey lines with customizable parameters
- **Line Planning**: Interactive line drawing with real-time elevation profiles

## Features

### Core Functionality
- **Multi-tab Interface**: Separate tabs for Calibration, Reference, and Line planning (left panel)
- **Dark Theme**: Qt GUI uses a dark theme; map (matplotlib) keeps default styling
- **GeoTIFF Support**: Load and visualize elevation data from GeoTIFF files
- **GMRT Download**: Optional download of GMRT bathymetry GeoTIFF when importing surveys (Calibration, Reference, Line tabs; configurable buffer)
- **Download GMRT GeoTIFF**: Button opens a separate "Download GMRT Grid" dialog; when split (topo/bathy) is used, SAT Planner loads the bathy grid
- **Dynamic Resolution**: Automatically adjust GeoTIFF resolution based on zoom level
- **Interactive Plotting**: Pan (middle mouse), zoom (scroll), no toolbar
- **Real-time Statistics**: Survey distances, times, and comprehensive statistics
- **Elevation Profiles**: View elevation and slope profiles for drawn lines
- **Activity Log**: Fixed on the right (320 px) below the map
- **Export**: CSV and Shapefile; statistics reports
- **Survey Import**: DDD, DMS, DMM, LNW, CSV, GeoJSON

### Calibration Survey Planning
- Draw pitch and roll calibration lines; generate heading lines from pitch line
- Turn time parameter; optional GMRT download on import (checkbox + buffer)
- Import from DDD/DMS/DMM/LNW, CSV, GeoJSON
- Validation when heading offset > 2× shallowest depth
- Export with detailed statistics

### Reference Survey Planning
- Parallel lines with line length, spacing, heading, speed, turn time
- Optional GMRT download on import; import from DDD/DMS/DMM/LNW, CSV, GeoJSON
- Survey time breakdown (main lines, crossline, transit, turn)
- Export with detailed statistics

### Line Planning
- Interactive line drawing; real-time elevation profiles; edit by dragging waypoints
- Optional GMRT download on import; import/export (DDD, DMS, DMM, LNW, CSV, GeoJSON)
- Survey statistics for drawn lines

### GeoTIFF Visualization
- Elevation/slope display; hillshade; dynamic resolution; multiple CRS
- Axis labels in degrees–decimal minutes (DDM)

### Download GMRT Grid Dialog
- **Download GMRT GeoTIFF** opens the "Download GMRT Grid" window (separate from main app)
- Set area (North/South/East/West or draw on map)
- Cell resolution: 100 m, 200 m, 400 m, or Custom (e.g. 50 m); default 100 m
- Optional split into bathymetry and topography; SAT Planner then loads the bathy grid
- Warning (orange) when estimated pixels > 16,000,000
- **Close GMRT Downloader** button at bottom
- Config: `~/.gmrtgrab_sat_planner_config.json`

## Requirements

- **Python**: 3.7+
- **Core**: PyQt6, matplotlib, numpy
- **Geospatial**: rasterio, pyproj, shapely, fiona
- **Optional**: Pillow (Imagery Basemap, NOAA ENC Charts), requests (GMRT download)

## Installation

### Pre-built Executable
- Download from [Releases](https://github.com/seamapper/SAT_Planner/releases): e.g. `SAT_Planner_v2026.08.exe` (Windows) or newer; see `sat_planner/constants.py` for version.

### From Source
1. Clone: `git clone https://github.com/seamapper/SAT_Planner.git && cd SAT_Planner`
2. Install: `pip install PyQt6 matplotlib numpy rasterio pyproj shapely fiona Pillow requests`
3. Run: `python SAT_Planner_PyQt.py`

## Building from Source

- **Windows**: `pip install pyinstaller` then run `build_exe.bat` (or use a PyInstaller spec); output in `dist/`.
- **macOS**: Use the provided macOS spec; icon `media/CCOM.icns` if available.

## Usage

1. **Load GeoTIFF** (optional): "Load GeoTIFF" or "Download GMRT GeoTIFF" (opens GMRT dialog; if split used, app loads bathy grid).
2. **Planning mode**: Calibration, Reference, or Line tab.
3. **Configure** parameters, generate or draw, view statistics, export.

### Map Controls (no toolbar)
- **Left click**: Add waypoint (drawing); pick points in Calibration/Reference/Line.
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
  - **`utils_geo.py`**, **`utils_ui.py`** – Helpers.
  - **`gmrt_dialog/`** – Embedded GMRT Download dialog (GeoTIFF-only, optional split).
  - **`mixins/`** – Basemap, GeoTIFF, Plotting, SurveyParsers, GMRTDownload, Reference, Calibration, LinePlanning, Profiles, MapInteraction, ExportImport, Config.

## Version History

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

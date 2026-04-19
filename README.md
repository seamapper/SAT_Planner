# SAT/QAT Planner

![Example Plot](media/SAT_Planner.jpg)

A comprehensive Shipboard Acceptance Testing (SAT) and Quality Assurance Testing (QAT) planning tool with GeoTIFF support, built with PyQt6.

## Overview

The SAT/QAT Planner is a desktop application designed for planning and visualizing multibeam testing operations. It supports four main planning modes:
- **Calibration Survey Planning**: Plan pitch, roll, and heading calibration lines
- **Accuracy Survey Planning**: Generate parallel survey lines with customizable parameters
- **Performance (Swath) Survey Planning**: Plan a four-heading swath performance test with optional RX-noise BIST legs relative to swell direction
- **Line Planning**: Interactive line drawing with real-time elevation profiles

## Features

### Core Functionality
- **Multi-tab Interface**: Separate tabs for Calibration, Accuracy, Performance, and Line planning (left panel)
- **Dark Theme**: Qt GUI always uses a dark theme; map (matplotlib) keeps default styling
- **GeoTIFF Support**: Load and visualize elevation data from GeoTIFF files
- **GMRT Download**: Optional download of GMRT bathymetry GeoTIFF when importing surveys (Calibration, Accuracy, Performance, Line tabs; configurable buffer)
- **Download Data**: Source selector in GeoTIFF Controls (default `Select Source`; current source `GMRT`) opens the "Download GMRT Grid" dialog immediately on selection and keeps the selected source after successful download
- **Dynamic Resolution**: Automatically adjust GeoTIFF resolution based on zoom level
- **Interactive Plotting**: Pan (middle mouse), zoom (scroll), and interact with survey plans on the map (no toolbar)
- **Real-time Statistics**: Calculate survey distances, times, and comprehensive statistics
- **Elevation Profiles**: View elevation and slope profiles for drawn lines
- **Activity Log**: Collapsible side panel below the map (expand/collapse strip)
- **Export Capabilities**: Export survey plans in CSV, Shapefile, GeoJSON, asciiplan, LNW, PNG, and companion text/statistics formats (varies by tab)
- **Survey Import**: Import calibration, accuracy, performance, and line plans from DDD, DMS, DMM, LNW, CSV, and GeoJSON (performance import uses an assignment dialog when geometry is ambiguous)
- **EEZ Overlay**: EEZ layer with opacity control (default 80%) and hover `GEONAME` tooltip lookup

### Calibration Survey Planning
- Draw pitch and roll calibration lines interactively
- Generate heading calibration lines from pitch line
- **Reverse Line Direction**: Flip start/end of any calibration line(s) (Pitch, Roll, Heading1, Heading2) via checkboxes
- Calculate heading line offset based on median depth
- Display pitch line depth statistics (shallowest, maximum, mean, median)
- Configure turn time for accurate time estimates
- Import calibration surveys (DDD/DMS/DMM/LNW, CSV, GeoJSON); **suggested line assignment** from file labels or geometry (Pitch = middle parallel line, Roll = non-parallel, Heading1/2 by file order); optional GMRT download after import (checkbox + buffer)
- GeoJSON import loads geometry from **`{name}.geojson`**; **`{name}_params.json`** in the same folder (if present) supplies **survey speed**, **turn time (min)**, **heading line offset**, and **export name** so you can edit those fields without changing the geometry file
- **Calibration Survey Info** dialog and *_info.txt with **Calibration Waypoints (DMM)** and **Calibration Waypoints (DDD)** sections (Pitch/Roll/Heading1/Heading2 start and end)
- Comprehensive statistics with survey time, transit time, and turn time breakdowns
- Validation warning when heading line offset exceeds 2x shallowest depth
- Export calibration survey plans with detailed statistics (shared DDD/DMM/DMS CSV and TXT, asciiplan, LNW via `sat_planner.export_utils`)

### Accuracy Survey Planning
- Generate parallel survey lines with customizable parameters
- Auto-regenerate plans when parameters change (with debounce)
- Configure line length, spacing, heading, speed, and turn time
- Import accuracy surveys (DDD/DMS/DMM/LNW, CSV, GeoJSON); **suggested crossline and reference line order** (crossline by orientation, reference lines in file order); optional GMRT download after import (checkbox + buffer)
- **Accuracy Survey Info** dialog and *_info.txt with **Accuracy Waypoints (DMM)** and **Accuracy Waypoints (DDD)** sections (L1S/L1E, L2S/L2E, …, CLS/CLE)
- Calculate comprehensive survey statistics with time breakdowns
- Survey time breakdown showing main lines, crossline, transit, and turn times
- Export accuracy survey plans with detailed statistics (shared DDD/DMM/DMS CSV and TXT, asciiplan, LNW via `sat_planner.export_utils`)

### Performance (Swath) Survey Planning
- **Purpose**: Support **swath performance** evaluation with **four legs** on headings **0°, 45°, 90°, and 135° relative to swell direction** (different aspects into/across/with/oblique to the seas). Each leg combines **swath collection** along the main segment (P1S–P4E) with **RX noise BIST** on the map extension from each line end when BIST time is non-zero.
- **Parameters**: Central lat/lon (**Pick Center from GeoTIFF** on Performance), **Swell Direction** (default **0°**), swath angle, sound velocity, test depth, pings, test speed, BIST time, turn time; **Line Length (m)** from speed and along-track collection time.
- **Plot Performance Lines**, **Zoom to Performance Lines**, **Remove Performance Lines**; **Show Performance Test Info** (pattern, legs, transits, times, waypoints DMM/DDD).
- **Auto-plot** (debounced) after editing test parameters or after performance pick-center when inputs are valid.
- **Profile** (Performance tab): line 1 swath + first BIST segment; **sienna** / **gold** to match the map.
- **Performance Import/Export**: same product family as Accuracy exports plus `{name}_performance_params.json`; optional **Download GMRT** on import; assignment dialog when import geometry is ambiguous. Default export basename: `Performance_<swath_angle>deg_<speed>_kts`.
- **Map markers**: **Perf Central Pt** vs **Acc Central Pt** (green accuracy center only when an accuracy plan is loaded).

### Line Planning
- Interactive line drawing with waypoint support
- **Reverse Line Direction**: Flip start and end of the line (one click)
- Real-time elevation profiles as you draw
- Edit existing lines by dragging waypoints
- Import/export line plans (DDD, DMS, DMM, LNW, CSV, GeoJSON; single polyline, no assignment dialog)
- Optional GMRT download after import (checkbox + buffer)
- **Survey Info** dialog and *_info.txt with **Line Plan Waypoints (DMM)** and **Line Plan Waypoints (DDD)** sections (WP1, WP2, …)
- Calculate survey statistics for drawn lines
- Export uses shared DDD/DMM/DMS CSV and TXT, asciiplan, LNW via `sat_planner.export_utils`

### GeoTIFF Visualization
- Display elevation data with color mapping
- Toggle between elevation and slope visualization
- Hillshade rendering for better terrain visualization
- Dynamic resolution loading for performance
- Support for various coordinate reference systems (CRS)
- Survey plan axis labels in degrees–decimal minutes (DDM)
- Contour interval and slope min/max entry updates are debounced while typing
- Survey legend is drawn above map overlays/layers
- EEZ overlay reloads on pan/zoom and supports paused-hover name lookup

### GeoJSON metadata (and sidecar JSON)
- **Calibration**: GeoJSON holds line geometry and `line_num` / `line_name` only. The FeatureCollection **`properties`** may include **`geotiff_path`** so a saved raster can be reopened on import. **Survey speed**, **turn time (min)**, **heading line offset**, and **export name** live in **`{export_name}_params.json`** next to the GeoJSON (not in the `.geojson`). On import, that sidecar is optional; if it is absent, the app uses defaults and may repopulate the heading offset from the pitch line after load.
- **Accuracy** and **Line** GeoJSON exports can include **`survey_speed`** and **`geotiff_path`** (collection and/or feature properties as applicable).
- **Performance** GeoJSON uses feature properties such as **`line_num`** (1–4 for swath legs, 11–14 for BIST extensions where applicable) and may include speed-related keys; full test parameters are also in **`{name}_performance_params.json`**.
- Whenever a saved **GeoTIFF** path is present in the imported files, SAT Planner tries to open it; if the file is missing, import continues with a warning.

## Requirements

### Python Version
- Python 3.7 or higher

### Core Dependencies
- PyQt6
- matplotlib
- numpy

### Geospatial Dependencies (Required for GeoTIFF support)
- rasterio
- pyproj
- shapely
- fiona

### Optional
- **Pillow (PIL)** – for Imagery Basemap and NOAA ENC Charts overlays
- **requests** – for GMRT bathymetry download (dialog and on import) – for GMRT bathymetry download when using “Download GMRT” on import

## Project structure

The application is organized as a package plus a launcher:

- **`SAT_Planner_PyQt.py`** – Entry point; creates the main window and runs the app (`python SAT_Planner_PyQt.py`).
- **`sat_planner/`** – Core package:
  - **`constants.py`** – Version, config path, geospatial library availability.
  - **`export_utils.py`** – Shared export writers: DDD/DMM/DMS CSV and TXT, SIS asciiplan, Hypack LNW; UTM zone from points (used by Calibration, Accuracy, Performance, Line planning).
  - **`performance_import_dialog.py`** – Assignment dialog for mapping imported segments to Performance swath lines 1–4 and optional BIST 1–4.
  - **`utils_geo.py`** – Coordinate helpers (e.g. decimal degrees to DDM).
  - **`utils_ui.py`** – UI helpers (message boxes, confirmations).
  - **`gmrt_dialog/`** – Embedded GMRT Download dialog (config, workers, map_widget, main_window); GeoTIFF-only output, optional split into topo/bathy.
  - **`mixins/`** – Feature mixins used by the main window:
    - **BasemapMixin** – Imagery basemap and NOAA ENC Charts overlays.
    - **GeoTIFFMixin** – Load/remove GeoTIFF, display mode, dynamic resolution, contours.
    - **PlottingMixin** – Survey plan plot, limits, colorbars, DDM axis labels.
    - **SurveyParsersMixin** – DDD/DMS/DMM/LNW parsers (lines and polylines), UTM zone dialog.
    - **GMRTDownloadMixin** – GMRT GridServer download and load GeoTIFF.
    - **ReferenceMixin** – Accuracy tab (reference/survey line planning), export/import, optional GMRT on import.
    - **CalibrationMixin** – Calibration tab, pitch/roll/heading lines, export/import, optional GMRT on import.
    - **LinePlanningMixin** – Line planning tab, draw/edit, profile, statistics, optional GMRT on import.
    - **PerformanceMixin** – Performance tab: ping/line-length math, pick center/depth, plot/info, import/export hooks, debounced auto-plot.
    - **ProfilesMixin** – Crossline, pitch, line-planning, and performance elevation profiles (performance: line 1 + BIST colors aligned with map).
    - **MapInteractionMixin** – Click, scroll, pan, zoom, pick center/pitch/roll.
    - **ExportImportMixin** – Save/load parameters, export accuracy survey files, export performance survey files.
    - **ConfigMixin** – Last-used directories, config load/save.

## Installation

### Option 1: Using Pre-built Executable

Download the latest executable from the [Releases](https://github.com/seamapper/SAT_Planner/releases) page:
- `SAT_Planner_v2026.14.exe` (Windows) or newer — version is in the filename (see `sat_planner/constants.py`).
- `SAT_Planner.app` (macOS) — if available

No installation required - just run the executable or app bundle.

### Option 2: From Source

1. Clone the repository:
```bash
git clone https://github.com/seamapper/SAT_Planner.git
cd SAT_Planner
```

2. Install dependencies:

**Using pip:**
```bash
pip install PyQt6 matplotlib numpy rasterio pyproj shapely fiona Pillow requests
```
(`requests` is used for GMRT bathymetry download.)

**Using conda (recommended for Windows and macOS):**
```bash
conda install -c conda-forge pyqt matplotlib numpy rasterio pyproj shapely fiona pillow
```

3. Run the application:
```bash
python SAT_Planner_PyQt.py
```

## Building from Source

### Building for Windows

1. Install PyInstaller:
```bash
pip install pyinstaller
```

2. Run the build script (edit `build_exe.bat` to set `PYTHON_PATH` if needed):
```bash
build_exe.bat
```

Or build manually: `pyinstaller SAT_Planner.spec`. The executable is created in `dist` as `SAT_Planner_v<version>.exe` (version from `sat_planner/constants.py`).

### Building for macOS

To build a macOS application (.app bundle), create a PyInstaller spec (one-file or one-folder) that includes the `sat_planner` package and `SAT_Planner_PyQt.py` as the entry point. If a macOS spec file is present in the repo (e.g. `Sat_Planner_macOS.spec`), run:

```bash
pip install pyinstaller
pyinstaller Sat_Planner_macOS.spec
```

Use `.icns` for the app icon (convert `media/CCOM.ico` with `sips -s format icns media/CCOM.ico --out media/CCOM.icns` if needed). The app bundle will be in `dist/`. Code signing and DMG packaging are optional for distribution.

## Usage

### Basic Workflow

1. **Load GeoTIFF** (optional): Click "Load GeoTIFF" to load elevation data, or use **Download Data** -> **GMRT** to open the GMRT dialog and download bathymetry (GeoTIFF-only; if you use split, the app loads the bathy grid)
2. **Enable Map Overlays** (optional): Toggle Imagery Basemap or NOAA ENC Charts checkboxes and adjust opacity sliders
3. **Select Planning Mode**: Choose Calibration, Accuracy, Performance, or Line tabs
4. **Configure Parameters**: Set survey parameters in the appropriate tab
5. **Generate/Plan**: Create survey lines based on parameters or draw interactively
6. **View Statistics / test info**: Open **Calibration Survey Info**, **Accuracy Survey Info**, **Show Performance Test Info** (Performance tab), or **Survey Info** (Line tab)
7. **Export**: Save survey plans using the Export buttons

### Calibration Survey Planning

1. Load a GeoTIFF (recommended)
2. Click "Draw a Pitch Line" and click start/end points on the map
3. Configure Turn Time (min) in Calibration Info - default: 5 minutes
4. Click "Add Heading Lines" to generate heading calibration lines (warning shown if offset > 2x shallowest depth)
5. Click "Draw a Roll Line" and click start/end points
6. View comprehensive statistics showing survey time, transit time, and turn time breakdowns
7. Export as needed (statistics file includes all details from dialog)

### Accuracy Survey Planning

1. Optionally load a GeoTIFF and pick a center point (view stays at current zoom)
2. Enter survey parameters:
   - Central Latitude/Longitude
   - Number of Lines
   - Line Length
   - Heading
   - Distance Between Lines
   - Survey Speed
   - Turn Time (min) - default: 5 minutes
3. The plan auto-regenerates as you change parameters
4. View comprehensive statistics with survey time, transit time, and turn time breakdowns
5. Export the plan

### Performance (Swath) Survey Planning

1. Load a GeoTIFF (needed for plotting and for depth at pick).
2. Open the **Performance** tab; set **Swell Direction** (default 0°) and other test parameters (speed, pings, BIST time, etc.).
3. Use **Pick Center from GeoTIFF** and click the map, or enter central latitude/longitude manually. With valid depth/speed/pings, the plan may **auto-plot**; otherwise click **Plot Performance Lines**.
4. Use **Show Performance Test Info** for a written pattern, distances, and times; use **Performance Import/Export** to share plans in the same family of formats as Accuracy.

### Line Planning

1. Load a GeoTIFF (recommended)
2. Click "Start Drawing Line"
3. Left-click to add waypoints
4. Right-click to finish the line
5. Edit by clicking "Edit Line Planning" and dragging waypoints
6. View elevation profile and statistics
7. Export the line plan

### Download GMRT Grid Dialog

Use **Download Data** -> **GMRT** in the GeoTIFF Control section to open the **Download GMRT Grid** window (separate from the main app). In the dialog you can:

- Set the area of interest (North/South/East/West or draw on the map)
- Choose **Cell Resolution**: 100 m, 200 m, 400 m, or **Custom** (e.g. 50 m); default is 100 m
- **Split Grid Into Bathymetry and Topography** is enabled by default; SAT Planner loads the bathymetry grid after download
- A warning (orange text) appears when estimated pixels exceed 16,000,000
- Download progress displays in-dialog (tile mode: `x of y`; single-grid mode: indeterminate "Downloading...")
- Use **"Close GMRT Downloader"** at the bottom to close the dialog

Downloads are always GeoTIFF. The dialog uses its own config: `~/.gmrtgrab_sat_planner_config.json`.

### Map Overlays

1. **Imagery Basemap**: Check the "Imagery Basemap" checkbox to overlay satellite imagery
   - Adjust opacity using the slider (0-100%)
   - The basemap updates automatically as you pan and zoom
   - Large-area/global views use wrap-aware alignment for correct placement

2. **NOAA ENC Charts**: Check the "NOAA ENC Charts" checkbox to overlay navigational charts
   - Adjust opacity using the slider (0-100%)
   - Charts are requested with equal cell size and properly reprojected for display
   - Charts update automatically as you pan and zoom

## Configuration

The application saves configuration in:
- `~/.cal_ref_planner_config.json` (user preferences)
- Last used directories
- Survey parameters

## Export Formats

- **CSV (DDD/DMM/DMS)**: Survey lines with coordinates and metadata (row format: line number, line name, point label, lat, lon)
- **TXT (DDD/DMM/DMS)**: Same coordinate formats as plain text
- **SIS asciiplan**: ASCII plan format for SIS
- **Hypack LNW**: LNW format for Hypack (UTM zone computed from survey points when needed)
- **Shapefile**: Geospatial vector format for GIS applications
- **Statistics / *_info.txt**: Text reports with distances, timing, survey details, and waypoint sections (DMM then DDD) for Calibration, Accuracy, and Line plans; Performance exports add a performance-oriented `*_info.txt` and **`{name}_performance_params.json`** (separate from Accuracy `{name}_params.json`)
- **Calibration `{name}_params.json`**: Sidecar to the calibration GeoJSON export—**`survey_speed`**, **`turn_time`**, **`line_offset`** (heading line offset), and **`export_name`**. Edit this file to tune those values and re-import the `.geojson` without regenerating geometry.

## Navigation

### Mouse Controls (no toolbar)
- **Left Click**: Add waypoint (in drawing mode); pick points in calibration/reference/line modes
- **Right Click**: Finish drawing line
- **Middle Mouse Button**: Pan the map
- **Scroll Wheel**: Zoom in/out

## Troubleshooting

### GeoTIFF Loading Issues

- Ensure the GeoTIFF file is valid and not corrupted
- Check that the file uses a supported CRS
- Large files may take time to load - be patient
- Try enabling "Dynamic Resolution" for better performance

### Import Errors

If you encounter `ModuleNotFoundError` for geospatial libraries:

1. **Install dependencies**: Run `pip install rasterio pyproj shapely fiona`
2. **Windows users**: Consider using conda or OSGeo4W for easier installation
3. **Check Python environment**: Ensure you're using the correct Python environment

### Application Runs Without GeoTIFF Support

The application will run with limited functionality if geospatial libraries aren't available:
- Basic survey planning still works
- GeoTIFF loading will be disabled
- A warning message will be displayed

### Performance Issues

- Large GeoTIFF files may cause slow loading
- Enable "Dynamic Resolution" to improve performance
- Consider using smaller tiles or downsampled data
- Close other applications to free up memory

## Version History

- **v2026.20**: README: document **calibration** export/import split—**`{name}_params.json`** holds survey speed, turn time, heading line offset, and export name; **calibration GeoJSON** is geometry plus line labels (and optional **`geotiff_path`** on the FeatureCollection). Clarified GeoJSON metadata behavior for Accuracy, Line, and Performance.
- **v2026.15**: **Performance (swath) survey planning** tab: four headings relative to swell, swath lines plus optional **RX noise BIST** extensions; **Plot Performance Lines**, zoom/remove, **Show Performance Test Info**; debounced **auto-plot** after parameter edits and after performance pick-center; profile for **line 1 + BIST** with map-matched colors; **Performance Import/Export** (same export family as Accuracy, assignment dialog on ambiguous import, optional GMRT on import); default swell direction **0°**; accuracy vs performance **central point** markers fixed so accuracy center shows only when an accuracy plan is present; UI labels (**Plot Performance Lines**, button layout). Added `performance_import_dialog.py` and **PerformanceMixin**; README updated for Performance workflow.
- **v2026.14**: Added GeoTIFF Controls **Download Data** source dropdown (`Select Source`, `GMRT`) that opens the source flow immediately and keeps selection after successful download. Added GMRT dialog download progress bar (`x of y` for tiled downloads, indeterminate for single download). Added EEZ pan/zoom-driven refresh, paused-hover EEZ `GEONAME` tooltip lookups (including without a loaded GeoTIFF), and default EEZ opacity 80%. Improved large-area/world-scale alignment for EEZ and Imagery Basemap overlays. Updated GMRT dialog so **Split Grid Into Bathymetry and Topography** is enabled by default. Fixed line-plan zoom GeoTIFF coverage refresh, line-plan profile refresh after import, and calibration Pitch Line Info refresh after survey import.
- **v2026.11**: UI and workflow updates. Accuracy tab naming in UI/docs (formerly Reference). Import/Export button labels updated (Import/Export Accuracy Survey, Import/Export Calibration Survey, Import/Export Line Survey) and reordered on tabs. "Download GMRT" defaults to unchecked and placement updated in import/export groups. Activity Log changed to collapsible side panel. GeoJSON export/import now includes `survey_speed` and saved `geotiff_path` for Calibration/Accuracy/Line; missing GeoTIFF path does not block import (warning + continue). Contour interval and slope min/max redraws are debounced while typing. Survey legend now draws above all overlays.
- **v2026.10**: Shared export utilities (`sat_planner/export_utils.py`): DDD/DMM/DMS CSV and TXT, SIS asciiplan, Hypack LNW; UTM zone from points. Calibration/Reference/Line planning use these writers. Calibration: **Reverse Line Direction** (Pitch/Roll/Heading1/Heading2); **import suggestion** from metadata (PLS/PLE, RLS/RLE, H1S/H1E, H2S/H2E) or geometry (Pitch = middle parallel, Roll = fourth line, Heading1/2 = outer two). Reference: **import suggestion** (crossline + reference lines). Line planning: **Reverse Line Direction**. Survey info dialogs: **Calibration Survey Info**, **Reference Survey Info**, **Survey Info** (Line); *_info.txt waypoints as **Calibration Waypoints (DMM/DDD)**, **Reference Waypoints (DMM/DDD)**, **Line Plan Waypoints (DMM/DDD)**. Exe build uses version from `sat_planner/constants.py`; output `dist/SAT_Planner_v2026.10.exe`.
- **v2026.09** (or later): Integrated GMRT Download dialog: "Download GMRT GeoTIFF" button opens a separate "Download GMRT Grid" window. Dialog is GeoTIFF-only (no output format selector); cell resolution 100/200/400 m or Custom (default 50 m), default preset 100 m. Large-area warning when estimated pixels > 16,000,000 (orange). When split (topo/bathy) is used, SAT Planner loads the bathy grid. "Close GMRT Downloader" button at bottom of dialog. Activity Log width 320 px.
- **v2026.08**: Dark theme for Qt GUI (Fusion + dark palette). Activity Log moved to right side below map (380 px wide). GMRT download option on Calibration, Reference, and Line import (checkbox + buffer). Line plan import from DDD/DMS/DMM/LNW/CSV/GeoJSON. Navigation toolbar removed; zoom (scroll) and pan (middle mouse) only. Survey parsers and GMRT download in dedicated mixins.
- **v2026.04**: Updated hover text to display coordinates in degrees and decimal minutes (DDM) format instead of decimal degrees. Changed default window height to 1110 pixels.
- **v2026.02**: Added Turn Time parameter to Calibration and Reference Info tabs. Enhanced statistics displays with Total Survey Time and Total Transit Time breakdowns. Fixed autozoom issue when picking center from GeoTIFF. Added validation warning when heading line offset exceeds 2x shallowest depth. Updated export functions to include comprehensive statistics matching dialog displays.
- **v2026.01**: Refactored into `sat_planner` package with mixins (Basemap, GeoTIFF, Plotting, Reference, Calibration, Line Planning, Profiles, Map Interaction, Export/Import, Config). Survey plan axes show DDM (degrees–decimal minutes) tick labels. Moved basemap/NOAA and geotiff/plotting helpers into mixins.
- **v2025.11**: Fixed Dynamic Resolution for toolbar zoom/pan operations, updated About this Program dialog
- **v2025.10**: Added Imagery Basemap and NOAA ENC Charts overlays with opacity controls, navigation toolbar at bottom of map, fixed Dynamic Resolution for toolbar zoom/pan, improved map visualization
- **v2025.09**: Made profile colors coordinate with survey plot
- **v2025.08**: Added About button to profile plot
- **v2025.07**: Added ability to plan all tests simultaneously, fixed profile plot updates
- **v2025.06**: Fixed labeling of waypoints in line planning tab, fixed preservation of lines when changing tabs
- **v2025.05**: Added ability to plan all tests simultaneously, fixed profile plot updates, improved dynamic resolution
- **v2025.04**: Converted to PyQt6
- **v2025.03**: Added line planning, import/export of lines
- **v2025.02**: Added metadata file saving/loading, contour interval synchronization
- **v2025.01**: Initial release

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](LICENSE) file for details.

Copyright (c) 2025, University of New Hampshire Center for Coastal and Ocean Mapping / Joint Hydrographic Center (UNH/CCOM-JHC)

## Contact

For questions or issues, please contact:
- Email: pjohnson@ccom.unh.edu
- Organization: UNH/CCOM-JHC

## GMRT

Bathymetry/topography data available through the Download GMRT Grid dialog and the GMRT download-on-import options are from the Global Multi-Resolution Topography (GMRT) synthesis:

Ryan, W. B. F., S.M. Carbotte, J. Coplan, S. O'Hara, A. Melkonian, R. Arko, R.A. Weissel, V. Ferrini, A. Goodwillie, F. Nitsche, J. Bonczkowski, and R. Zemsky (2009), Global Multi-Resolution Topography (GMRT) synthesis data set, *Geochem. Geophys. Geosyst.*, 10, Q03014, doi:[10.1029/2008GC002332](https://doi.org/10.1029/2008GC002332).

Data doi: [10.1594/IEDA.100001](https://doi.org/10.1594/IEDA.100001).

## Acknowledgments

Developed at the University of New Hampshire, Center for Coastal and Ocean Mapping - Joint Hydrographic Center (UNH/CCOM-JHC) under grant NA25NOSX400C0001-T1-01 from the National Oceanic and Atmospheric Administration (NOAA).


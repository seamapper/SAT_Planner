# SAT/QAT Planner

A comprehensive survey planning application for multibeam survey operations with GeoTIFF support, built with PyQt6.

## Overview

The SAT/QAT Planner is a desktop application designed for planning and visualizing multibeam survey operations. It supports three main planning modes:
- **Calibration Survey Planning**: Plan pitch, roll, and heading calibration lines
- **Reference Survey Planning**: Generate parallel survey lines with customizable parameters
- **Line Planning**: Interactive line drawing with real-time elevation profiles

## Features

### Core Functionality
- **Multi-tab Interface**: Separate tabs for Calibration, Reference, and Line planning
- **GeoTIFF Support**: Load and visualize elevation data from GeoTIFF files
- **Dynamic Resolution**: Automatically adjust GeoTIFF resolution based on zoom level
- **Interactive Plotting**: Pan, zoom, and interact with survey plans on the map
- **Real-time Statistics**: Calculate survey distances, times, and comprehensive statistics
- **Elevation Profiles**: View elevation and slope profiles for drawn lines
- **Export Capabilities**: Export survey plans in CSV and Shapefile formats

### Calibration Survey Planning
- Draw pitch and roll calibration lines interactively
- Generate heading calibration lines from pitch line
- Calculate heading line offset based on median depth
- Display pitch line depth statistics (shallowest, maximum, mean, median)
- Export calibration survey plans

### Reference Survey Planning
- Generate parallel survey lines with customizable parameters
- Auto-regenerate plans when parameters change (with debounce)
- Configure line length, spacing, heading, and speed
- Calculate comprehensive survey statistics
- Export reference survey plans

### Line Planning
- Interactive line drawing with waypoint support
- Real-time elevation profiles as you draw
- Edit existing lines by dragging waypoints
- Import/export line plans
- Calculate survey statistics for drawn lines

### GeoTIFF Visualization
- Display elevation data with color mapping
- Toggle between elevation and slope visualization
- Hillshade rendering for better terrain visualization
- Dynamic resolution loading for performance
- Support for various coordinate reference systems (CRS)

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

## Installation

### Option 1: Using Pre-built Executable

Download the latest executable from the [Releases](https://github.com/seamapper/SAT_Planner/releases) page:
- `Sat_Planner_v2025.05.exe` (Windows)

No installation required - just run the executable.

### Option 2: From Source

1. Clone the repository:
```bash
git clone https://github.com/seamapper/SAT_Planner.git
cd SAT_Planner
```

2. Install dependencies:

**Using pip:**
```bash
pip install PyQt6 matplotlib numpy rasterio pyproj shapely fiona
```

**Using conda (recommended for Windows):**
```bash
conda install -c conda-forge pyqt matplotlib numpy rasterio pyproj shapely fiona
```

3. Run the application:
```bash
python SAT_Planner_PyQt.py
```

## Building from Source

To build your own executable:

1. Install PyInstaller:
```bash
pip install pyinstaller
```

2. Run the build script:
```bash
build_exe.bat
```

Or manually:
```bash
pyinstaller Sat_Planner_v2025.05.spec
```

The executable will be created in the `dist` folder.

## Usage

### Basic Workflow

1. **Load GeoTIFF** (optional): Click "Load GeoTIFF" to load elevation data for visualization
2. **Select Planning Mode**: Choose between Calibration, Reference, or Line tabs
3. **Configure Parameters**: Set survey parameters in the appropriate tab
4. **Generate/Plan**: Create survey lines based on parameters or draw interactively
5. **View Statistics**: Click "Show [Type] Test Info" to view survey statistics
6. **Export**: Save survey plans using the Export buttons

### Calibration Survey Planning

1. Load a GeoTIFF (recommended)
2. Click "Draw a Pitch Line" and click start/end points on the map
3. Click "Add Heading Lines" to generate heading calibration lines
4. Click "Draw a Roll Line" and click start/end points
5. View statistics and export as needed

### Reference Survey Planning

1. Optionally load a GeoTIFF and pick a center point
2. Enter survey parameters:
   - Central Latitude/Longitude
   - Number of Lines
   - Line Length
   - Heading
   - Distance Between Lines
   - Survey Speed
3. The plan auto-regenerates as you change parameters
4. View statistics and export the plan

### Line Planning

1. Load a GeoTIFF (recommended)
2. Click "Start Drawing Line"
3. Left-click to add waypoints
4. Right-click to finish the line
5. Edit by clicking "Edit Line Planning" and dragging waypoints
6. View elevation profile and statistics
7. Export the line plan

## Configuration

The application saves configuration in:
- `~/.cal_ref_planner_config.json` (user preferences)
- Last used directories
- Survey parameters

## Export Formats

- **CSV**: Survey lines with coordinates and metadata
- **Shapefile**: Geospatial vector format for GIS applications
- **Statistics Reports**: Text reports with distances, timing, and survey details

## Keyboard Shortcuts

- **Left Click**: Add waypoint (in drawing mode)
- **Right Click**: Finish drawing line
- **Pan**: Click and drag on the plot
- **Zoom**: Use mouse wheel or zoom buttons

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

- **v2025.05**: Added ability to plan all tests simultaneously, fixed profile plot updates, improved dynamic resolution
- **v2025.04**: Converted to PyQt6
- **v2025.03**: Added line planning, import/export of lines
- **v2025.02**: Added metadata file saving/loading, contour interval synchronization
- **v2025.01**: Initial release

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

[Add your license information here]

## Contact

For questions or issues, please contact:
- Email: pjohnson@ccom.unh.edu
- Organization: UNH/CCOM-JHC

## Acknowledgments

Developed at the University of New Hampshire Center for Coastal and Ocean Mapping / Joint Hydrographic Center (UNH/CCOM-JHC).


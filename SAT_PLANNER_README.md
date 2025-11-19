# SAT QAT Planner

A comprehensive survey planning application with GeoTIFF support for multibeam survey operations.

## Features

- **Survey Plan Generation**: Create and visualize survey lines with customizable parameters
- **GeoTIFF Support**: Load and display elevation data from GeoTIFF files
- **Line Planning**: Interactive line drawing with elevation profiles
- **Calibration Survey Planning**: Plan pitch, roll, and heading calibration lines
- **Export Capabilities**: Export survey plans in various formats
- **Real-time Statistics**: Calculate survey distances, times, and statistics

## Installation

### Prerequisites

The application requires Python 3.7+ and the following core dependencies:
- tkinter (usually included with Python)
- matplotlib
- numpy

### Geospatial Dependencies (Optional but Recommended)

For full GeoTIFF and Shapefile support, install the geospatial dependencies:

#### Option 1: Automatic Installation
```bash
# Run the installation script
python install_geospatial_deps.py

# Or on Windows, use the batch file
install_geospatial_deps.bat
```

#### Option 2: Manual Installation
```bash
# Using pip
pip install rasterio pyproj shapely fiona

# Using conda (recommended for Windows)
conda install -c conda-forge gdal rasterio pyproj shapely fiona
```

#### Option 3: Windows-specific
If you encounter issues on Windows, consider using:
- OSGeo4W distribution
- Conda-forge channel
- Pre-compiled wheels from Christoph Gohlke's repository

## Usage

### Running the Application

```bash
python multibeam_tools/apps/sat_qat_planner_01.py
```

### Basic Workflow

1. **Load GeoTIFF** (optional): Load elevation data for visualization
2. **Set Parameters**: Configure survey parameters in the tabs
3. **Generate Survey Plan**: Create survey lines based on parameters
4. **Interactive Planning**: Use line planning mode for custom routes
5. **Export Results**: Save survey plans and statistics

### Key Features

#### Survey Plan Generation
- Configure line length, spacing, and speed
- Generate parallel survey lines
- Calculate survey statistics and timing

#### Line Planning Mode
- Interactive line drawing on the map
- Real-time elevation profiles
- Export drawn lines for survey execution

#### Calibration Survey Planning
- Plan pitch and roll calibration lines
- Generate heading calibration lines
- Calculate comprehensive calibration statistics

#### GeoTIFF Visualization
- Display elevation data with color mapping
- Toggle between elevation and slope visualization
- Hillshade rendering for better terrain visualization

## Troubleshooting

### Import Errors

If you encounter `ModuleNotFoundError` for geospatial libraries:

1. **Install dependencies**: Run `python install_geospatial_deps.py`
2. **Check Python environment**: Ensure you're using the correct Python environment
3. **Windows users**: Consider using conda or OSGeo4W for easier installation

### Application Runs Without GeoTIFF Support

The application will run with limited functionality if geospatial libraries aren't available:
- Basic survey planning still works
- GeoTIFF loading will be disabled
- A warning message will be displayed

### Performance Issues

- Large GeoTIFF files may cause slow loading
- Consider using smaller tiles or downsampled data
- Close other applications to free up memory

## File Structure

```
multibeam_tools/apps/
├── sat_qat_planner_01.py    # Main application
├── install_geospatial_deps.py  # Dependency installer
└── install_geospatial_deps.bat # Windows batch installer
```

## Configuration

The application saves configuration in:
- `~/.cal_ref_planner_config.json` (user preferences)
- Last used directories
- Survey parameters

## Export Formats

- **Survey Lines**: CSV format with coordinates
- **Statistics**: Text reports with distances and timing
- **Profiles**: Elevation data along survey lines
- **Calibration Plans**: Comprehensive calibration survey data

## Technical Details

### Dependencies
- **Core**: tkinter, matplotlib, numpy
- **Geospatial**: rasterio, pyproj, shapely, fiona
- **Optional**: GDAL (for advanced GeoTIFF features)

### Architecture
- **GUI**: Tkinter-based interface with matplotlib integration
- **Data Processing**: NumPy for numerical operations
- **Geospatial**: Rasterio for GeoTIFF handling, PyProj for coordinate transformations
- **Visualization**: Matplotlib for plotting and interactive features

## Recent Fixes

### Import Error Resolution
- **Issue**: `ModuleNotFoundError: No module named 'rasterio'` on startup
- **Solution**: Moved all geospatial imports into conditional import blocks
- **Result**: Application now starts successfully even without geospatial libraries
- **Benefit**: Graceful degradation - runs with limited functionality when libraries are missing

### Installation Improvements
- Created automated installation script for geospatial dependencies
- Added Windows batch file for easy installation
- Provided multiple installation options for different environments

## Support

For issues or questions:
1. Check the troubleshooting section above
2. Ensure all dependencies are properly installed
3. Verify Python environment and version compatibility
4. Check that GeoTIFF files are valid and accessible

## License

This application is part of the MultibeamToolsAI project. See the main project LICENSE file for details. 
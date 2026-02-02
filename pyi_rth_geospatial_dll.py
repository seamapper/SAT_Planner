# Runtime hook for PyInstaller one-file exe on Windows.
# Add the extraction dir to the DLL search path so rasterio/fiona (GDAL) can load their native DLLs.
# Must run before any geospatial imports.
import os
import sys

if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
    try:
        os.add_dll_directory(sys._MEIPASS)
    except (AttributeError, OSError):
        pass  # add_dll_directory is 3.8+; ignore on older Python or if path invalid

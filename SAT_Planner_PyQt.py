"""
UNH/CCOM-JHC Shipboard Acceptance Testing (SAT) and Quality Assurance Testing (QAT) Planner.
Launcher: runs the SAT Planner application.

Program by Paul Johnson, pjohnson@ccom.unh.edu
Center for Coastal and Ocean Mapping/Joint Hydrographic Center, University of New Hampshire

This program was developed at UNH/CCOM-JHC under the grant NA20NOS4000196 from NOAA.

Copyright (c) 2025, UNH/CCOM-JHC. BSD 3-Clause License (see LICENSE in repo root).
"""
import sys
import os

# Fledermaus ships DLLs that collide with conda rasterio/GDAL. Drop any
# Fledermaus folder from this process only so the system PATH stays unchanged.
def _is_fledermaus_path(entry):
    parts = os.path.normcase(os.path.normpath(entry)).split(os.sep)
    return "fledermaus" in parts

os.environ["PATH"] = os.pathsep.join(
    entry for entry in os.environ.get("PATH", "").split(os.pathsep)
    if entry and not _is_fledermaus_path(entry)
)

# Import geospatial libs here so PyInstaller (one-file exe) traces and bundles them.
# They are used via sat_planner.constants; without this, the exe misses them after refactoring.
try:
    import rasterio
    import pyproj
    import shapely
    import fiona
except ImportError:
    pass
from PyQt6.QtWidgets import QApplication
from PyQt6.QtGui import QIcon
# Qt WebEngine must be imported before the QApplication object is created.
import PyQt6.QtWebEngineWidgets  # noqa: F401
from sat_planner import SurveyPlanApp
from sat_planner.utils_ui import apply_dark_theme, media_path

if __name__ == "__main__":
    app = QApplication(sys.argv)
    apply_dark_theme(app)
    icon_path = media_path("CCOM.ico")
    if os.path.exists(icon_path):
        app.setWindowIcon(QIcon(icon_path))
    window = SurveyPlanApp()
    window.show()
    sys.exit(app.exec())

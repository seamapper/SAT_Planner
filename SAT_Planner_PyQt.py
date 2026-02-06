"""
UNH/CCOM-JHC Shipboard Acceptance Testing (SAT) and Quality Assurance Testing (QAT) Planner.
Launcher: runs the SAT Planner application.

Program by Paul Johnson, pjohnson@ccom.unh.edu
Center for Coastal and Ocean Mapping/Joint Hydrographic Center, University of New Hampshire

This program was developed at UNH/CCOM-JHC under the grant NA20NOS4000196 from NOAA.

Copyright (c) 2025, UNH/CCOM-JHC. BSD 3-Clause License (see LICENSE in repo root).
"""
import sys
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
from sat_planner import SurveyPlanApp

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = SurveyPlanApp()
    window.show()
    sys.exit(app.exec())

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


def _configure_tls_ca_bundle():
    """Point HTTPS clients at a CA bundle that is actually on disk.

    Conda-forge certifi does not ship cacert.pem. It looks under
    sys.base_prefix, which is the conda install when running from source and
    the unpacked app folder when frozen. requests then raises "Could not find
    a suitable TLS CA certificate bundle" if that path is missing.
    """
    roots = []
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        roots.append(sys._MEIPASS)
    roots.append(getattr(sys, "base_prefix", ""))
    roots.append(getattr(sys, "prefix", ""))
    candidates = []
    for key in ("SSL_CERT_FILE", "REQUESTS_CA_BUNDLE", "CURL_CA_BUNDLE"):
        value = os.environ.get(key)
        if value:
            candidates.append(value)
    for root in roots:
        if not root:
            continue
        candidates.append(os.path.join(root, "Library", "ssl", "cacert.pem"))
        candidates.append(os.path.join(root, "ssl", "cacert.pem"))
    chosen = next((path for path in candidates if path and os.path.isfile(path)), None)
    if not chosen:
        return
    for key in ("SSL_CERT_FILE", "REQUESTS_CA_BUNDLE", "CURL_CA_BUNDLE"):
        current = os.environ.get(key)
        if not current or not os.path.isfile(current):
            os.environ[key] = chosen


_configure_tls_ca_bundle()

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
# Qt WebEngine must be imported before the QApplication object is created.
import PyQt6.QtWebEngineWidgets  # noqa: F401
from sat_planner import SurveyPlanApp
from sat_planner.utils_ui import application_icon, apply_dark_theme

if __name__ == "__main__":
    app = QApplication(sys.argv)
    apply_dark_theme(app)
    icon = application_icon()
    if not icon.isNull():
        app.setWindowIcon(icon)
    window = SurveyPlanApp()
    window.show()
    sys.exit(app.exec())

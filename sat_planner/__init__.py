"""
SAT Planner package.
UNH/CCOM-JHC Shipboard Acceptance Testing (SAT) and Quality Assurance Testing (QAT) Planner.
"""
from .constants import __version__, CONFIG_FILENAME, GEOSPATIAL_LIBS_AVAILABLE
from .utils_geo import decimal_degrees_to_ddm

__all__ = [
    "__version__",
    "CONFIG_FILENAME",
    "GEOSPATIAL_LIBS_AVAILABLE",
    "decimal_degrees_to_ddm",
]

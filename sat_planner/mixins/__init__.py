"""
Mixins for SurveyPlanApp. Each mixin holds a logical group of methods.
SurveyPlanApp will inherit from these (e.g. GeoTIFFMixin, PlottingMixin, ReferenceMixin, ...).
"""
from .geotiff_mixin import GeoTIFFMixin
from .plotting_mixin import PlottingMixin
from .reference_mixin import ReferenceMixin
from .calibration_mixin import CalibrationMixin
from .line_planning_mixin import LinePlanningMixin
from .profiles_mixin import ProfilesMixin
from .map_interaction_mixin import MapInteractionMixin

__all__ = ["GeoTIFFMixin", "PlottingMixin", "ReferenceMixin", "CalibrationMixin", "LinePlanningMixin", "ProfilesMixin", "MapInteractionMixin"]

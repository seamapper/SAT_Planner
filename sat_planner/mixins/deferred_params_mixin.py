"""
Mixin: commit-on-Enter/blur for parameter line edits with pending (amber) styling.
"""
from __future__ import annotations

from typing import Callable, Dict, List, Optional

from PyQt6.QtWidgets import QLineEdit

from ..deferred_params import (
    bind_deferred_line_edit,
    commit_all_deferred_line_edits,
    commit_deferred_line_edit,
    deferred_mark_applied,
    deferred_set_line_edit,
)


class DeferredParamsMixin:
    """Track applied vs current values for bound QLineEdit parameter fields."""

    def _init_deferred_params(self) -> None:
        if getattr(self, "_deferred_params_initialized", False):
            return
        self._deferred_params_initialized = True
        self._deferred_applied_values: Dict[int, str] = {}
        self._deferred_bound_widgets: List[QLineEdit] = []
        self._deferred_commit_handlers: Dict[int, Callable[[QLineEdit], None]] = {}

    def _bind_deferred_param(
        self,
        widget: QLineEdit,
        on_commit: Optional[Callable[[QLineEdit], None]] = None,
    ) -> None:
        self._init_deferred_params()
        self._deferred_commit_handlers[id(widget)] = on_commit or (lambda _w: None)
        bind_deferred_line_edit(
            widget,
            self._deferred_applied_values,
            on_commit=on_commit,
            bound_registry=self._deferred_bound_widgets,
        )
        self._deferred_mark_applied(widget)

    def _deferred_set_line_edit(
        self,
        widget: QLineEdit,
        text: str,
        *,
        mark_applied: bool = True,
    ) -> None:
        self._init_deferred_params()
        deferred_set_line_edit(
            widget,
            text,
            self._deferred_applied_values,
            mark_applied=mark_applied,
        )

    def _deferred_mark_applied(self, widget: QLineEdit) -> None:
        self._init_deferred_params()
        deferred_mark_applied(widget, self._deferred_applied_values)

    def _commit_deferred_param(self, widget: QLineEdit) -> None:
        self._init_deferred_params()
        handler = self._deferred_commit_handlers.get(id(widget))
        commit_deferred_line_edit(widget, self._deferred_applied_values, on_commit=handler)

    def _commit_all_deferred_line_edits(self) -> None:
        self._init_deferred_params()
        commit_all_deferred_line_edits(
            self._deferred_bound_widgets,
            self._deferred_applied_values,
            self._deferred_commit_handlers,
        )

    def _deferred_sync_all_bound_params(self) -> None:
        """After import/programmatic fill, treat all bound fields as applied (no pending style)."""
        self._init_deferred_params()
        for widget in self._deferred_bound_widgets:
            self._deferred_mark_applied(widget)

    def _apply_performance_param_commit(self, _widget=None) -> None:
        if hasattr(self, "_update_performance_ping_time"):
            self._update_performance_ping_time()
        if hasattr(self, "_update_performance_total_test_time"):
            self._update_performance_total_test_time()
        if hasattr(self, "_update_performance_line_length"):
            self._update_performance_line_length()
        if hasattr(self, "_update_performance_export_name"):
            self._update_performance_export_name()
        if hasattr(self, "_run_autoplot_performance_test_lines"):
            self._run_autoplot_performance_test_lines()

    def _apply_cal_survey_info_commit(self, _widget=None) -> None:
        if hasattr(self, "_apply_cal_lead_in_change"):
            self._apply_cal_lead_in_change()
        elif hasattr(self, "_update_cal_line_times"):
            self._update_cal_line_times()

    def _setup_deferred_parameter_bindings(self) -> None:
        """Wire line-edit parameters: amber highlight while typing, apply on Enter / blur."""
        self._init_deferred_params()

        accuracy_plan = getattr(self, "_commit_accuracy_plan_field", None)
        if accuracy_plan is not None:
            for name in (
                "central_lat_entry",
                "central_lon_entry",
                "line_length_entry",
                "heading_entry",
                "dist_between_lines_entry",
                "num_lines_entry",
                "bisect_lead_entry",
                "survey_speed_entry",
                "crossline_passes_entry",
                "ref_turn_time_entry",
            ):
                w = getattr(self, name, None)
                if w is not None:
                    self._bind_deferred_param(w, accuracy_plan)

            for name in ("multiplier_entry_len", "multiplier_entry_dist"):
                w = getattr(self, name, None)
                if w is not None:
                    fn = (
                        self._on_line_length_multiplier_changed
                        if name == "multiplier_entry_len"
                        else self._on_separation_multiplier_changed
                    )
                    self._bind_deferred_param(w, lambda _widget, f=fn: f())

            export_entry = getattr(self, "export_name_entry", None)
            if export_entry is not None and hasattr(self, "_update_export_name"):
                self._bind_deferred_param(export_entry, lambda _w: self._update_export_name())

        for name, handler in (
            ("geotiff_nan_entry", "_on_geotiff_nan_value_changed"),
            ("contour_interval_entry", "_on_contour_interval_changed"),
            ("slope_overlay_min_entry", "_on_slope_overlay_min_changed"),
            ("slope_overlay_max_entry", "_on_slope_overlay_max_changed"),
        ):
            w = getattr(self, name, None)
            fn = getattr(self, handler, None)
            if w is not None and fn is not None:
                self._bind_deferred_param(w, lambda _widget, f=fn: f())

        cal_lead = getattr(self, "cal_lead_in_entry", None)
        if cal_lead is not None:
            self._bind_deferred_param(cal_lead, self._apply_cal_survey_info_commit)

        cal_offset = getattr(self, "cal_line_offset_entry", None)
        if cal_offset is not None and hasattr(self, "_on_cal_line_offset_user_edited"):
            self._bind_deferred_param(cal_offset, lambda _w: self._on_cal_line_offset_user_edited())

        for name in ("cal_survey_speed_entry", "cal_turn_time_entry"):
            w = getattr(self, name, None)
            if w is not None:
                self._bind_deferred_param(w, self._apply_cal_survey_info_commit)

        backscatter_handlers = (
            ("backscatter_depth_min_m_entry", "_on_backscatter_depth_text_changed"),
            ("backscatter_depth_max_m_entry", "_on_backscatter_depth_text_changed"),
            ("backscatter_min_area_m2_entry", "_on_backscatter_min_area_m2_text_changed"),
            ("backscatter_min_width_m_entry", "_on_backscatter_extent_m_text_changed"),
            ("backscatter_min_height_m_entry", "_on_backscatter_extent_m_text_changed"),
            ("backscatter_nan_entry", "_on_backscatter_nan_value_changed"),
            ("backscatter_percent_clip_min_entry", "_on_backscatter_percent_clip_changed"),
            ("backscatter_percent_clip_max_entry", "_on_backscatter_percent_clip_changed"),
            ("backscatter_lead_in_entry", "_on_backscatter_line_info_text_changed"),
            ("backscatter_survey_speed_entry", "_on_backscatter_line_info_text_changed"),
            ("backscatter_swath_ang_entry", "_on_backscatter_line_info_text_changed"),
            ("backscatter_sv_entry", "_on_backscatter_line_info_text_changed"),
            ("backscatter_box_width_entry", "_on_backscatter_line_info_text_changed"),
        )
        for name, handler_name in backscatter_handlers:
            w = getattr(self, name, None)
            fn = getattr(self, handler_name, None)
            if w is not None and fn is not None:
                self._bind_deferred_param(w, lambda _widget, f=fn: f())

        perf_commit = self._apply_performance_param_commit
        for name in (
            "performance_central_lat_entry",
            "performance_central_lon_entry",
            "performance_test_depth_entry",
            "performance_swath_angle_entry",
            "performance_sound_velocity_entry",
            "performance_swell_direction_entry",
            "performance_num_pings_entry",
            "performance_bist_time_entry",
            "performance_test_speed_entry",
            "performance_turn_time_entry",
            "performance_export_name_entry",
        ):
            w = getattr(self, name, None)
            if w is not None:
                self._bind_deferred_param(w, lambda _widget, f=perf_commit: f())

        cal_export = getattr(self, "cal_export_name_entry", None)
        if cal_export is not None:
            self._bind_deferred_param(cal_export, lambda _w: None)

        adcp_commit = getattr(self, "_apply_adcp_params_commit", None)
        if adcp_commit is not None:
            for name in (
                "adcp_circle_diameter_entry",
                "adcp_survey_speed_entry",
                "adcp_turn_time_entry",
            ):
                w = getattr(self, name, None)
                if w is not None:
                    self._bind_deferred_param(w, lambda _widget, f=adcp_commit: f())

        adcp_export = getattr(self, "adcp_export_name_entry", None)
        if adcp_export is not None:
            self._bind_deferred_param(adcp_export, lambda _w: None)

        if hasattr(self, "offset_direction_combo") and hasattr(self, "_apply_accuracy_survey_plan_now"):
            self.offset_direction_combo.currentTextChanged.connect(
                lambda _text: self._apply_accuracy_survey_plan_now()
            )

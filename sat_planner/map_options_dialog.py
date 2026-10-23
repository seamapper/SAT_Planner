"""Non-modal Map Options dialog for map overlays and GeoTIFF display options."""

from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QGroupBox,
    QLabel,
    QSizePolicy,
    QPushButton,
)


class MapOptionsDialog(QDialog):
    """Hosts map layer and GeoTIFF overlay controls moved off the main layout."""

    def __init__(
        self,
        parent,
        imagery_basemap_checkbox,
        noaa_charts_checkbox,
        noaa_charts_opacity_label,
        noaa_charts_opacity_slider,
        eez_checkbox,
        eez_opacity_label,
        eez_opacity_slider,
        add_shapefile_btn,
        show_contours_checkbox,
        contour_interval_entry,
        show_slope_overlay_checkbox,
        slope_overlay_band_widgets,
        slope_overlay_opacity_label,
        slope_overlay_opacity_slider,
        show_depth_overlay_checkbox,
        depth_overlay_band_widgets,
        depth_overlay_opacity_label,
        depth_overlay_opacity_slider,
        elevation_slope_combo,
        shaded_relief_cmap_combo,
        show_map_colorbar_checkbox,
        show_topography_checkbox,
        show_map_legend_checkbox,
        dyn_vert_exag_btn,
        dynamic_resolution_btn,
    ):
        super().__init__(parent)
        self.setWindowTitle("Map Options")
        self.setModal(False)
        self.resize(420, 620)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)

        layers_group = QGroupBox("Map Overlay Options")
        layers_layout = QVBoxLayout(layers_group)
        layers_layout.setSpacing(6)
        layers_layout.addWidget(imagery_basemap_checkbox)

        noaa_row = QHBoxLayout()
        noaa_row.addWidget(noaa_charts_checkbox)
        noaa_row.addWidget(noaa_charts_opacity_label)
        noaa_charts_opacity_slider.setMaximumWidth(120)
        noaa_row.addWidget(noaa_charts_opacity_slider, 1)
        layers_layout.addLayout(noaa_row)

        eez_row = QHBoxLayout()
        eez_row.addWidget(eez_checkbox)
        eez_row.addWidget(eez_opacity_label)
        eez_opacity_slider.setMaximumWidth(120)
        eez_row.addWidget(eez_opacity_slider, 1)
        layers_layout.addLayout(eez_row)

        layers_layout.addWidget(add_shapefile_btn)
        layout.addWidget(layers_group)

        geotiff_group = QGroupBox("GeoTIFF Overlay Options")
        geotiff_layout = QVBoxLayout(geotiff_group)
        geotiff_layout.setSpacing(6)

        half_width = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

        map_display_row = QHBoxLayout()
        map_display_row.addWidget(QLabel("Map Display:"))
        elevation_slope_combo.setSizePolicy(half_width)
        map_display_row.addWidget(elevation_slope_combo, 1)
        geotiff_layout.addLayout(map_display_row)

        cmap_row = QHBoxLayout()
        cmap_row.addWidget(QLabel("Color Map:"))
        shaded_relief_cmap_combo.setSizePolicy(half_width)
        cmap_row.addWidget(shaded_relief_cmap_combo, 1)
        geotiff_layout.addLayout(cmap_row)
        colorbar_row = QHBoxLayout()
        colorbar_row.addWidget(show_map_colorbar_checkbox)
        colorbar_row.addWidget(show_topography_checkbox)
        colorbar_row.addWidget(show_map_legend_checkbox)
        colorbar_row.addStretch()
        geotiff_layout.addLayout(colorbar_row)

        display_row = QHBoxLayout()
        dyn_vert_exag_btn.setSizePolicy(half_width)
        dynamic_resolution_btn.setSizePolicy(half_width)
        display_row.addWidget(dyn_vert_exag_btn, 1)
        display_row.addWidget(dynamic_resolution_btn, 1)
        geotiff_layout.addLayout(display_row)

        contours_row = QHBoxLayout()
        contours_row.addWidget(show_contours_checkbox)
        contour_interval_entry.setMaximumWidth(80)
        contours_row.addWidget(contour_interval_entry)
        contours_row.addStretch()
        geotiff_layout.addLayout(contours_row)

        slopes_group = QGroupBox("Slopes")
        slopes_layout = QVBoxLayout(slopes_group)
        slopes_layout.setSpacing(6)
        slopes_layout.addWidget(show_slope_overlay_checkbox)
        for band_idx, (min_entry, max_entry, color_btn) in enumerate(slope_overlay_band_widgets):
            slopes_row = QHBoxLayout()
            slopes_row.addWidget(QLabel(f"Range {band_idx + 1}"))
            min_entry.setMaximumWidth(60)
            slopes_row.addWidget(min_entry)
            slopes_row.addWidget(QLabel("To"))
            max_entry.setMaximumWidth(60)
            slopes_row.addWidget(max_entry)
            slopes_row.addWidget(color_btn)
            slopes_row.addStretch()
            slopes_layout.addLayout(slopes_row)

        opacity_row = QHBoxLayout()
        opacity_row.addWidget(slope_overlay_opacity_label)
        slope_overlay_opacity_slider.setMaximumWidth(160)
        opacity_row.addWidget(slope_overlay_opacity_slider, 1)
        slopes_layout.addLayout(opacity_row)
        geotiff_layout.addWidget(slopes_group)

        depths_group = QGroupBox("Depths")
        depths_layout = QVBoxLayout(depths_group)
        depths_layout.setSpacing(6)
        depths_layout.addWidget(show_depth_overlay_checkbox)
        for band_idx, (min_entry, max_entry, color_btn) in enumerate(depth_overlay_band_widgets):
            depths_row = QHBoxLayout()
            depths_row.addWidget(QLabel(f"Range {band_idx + 1}"))
            min_entry.setMaximumWidth(60)
            depths_row.addWidget(min_entry)
            depths_row.addWidget(QLabel("To"))
            max_entry.setMaximumWidth(60)
            depths_row.addWidget(max_entry)
            depths_row.addWidget(color_btn)
            depths_row.addStretch()
            depths_layout.addLayout(depths_row)

        depth_opacity_row = QHBoxLayout()
        depth_opacity_row.addWidget(depth_overlay_opacity_label)
        depth_overlay_opacity_slider.setMaximumWidth(160)
        depth_opacity_row.addWidget(depth_overlay_opacity_slider, 1)
        depths_layout.addLayout(depth_opacity_row)
        geotiff_layout.addWidget(depths_group)

        layout.addWidget(geotiff_group)

        layout.addStretch()

        # Prevent Enter in text fields from activating Add Shapefile / other buttons.
        self._disable_dialog_button_defaults()

    def _disable_dialog_button_defaults(self):
        """Make sure no button becomes the dialog default on Return."""
        for btn in self.findChildren(QPushButton):
            btn.setAutoDefault(False)
            btn.setDefault(False)

    def showEvent(self, event):
        super().showEvent(event)
        # Buttons can regain auto-default behavior when reparented/shown.
        self._disable_dialog_button_defaults()

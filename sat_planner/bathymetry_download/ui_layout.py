"""Main window UI construction for the Bathymetry Downloader."""

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QProgressBar,
    QPushButton,
    QSizePolicy,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from .ui_widgets import ClickableLabel


def build_main_window_ui(window):
    """Build widgets and layouts on ``window`` (a MainWindow instance)."""
    window.setWindowTitle("Download Bathymetry")
    window.resize(1200, 800)

    main_layout = QHBoxLayout(window)

    left_container = _build_left_panel(window)
    right_panel = _build_right_panel(window)
    main_layout.addWidget(left_container)
    main_layout.addWidget(right_panel)
    main_layout.setStretch(0, 1)
    main_layout.setStretch(1, 0)

    window._update_output_options_visibility()
    window._update_attribution()


def _build_left_panel(window):
    """Map group, black canvas, and attribution strip."""
    window.map_group = QGroupBox("Map")
    window.map_group.setObjectName("Map")
    map_layout = QVBoxLayout()

    map_controls = QHBoxLayout()

    window.legend_checkbox = QCheckBox("Legend")
    window.legend_checkbox.setChecked(True)
    window.legend_checkbox.stateChanged.connect(window.on_legend_toggled)
    map_controls.addWidget(window.legend_checkbox)

    window.aoi_checkbox = QCheckBox("AoI")
    window.aoi_checkbox.setChecked(True)
    window.aoi_checkbox.stateChanged.connect(window.on_aoi_toggled)
    map_controls.addWidget(window.aoi_checkbox)

    window.gmrt_mask_checkbox = QCheckBox("GMRT Hi-Res Mask")
    window.gmrt_mask_checkbox.setChecked(False)
    window.gmrt_mask_checkbox.setToolTip(
        "Highlight high-resolution GMRT coverage in the map preview"
    )
    window.gmrt_mask_checkbox.stateChanged.connect(window.on_gmrt_mask_toggled)
    window.gmrt_mask_checkbox.setVisible(False)
    map_controls.addWidget(window.gmrt_mask_checkbox)

    window.zoom_back_btn = QPushButton("Zoom Prev")
    window.zoom_back_btn.clicked.connect(window.zoom_back)
    window.zoom_back_btn.setEnabled(False)
    window.zoom_next_btn = QPushButton("Zoom Next")
    window.zoom_next_btn.clicked.connect(window.zoom_next)
    window.zoom_next_btn.setEnabled(False)
    window.fit_extent_btn = QPushButton("Zoom to Full Extent")
    window.fit_extent_btn.clicked.connect(window.fit_to_extent)
    window.clear_selection_btn = QPushButton("Clear Selection")
    window.clear_selection_btn.clicked.connect(window.clear_selection)
    window.refresh_map_btn = QPushButton("Refresh Map")
    window.refresh_map_btn.clicked.connect(window.refresh_map)
    map_controls.addWidget(window.zoom_back_btn)
    map_controls.addWidget(window.zoom_next_btn)
    map_controls.addWidget(window.fit_extent_btn)
    map_controls.addWidget(window.clear_selection_btn)
    map_controls.addWidget(window.refresh_map_btn)
    map_controls.addStretch()
    map_layout.addLayout(map_controls)

    window.map_canvas = QWidget()
    window.map_canvas.setObjectName("MapCanvas")
    window.map_canvas.setStyleSheet("QWidget#MapCanvas { background-color: #000000; }")
    window.map_canvas.setMinimumSize(600, 400)
    map_canvas_layout = QVBoxLayout(window.map_canvas)
    map_canvas_layout.setContentsMargins(0, 0, 0, 0)
    map_canvas_layout.setSpacing(0)

    window.map_widget = None
    window.loading_label = QLabel("Loading service info...")
    window.loading_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
    window.loading_label.setStyleSheet("background-color: #000000; color: #cccccc;")
    map_canvas_layout.addWidget(window.loading_label)
    map_layout.addWidget(window.map_canvas)
    window.map_group.setLayout(map_layout)

    left_container = QWidget()
    left_container_layout = QVBoxLayout(left_container)
    left_container_layout.setContentsMargins(0, 0, 0, 0)
    left_container_layout.addWidget(window.map_group)

    attribution_group = QGroupBox("Data Set Attribution")
    attribution_group.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
    attribution_group.setMinimumHeight(52)
    attribution_group.setMaximumHeight(64)
    attribution_layout = QVBoxLayout()
    attribution_layout.setContentsMargins(2, 4, 2, 4)
    attribution_layout.setSpacing(0)
    window.attribution_label = ClickableLabel()
    window.attribution_label.setWordWrap(True)
    window.attribution_label.setStyleSheet(
        "color: orange; text-decoration: underline; cursor: pointer; "
        "padding: 0px; margin: 0px; border: none;"
    )
    window.attribution_label.setContentsMargins(0, 0, 0, 0)
    window.attribution_label.clicked.connect(window._open_attribution_url)
    window._current_attribution_url = None
    window._legend_was_on_before_aoi_off = False
    attribution_layout.addWidget(window.attribution_label)
    attribution_group.setLayout(attribution_layout)
    left_container_layout.addWidget(attribution_group)

    return left_container


def _build_right_panel(window):
    """Controls, selection fields, output options, progress, and activity log."""
    right_panel = QWidget()
    right_panel.setFixedWidth(480)
    right_layout = QVBoxLayout(right_panel)

    data_source_group = QGroupBox("Data Source")
    data_source_layout = QVBoxLayout()
    window.data_source_combo = QComboBox()
    window.data_source_combo.addItems(list(window.ui_source_names))
    window.data_source_combo.setCurrentText(window.current_data_source)
    window.data_source_combo.currentTextChanged.connect(window.on_data_source_changed)
    data_source_layout.addWidget(window.data_source_combo)
    data_source_group.setLayout(data_source_layout)
    right_layout.addWidget(data_source_group)

    selection_group = QGroupBox("Selected Area")
    selection_main_layout = QGridLayout()
    window.west_edit = QLineEdit()
    window.west_edit.setPlaceholderText("West")
    window.south_edit = QLineEdit()
    window.south_edit.setPlaceholderText("South")
    window.east_edit = QLineEdit()
    window.east_edit.setPlaceholderText("East")
    window.north_edit = QLineEdit()
    window.north_edit.setPlaceholderText("North")

    window._coordinate_fields = (
        window.west_edit,
        window.south_edit,
        window.east_edit,
        window.north_edit,
    )
    for field in window._coordinate_fields:
        field.returnPressed.connect(window.on_geographic_changed)
        field.installEventFilter(window)

    north_layout = QHBoxLayout()
    north_layout.addWidget(QLabel("North:"))
    north_layout.addWidget(window.north_edit)
    selection_main_layout.addLayout(north_layout, 0, 1)

    west_layout = QHBoxLayout()
    west_layout.addWidget(QLabel("West:"))
    west_layout.addWidget(window.west_edit)
    selection_main_layout.addLayout(west_layout, 1, 0)

    east_layout = QHBoxLayout()
    east_layout.addWidget(QLabel("East:"))
    east_layout.addWidget(window.east_edit)
    selection_main_layout.addLayout(east_layout, 1, 2)

    south_layout = QHBoxLayout()
    south_layout.addWidget(QLabel("South:"))
    south_layout.addWidget(window.south_edit)
    selection_main_layout.addLayout(south_layout, 2, 1)

    selection_group.setLayout(selection_main_layout)
    right_layout.addWidget(selection_group)

    output_group = QGroupBox("Output Options")
    output_layout = QVBoxLayout()

    window.output_data_types_group = QGroupBox("Output Grid Data Types")
    output_data_types_layout = QVBoxLayout()
    window.download_mode_container = QWidget()
    download_mode_layout = QGridLayout(window.download_mode_container)
    window.check_combined = QCheckBox("Combined Bathymetry && Land")
    window.check_combined.setChecked(True)
    window.check_combined.setToolTip("Native grid with bathymetry and elevation")
    download_mode_layout.addWidget(window.check_combined, 0, 0)
    window.check_direct_measurements_only = QCheckBox("Direct Measurements")
    window.check_direct_measurements_only.setToolTip(
        "Only cells where TID is 10–20 (direct measurements)"
    )
    download_mode_layout.addWidget(window.check_direct_measurements_only, 0, 1)
    window.check_direct_unknown_measurements_only = QCheckBox(
        "Direct && Unknown Measurement"
    )
    window.check_direct_unknown_measurements_only.setToolTip(
        "Only cells where TID is 10–20, 44, or 70 (direct and unknown measurements)"
    )
    download_mode_layout.addWidget(window.check_direct_unknown_measurements_only, 1, 1)
    window.check_bathymetry_only = QCheckBox("Bathymetry")
    window.check_bathymetry_only.setToolTip("Only cells where TID is not 0 (water)")
    download_mode_layout.addWidget(window.check_bathymetry_only, 1, 0)
    window.check_land_only = QCheckBox("Land")
    window.check_land_only.setToolTip("Only cells where TID is 0 (land)")
    download_mode_layout.addWidget(window.check_land_only, 2, 0)
    for cb in (
        window.check_combined,
        window.check_bathymetry_only,
        window.check_land_only,
        window.check_direct_measurements_only,
        window.check_direct_unknown_measurements_only,
    ):
        cb.toggled.connect(window.check_and_update_download_button)
    window.check_direct_measurements_only.toggled.connect(
        lambda checked: window.log_message(
            "Only extracting bathymetry values with associated TID values from 10 to 17",
            bold=True,
        )
        if checked
        else None
    )
    window.check_direct_unknown_measurements_only.toggled.connect(
        lambda checked: window.log_message(
            "Only extracting bathymetry values with TID 10 to 17, 44 and 70",
            bold=True,
        )
        if checked
        else None
    )
    output_data_types_layout.addWidget(window.download_mode_container)
    window.output_data_types_group.setLayout(output_data_types_layout)
    output_layout.addWidget(window.output_data_types_group)

    window.cell_size_meters_container = QWidget()
    cell_size_m_layout = QHBoxLayout(window.cell_size_meters_container)
    cell_size_m_layout.setContentsMargins(0, 0, 0, 0)
    window.cell_size_label = QLabel("Cell Size (m):")
    cell_size_m_layout.addWidget(window.cell_size_label)
    window.cell_size_combo = QComboBox()
    window.cell_size_combo.setMinimumWidth(100)
    window.cell_size_combo.currentTextChanged.connect(window.on_cell_size_changed)
    cell_size_m_layout.addWidget(window.cell_size_combo)
    cell_size_m_layout.addStretch()
    output_layout.addWidget(window.cell_size_meters_container)

    window.cell_size_degrees_container = QWidget()
    cell_size_deg_layout = QHBoxLayout(window.cell_size_degrees_container)
    cell_size_deg_layout.setContentsMargins(0, 0, 0, 0)
    cell_size_deg_layout.addWidget(QLabel("Cell Size (deg):"))
    window.cell_size_degrees_edit = QLineEdit()
    window.cell_size_degrees_edit.setPlaceholderText("Service native resolution")
    cell_size_deg_layout.addWidget(window.cell_size_degrees_edit)
    output_layout.addWidget(window.cell_size_degrees_container)
    window.cell_size_degrees_edit.editingFinished.connect(
        window._on_cell_size_degrees_changed
    )

    window.bathymetry_only_notice_label = QLabel(
        "Download is bathymetry only; the land layer is for map display reference."
    )
    window.bathymetry_only_notice_label.setWordWrap(True)
    output_layout.addWidget(window.bathymetry_only_notice_label)

    window.pixel_count_label = QLabel("Pixels: --")
    window.pixel_count_label.setStyleSheet("font-weight: bold; padding: 5px;")
    output_layout.addWidget(window.pixel_count_label)

    output_group.setLayout(output_layout)
    right_layout.addWidget(output_group)

    output_dir_btn = QPushButton("Select Output Directory")
    output_dir_btn.clicked.connect(window.select_output_directory)
    right_layout.addWidget(output_dir_btn)

    output_dir_layout = QHBoxLayout()
    output_dir_layout.addWidget(QLabel("Directory:"))
    window.output_dir_edit = QLineEdit()
    window.output_dir_edit.setPlaceholderText("Not set")
    window.output_dir_edit.setReadOnly(True)
    output_dir_layout.addWidget(window.output_dir_edit, stretch=1)
    right_layout.addLayout(output_dir_layout)

    button_layout = QHBoxLayout()
    window.download_btn = QPushButton("Download Selected Area")
    window.download_btn.clicked.connect(window.start_download)
    window.download_btn.setEnabled(False)
    button_layout.addWidget(window.download_btn)
    right_layout.addLayout(button_layout)

    window.tile_download_checkbox = QCheckBox("Tile Download")
    window.tile_download_checkbox.setChecked(True)
    right_layout.addWidget(window.tile_download_checkbox)

    window.split_topo_depths_checkbox = QCheckBox("Split Topo/Depths after download")
    window.split_topo_depths_checkbox.setChecked(True)
    window.split_topo_depths_checkbox.setToolTip(
        "When on, split the downloaded grid and load only the bathymetry file into SAT Planner."
    )
    right_layout.addWidget(window.split_topo_depths_checkbox)

    progress_group = QGroupBox("Progress")
    progress_layout = QVBoxLayout()
    window.progress_bar = QProgressBar()
    window.progress_bar.setValue(0)
    progress_layout.addWidget(window.progress_bar)
    window.status_label = QLabel("Ready")
    progress_layout.addWidget(window.status_label)
    progress_group.setLayout(progress_layout)
    right_layout.addWidget(progress_group)

    log_group = QGroupBox("Activity Log")
    log_layout = QVBoxLayout()
    window.log_text = QTextEdit()
    window.log_text.setReadOnly(True)
    log_layout.addWidget(window.log_text)
    log_group.setLayout(log_layout)
    right_layout.addWidget(log_group, 1)

    return right_panel

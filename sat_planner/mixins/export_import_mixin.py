"""
Export/import: save/load survey parameters (JSON), export survey files (CSV, Shapefile, GeoJSON, etc.).
_save_survey_parameters, _load_survey_parameters, _load_survey_parameters_dialog, _export_survey_files.
"""
import csv
import datetime
import json
import os

from PyQt6.QtWidgets import QFileDialog

from sat_planner import GEOSPATIAL_LIBS_AVAILABLE, decimal_degrees_to_ddm
from sat_planner.constants import LineString, fiona, pyproj
from sat_planner import export_utils

try:
    from shapely.geometry import mapping
except ImportError:
    mapping = None


class ExportImportMixin:
    """Mixin for save/load survey parameters and export survey files."""

    def _export_type_enabled(self, key):
        defaults = (
            self._default_export_type_options()
            if hasattr(self, "_default_export_type_options")
            else {}
        )
        options = getattr(self, "export_type_options", {}) or {}
        if key in options:
            return bool(options[key])
        return bool(defaults.get(key, True))

    def _save_survey_parameters(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot save parameters.")
            return

        is_valid, values = self._validate_inputs()
        if not is_valid:
            return

        save_dir = QFileDialog.getExistingDirectory(self, "Select Directory to Save Survey Parameters", self.last_survey_params_dir)
        if not save_dir:
            return
        self.last_survey_params_dir = save_dir
        self._save_last_survey_params_dir()

        # Use the export name from the form, with .json extension
        export_name = values['export_name']
        if not export_name:
            # Fallback to default naming if export name is empty
            dist_between_lines = int(values['dist_between_lines'])
            heading_int = int(values['heading'])
            export_name = f"Accuracy_{dist_between_lines}m_{heading_int}deg"
        filename = f"{export_name}_params.json"
        file_path = os.path.join(save_dir, filename)

        try:
            # Create parameters dictionary
            params = {
                'central_lat': float(self.central_lat_entry.text()),
                'central_lon': float(self.central_lon_entry.text()),
                'central_point_depth_m': float(getattr(self, '_depth_at_picked_point', 0.0)) if getattr(self, '_depth_at_picked_point', None) is not None else None,
                'line_length': float(self.line_length_entry.text()),
                'heading': float(self.heading_entry.text()),
                'dist_between_lines': float(self.dist_between_lines_entry.text()),
                'num_lines': int(self.num_lines_entry.text()),
                'bisect_lead': float(self.bisect_lead_entry.text()),
                'survey_speed': float(self.survey_speed_entry.text()),
                'export_name': self.export_name_entry.text().strip(),
                'geotiff_path': (self.current_geotiff_path if hasattr(self, 'current_geotiff_path') and self.current_geotiff_path else None),
                'visualization_shapefile_paths': list(getattr(self, 'visualization_shapefile_paths', []) or []),
                'offset_direction': self.offset_direction_var,
                'line_length_multiplier': self.line_length_multiplier,
                    'dist_between_lines_multiplier': self.dist_between_lines_multiplier
            }

            # Save to JSON file
            with open(file_path, 'w') as f:
                json.dump(params, f, indent=4)

            self.set_ref_info_text(f"Survey parameters saved successfully to:\n{file_path}", append=False)

        except Exception as e:
            self._show_message("error","Save Error", f"Failed to save survey parameters: {e}")

    def _load_survey_parameters(self, file_path):
        """Load survey parameters from a JSON file."""
        try:
            with open(file_path, 'r') as f:
                params = json.load(f)

            # Update all input fields
            self.central_lat_entry.clear()
            self.central_lat_entry.setText(str(params['central_lat']))

            self.central_lon_entry.clear()
            self.central_lon_entry.setText(str(params['central_lon']))

            self.line_length_entry.setText(str(params['line_length']))

            self.heading_entry.setText(str(params['heading']))

            self.dist_between_lines_entry.setText(str(params['dist_between_lines']))

            self.num_lines_entry.setText(str(params['num_lines']))

            self.bisect_lead_entry.setText(str(params['bisect_lead']))

            self.survey_speed_entry.setText(str(params.get('survey_speed', '')))

            self.export_name_entry.clear()
            self.export_name_entry.setText(params['export_name'])

            self.offset_direction_var.set(params['offset_direction'])

            self.line_length_multiplier = float(params['line_length_multiplier'])
            self._update_multiplier_label_len(self.line_length_multiplier)

            self.dist_between_lines_multiplier = float(params['dist_between_lines_multiplier'])
            self._update_multiplier_label_dist(self.dist_between_lines_multiplier)

            # Regenerate the plot with loaded parameters
            self._generate_and_plot()

            self.set_ref_info_text(f"Survey parameters loaded from: {file_path}", append=False)

        except Exception as e:
            self._show_message("error","Load Error", f"Failed to load survey parameters: {e}")

    def _load_survey_parameters_dialog(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot load parameters.")
            return

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Load Survey Parameters",
            self.last_survey_params_dir,
            "JSON files (*.json);;All files (*.*)"
        )
        if file_path:
            self.last_survey_params_dir = os.path.dirname(file_path)
            self._save_last_survey_params_dir()
            self._load_survey_parameters(file_path)

    def _export_survey_files(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot export survey files.")
            return

        is_valid, values = self._validate_inputs()
        if not is_valid:
            return

        export_name = values['export_name']

        if not self.survey_lines_data and not self.cross_line_data:
            self._show_message("warning","No Data", "No survey lines to export. Generate them first.")
            return

        export_dir = QFileDialog.getExistingDirectory(self, "Select Export Directory", self.last_export_dir)
        if not export_dir:
            return
        self.last_export_dir = export_dir
        self._save_last_export_dir()

        try:
            profile_csv_path = None
            profile_png_paths = []
            export_shapefile = self._export_type_enabled("esri_shapefile")
            export_sis = self._export_type_enabled("sis_asciiplan")
            export_gpx = self._export_type_enabled("gpx")
            export_text_csv = self._export_type_enabled("text_csv")
            export_text_txt = self._export_type_enabled("text_txt")
            export_hypack = self._export_type_enabled("hypack_lnw")
            export_map_png = self._export_type_enabled("map_png")
            export_profiles_png = self._export_type_enabled("profiles_png")
            if export_shapefile and mapping is None:
                raise ImportError("shapely.geometry.mapping is required for shapefile export")

            # --- Build common rows (line_num, line_name, point_label, lat, lon) for reference survey ---
            ref_rows = []
            for i, line in enumerate(self.survey_lines_data):
                if i % 2 == 0:
                    start, end = line[0], line[1]
                else:
                    start, end = line[1], line[0]
                line_name = f'ReferenceLine{i + 1}'
                ref_rows.append((i + 1, line_name, f'L{i+1}S', start[0], start[1]))
                ref_rows.append((i + 1, line_name, f'L{i+1}E', end[0], end[1]))
            if self.cross_line_data:
                ref_rows.append((0, 'Crossline', 'CLS', self.cross_line_data[0][0], self.cross_line_data[0][1]))
                ref_rows.append((0, 'Crossline', 'CLE', self.cross_line_data[1][0], self.cross_line_data[1][1]))

            csv_file_path = os.path.join(export_dir, f"{export_name}_DDD.csv")
            ddm_file_path = os.path.join(export_dir, f"{export_name}_DMM.csv")
            dms_file_path = os.path.join(export_dir, f"{export_name}_DMS.csv")
            ddm_txt_file_path = os.path.join(export_dir, f"{export_name}_DMM.txt")
            dms_txt_file_path = os.path.join(export_dir, f"{export_name}_DMS.txt")
            txt_file_path = os.path.join(export_dir, f"{export_name}_DDD.txt")
            if export_text_csv:
                export_utils.write_ddd_csv(csv_file_path, ref_rows, newline='')
                export_utils.write_dmm_csv(ddm_file_path, ref_rows)
                export_utils.write_dms_csv(dms_file_path, ref_rows)
            if export_text_txt:
                export_utils.write_dmm_txt(ddm_txt_file_path, ref_rows)
                export_utils.write_dms_txt(dms_txt_file_path, ref_rows)
                export_utils.write_ddd_txt(txt_file_path, ref_rows)

            # --- Export to ESRI Shapefile (.shp) ---
            shapefile_path = os.path.join(export_dir, f"{export_name}.shp")
            if export_shapefile:
                schema = {
                    'geometry': 'LineString',
                    'properties': {'line_num': 'int', 'line_name': 'str'},
                }
                crs_epsg = 'EPSG:4326'  # WGS 84
                features = []
                # Add main survey lines (names match CSV / GPX: ReferenceLine1, …)
                for i, line_coords in enumerate(self.survey_lines_data):
                    shapely_line = LineString([(p[1], p[0]) for p in line_coords])
                    features.append({
                        'geometry': mapping(shapely_line),
                        'properties': {'line_num': i + 1, 'line_name': f'ReferenceLine{i + 1}'},
                    })
                # Add crossline
                if self.cross_line_data:
                    shapely_cross_line = LineString([(p[1], p[0]) for p in self.cross_line_data])
                    features.append({
                        'geometry': mapping(shapely_cross_line),
                        'properties': {'line_num': 0, 'line_name': 'Crossline'},
                    })
                with fiona.open(shapefile_path, 'w', driver='ESRI Shapefile', crs=crs_epsg, schema=schema) as collection:
                    collection.writerecords(features)

            # --- Export to GeoJSON ---
            geojson_file_path = os.path.join(export_dir, f"{export_name}.geojson")
            try:
                ref_export_speed = float(self.survey_speed_entry.text()) if hasattr(self, 'survey_speed_entry') and self.survey_speed_entry.text() else 8.0
            except (ValueError, TypeError):
                ref_export_speed = 8.0
            geojson_features = []
            # Main survey lines
            for i, line in enumerate(self.survey_lines_data):
                geojson_features.append({
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": [[line[0][1], line[0][0]], [line[1][1], line[1][0]]]
                    },
                    "properties": {
                        "line_num": i + 1,
                        "survey_speed": ref_export_speed,
                        "points": [
                            {"point_num": 1, "lat": line[0][0], "lon": line[0][1]},
                            {"point_num": 2, "lat": line[1][0], "lon": line[1][1]}
                        ]
                    }
                })
            # Crossline
            if self.cross_line_data:
                geojson_features.append({
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": [[self.cross_line_data[0][1], self.cross_line_data[0][0]], [self.cross_line_data[1][1], self.cross_line_data[1][0]]]
                    },
                    "properties": {
                        "line_num": 0,
                        "survey_speed": ref_export_speed,
                        "points": [
                            {"point_num": 1, "lat": self.cross_line_data[0][0], "lon": self.cross_line_data[0][1]},
                            {"point_num": 2, "lat": self.cross_line_data[1][0], "lon": self.cross_line_data[1][1]}
                        ]
                    }
                })
            geojson_collection = {
                "type": "FeatureCollection",
                "properties": {},
                "features": geojson_features
            }
            with open(geojson_file_path, 'w') as f:
                json.dump(geojson_collection, f, indent=2)

            # --- Export to Hypack LNW format (LIN/PTS/UTM), filename includes UTM zone ---
            lnw_file_path = None
            if export_hypack:
                lnw_lines = [(f"LINE{i+1:03d}", [line[0], line[1]]) for i, line in enumerate(self.survey_lines_data)]
                if self.cross_line_data:
                    lnw_lines.append(("CROSS", [self.cross_line_data[0], self.cross_line_data[1]]))
                if lnw_lines:
                    all_pts = [p for _name, pts in lnw_lines for p in pts]
                    zone, hem = export_utils.compute_utm_zone_from_points(all_pts)
                    utm_suffix = f"_UTM{zone}{'N' if hem == 'North' else 'S'}"
                    lnw_file_path = os.path.join(export_dir, f"{export_name}{utm_suffix}.lnw")
                    if not export_utils.write_lnw(lnw_file_path, lnw_lines):
                        lnw_file_path = None

            # --- Export to Kongsberg SIS ASCII Plan format ---
            sis_file_path = os.path.join(export_dir, f"{export_name}.asciiplan")
            ref_ascii_lines = []
            if self.cross_line_data:
                ref_ascii_lines.append(('Crossline', [self.cross_line_data[0], self.cross_line_data[1]]))
            for i, line in enumerate(self.survey_lines_data):
                ref_ascii_lines.append((f'Reference{i + 1}', [line[0], line[1]]))
            if export_sis:
                export_utils.write_asciiplan(sis_file_path, ref_ascii_lines)
            gpx_file_path = os.path.join(export_dir, f"{export_name}.gpx")
            ref_gpx_lines = []
            if self.cross_line_data:
                ref_gpx_lines.append(("Crossline", [self.cross_line_data[0], self.cross_line_data[1]]))
            for i, line in enumerate(self.survey_lines_data):
                if i % 2 == 0:
                    start, end = line[0], line[1]
                else:
                    start, end = line[1], line[0]
                ref_gpx_lines.append((f"ReferenceLine{i + 1}", [start, end]))
            gpx_written = False
            gpx_per_test_names = []
            if export_gpx:
                gpx_written = export_utils.write_gpx(gpx_file_path, ref_gpx_lines, creator="SAT Planner Accuracy")
                acc_gpx_tests = []
                if self.cross_line_data:
                    acc_gpx_tests.append(
                        ("Crossline", "Crossline", [self.cross_line_data[0], self.cross_line_data[1]])
                    )
                for i, line in enumerate(self.survey_lines_data):
                    if i % 2 == 0:
                        start, end = line[0], line[1]
                    else:
                        start, end = line[1], line[0]
                    suf = f"Line{i + 1:02d}"
                    acc_gpx_tests.append((suf, f"ReferenceLine{i + 1}", [start, end]))
                gpx_per_test_names = export_utils.write_gpx_per_test_files(
                    export_dir, export_name, acc_gpx_tests, creator="SAT Planner Accuracy"
                )

            # --- Export Accuracy Survey Information (same content as Accuracy Survey Info dialog) ---
            stats_file_path = os.path.join(export_dir, f"{export_name}_info.txt")
            info_text = None
            if hasattr(self, "_build_reference_planning_info_text"):
                info_text = self._build_reference_planning_info_text(
                    export_name=export_name,
                    export_date=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                )
            if not info_text:
                info_text = "ACCURACY SURVEY INFORMATION\nUnable to generate survey info text.\n"

            def _normalize_degree_symbols_for_export(text):
                normalized_lines = []
                for line in str(text).splitlines():
                    if line.endswith("°"):
                        line = line[:-1]
                    line = line.replace("°", " ")
                    normalized_lines.append(line)
                return "\n".join(normalized_lines) + "\n"

            with open(stats_file_path, 'w', encoding='utf-8') as f:
                f.write(_normalize_degree_symbols_for_export(info_text))
            total_survey_time = self._calculate_total_survey_time()
            include_crossline = (
                total_survey_time.get('num_crossline_passes', 0) > 0
                and total_survey_time.get('crossline_total_distance_m', 0) > 0
            )

            # --- Export Survey Plan as PNG ---
            map_png_path = os.path.join(export_dir, f"{export_name}_map.png")
            if export_map_png:
                if hasattr(self, "_hide_map_hover_tooltip_for_export"):
                    self._hide_map_hover_tooltip_for_export()
                self.figure.savefig(map_png_path, dpi=300, bbox_inches='tight', facecolor='white')
            
            # --- Export Profile Plots as PNG (crossline + each main line) ---
            if export_profiles_png and hasattr(self, 'profile_fig') and self.profile_fig is not None and hasattr(self, "_draw_segment_profile"):
                combo = getattr(self, "ref_profile_select_combo", None)
                previous_profile_selection = combo.currentText().strip() if combo is not None and combo.count() > 0 else "Crossline"

                # Crossline profile first (when crossline is enabled/present in the active plan)
                if include_crossline and self.cross_line_data and len(self.cross_line_data) == 2:
                    self._draw_segment_profile(self.cross_line_data, "Crossline Elevation Profile", "darkorchid")
                    crossline_profile_png_path = os.path.join(export_dir, f"{export_name}_profile_crossline.png")
                    self.profile_fig.savefig(crossline_profile_png_path, dpi=300, bbox_inches='tight', facecolor='white')
                    profile_png_paths.append(crossline_profile_png_path)

                # Main line profiles
                for i, line_data in enumerate(self.survey_lines_data or []):
                    if not line_data or len(line_data) != 2:
                        continue
                    line_num = i + 1
                    self._draw_segment_profile(line_data, f"Main Line {line_num} Elevation Profile", "blue")
                    line_profile_png_path = os.path.join(
                        export_dir,
                        f"{export_name}_profile_main_line_{line_num:02d}.png"
                    )
                    self.profile_fig.savefig(line_profile_png_path, dpi=300, bbox_inches='tight', facecolor='white')
                    profile_png_paths.append(line_profile_png_path)

                # Restore previously selected profile for user continuity
                if combo is not None:
                    combo.blockSignals(True)
                    if combo.findText(previous_profile_selection) >= 0:
                        combo.setCurrentText(previous_profile_selection)
                    combo.blockSignals(False)
                if hasattr(self, "_draw_current_profile"):
                    self._draw_current_profile()

            # --- Export crossline elevation profile as CSV (Distance m, Elevation m, Slope deg) ---
            if export_text_csv and include_crossline and self.cross_line_data and len(self.cross_line_data) == 2:
                cl_a, cl_b = self.cross_line_data
                prof = self._profile_arrays_along_segment_endpoints(cl_a[0], cl_a[1], cl_b[0], cl_b[1])
                if prof[0] is not None:
                    profile_csv_path = os.path.join(export_dir, f"{export_name}_profile.csv")
                    export_utils.write_profile_csv(profile_csv_path, prof[0], prof[1], prof[2])

            # --- Export parameters metadata as JSON ---
            json_metadata_path = os.path.join(export_dir, f"{export_name}_params.json")
            try:
                # Get current parameter values (even if not all fields are filled, use what we have)
                params = {}
                try:
                    params['central_lat'] = float(self.central_lat_entry.text())
                except:
                    # Calculate from survey lines if not available
                    if self.survey_lines_data:
                        all_lats = [p[0] for line in self.survey_lines_data for p in line]
                        params['central_lat'] = (min(all_lats) + max(all_lats)) / 2.0
                    else:
                        params['central_lat'] = None
                
                try:
                    params['central_lon'] = float(self.central_lon_entry.text())
                except:
                    if self.survey_lines_data:
                        all_lons = [p[1] for line in self.survey_lines_data for p in line]
                        params['central_lon'] = (min(all_lons) + max(all_lons)) / 2.0
                    else:
                        params['central_lon'] = None
                
                try:
                    params['line_length'] = float(self.line_length_entry.text())
                except:
                    params['line_length'] = None
                
                try:
                    params['heading'] = float(self.heading_entry.text())
                except:
                    params['heading'] = None
                
                try:
                    params['dist_between_lines'] = float(self.dist_between_lines_entry.text())
                except:
                    params['dist_between_lines'] = None
                
                try:
                    params['num_lines'] = int(self.num_lines_entry.text())
                except:
                    params['num_lines'] = len(self.survey_lines_data) if self.survey_lines_data else None
                
                try:
                    params['bisect_lead'] = float(self.bisect_lead_entry.text())
                except:
                    params['bisect_lead'] = None
                
                try:
                    params['survey_speed'] = float(self.survey_speed_entry.text())
                except:
                    params['survey_speed'] = 8.0  # Default

                params['geotiff_path'] = (
                    self.current_geotiff_path
                    if hasattr(self, 'current_geotiff_path') and self.current_geotiff_path
                    else None
                )
                params['central_point_depth_m'] = (
                    float(self._depth_at_picked_point)
                    if hasattr(self, '_depth_at_picked_point') and self._depth_at_picked_point is not None
                    else None
                )
                params['visualization_shapefile_paths'] = list(getattr(self, 'visualization_shapefile_paths', []) or [])
                
                try:
                    params['crossline_passes'] = int(self.crossline_passes_entry.text())
                except:
                    params['crossline_passes'] = 2  # Default
                
                try:
                    params['export_name'] = self.export_name_entry.text().strip()
                except:
                    params['export_name'] = export_name
                
                try:
                    params['offset_direction'] = self.offset_direction_var
                except:
                    params['offset_direction'] = 'North'  # Default
                
                try:
                    params['line_length_multiplier'] = self.line_length_multiplier
                except:
                    params['line_length_multiplier'] = 8.0  # Default
                
                try:
                    params['dist_between_lines_multiplier'] = self.dist_between_lines_multiplier
                except:
                    params['dist_between_lines_multiplier'] = 1.0  # Default
                
                # Save metadata
                with open(json_metadata_path, 'w', encoding='utf-8') as f:
                    json.dump(params, f, indent=2)
            except Exception as e:
                # If metadata export fails, continue without it
                print(f"Warning: Could not export metadata: {e}")
            
            # Update success message to include stats file, text file, PNG files, and metadata
            success_files = [f"- {os.path.basename(geojson_file_path)}"]
            if export_text_csv:
                success_files.extend(
                    [
                        f"- {os.path.basename(csv_file_path)}",
                        f"- {os.path.basename(ddm_file_path)}",
                        f"- {os.path.basename(dms_file_path)}",
                    ]
                )
            if export_text_txt:
                success_files.extend(
                    [
                        f"- {os.path.basename(txt_file_path)}",
                        f"- {os.path.basename(ddm_txt_file_path)}",
                        f"- {os.path.basename(dms_txt_file_path)}",
                    ]
                )
            if export_shapefile:
                success_files.append(f"- {os.path.basename(shapefile_path)} (and associated files)")
            if lnw_file_path:
                success_files.append(f"- {os.path.basename(lnw_file_path)}")
            if gpx_written:
                success_files.append(f"- {os.path.basename(gpx_file_path)}")
            for bn in gpx_per_test_names:
                success_files.append(f"- {bn}")
            if export_sis:
                success_files.append(f"- {os.path.basename(sis_file_path)}")
            success_files.append(f"- {os.path.basename(stats_file_path)}")
            if export_map_png:
                success_files.append(f"- {os.path.basename(map_png_path)}")
            
            # Add metadata JSON if it was created
            try:
                if os.path.exists(json_metadata_path):
                    success_files.append(f"- {os.path.basename(json_metadata_path)}")
            except:
                pass
            
            # Add profile PNGs if they were created
            for profile_png_path in profile_png_paths:
                success_files.append(f"- {os.path.basename(profile_png_path)}")
            if profile_csv_path and os.path.isfile(profile_csv_path):
                success_files.append(f"- {os.path.basename(profile_csv_path)}")

            self.set_ref_info_text(
                f"Survey exported successfully to:\n" + "\n".join(success_files) + 
                f"\nin directory: {export_dir}", append=False)

        except Exception as e:
            self._show_message("error","Export Error", f"Failed to export survey files: {e}")

    def _export_performance_survey_files(self):
        """Export performance swath lines and BIST segments (same products as Accuracy export)."""
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning", "Disabled Feature", "Geospatial libraries not loaded. Cannot export.")
            return
        export_shapefile = self._export_type_enabled("esri_shapefile")
        export_sis = self._export_type_enabled("sis_asciiplan")
        export_gpx = self._export_type_enabled("gpx")
        export_text_csv = self._export_type_enabled("text_csv")
        export_text_txt = self._export_type_enabled("text_txt")
        export_hypack = self._export_type_enabled("hypack_lnw")
        export_map_png = self._export_type_enabled("map_png")
        export_profiles_png = self._export_type_enabled("profiles_png")
        if export_shapefile and mapping is None:
            self._show_message("warning", "Export Error", "Shapely is required for shapefile export.")
            return

        lines = getattr(self, "performance_test_lines_data", None) or []
        if len(lines) != 4:
            self._show_message("warning", "No Data", "Plot four performance test lines first (Plot Performance Lines).")
            return

        export_name = ""
        if hasattr(self, "performance_export_name_entry"):
            export_name = self.performance_export_name_entry.text().strip()
        if not export_name:
            export_name = "Performance_export"
        bad = '<>:"/\\|?*'
        for c in bad:
            export_name = export_name.replace(c, "_")
        export_name = export_name.strip().strip(".")

        export_dir = QFileDialog.getExistingDirectory(self, "Select Export Directory", self.last_export_dir)
        if not export_dir:
            return
        self.last_export_dir = export_dir
        self._save_last_export_dir()

        bist_segs = getattr(self, "performance_bist_segments_data", None) or []
        has_bist = len(bist_segs) == 4

        try:
            speed_kts = 8.0
            if hasattr(self, "performance_test_speed_entry"):
                try:
                    speed_kts = float(self.performance_test_speed_entry.text().strip())
                except (ValueError, TypeError):
                    speed_kts = 8.0

            perf_rows = []
            for i, line in enumerate(lines):
                n = i + 1
                start, end = line[0], line[1]
                lname = f"PerformanceLine{n}"
                perf_rows.append((n, lname, f"P{n}S", start[0], start[1]))
                perf_rows.append((n, lname, f"P{n}E", end[0], end[1]))
            if has_bist:
                for i, seg in enumerate(bist_segs):
                    bn = i + 1
                    line_num = 10 + bn
                    lname = f"BISTLine{bn}"
                    perf_rows.append((line_num, lname, f"B{bn}S", seg[0][0], seg[0][1]))
                    perf_rows.append((line_num, lname, f"B{bn}E", seg[1][0], seg[1][1]))

            csv_file_path = os.path.join(export_dir, f"{export_name}_DDD.csv")
            ddm_file_path = os.path.join(export_dir, f"{export_name}_DMM.csv")
            dms_file_path = os.path.join(export_dir, f"{export_name}_DMS.csv")
            ddm_txt_file_path = os.path.join(export_dir, f"{export_name}_DMM.txt")
            dms_txt_file_path = os.path.join(export_dir, f"{export_name}_DMS.txt")
            txt_file_path = os.path.join(export_dir, f"{export_name}_DDD.txt")
            if export_text_csv:
                export_utils.write_ddd_csv(csv_file_path, perf_rows, newline="")
                export_utils.write_dmm_csv(ddm_file_path, perf_rows)
                export_utils.write_dms_csv(dms_file_path, perf_rows)
            if export_text_txt:
                export_utils.write_dmm_txt(ddm_txt_file_path, perf_rows)
                export_utils.write_dms_txt(dms_txt_file_path, perf_rows)
                export_utils.write_ddd_txt(txt_file_path, perf_rows)

            shapefile_path = os.path.join(export_dir, f"{export_name}.shp")
            if export_shapefile:
                schema = {"geometry": "LineString", "properties": {"line_num": "int", "line_name": "str"}}
                crs_epsg = "EPSG:4326"
                features = []
                for i, line_coords in enumerate(lines):
                    shapely_line = LineString([(p[1], p[0]) for p in line_coords])
                    n = i + 1
                    features.append({
                        "geometry": mapping(shapely_line),
                        "properties": {"line_num": n, "line_name": f"PerformanceLine{n}"},
                    })
                if has_bist:
                    for i, seg in enumerate(bist_segs):
                        shapely_line = LineString([(p[1], p[0]) for p in seg])
                        bn = i + 1
                        features.append({
                            "geometry": mapping(shapely_line),
                            "properties": {"line_num": 10 + bn, "line_name": f"BISTLine{bn}"},
                        })
                with fiona.open(shapefile_path, "w", driver="ESRI Shapefile", crs=crs_epsg, schema=schema) as collection:
                    collection.writerecords(features)

            geojson_features = []
            geotiff_path = self.current_geotiff_path if hasattr(self, "current_geotiff_path") and self.current_geotiff_path else None
            for i, line in enumerate(lines):
                geojson_features.append(
                    {
                        "type": "Feature",
                        "geometry": {
                            "type": "LineString",
                            "coordinates": [[line[0][1], line[0][0]], [line[1][1], line[1][0]]],
                        },
                        "properties": {
                            "line_num": i + 1,
                            "performance_segment": "swath",
                            "survey_speed": speed_kts,
                            "test_speed_kts": speed_kts,
                            "geotiff_path": geotiff_path,
                            "points": [
                                {"point_num": 1, "lat": line[0][0], "lon": line[0][1]},
                                {"point_num": 2, "lat": line[1][0], "lon": line[1][1]},
                            ],
                        },
                    }
                )
            if has_bist:
                for i, seg in enumerate(bist_segs):
                    geojson_features.append(
                        {
                            "type": "Feature",
                            "geometry": {
                                "type": "LineString",
                                "coordinates": [[seg[0][1], seg[0][0]], [seg[1][1], seg[1][0]]],
                            },
                            "properties": {
                                "line_num": 11 + i,
                                "performance_segment": "bist",
                                "survey_speed": speed_kts,
                                "test_speed_kts": speed_kts,
                                "geotiff_path": geotiff_path,
                                "points": [
                                    {"point_num": 1, "lat": seg[0][0], "lon": seg[0][1]},
                                    {"point_num": 2, "lat": seg[1][0], "lon": seg[1][1]},
                                ],
                            },
                        }
                    )
            geojson_collection = {
                "type": "FeatureCollection",
                "properties": {"geotiff_path": geotiff_path, "plan_type": "performance"},
                "features": geojson_features,
            }
            geojson_file_path = os.path.join(export_dir, f"{export_name}.geojson")
            with open(geojson_file_path, "w", encoding="utf-8") as f:
                json.dump(geojson_collection, f, indent=2)

            lnw_lines = [(f"PLN{i + 1:03d}", [line[0], line[1]]) for i, line in enumerate(lines)]
            if has_bist:
                for i, seg in enumerate(bist_segs):
                    lnw_lines.append((f"BIST{i + 1:03d}", [seg[0], seg[1]]))
            lnw_file_path = None
            if export_hypack and lnw_lines:
                all_pts = [p for _name, pts in lnw_lines for p in pts]
                zone, hem = export_utils.compute_utm_zone_from_points(all_pts)
                utm_suffix = f"_UTM{zone}{'N' if hem == 'North' else 'S'}"
                lnw_file_path = os.path.join(export_dir, f"{export_name}{utm_suffix}.lnw")
                if not export_utils.write_lnw(lnw_file_path, lnw_lines):
                    lnw_file_path = None

            sis_file_path = os.path.join(export_dir, f"{export_name}.asciiplan")
            perf_ascii_lines = [(f"Performance{i + 1}", [line[0], line[1]]) for i, line in enumerate(lines)]
            if has_bist:
                for i, seg in enumerate(bist_segs):
                    perf_ascii_lines.append((f"BIST{i + 1}", [seg[0], seg[1]]))
            if export_sis:
                export_utils.write_asciiplan(sis_file_path, perf_ascii_lines)

            gpx_file_path = os.path.join(export_dir, f"{export_name}.gpx")
            perf_gpx_lines = [(f"Performance{i + 1}", [line[0], line[1]]) for i, line in enumerate(lines)]
            if has_bist:
                for i, seg in enumerate(bist_segs):
                    perf_gpx_lines.append((f"BIST{i + 1}", [seg[0], seg[1]]))
            gpx_written = False
            gpx_per_test_names = []
            if export_gpx:
                gpx_written = export_utils.write_gpx(
                    gpx_file_path, perf_gpx_lines, creator="SAT Planner Performance"
                )
                perf_gpx_tests = [
                    (f"Performance{i + 1}", f"Performance{i + 1}", [line[0], line[1]])
                    for i, line in enumerate(lines)
                ]
                if has_bist:
                    for i, seg in enumerate(bist_segs):
                        perf_gpx_tests.append((f"BIST{i + 1}", f"BIST{i + 1}", [seg[0], seg[1]]))
                gpx_per_test_names = export_utils.write_gpx_per_test_files(
                    export_dir, export_name, perf_gpx_tests, creator="SAT Planner Performance"
                )

            stats_file_path = os.path.join(export_dir, f"{export_name}_info.txt")
            try:
                turn_min = 10.0
                if hasattr(self, "performance_turn_time_entry"):
                    raw_turn = self.performance_turn_time_entry.text().strip()
                    if raw_turn:
                        turn_min = float(raw_turn)
                if turn_min < 0:
                    turn_min = 0.0
            except Exception:
                turn_min = 10.0

            try:
                perf_info_text = self._build_performance_test_info_text(
                    lines, bist_segs if has_bist else [], float(speed_kts), float(turn_min)
                )
            except Exception as e:
                perf_info_text = (
                    "PERFORMANCE TEST SUMMARY\n"
                    + "=" * 24
                    + "\n"
                    + f"Could not build detailed performance info: {e}\n"
                )

            with open(stats_file_path, "w", encoding="utf-8") as f:
                f.write("PERFORMANCE SURVEY INFORMATION\n")
                f.write("=" * 50 + "\n\n")
                f.write(f"Performance Survey: {export_name}\n")
                f.write(f"Export Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                f.write(perf_info_text)

            json_metadata_path = os.path.join(export_dir, f"{export_name}_performance_params.json")
            pmeta = {
                "export_name": export_name,
                "plan_type": "performance",
                "geotiff_path": geotiff_path,
                "visualization_shapefile_paths": list(getattr(self, 'visualization_shapefile_paths', []) or []),
                "test_speed_kts": speed_kts,
            }
            for key, attr in (
                ("central_lat", "performance_central_lat_entry"),
                ("central_lon", "performance_central_lon_entry"),
                ("test_depth_m", "performance_test_depth_entry"),
                ("swell_direction_deg", "performance_swell_direction_entry"),
                ("swath_angle_deg", "performance_swath_angle_entry"),
                ("bist_time_min", "performance_bist_time_entry"),
                ("turn_time_min", "performance_turn_time_entry"),
                ("num_pings", "performance_num_pings_entry"),
                ("sound_velocity", "performance_sound_velocity_entry"),
            ):
                w = getattr(self, attr, None)
                if w is not None and hasattr(w, "text"):
                    try:
                        pmeta[key] = w.text().strip()
                    except Exception:
                        pmeta[key] = None
                else:
                    pmeta[key] = None
            try:
                ll = getattr(self, "performance_line_length_m_entry", None)
                if ll is not None and hasattr(ll, "text"):
                    pmeta["line_length_m_label"] = ll.text().strip()
            except Exception:
                pass
            with open(json_metadata_path, "w", encoding="utf-8") as f:
                json.dump(pmeta, f, indent=2)

            map_png_path = os.path.join(export_dir, f"{export_name}_map.png")
            if export_map_png:
                if hasattr(self, "_hide_map_hover_tooltip_for_export"):
                    self._hide_map_hover_tooltip_for_export()
                self.figure.savefig(map_png_path, dpi=300, bbox_inches="tight", facecolor="white")

            profile_png_path = os.path.join(export_dir, f"{export_name}_profile.png")
            if export_profiles_png and hasattr(self, "profile_fig") and self.profile_fig is not None:
                self.profile_fig.savefig(profile_png_path, dpi=300, bbox_inches="tight", facecolor="white")

            msg_lines = [f"- {os.path.basename(geojson_file_path)}"]
            if export_text_csv:
                msg_lines.extend(
                    [
                        f"- {os.path.basename(csv_file_path)}",
                        f"- {os.path.basename(ddm_file_path)}",
                        f"- {os.path.basename(dms_file_path)}",
                    ]
                )
            if export_text_txt:
                msg_lines.extend(
                    [
                        f"- {os.path.basename(txt_file_path)}",
                        f"- {os.path.basename(ddm_txt_file_path)}",
                        f"- {os.path.basename(dms_txt_file_path)}",
                    ]
                )
            if export_shapefile:
                msg_lines.append(f"- {os.path.basename(shapefile_path)} (and sidecars)")
            if gpx_written:
                msg_lines.append(f"- {os.path.basename(gpx_file_path)}")
            for bn in gpx_per_test_names:
                msg_lines.append(f"- {bn}")
            if lnw_file_path:
                msg_lines.append(f"- {os.path.basename(lnw_file_path)}")
            if export_sis:
                msg_lines.append(f"- {os.path.basename(sis_file_path)}")
            msg_lines.append(f"- {os.path.basename(stats_file_path)}")
            if export_map_png:
                msg_lines.append(f"- {os.path.basename(map_png_path)}")
            msg_lines.append(f"- {os.path.basename(json_metadata_path)}")
            if export_profiles_png and hasattr(self, "profile_fig") and self.profile_fig is not None:
                msg_lines.append(f"- {os.path.basename(profile_png_path)}")
            log_msg = "Performance survey exported successfully to:\n" + "\n".join(msg_lines) + f"\nin directory: {export_dir}"
            if hasattr(self, "set_performance_activity_text"):
                self.set_performance_activity_text(log_msg, append=False)
            else:
                self.set_ref_info_text(log_msg, append=False)

        except Exception as e:
            self._show_message("error", "Export Error", f"Failed to export performance survey: {e}")

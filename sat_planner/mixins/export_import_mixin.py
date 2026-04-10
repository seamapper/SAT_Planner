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
                'line_length': float(self.line_length_entry.text()),
                'heading': float(self.heading_entry.text()),
                'dist_between_lines': float(self.dist_between_lines_entry.text()),
                'num_lines': int(self.num_lines_entry.text()),
                'bisect_lead': float(self.bisect_lead_entry.text()),
                'survey_speed': float(self.survey_speed_entry.text()),
                'export_name': self.export_name_entry.text().strip(),
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

            self.line_length_multiplier.set(params['line_length_multiplier'])
            self._update_multiplier_label_len(params['line_length_multiplier'])

            self.dist_between_lines_multiplier.set(params['dist_between_lines_multiplier'])
            self._update_multiplier_label_dist(params['dist_between_lines_multiplier'])

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
            if mapping is None:
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
            export_utils.write_ddd_csv(csv_file_path, ref_rows, newline='')
            ddm_file_path = os.path.join(export_dir, f"{export_name}_DMM.csv")
            export_utils.write_dmm_csv(ddm_file_path, ref_rows)
            dms_file_path = os.path.join(export_dir, f"{export_name}_DMS.csv")
            export_utils.write_dms_csv(dms_file_path, ref_rows)
            ddm_txt_file_path = os.path.join(export_dir, f"{export_name}_DMM.txt")
            export_utils.write_dmm_txt(ddm_txt_file_path, ref_rows)
            dms_txt_file_path = os.path.join(export_dir, f"{export_name}_DMS.txt")
            export_utils.write_dms_txt(dms_txt_file_path, ref_rows)
            txt_file_path = os.path.join(export_dir, f"{export_name}_DDD.txt")
            export_utils.write_ddd_txt(txt_file_path, ref_rows)

            # --- Export to ESRI Shapefile (.shp) ---
            schema = {
                'geometry': 'LineString',
                'properties': {'line_num': 'int'},
            }
            crs_epsg = 'EPSG:4326'  # WGS 84
            features = []
            # Add main survey lines
            for i, line_coords in enumerate(self.survey_lines_data):
                shapely_line = LineString([(p[1], p[0]) for p in line_coords])
                features.append({
                    'geometry': mapping(shapely_line),
                    'properties': {'line_num': i + 1},
                })
            # Add crossline
            if self.cross_line_data:
                shapely_cross_line = LineString([(p[1], p[0]) for p in self.cross_line_data])
                features.append({
                    'geometry': mapping(shapely_cross_line),
                    'properties': {'line_num': 0},
                })
            shapefile_path = os.path.join(export_dir, f"{export_name}.shp")
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
                        "geotiff_path": (self.current_geotiff_path if hasattr(self, 'current_geotiff_path') and self.current_geotiff_path else None),
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
                        "geotiff_path": (self.current_geotiff_path if hasattr(self, 'current_geotiff_path') and self.current_geotiff_path else None),
                        "survey_speed": ref_export_speed,
                        "points": [
                            {"point_num": 1, "lat": self.cross_line_data[0][0], "lon": self.cross_line_data[0][1]},
                            {"point_num": 2, "lat": self.cross_line_data[1][0], "lon": self.cross_line_data[1][1]}
                        ]
                    }
                })
            geojson_collection = {
                "type": "FeatureCollection",
                "properties": {
                    "geotiff_path": (self.current_geotiff_path if hasattr(self, 'current_geotiff_path') and self.current_geotiff_path else None)
                },
                "features": geojson_features
            }
            with open(geojson_file_path, 'w') as f:
                json.dump(geojson_collection, f, indent=2)

            self.set_ref_info_text(
                f"Survey exported successfully to:\n"
                f"- {os.path.basename(csv_file_path)}\n"
                f"- {os.path.basename(shapefile_path)} (and associated files)\n"
                f"- {os.path.basename(geojson_file_path)}\n"
                f"in directory: {export_dir}", append=False)

            # --- Export to Hypack LNW format (LIN/PTS/UTM), filename includes UTM zone ---
            lnw_file_path = None
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

            ref_msg = (f"Survey exported successfully to:\n"
                f"- {os.path.basename(csv_file_path)}\n"
                f"- {os.path.basename(shapefile_path)} (and associated files)\n"
                f"- {os.path.basename(geojson_file_path)}\n")
            if lnw_file_path:
                ref_msg += f"- {os.path.basename(lnw_file_path)}\n"
            ref_msg += f"in directory: {export_dir}"
            self.set_ref_info_text(ref_msg, append=False)

            # --- Export to Kongsberg SIS ASCII Plan format ---
            sis_file_path = os.path.join(export_dir, f"{export_name}.asciiplan")
            ref_ascii_lines = []
            if self.cross_line_data:
                ref_ascii_lines.append(('Crossline', [self.cross_line_data[0], self.cross_line_data[1]]))
            for i, line in enumerate(self.survey_lines_data):
                ref_ascii_lines.append((f'Reference{i + 1}', [line[0], line[1]]))
            export_utils.write_asciiplan(sis_file_path, ref_ascii_lines)

            # Update success message to include SIS
            ref_msg2 = (f"Survey exported successfully to:\n"
                f"- {os.path.basename(csv_file_path)}\n"
                f"- {os.path.basename(shapefile_path)} (and associated files)\n"
                f"- {os.path.basename(geojson_file_path)}\n")
            if lnw_file_path:
                ref_msg2 += f"- {os.path.basename(lnw_file_path)}\n"
            ref_msg2 += f"- {os.path.basename(sis_file_path)}\n"
            ref_msg2 += f"in directory: {export_dir}"
            self.set_ref_info_text(ref_msg2, append=False)

            # --- Export Comprehensive Survey Statistics ---
            stats_file_path = os.path.join(export_dir, f"{export_name}_info.txt")
            total_survey_time = self._calculate_total_survey_time()
            
            with open(stats_file_path, 'w') as f:
                f.write("COMPREHENSIVE SURVEY STATISTICS\n")
                f.write("=" * 50 + "\n\n")
                f.write(f"Survey Plan: {export_name}\n")
                f.write(f"Export Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                f.write("SURVEY INPUT PARAMETERS\n")
                f.write("-" * 30 + "\n")
                for key, value in values.items():
                    f.write(f"{key.replace('_', ' ').title()}: {value}\n")
                f.write("\n")
                
                f.write("SURVEY DISTANCE BREAKDOWN\n")
                f.write("-" * 30 + "\n")
                # Calculate single line length and heading
                num_main_lines = len(self.survey_lines_data) if self.survey_lines_data else 0
                if num_main_lines > 0:
                    # Calculate heading for the first main line
                    try:
                        geod = pyproj.Geod(ellps="WGS84")
                        first_line = self.survey_lines_data[0]
                        lat1, lon1 = first_line[0]
                        lat2, lon2 = first_line[1]
                        fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)
                        
                        heading = fwd_az % 360
                        reciprocal_heading = back_az % 360
                        
                        f.write(f"Heading: {heading:.1f}°\n")
                        f.write(f"Reciprocal Heading: {reciprocal_heading:.1f}°\n")
                    except Exception:
                        f.write("Heading: Unable to calculate\n")
                    
                    single_line_length_m = total_survey_time['main_lines_distance_m'] / num_main_lines
                    single_line_length_km = single_line_length_m / 1000
                    single_line_length_nm = single_line_length_m / 1852
                    f.write(f"Length of Main Line: {single_line_length_m:.1f} m ({single_line_length_km:.3f} km, {single_line_length_nm:.3f} nm)\n")
                f.write(f"Main Lines Total Distance: {total_survey_time['main_lines_distance_m']:.1f} m\n")
                f.write(f"Main Lines Total Distance: {total_survey_time['main_lines_distance_km']:.3f} km\n")
                f.write(f"Main Lines Total Distance: {total_survey_time['main_lines_distance_nm']:.3f} nm\n")
                f.write(f"Travel Between Lines: {total_survey_time['travel_between_lines_distance_m']:.1f} m\n")
                f.write(f"Travel Between Lines: {total_survey_time['travel_between_lines_distance_km']:.3f} km\n")
                f.write(f"Travel Between Lines: {total_survey_time['travel_between_lines_distance_nm']:.3f} nm\n")
                f.write(f"Travel to Crossline: {total_survey_time['travel_to_crossline_distance_m']:.1f} m\n")
                f.write(f"Travel to Crossline: {total_survey_time['travel_to_crossline_distance_km']:.3f} km\n")
                f.write(f"Travel to Crossline: {total_survey_time['travel_to_crossline_distance_nm']:.3f} nm\n")
                f.write(f"Crossline (per pass): {total_survey_time['crossline_single_pass_distance_m']:.1f} m\n")
                f.write(f"Crossline (per pass): {total_survey_time['crossline_single_pass_distance_km']:.3f} km\n")
                f.write(f"Crossline (per pass): {total_survey_time['crossline_single_pass_distance_nm']:.3f} nm\n")
                f.write(f"Crossline Passes: {total_survey_time['num_crossline_passes']}\n")
                
                # Add crossline heading information
                if self.cross_line_data and len(self.cross_line_data) >= 2:
                    try:
                        if pyproj is not None:
                            geod = pyproj.Geod(ellps="WGS84")
                            lat1, lon1 = self.cross_line_data[0]
                            lat2, lon2 = self.cross_line_data[1]
                            fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)
                            
                            crossline_heading = fwd_az % 360
                            crossline_reciprocal_heading = back_az % 360
                            
                            f.write(f"Crossline Heading: {crossline_heading:.1f}°\n")
                            f.write(f"Crossline Reciprocal Heading: {crossline_reciprocal_heading:.1f}°\n")
                        else:
                            f.write("Crossline Heading: pyproj not available\n")
                    except Exception:
                        f.write("Crossline Heading: Unable to calculate\n")
                
                f.write(f"Crossline Total Distance: {total_survey_time['crossline_total_distance_m']:.1f} m\n")
                f.write(f"Crossline Total Distance: {total_survey_time['crossline_total_distance_km']:.3f} km\n")
                f.write(f"Crossline Total Distance: {total_survey_time['crossline_total_distance_nm']:.3f} nm\n\n")
                
                f.write("TOTAL SURVEY DISTANCE\n")
                f.write("-" * 25 + "\n")
                f.write(f"Total Survey Distance: {total_survey_time['total_distance_m']:.1f} m\n")
                f.write(f"Total Survey Distance: {total_survey_time['total_distance_km']:.3f} km\n")
                f.write(f"Total Survey Distance: {total_survey_time['total_distance_nm']:.3f} nm\n\n")
                
                f.write("SURVEY TIME BREAKDOWN\n")
                f.write("-" * 25 + "\n")
                f.write(f"Main Lines Survey: {total_survey_time['main_lines_minutes']:.1f} min\n")
                f.write(f"Crossline Survey: {total_survey_time['crossline_minutes']:.1f} min\n")
                f.write(f"Travel Between Lines: {total_survey_time['travel_minutes']:.1f} min\n")
                f.write(f"Travel to Crossline: {total_survey_time['travel_to_crossline_minutes']:.1f} min\n")
                if 'total_turn_time_min' in total_survey_time:
                    num_turns = (total_survey_time.get('num_turns_between_main_lines', 0) + 
                               total_survey_time.get('num_turns_to_crossline', 0) + 
                               total_survey_time.get('num_turns_between_crossline_passes', 0))
                    turn_time_per_turn = total_survey_time.get('turn_time_per_turn_min', 5.0)
                    f.write(f"Total Turn Time: {total_survey_time['total_turn_time_min']:.1f} min ({num_turns} turns × {turn_time_per_turn:.1f} min)\n\n")
                else:
                    f.write("\n")
                
                # Calculate total survey time and transit time
                crossline_survey_time_only = total_survey_time['crossline_minutes']
                if 'num_turns_between_crossline_passes' in total_survey_time and total_survey_time['num_turns_between_crossline_passes'] > 0:
                    crossline_survey_time_only -= total_survey_time.get('turn_time_per_turn_min', 5.0) * total_survey_time['num_turns_between_crossline_passes']
                total_survey_time_only = total_survey_time['main_lines_minutes'] + crossline_survey_time_only
                
                travel_between_lines_time_only = total_survey_time['travel_minutes']
                if 'num_turns_between_main_lines' in total_survey_time and total_survey_time['num_turns_between_main_lines'] > 0:
                    travel_between_lines_time_only -= total_survey_time.get('turn_time_per_turn_min', 5.0) * total_survey_time['num_turns_between_main_lines']
                
                travel_to_crossline_time_only = total_survey_time['travel_to_crossline_minutes']
                if 'num_turns_to_crossline' in total_survey_time and total_survey_time['num_turns_to_crossline'] > 0:
                    travel_to_crossline_time_only -= total_survey_time.get('turn_time_per_turn_min', 5.0) * total_survey_time['num_turns_to_crossline']
                
                total_transit_time_only = travel_between_lines_time_only + travel_to_crossline_time_only
                
                f.write("TOTAL SURVEY TIME\n")
                f.write("-" * 20 + "\n")
                f.write(f"Total Survey Time: {total_survey_time_only:.1f} min ({total_survey_time_only / 60.0:.2f} hr)\n")
                f.write(f"Total Transit Time: {total_transit_time_only:.1f} min ({total_transit_time_only / 60.0:.2f} hr)\n")
                if 'total_turn_time_min' in total_survey_time:
                    f.write(f"Total Turn Time: {total_survey_time['total_turn_time_min']:.1f} min\n")
                f.write(f"Total Time: {total_survey_time['total_minutes']:.1f} min ({total_survey_time['total_hours']:.2f} hr)\n\n")
                
                f.write("SURVEY LINE DETAILS\n")
                f.write("-" * 20 + "\n")
                f.write(f"Number of Main Survey Lines: {len(self.survey_lines_data)}\n")
                if self.cross_line_data:
                    f.write("Crossline: Present\n")
                else:
                    f.write("Crossline: Not present\n")
                f.write(f"Survey Speed: {values.get('survey_speed', '8.0')} knots\n")
                if 'turn_time_per_turn_min' in total_survey_time:
                    f.write(f"Turn Time (per turn): {total_survey_time['turn_time_per_turn_min']:.1f} min\n")
                f.write(f"Crossline Passes: {total_survey_time['num_crossline_passes']}\n\n")
                
                f.write("SURVEY PATTERN\n")
                f.write("-" * 15 + "\n")
                f.write("Main lines are surveyed in a zigzag pattern:\n")
                f.write("- Even-numbered lines (2, 4, 6...): Survey from start to end\n")
                f.write("- Odd-numbered lines (1, 3, 5...): Survey from end to start (flipped)\n")
                f.write("- This minimizes travel distance between lines\n")
                f.write("- After completing all main lines, travel to crossline start\n")
                f.write("- Complete crossline survey with specified number of passes\n\n")
                
                # Accuracy Waypoints (DMM)
                dmm_heading = "Accuracy Waypoints (DMM)"
                f.write(f"\n{dmm_heading}\n")
                f.write("-" * len(dmm_heading) + "\n")
                if self.survey_lines_data:
                    for i, line in enumerate(self.survey_lines_data):
                        if i % 2 == 0:
                            start, end = line[0], line[1]
                        else:
                            start, end = line[1], line[0]
                        start_label = f'L{i+1}S'
                        end_label = f'L{i+1}E'
                        start_lat_ddm = decimal_degrees_to_ddm(start[0], is_latitude=True)
                        start_lon_ddm = decimal_degrees_to_ddm(start[1], is_latitude=False)
                        end_lat_ddm = decimal_degrees_to_ddm(end[0], is_latitude=True)
                        end_lon_ddm = decimal_degrees_to_ddm(end[1], is_latitude=False)
                        f.write(f"{start_label}: {start_lat_ddm}, {start_lon_ddm}\n")
                        f.write(f"{end_label}: {end_lat_ddm}, {end_lon_ddm}\n")
                if self.cross_line_data:
                    cls_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][0], is_latitude=True)
                    cls_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][1], is_latitude=False)
                    cle_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][0], is_latitude=True)
                    cle_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][1], is_latitude=False)
                    f.write(f"CLS: {cls_lat_ddm}, {cls_lon_ddm}\n")
                    f.write(f"CLE: {cle_lat_ddm}, {cle_lon_ddm}\n")

                # Accuracy Waypoints (DDD)
                ddd_heading = "Accuracy Waypoints (DDD)"
                f.write(f"\n{ddd_heading}\n")
                f.write("-" * len(ddd_heading) + "\n")
                if self.survey_lines_data:
                    for i, line in enumerate(self.survey_lines_data):
                        if i % 2 == 0:
                            start, end = line[0], line[1]
                        else:
                            start, end = line[1], line[0]
                        start_label = f'L{i+1}S'
                        end_label = f'L{i+1}E'
                        f.write(f"{start_label}: {start[0]:.6f}, {start[1]:.6f}\n")
                        f.write(f"{end_label}: {end[0]:.6f}, {end[1]:.6f}\n")
                if self.cross_line_data:
                    f.write(f"CLS: {self.cross_line_data[0][0]:.6f}, {self.cross_line_data[0][1]:.6f}\n")
                    f.write(f"CLE: {self.cross_line_data[1][0]:.6f}, {self.cross_line_data[1][1]:.6f}\n")
                f.write("\n")
                
                f.write("NOTES\n")
                f.write("-" * 5 + "\n")
                f.write("- All distances calculated using geodetic (spherical) geometry\n")
                f.write("- Times calculated based on specified survey speed\n")
                f.write("- Crossline lead-in/out distance extends beyond main survey area\n")
                f.write("- Survey pattern optimized for efficiency\n")

            # --- Export Survey Plan as PNG ---
            map_png_path = os.path.join(export_dir, f"{export_name}_map.png")
            if hasattr(self, "_hide_map_hover_tooltip_for_export"):
                self._hide_map_hover_tooltip_for_export()
            self.figure.savefig(map_png_path, dpi=300, bbox_inches='tight', facecolor='white')
            
            # --- Export Profile Plot as PNG (if it exists) ---
            profile_png_path = os.path.join(export_dir, f"{export_name}_profile.png")
            if hasattr(self, 'profile_fig') and self.profile_fig is not None:
                self.profile_fig.savefig(profile_png_path, dpi=300, bbox_inches='tight', facecolor='white')

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
            success_files = [
                f"- {os.path.basename(csv_file_path)}",
                f"- {os.path.basename(shapefile_path)} (and associated files)",
                f"- {os.path.basename(geojson_file_path)}",
            ]
            if lnw_file_path:
                success_files.append(f"- {os.path.basename(lnw_file_path)}")
            success_files.extend([
                f"- {os.path.basename(sis_file_path)}",
                f"- {os.path.basename(txt_file_path)}",
                f"- {os.path.basename(ddm_file_path)}",
                f"- {os.path.basename(dms_file_path)}",
                f"- {os.path.basename(ddm_txt_file_path)}",
                f"- {os.path.basename(dms_txt_file_path)}",
                f"- {os.path.basename(stats_file_path)}",
                f"- {os.path.basename(map_png_path)}"
            ])
            
            # Add metadata JSON if it was created
            try:
                if os.path.exists(json_metadata_path):
                    success_files.append(f"- {os.path.basename(json_metadata_path)}")
            except:
                pass
            
            # Add profile PNG if it was created
            if hasattr(self, 'profile_fig') and self.profile_fig is not None:
                success_files.append(f"- {os.path.basename(profile_png_path)}")
            
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
        if mapping is None:
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
            export_utils.write_ddd_csv(csv_file_path, perf_rows, newline="")
            ddm_file_path = os.path.join(export_dir, f"{export_name}_DMM.csv")
            export_utils.write_dmm_csv(ddm_file_path, perf_rows)
            dms_file_path = os.path.join(export_dir, f"{export_name}_DMS.csv")
            export_utils.write_dms_csv(dms_file_path, perf_rows)
            ddm_txt_file_path = os.path.join(export_dir, f"{export_name}_DMM.txt")
            export_utils.write_dmm_txt(ddm_txt_file_path, perf_rows)
            dms_txt_file_path = os.path.join(export_dir, f"{export_name}_DMS.txt")
            export_utils.write_dms_txt(dms_txt_file_path, perf_rows)
            txt_file_path = os.path.join(export_dir, f"{export_name}_DDD.txt")
            export_utils.write_ddd_txt(txt_file_path, perf_rows)

            schema = {"geometry": "LineString", "properties": {"line_num": "int"}}
            crs_epsg = "EPSG:4326"
            features = []
            for i, line_coords in enumerate(lines):
                shapely_line = LineString([(p[1], p[0]) for p in line_coords])
                features.append({"geometry": mapping(shapely_line), "properties": {"line_num": i + 1}})
            if has_bist:
                for i, seg in enumerate(bist_segs):
                    shapely_line = LineString([(p[1], p[0]) for p in seg])
                    features.append({"geometry": mapping(shapely_line), "properties": {"line_num": 11 + i}})
            shapefile_path = os.path.join(export_dir, f"{export_name}.shp")
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
            if lnw_lines:
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
            export_utils.write_asciiplan(sis_file_path, perf_ascii_lines)

            stats_file_path = os.path.join(export_dir, f"{export_name}_info.txt")
            with open(stats_file_path, "w", encoding="utf-8") as f:
                f.write("PERFORMANCE SURVEY EXPORT\n")
                f.write("=" * 40 + "\n\n")
                f.write(f"Plan name: {export_name}\n")
                f.write(f"Export date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                f.write("TAB PARAMETERS (at export)\n")
                f.write("-" * 28 + "\n")
                for label, attr in (
                    ("Central latitude", "performance_central_lat_entry"),
                    ("Central longitude", "performance_central_lon_entry"),
                    ("Swell direction (deg)", "performance_swell_direction_entry"),
                    ("Swath angle (deg)", "performance_swath_angle_entry"),
                    ("Test speed (kts)", "performance_test_speed_entry"),
                    ("BIST time (min)", "performance_bist_time_entry"),
                    ("Turn time (min)", "performance_turn_time_entry"),
                ):
                    w = getattr(self, attr, None)
                    if w is not None and hasattr(w, "text"):
                        try:
                            f.write(f"{label}: {w.text()}\n")
                        except Exception:
                            f.write(f"{label}: (n/a)\n")
                f.write("\nSEGMENTS\n")
                f.write("-" * 10 + "\n")
                f.write("Four performance swath lines (P1S–P4E).\n")
                if has_bist:
                    f.write("Four BIST collection segments (B1S–B4E).\n")
                else:
                    f.write("No BIST segments in this export.\n")
                f.write("\n")
                dmm_h = "Waypoints (DMM)"
                f.write(f"{dmm_h}\n")
                f.write("-" * len(dmm_h) + "\n")
                for i, line in enumerate(lines):
                    n = i + 1
                    s, e = line[0], line[1]
                    f.write(
                        f"P{n}S: {decimal_degrees_to_ddm(s[0], is_latitude=True)}, "
                        f"{decimal_degrees_to_ddm(s[1], is_latitude=False)}\n"
                    )
                    f.write(
                        f"P{n}E: {decimal_degrees_to_ddm(e[0], is_latitude=True)}, "
                        f"{decimal_degrees_to_ddm(e[1], is_latitude=False)}\n"
                    )
                if has_bist:
                    for i, seg in enumerate(bist_segs):
                        n = i + 1
                        b0, b1 = seg[0], seg[1]
                        f.write(
                            f"B{n}S: {decimal_degrees_to_ddm(b0[0], is_latitude=True)}, "
                            f"{decimal_degrees_to_ddm(b0[1], is_latitude=False)}\n"
                        )
                        f.write(
                            f"B{n}E: {decimal_degrees_to_ddm(b1[0], is_latitude=True)}, "
                            f"{decimal_degrees_to_ddm(b1[1], is_latitude=False)}\n"
                        )

            json_metadata_path = os.path.join(export_dir, f"{export_name}_performance_params.json")
            pmeta = {
                "export_name": export_name,
                "plan_type": "performance",
                "geotiff_path": geotiff_path,
                "test_speed_kts": speed_kts,
            }
            for key, attr in (
                ("central_lat", "performance_central_lat_entry"),
                ("central_lon", "performance_central_lon_entry"),
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
            if hasattr(self, "_hide_map_hover_tooltip_for_export"):
                self._hide_map_hover_tooltip_for_export()
            self.figure.savefig(map_png_path, dpi=300, bbox_inches="tight", facecolor="white")

            profile_png_path = os.path.join(export_dir, f"{export_name}_profile.png")
            if hasattr(self, "profile_fig") and self.profile_fig is not None:
                self.profile_fig.savefig(profile_png_path, dpi=300, bbox_inches="tight", facecolor="white")

            msg_lines = [
                f"- {os.path.basename(csv_file_path)}",
                f"- {os.path.basename(shapefile_path)} (and sidecars)",
                f"- {os.path.basename(geojson_file_path)}",
            ]
            if lnw_file_path:
                msg_lines.append(f"- {os.path.basename(lnw_file_path)}")
            msg_lines.extend(
                [
                    f"- {os.path.basename(sis_file_path)}",
                    f"- {os.path.basename(txt_file_path)}",
                    f"- {os.path.basename(ddm_file_path)}",
                    f"- {os.path.basename(dms_file_path)}",
                    f"- {os.path.basename(stats_file_path)}",
                    f"- {os.path.basename(map_png_path)}",
                    f"- {os.path.basename(json_metadata_path)}",
                ]
            )
            if hasattr(self, "profile_fig") and self.profile_fig is not None:
                msg_lines.append(f"- {os.path.basename(profile_png_path)}")
            log_msg = "Performance survey exported successfully to:\n" + "\n".join(msg_lines) + f"\nin directory: {export_dir}"
            if hasattr(self, "set_performance_activity_text"):
                self.set_performance_activity_text(log_msg, append=False)
            else:
                self.set_ref_info_text(log_msg, append=False)

        except Exception as e:
            self._show_message("error", "Export Error", f"Failed to export performance survey: {e}")

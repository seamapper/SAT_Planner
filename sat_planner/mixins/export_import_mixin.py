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
            export_name = f"Reference_{dist_between_lines}m_{heading_int}deg"
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

            # --- Export to CSV ---
            # CSV files for decimal degrees always use _DD suffix
            csv_file_path = os.path.join(export_dir, f"{export_name}_DD.csv")
            with open(csv_file_path, 'w', newline='') as csvfile:
                csv_writer = csv.writer(csvfile)
                csv_writer.writerow(['Line Number', 'Point Label', 'Latitude', 'Longitude'])
                # Main survey lines
                for i, line in enumerate(self.survey_lines_data):
                    if i % 2 == 0:  # Even lines - normal order
                        start, end = line[0], line[1]
                    else:  # Odd lines - flipped order
                        start, end = line[1], line[0]
                    start_label = f'L{i+1}S'
                    end_label = f'L{i+1}E'
                    csv_writer.writerow([i + 1, start_label, start[0], start[1]])
                    csv_writer.writerow([i + 1, end_label, end[0], end[1]])
                # Crossline
                if self.cross_line_data:
                    csv_writer.writerow([0, 'CLS', self.cross_line_data[0][0], self.cross_line_data[0][1]])
                    csv_writer.writerow([0, 'CLE', self.cross_line_data[1][0], self.cross_line_data[1][1]])
            
            # --- Export to DDM format (Decimal Minutes) ---
            ddm_file_path = os.path.join(export_dir, f"{export_name}_DM.csv")
            with open(ddm_file_path, 'w', newline='', encoding='utf-8') as ddmfile:
                ddm_writer = csv.writer(ddmfile)
                ddm_writer.writerow(['Line Number', 'Point Label', 'Latitude with Decimal Minutes', 'Longitude with Decimal Minutes'])
                # Main survey lines
                for i, line in enumerate(self.survey_lines_data):
                    if i % 2 == 0:  # Even lines - normal order
                        start, end = line[0], line[1]
                    else:  # Odd lines - flipped order
                        start, end = line[1], line[0]
                    start_label = f'L{i+1}S'
                    end_label = f'L{i+1}E'
                    start_lat_ddm = decimal_degrees_to_ddm(start[0], is_latitude=True)
                    start_lon_ddm = decimal_degrees_to_ddm(start[1], is_latitude=False)
                    end_lat_ddm = decimal_degrees_to_ddm(end[0], is_latitude=True)
                    end_lon_ddm = decimal_degrees_to_ddm(end[1], is_latitude=False)
                    ddm_writer.writerow([i + 1, start_label, start_lat_ddm, start_lon_ddm])
                    ddm_writer.writerow([i + 1, end_label, end_lat_ddm, end_lon_ddm])
                # Crossline
                if self.cross_line_data:
                    cls_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][0], is_latitude=True)
                    cls_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][1], is_latitude=False)
                    cle_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][0], is_latitude=True)
                    cle_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][1], is_latitude=False)
                    ddm_writer.writerow([0, 'CLS', cls_lat_ddm, cls_lon_ddm])
                    ddm_writer.writerow([0, 'CLE', cle_lat_ddm, cle_lon_ddm])
            
            # --- Export to DDM text format (Decimal Minutes) ---
            ddm_txt_file_path = os.path.join(export_dir, f"{export_name}_DM.txt")
            with open(ddm_txt_file_path, 'w', encoding='utf-8') as ddm_txt_file:
                # Main survey lines
                for i, line in enumerate(self.survey_lines_data):
                    if i % 2 == 0:  # Even lines - normal order
                        start, end = line[0], line[1]
                    else:  # Odd lines - flipped order
                        start, end = line[1], line[0]
                    start_label = f'L{i+1}S'
                    end_label = f'L{i+1}E'
                    start_lat_ddm = decimal_degrees_to_ddm(start[0], is_latitude=True)
                    start_lon_ddm = decimal_degrees_to_ddm(start[1], is_latitude=False)
                    end_lat_ddm = decimal_degrees_to_ddm(end[0], is_latitude=True)
                    end_lon_ddm = decimal_degrees_to_ddm(end[1], is_latitude=False)
                    ddm_txt_file.write(f"{start_label}, {start_lat_ddm}, {start_lon_ddm}\n")
                    ddm_txt_file.write(f"{end_label}, {end_lat_ddm}, {end_lon_ddm}\n")
                # Crossline
                if self.cross_line_data:
                    cls_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][0], is_latitude=True)
                    cls_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][1], is_latitude=False)
                    cle_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][0], is_latitude=True)
                    cle_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][1], is_latitude=False)
                    ddm_txt_file.write(f"CLS, {cls_lat_ddm}, {cls_lon_ddm}\n")
                    ddm_txt_file.write(f"CLE, {cle_lat_ddm}, {cle_lon_ddm}\n")

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
                        "points": [
                            {"point_num": 1, "lat": self.cross_line_data[0][0], "lon": self.cross_line_data[0][1]},
                            {"point_num": 2, "lat": self.cross_line_data[1][0], "lon": self.cross_line_data[1][1]}
                        ]
                    }
                })
            geojson_collection = {
                "type": "FeatureCollection",
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

            # --- Export to Hypack LNW format ---
            lnw_file_path = os.path.join(export_dir, f"{export_name}.lnw")
            with open(lnw_file_path, 'w') as f:
                f.write("LNW 1.0\n")
                # Main survey lines
                for i, line in enumerate(self.survey_lines_data):
                    line_num = i + 1
                    # Start point
                    lat1, lon1 = line[0]
                    lat2, lon2 = line[1]
                    # Get depth from GeoTIFF if available
                    depth1 = self._get_depth_at_point(lat1, lon1) if self.geotiff_data_array is not None else 0.0
                    depth2 = self._get_depth_at_point(lat2, lon2) if self.geotiff_data_array is not None else 0.0
                    # Get speed from survey speed entry
                    try:
                        speed = float(self.survey_speed_entry.text()) if self.survey_speed_entry.text() else 8.0
                    except:
                        speed = 8.0
                    f.write(f"LINE{line_num:03d}_001,{lat1:.6f},{lon1:.6f},{abs(depth1):.1f},{speed:.1f},50.0,{line_num},1\n")
                    f.write(f"LINE{line_num:03d}_002,{lat2:.6f},{lon2:.6f},{abs(depth2):.1f},{speed:.1f},50.0,{line_num},1\n")
                # Crossline
                if self.cross_line_data:
                    lat1, lon1 = self.cross_line_data[0]
                    lat2, lon2 = self.cross_line_data[1]
                    depth1 = self._get_depth_at_point(lat1, lon1) if self.geotiff_data_array is not None else 0.0
                    depth2 = self._get_depth_at_point(lat2, lon2) if self.geotiff_data_array is not None else 0.0
                    try:
                        speed = float(self.survey_speed_entry.text()) if self.survey_speed_entry.text() else 8.0
                    except:
                        speed = 8.0
                    f.write(f"CROSS_001,{lat1:.6f},{lon1:.6f},{abs(depth1):.1f},{speed:.1f},50.0,0,2\n")
                    f.write(f"CROSS_002,{lat2:.6f},{lon2:.6f},{abs(depth2):.1f},{speed:.1f},50.0,0,2\n")

            self.set_ref_info_text(
                f"Survey exported successfully to:\n"
                f"- {os.path.basename(csv_file_path)}\n"
                f"- {os.path.basename(shapefile_path)} (and associated files)\n"
                f"- {os.path.basename(geojson_file_path)}\n"
                f"- {os.path.basename(lnw_file_path)}\n"
                f"in directory: {export_dir}", append=False)

            # --- Export to SIS ASCII Plan format ---
            sis_file_path = os.path.join(export_dir, f"{export_name}.asciiplan")
            with open(sis_file_path, 'w') as f:
                f.write("SIS ASCII Plan\n")
                # Main survey lines
                for i, line in enumerate(self.survey_lines_data):
                    line_num = i + 1
                    # Start point
                    lat1, lon1 = line[0]
                    lat2, lon2 = line[1]
                    # Get depth from GeoTIFF if available
                    depth1 = self._get_depth_at_point(lat1, lon1) if self.geotiff_data_array is not None else 0.0
                    depth2 = self._get_depth_at_point(lat2, lon2) if self.geotiff_data_array is not None else 0.0
                    # Get speed from survey speed entry
                    try:
                        speed = float(self.survey_speed_entry.text()) if self.survey_speed_entry.text() else 8.0
                    except:
                        speed = 8.0
                    f.write(f"LINE{line_num:03d}_001, {lat1:.6f}, {lon1:.6f}, {abs(depth1):.1f}, {speed:.1f}, {line_num}, {line_num}\n")
                    f.write(f"LINE{line_num:03d}_002, {lat2:.6f}, {lon2:.6f}, {abs(depth2):.1f}, {speed:.1f}, {line_num}, {line_num}\n")
                
                # Crossline if present
                if self.cross_line_data:
                    lat1, lon1 = self.cross_line_data[0]
                    lat2, lon2 = self.cross_line_data[1]
                    depth1 = self._get_depth_at_point(lat1, lon1) if self.geotiff_data_array is not None else 0.0
                    depth2 = self._get_depth_at_point(lat2, lon2) if self.geotiff_data_array is not None else 0.0
                    try:
                        speed = float(self.survey_speed_entry.text()) if self.survey_speed_entry.text() else 8.0
                    except:
                        speed = 8.0
                    crossline_num = len(self.survey_lines_data) + 1
                    f.write(f"CROSSLINE_001, {lat1:.6f}, {lon1:.6f}, {abs(depth1):.1f}, {speed:.1f}, {crossline_num}, {crossline_num}\n")
                    f.write(f"CROSSLINE_002, {lat2:.6f}, {lon2:.6f}, {abs(depth2):.1f}, {speed:.1f}, {crossline_num}, {crossline_num}\n")

            # Update success message to include SIS
            self.set_ref_info_text(
                f"Survey exported successfully to:\n"
                f"- {os.path.basename(csv_file_path)}\n"
                f"- {os.path.basename(shapefile_path)} (and associated files)\n"
                f"- {os.path.basename(geojson_file_path)}\n"
                f"- {os.path.basename(lnw_file_path)}\n"
                f"- {os.path.basename(sis_file_path)}\n"
                f"in directory: {export_dir}", append=False)

            # --- Export Comprehensive Survey Statistics ---
            stats_file_path = os.path.join(export_dir, f"{export_name}_stats.txt")
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
                f.write(f"Travel to Crossline: {total_survey_time['travel_to_crossline_minutes']:.1f} min\n\n")
                
                f.write("TOTAL SURVEY TIME\n")
                f.write("-" * 20 + "\n")
                f.write(f"Total Survey Time: {total_survey_time['total_minutes']:.1f} min\n")
                f.write(f"Total Survey Time: {total_survey_time['total_hours']:.2f} hr\n\n")
                
                f.write("SURVEY LINE DETAILS\n")
                f.write("-" * 20 + "\n")
                f.write(f"Number of Main Survey Lines: {len(self.survey_lines_data)}\n")
                if self.cross_line_data:
                    f.write("Crossline: Present\n")
                else:
                    f.write("Crossline: Not present\n")
                f.write(f"Survey Speed: {values.get('survey_speed', '8.0')} knots\n")
                f.write(f"Crossline Passes: {total_survey_time['num_crossline_passes']}\n\n")
                
                f.write("SURVEY PATTERN\n")
                f.write("-" * 15 + "\n")
                f.write("Main lines are surveyed in a zigzag pattern:\n")
                f.write("- Even-numbered lines (2, 4, 6...): Survey from start to end\n")
                f.write("- Odd-numbered lines (1, 3, 5...): Survey from end to start (flipped)\n")
                f.write("- This minimizes travel distance between lines\n")
                f.write("- After completing all main lines, travel to crossline start\n")
                f.write("- Complete crossline survey with specified number of passes\n\n")
                
                # Waypoints with labels
                f.write("WAYPOINTS\n")
                f.write("-" * 10 + "\n")
                if self.survey_lines_data:
                    for i, line in enumerate(self.survey_lines_data):
                        # Determine start and end based on zigzag pattern
                        if i % 2 == 0:  # Even lines - normal order
                            start, end = line[0], line[1]
                        else:  # Odd lines - flipped order
                            start, end = line[1], line[0]
                        start_label = f'L{i+1}S'
                        end_label = f'L{i+1}E'
                        start_lat_ddm = decimal_degrees_to_ddm(start[0], is_latitude=True)
                        start_lon_ddm = decimal_degrees_to_ddm(start[1], is_latitude=False)
                        end_lat_ddm = decimal_degrees_to_ddm(end[0], is_latitude=True)
                        end_lon_ddm = decimal_degrees_to_ddm(end[1], is_latitude=False)
                        f.write(f"{start_label}: {start_lat_ddm}, {start_lon_ddm} ({start[0]:.6f}, {start[1]:.6f})\n")
                        f.write(f"{end_label}: {end_lat_ddm}, {end_lon_ddm} ({end[0]:.6f}, {end[1]:.6f})\n")
                if self.cross_line_data:
                    cls_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][0], is_latitude=True)
                    cls_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][1], is_latitude=False)
                    cle_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][0], is_latitude=True)
                    cle_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][1], is_latitude=False)
                    f.write(f"CLS: {cls_lat_ddm}, {cls_lon_ddm} ({self.cross_line_data[0][0]:.6f}, {self.cross_line_data[0][1]:.6f})\n")
                    f.write(f"CLE: {cle_lat_ddm}, {cle_lon_ddm} ({self.cross_line_data[1][0]:.6f}, {self.cross_line_data[1][1]:.6f})\n")
                f.write("\n")
                
                f.write("NOTES\n")
                f.write("-" * 5 + "\n")
                f.write("- All distances calculated using geodetic (spherical) geometry\n")
                f.write("- Times calculated based on specified survey speed\n")
                f.write("- Crossline lead-in/out distance extends beyond main survey area\n")
                f.write("- Survey pattern optimized for efficiency\n")

            # --- Export to Text format (point label, latitude, longitude) ---
            txt_file_path = os.path.join(export_dir, f"{export_name}_DD.txt")
            with open(txt_file_path, 'w') as f:
                # Main survey lines
                for i, line in enumerate(self.survey_lines_data):
                    if i % 2 == 0:  # Even lines - normal order
                        start, end = line[0], line[1]
                    else:  # Odd lines - flipped order
                        start, end = line[1], line[0]
                    start_label = f'L{i+1}S'
                    end_label = f'L{i+1}E'
                    f.write(f"{start_label} {start[0]:.6f} {start[1]:.6f}\n")
                    f.write(f"{end_label} {end[0]:.6f} {end[1]:.6f}\n")
                # Crossline
                if self.cross_line_data:
                    f.write(f"CLS {self.cross_line_data[0][0]:.6f} {self.cross_line_data[0][1]:.6f}\n")
                    f.write(f"CLE {self.cross_line_data[1][0]:.6f} {self.cross_line_data[1][1]:.6f}\n")

            # --- Export Survey Plan as PNG ---
            survey_plan_png_path = os.path.join(export_dir, f"{export_name}_survey_plan.png")
            self.figure.savefig(survey_plan_png_path, dpi=300, bbox_inches='tight', facecolor='white')
            
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
                f"- {os.path.basename(lnw_file_path)}",
                f"- {os.path.basename(sis_file_path)}",
                f"- {os.path.basename(txt_file_path)}",
                f"- {os.path.basename(ddm_file_path)}",
                f"- {os.path.basename(ddm_txt_file_path)}",
                f"- {os.path.basename(stats_file_path)}",
                f"- {os.path.basename(survey_plan_png_path)}"
            ]
            
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

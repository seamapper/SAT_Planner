"""
Reference tab: survey lines validation, export/import, info text, reset, statistics.
"""
import os
import csv
import json

from PyQt6.QtWidgets import QFileDialog
from PyQt6.QtGui import QTextCursor

from sat_planner.constants import GEOSPATIAL_LIBS_AVAILABLE, pyproj, LineString, fiona
from sat_planner import decimal_degrees_to_ddm


class ReferenceMixin:
    """Mixin for reference planning: _validate_inputs, _export_survey_data, _import_survey_files,
    set_ref_info_text, _reset_reference_tab, _calculate_total_survey_time,
    _calculate_reference_survey_statistics, _show_reference_planning_info."""

    def _validate_inputs(self):
        """Validates all input fields and returns a tuple (is_valid, values_dict)."""
        values = {}
        try:
            # Central Latitude
            central_lat_str = self.central_lat_entry.text().strip()
            if not central_lat_str:
                self._show_message("error","Input Error", "Central Latitude cannot be empty.")
                return False, {}
            values['central_lat'] = float(central_lat_str)
            if not (-90 <= values['central_lat'] <= 90):
                self._show_message("error","Input Error", "Central Latitude must be between -90 and 90.")
                return False, {}

            # Central Longitude
            central_lon_str = self.central_lon_entry.text().strip()
            if not central_lon_str:
                self._show_message("error","Input Error", "Central Longitude cannot be empty.")
                return False, {}
            values['central_lon'] = float(central_lon_str)
            if not (-180 <= values['central_lon'] <= 180):
                self._show_message("error","Input Error", "Central Longitude must be between -180 and 180.")
                return False, {}

            # Line Length
            line_length_str = self.line_length_entry.text().strip()
            if not line_length_str:
                self._show_message("error","Input Error", "Line Length cannot be empty.")
                return False, {}
            values['line_length'] = float(line_length_str)
            if not (values['line_length'] > 0):
                self._show_message("error","Input Error", "Line Length must be positive.")
                return False, {}

            # Heading
            heading_str = self.heading_entry.text().strip()
            if not heading_str:
                self._show_message("error","Input Error", "Heading cannot be empty.")
                return False, {}
            values['heading'] = float(heading_str)
            if not (0 <= values['heading'] <= 360):
                self._show_message("error","Input Error", "Heading must be between 0 and 360 degrees.")
                return False, {}

            # Distance Between Lines
            dist_between_lines_str = self.dist_between_lines_entry.text().strip()
            if not dist_between_lines_str:
                self._show_message("error","Input Error", "Distance Between Lines cannot be empty.")
                return False, {}
            values['dist_between_lines'] = float(dist_between_lines_str)
            if not (values['dist_between_lines'] > 0):
                self._show_message("error","Input Error", "Distance Between Lines must be positive.")
                return False, {}

            # Number of Lines
            num_lines_str = self.num_lines_entry.text().strip()
            if not num_lines_str:
                self._show_message("error","Input Error", "Number of Lines cannot be empty.")
                return False, {}
            values['num_lines'] = int(num_lines_str)
            if not (values['num_lines'] > 0):
                self._show_message("error","Input Error", "Number of Lines must be a positive integer.")
                return False, {}

            # Crossline Lead-in/out
            bisect_lead_str = self.bisect_lead_entry.text().strip()
            if not bisect_lead_str:
                self._show_message("error","Input Error", "Crossline Lead-in/out cannot be empty.")
                return False, {}
            values['bisect_lead'] = float(bisect_lead_str)
            if not (values['bisect_lead'] >= 0):
                self._show_message("error","Input Error", "Crossline Lead-in/out cannot be negative.")
                return False, {}

            # Survey Speed
            survey_speed_str = self.survey_speed_entry.text().strip()
            if not survey_speed_str:
                self._show_message("error","Input Error", "Survey Speed cannot be empty.")
                return False, {}
            values['survey_speed'] = float(survey_speed_str)
            if not (values['survey_speed'] > 0):
                self._show_message("error","Input Error", "Survey Speed must be positive.")
                return False, {}

            # Export Name
            values['export_name'] = self.export_name_entry.text().strip()
            if not values['export_name']:
                self._show_message("error","Input Error", "Export Name cannot be empty.")
                return False, {}

            values['offset_direction'] = self.offset_direction_var
            values['line_length_multiplier'] = self.line_length_multiplier  # Get slider value
            values['dist_between_lines_multiplier'] = self.dist_between_lines_multiplier  # Get new slider value

        except ValueError as e:
            self._show_message("error","Input Error", f"Invalid numeric input: {str(e)}")
            return False, {}
        except Exception as e:
            self._show_message("error","Error", f"An unexpected error occurred during input validation: {str(e)}")
            return False, {}

        return True, values

    def _export_survey_data(self):
        if not GEOSPATIAL_LIBS_AVAILABLE:
            self._show_message("warning","Disabled Feature", "Geospatial libraries not loaded. Cannot export data.")
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

            # --- Export Input Parameters to TXT ---
            params_file_path = os.path.join(export_dir, f"{export_name}_params_DD.txt")
            with open(params_file_path, 'w') as f:
                f.write("Survey Plan Input Parameters:\n\n")
                for key, value in values.items():
                    f.write(f"{key.replace('_', ' ').title()}: {value}\n")

            # --- Export to ESRI Shapefile (.shp) ---
            if LineString is not None and fiona is not None:
                schema = {
                    'geometry': 'LineString',
                    'properties': {'line_num': 'int'},
                }
                crs_epsg = 'EPSG:4326'  # WGS 84

                features = []
                # Add main survey lines
                for i, line_coords in enumerate(self.survey_lines_data):
                    # Shapely expects (lon, lat) order
                    shapely_line = LineString([(p[1], p[0]) for p in line_coords])
                    features.append({
                        'geometry': shapely_line.__geo_interface__,
                        'properties': {'line_num': i + 1},
                    })

                # Add bisecting line
                if self.cross_line_data:
                    shapely_cross_line = LineString([(p[1], p[0]) for p in self.cross_line_data])
                    features.append({
                        'geometry': shapely_cross_line.__geo_interface__,
                        'properties': {'line_num': 0},  # Using 0 for crossline
                    })

                shapefile_path = os.path.join(export_dir, f"{export_name}.shp")
                with fiona.open(shapefile_path, 'w', driver='ESRI Shapefile', crs=crs_epsg, schema=schema) as collection:
                    collection.writerecords(features)
                self.set_ref_info_text(
                    f"Data exported successfully to:\n"
                    f"- {os.path.basename(csv_file_path)}\n"
                    f"- {os.path.basename(params_file_path)}\n"
                    f"- {os.path.basename(shapefile_path)} (and associated files)\n"
                    f"in directory: {export_dir}", append=False)
            else:
                self.set_ref_info_text(
                    f"Data exported successfully to:\n"
                    f"- {os.path.basename(csv_file_path)}\n"
                    f"- {os.path.basename(params_file_path)}\n"
                    f"(Shapefile export requires fiona/shapely)\n"
                    f"in directory: {export_dir}", append=False)

        except Exception as e:
            self._show_message("error","Export Error", f"Failed to export data: {e}")

    def _import_survey_files(self):
        """Import reference planning survey lines from CSV or GeoJSON file."""
        # Open file dialog to select import file
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Survey File to Import",
            self.last_ref_import_dir,
            "Decimal Degree CSV files (*_DD.csv);;CSV files (*.csv);;GeoJSON files (*.geojson);;JSON files (*.json);;All files (*.*)"
        )

        if not file_path:
            return

        # Save the directory for next time
        import_dir = os.path.dirname(file_path)
        if import_dir and os.path.isdir(import_dir):
            self.last_ref_import_dir = import_dir
            self._save_last_ref_import_dir()

        try:
            # Clear existing lines
            self.survey_lines_data = []
            self.cross_line_data = []

            # Determine file type and import accordingly
            file_ext = os.path.splitext(file_path)[1].lower()

            if file_ext == '.csv':
                # Import from CSV
                with open(file_path, 'r', encoding='utf-8') as csvfile:
                    csv_reader = csv.DictReader(csvfile)

                    # Dictionary to collect points for each line
                    lines_data = {}
                    crossline_points = {}

                    for row in csv_reader:
                        try:
                            line_num = int(row.get('Line Number', -1))
                            point_label = row.get('Point Label', '').strip()
                            lat = float(row.get('Latitude', 0))
                            lon = float(row.get('Longitude', 0))
                        except (ValueError, TypeError):
                            continue

                        # Handle crossline (line_num == 0)
                        if line_num == 0:
                            if point_label == 'CLS':
                                crossline_points['start'] = (lat, lon)
                            elif point_label == 'CLE':
                                crossline_points['end'] = (lat, lon)
                        else:
                            # Main survey line
                            if line_num not in lines_data:
                                lines_data[line_num] = {}
                            lines_data[line_num][point_label] = (lat, lon)

                    # Process main survey lines
                    for line_num in sorted(lines_data.keys()):
                        points = lines_data[line_num]
                        # Look for L{n}S and L{n}E labels
                        start_label = f'L{line_num}S'
                        end_label = f'L{line_num}E'
                        if start_label in points and end_label in points:
                            self.survey_lines_data.append([points[start_label], points[end_label]])

                    # Process crossline
                    if 'start' in crossline_points and 'end' in crossline_points:
                        self.cross_line_data = [crossline_points['start'], crossline_points['end']]

            elif file_ext in ['.geojson', '.json']:
                # Import from GeoJSON
                with open(file_path, 'r', encoding='utf-8') as f:
                    geojson_data = json.load(f)

                # Ensure geojson_data is a dictionary
                if not isinstance(geojson_data, dict):
                    self._show_message("error","Import Error", "Invalid GeoJSON format: root must be an object")
                    return

                # Handle FeatureCollection
                if geojson_data.get('type') == 'FeatureCollection':
                    features = geojson_data.get('features', [])
                    if not isinstance(features, list):
                        features = []
                elif geojson_data.get('type') == 'Feature':
                    features = [geojson_data]
                else:
                    features = []

                # Sort features by line_num to maintain order
                features_with_num = []
                for feature in features:
                    if not isinstance(feature, dict) or feature.get('type') != 'Feature':
                        continue
                    geometry = feature.get('geometry', {})
                    if not isinstance(geometry, dict) or geometry.get('type') != 'LineString':
                        continue
                    properties = feature.get('properties', {})
                    if not isinstance(properties, dict):
                        properties = {}
                    line_num = properties.get('line_num', 0)
                    features_with_num.append((line_num, feature))

                # Sort by line_num
                features_with_num.sort(key=lambda x: x[0])

                for line_num, feature in features_with_num:
                    if not isinstance(feature, dict):
                        continue
                    geometry = feature.get('geometry', {})
                    if not isinstance(geometry, dict):
                        continue
                    coordinates = geometry.get('coordinates', [])
                    if not isinstance(coordinates, list) or len(coordinates) < 2:
                        continue

                    # GeoJSON uses [lon, lat] format, we need [lat, lon]
                    point1 = (coordinates[0][1], coordinates[0][0])
                    point2 = (coordinates[-1][1], coordinates[-1][0])

                    if line_num == 0:
                        # Crossline
                        self.cross_line_data = [point1, point2]
                    else:
                        # Main survey line
                        self.survey_lines_data.append([point1, point2])

            else:
                self._show_message("error","Import Error", f"Unsupported file format: {file_ext}")
                return

            if not self.survey_lines_data and not self.cross_line_data:
                self._show_message("warning","Import Warning", "No valid survey lines found in the selected file.")
                return

            # Try to load parameters metadata JSON file
            base_name = os.path.splitext(os.path.basename(file_path))[0]
            dir_name = os.path.dirname(file_path)
            metadata_path = os.path.join(dir_name, f"{base_name}_params.json")
            params = None

            if os.path.exists(metadata_path):
                try:
                    with open(metadata_path, 'r', encoding='utf-8') as f:
                        params = json.load(f)
                except Exception as e:
                    print(f"Warning: Could not load metadata file: {e}")

            # Populate parameter fields from metadata or calculate from lines
            try:
                if params:
                    # Use metadata if available
                    if params.get('central_lat') is not None:
                        self.central_lat_entry.clear()
                        self.central_lat_entry.setText(str(params['central_lat']))
                        if hasattr(self, 'pick_center_btn'):
                            self.pick_center_btn.setStyleSheet("")

                    if params.get('central_lon') is not None:
                        self.central_lon_entry.clear()
                        self.central_lon_entry.setText(str(params['central_lon']))

                    if params.get('line_length') is not None:
                        self.line_length_entry.setText(str(params['line_length']))

                    if params.get('heading') is not None:
                        self.heading_entry.setText(str(params['heading']))

                    if params.get('dist_between_lines') is not None:
                        self.dist_between_lines_entry.setText(str(params['dist_between_lines']))

                    if params.get('num_lines') is not None:
                        self.num_lines_entry.setText(str(params['num_lines']))

                    if params.get('bisect_lead') is not None:
                        self.bisect_lead_entry.setText(str(params['bisect_lead']))

                    if params.get('survey_speed') is not None:
                        self.survey_speed_entry.setText(str(params['survey_speed']))

                    if params.get('crossline_passes') is not None:
                        self.crossline_passes_entry.setText(str(params['crossline_passes']))

                    if params.get('export_name'):
                        self.export_name_entry.clear()
                        self.export_name_entry.setText(params['export_name'])

                    if params.get('offset_direction'):
                        direction = params['offset_direction']
                        if hasattr(self, 'offset_direction_combo'):
                            idx = self.offset_direction_combo.findText(direction)
                            if idx >= 0:
                                self.offset_direction_combo.setCurrentIndex(idx)
                        self.offset_direction_var = direction

                    if params.get('line_length_multiplier') is not None:
                        self.line_length_multiplier.set(params['line_length_multiplier'])
                        self._update_multiplier_label_len(params['line_length_multiplier'])

                    if params.get('dist_between_lines_multiplier') is not None:
                        self.dist_between_lines_multiplier.set(params['dist_between_lines_multiplier'])
                        self._update_multiplier_label_dist(params['dist_between_lines_multiplier'])
                else:
                    # Calculate parameters from imported lines
                    if pyproj is not None:
                        geod = pyproj.Geod(ellps="WGS84")

                        # Calculate central lat/lon from all points
                        all_points = []
                        for line in self.survey_lines_data:
                            all_points.extend(line)
                        if self.cross_line_data:
                            all_points.extend(self.cross_line_data)

                        if all_points:
                            all_lats = [p[0] for p in all_points]
                            all_lons = [p[1] for p in all_points]
                            central_lat = (min(all_lats) + max(all_lats)) / 2.0
                            central_lon = (min(all_lons) + max(all_lons)) / 2.0

                            self.central_lat_entry.clear()
                            self.central_lat_entry.setText(f"{central_lat:.6f}")
                            if hasattr(self, 'pick_center_btn'):
                                self.pick_center_btn.setStyleSheet("")

                            self.central_lon_entry.clear()
                            self.central_lon_entry.setText(f"{central_lon:.6f}")

                        # Calculate heading from first line
                        if len(self.survey_lines_data) > 0:
                            first_line = self.survey_lines_data[0]
                            try:
                                lat1, lon1 = first_line[0]
                                lat2, lon2 = first_line[1]
                                fwd_az, back_az, dist = geod.inv(lon1, lat1, lon2, lat2)
                                heading = fwd_az % 360

                                self.heading_entry.clear()
                                self.heading_entry.setText(f"{heading:.1f}")

                                line_length = dist
                                self.line_length_entry.clear()
                                self.line_length_entry.setText(f"{line_length:.1f}")
                            except Exception:
                                pass

                        # Count number of lines
                        if self.survey_lines_data:
                            self.num_lines_entry.clear()
                            self.num_lines_entry.setText(str(len(self.survey_lines_data)))

                        # Calculate distance between lines (approximate)
                        if len(self.survey_lines_data) > 1:
                            try:
                                line1_mid = ((self.survey_lines_data[0][0][0] + self.survey_lines_data[0][1][0]) / 2,
                                           (self.survey_lines_data[0][0][1] + self.survey_lines_data[0][1][1]) / 2)
                                line2_mid = ((self.survey_lines_data[1][0][0] + self.survey_lines_data[1][1][0]) / 2,
                                           (self.survey_lines_data[1][0][1] + self.survey_lines_data[1][1][1]) / 2)
                                _, _, dist = geod.inv(line1_mid[1], line1_mid[0], line2_mid[1], line2_mid[0])
                                self.dist_between_lines_entry.clear()
                                self.dist_between_lines_entry.setText(f"{dist:.1f}")
                            except Exception:
                                pass

                    # Update export name if not set
                    if not self.export_name_entry.text().strip():
                        self.export_name_entry.clear()
                        self.export_name_entry.setText(base_name)
            except Exception as e:
                print(f"Warning: Error populating parameter fields: {e}")

            # Update plot
            self._plot_survey_plan(preserve_view_limits=True)

            # Show success message
            imported_items = []
            if self.survey_lines_data:
                imported_items.append(f"{len(self.survey_lines_data)} survey line(s)")
            if self.cross_line_data:
                imported_items.append("crossline")

            if imported_items:
                msg = f"Successfully imported: {', '.join(imported_items)}"
                if params:
                    msg += " (parameters loaded from metadata)"
                else:
                    msg += " (parameters calculated from lines)"
                self.set_ref_info_text(msg)
            else:
                self._show_message("warning","Import Warning", "No valid survey lines found in the selected file.")

        except Exception as e:
            self._show_message("error","Import Error", f"Failed to import survey files: {e}")
            import traceback
            traceback.print_exc()

    def set_ref_info_text(self, message, append=False):
        """Add a message to the Activity Log with [Reference] prefix. Maintains up to 200 lines of history."""
        if not hasattr(self, 'activity_log_text'):
            return
        prefixed_message = f"[Reference] {message}"
        self.activity_log_text.setReadOnly(False)
        if append:
            # Prepend new message (newest at top) - only when explicitly requested
            current_text = self.activity_log_text.toPlainText().rstrip('\n')
            new_text = prefixed_message + "\n" + current_text if current_text else prefixed_message + "\n"
            # Split into lines and limit to 200
            lines = new_text.splitlines()
            if len(lines) > 200:
                lines = lines[:200]
            self.activity_log_text.setPlainText("\n".join(lines) + "\n")
            # Move cursor to top
            cursor = self.activity_log_text.textCursor()
            cursor.movePosition(QTextCursor.MoveOperation.Start)
            self.activity_log_text.setTextCursor(cursor)
        else:
            # Append to end (default behavior - newest at bottom)
            self.activity_log_text.append(prefixed_message)
            # Maintain up to 200 lines by removing oldest from top
            num_lines = self.activity_log_text.document().blockCount()
            if num_lines > 200:
                cursor = self.activity_log_text.textCursor()
                cursor.movePosition(QTextCursor.MoveOperation.Start)
                cursor.movePosition(QTextCursor.MoveOperation.Down, QTextCursor.MoveMode.MoveAnchor, num_lines - 200)
                cursor.movePosition(QTextCursor.MoveOperation.Start, QTextCursor.MoveMode.KeepAnchor)
                cursor.removeSelectedText()
            # Move cursor to end and scroll to show newest message
            cursor = self.activity_log_text.textCursor()
            cursor.movePosition(QTextCursor.MoveOperation.End)
            self.activity_log_text.setTextCursor(cursor)
        self.activity_log_text.ensureCursorVisible()
        self.activity_log_text.setReadOnly(True)

    def _reset_reference_tab(self):
        """Reset Reference tab data and entry fields to defaults."""
        # Clear line data
        self.survey_lines_data = []
        self.cross_line_data = []
        self.central_point_coords = (None, None)

        # Reset entry fields to defaults
        if hasattr(self, 'central_lat_entry'):
            self.central_lat_entry.clear()
        if hasattr(self, 'central_lon_entry'):
            self.central_lon_entry.clear()
        if hasattr(self, 'line_length_entry'):
            self.line_length_entry.clear()
        if hasattr(self, 'heading_entry'):
            self.heading_entry.setText("0")
        if hasattr(self, 'dist_between_lines_entry'):
            self.dist_between_lines_entry.clear()
        if hasattr(self, 'num_lines_entry'):
            self.num_lines_entry.setText("5")
        if hasattr(self, 'bisect_lead_entry'):
            self.bisect_lead_entry.setText("100")
        if hasattr(self, 'survey_speed_entry'):
            self.survey_speed_entry.setText("8")
        if hasattr(self, 'crossline_passes_entry'):
            self.crossline_passes_entry.setText("2")
        if hasattr(self, 'export_name_entry'):
            self.export_name_entry.clear()
            try:
                dist = int(float(self.dist_between_lines_entry.text())) if self.dist_between_lines_entry.text() else 0
                heading = int(float(self.heading_entry.text())) if self.heading_entry.text() else 0
                export_name = f"Reference_{dist}m_{heading}deg"
            except Exception:
                export_name = "Reference_0m_0deg"
            self.export_name_entry.setText(export_name)
        if hasattr(self, 'offset_direction_var'):
            self.offset_direction_var = "North"
            if hasattr(self, 'offset_direction_combo'):
                self.offset_direction_combo.setCurrentIndex(0)
        if hasattr(self, 'line_length_multiplier'):
            if hasattr(self.line_length_multiplier, 'set'):
                self.line_length_multiplier.set(8.0)
            else:
                self.line_length_multiplier = 8.0
            if hasattr(self, 'multiplier_slider_len'):
                self.multiplier_slider_len.setValue(int(8.0 * 10))
            if hasattr(self, 'multiplier_label_len'):
                self.multiplier_label_len.setText("8.0")
        if hasattr(self, 'dist_between_lines_multiplier'):
            if hasattr(self.dist_between_lines_multiplier, 'set'):
                self.dist_between_lines_multiplier.set(1.0)
            else:
                self.dist_between_lines_multiplier = 1.0
            if hasattr(self, 'multiplier_slider_dist'):
                self.multiplier_slider_dist.setValue(int(1.0 * 10))
            if hasattr(self, 'multiplier_label_dist'):
                self.multiplier_label_dist.setText("1.0")

        # Redraw plot (without clearing GeoTIFF)
        self._plot_survey_plan(preserve_view_limits=True)

    def _calculate_total_survey_time(self):
        """Calculate total survey time including travel between lines and crossline."""
        try:
            geod = pyproj.Geod(ellps="WGS84")
            speed_knots = float(self.survey_speed_entry.text()) if self.survey_speed_entry.text() else 8.0
            speed_m_per_h = speed_knots * 1852

            # Initialize distance calculations
            main_lines_total_distance = 0
            travel_between_lines_total_distance = 0
            travel_to_crossline_distance = 0
            crossline_single_pass_distance = 0

            # Calculate main lines survey time and distance (all main lines)
            main_lines_survey_time_minutes = 0
            for line in self.survey_lines_data:
                lat1, lon1 = line[0]
                lat2, lon2 = line[1]
                _, _, line_length = geod.inv(lon1, lat1, lon2, lat2)
                main_lines_total_distance += line_length
                line_time_hours = line_length / speed_m_per_h if speed_m_per_h > 0 else 0
                main_lines_survey_time_minutes += line_time_hours * 60

            # Calculate travel time and distance between main lines (zigzag pattern)
            travel_between_lines_minutes = 0
            if len(self.survey_lines_data) > 1:
                for i in range(1, len(self.survey_lines_data)):
                    if (i-1) % 2 == 0:
                        end_lat, end_lon = self.survey_lines_data[i-1][1]
                    else:
                        end_lat, end_lon = self.survey_lines_data[i-1][0]

                    if i % 2 == 0:
                        start_lat, start_lon = self.survey_lines_data[i][0]
                    else:
                        start_lat, start_lon = self.survey_lines_data[i][1]

                    _, _, travel_distance = geod.inv(end_lon, end_lat, start_lon, start_lat)
                    travel_between_lines_total_distance += travel_distance
                    travel_time_hours = travel_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                    travel_between_lines_minutes += travel_time_hours * 60

            # Calculate travel time and distance from last main line to crossline
            travel_to_crossline_minutes = 0
            if self.cross_line_data and self.survey_lines_data:
                last_line_index = len(self.survey_lines_data) - 1
                if last_line_index % 2 == 0:
                    last_end_lat, last_end_lon = self.survey_lines_data[last_line_index][1]
                else:
                    last_end_lat, last_end_lon = self.survey_lines_data[last_line_index][0]

                crossline_start_lat, crossline_start_lon = self.cross_line_data[0]
                crossline_end_lat, crossline_end_lon = self.cross_line_data[1]

                _, _, dist_to_crossline_start = geod.inv(last_end_lon, last_end_lat, crossline_start_lon, crossline_start_lat)
                _, _, dist_to_crossline_end = geod.inv(last_end_lon, last_end_lat, crossline_end_lon, crossline_end_lat)

                if dist_to_crossline_start <= dist_to_crossline_end:
                    travel_distance = dist_to_crossline_start
                else:
                    travel_distance = dist_to_crossline_end

                travel_to_crossline_distance = travel_distance
                travel_time_hours = travel_distance / speed_m_per_h if speed_m_per_h > 0 else 0
                travel_to_crossline_minutes = travel_time_hours * 60

            # Calculate crossline survey time and distance
            crossline_survey_minutes = 0
            num_passes = 1
            if self.cross_line_data:
                lat1, lon1 = self.cross_line_data[0]
                lat2, lon2 = self.cross_line_data[1]
                _, _, crossline_length = geod.inv(lon1, lat1, lon2, lat2)
                crossline_single_pass_distance = crossline_length
                crossline_time_hours = crossline_length / speed_m_per_h if speed_m_per_h > 0 else 0
                crossline_survey_minutes = crossline_time_hours * 60

                try:
                    num_passes = int(self.crossline_passes_entry.text()) if self.crossline_passes_entry.text() else 2
                    crossline_survey_minutes *= num_passes
                except (ValueError, AttributeError):
                    num_passes = 2
                    crossline_survey_minutes *= num_passes

            # Calculate total time and distance
            total_minutes = (main_lines_survey_time_minutes +
                           travel_between_lines_minutes +
                           travel_to_crossline_minutes +
                           crossline_survey_minutes)
            total_hours = total_minutes / 60

            total_distance = (main_lines_total_distance +
                            travel_between_lines_total_distance +
                            travel_to_crossline_distance +
                            (crossline_single_pass_distance * num_passes))

            return {
                'total_minutes': total_minutes,
                'total_hours': total_hours,
                'main_lines_minutes': main_lines_survey_time_minutes,
                'crossline_minutes': crossline_survey_minutes,
                'travel_minutes': travel_between_lines_minutes,
                'travel_to_crossline_minutes': travel_to_crossline_minutes,
                'total_distance_m': total_distance,
                'total_distance_km': total_distance / 1000.0,
                'total_distance_nm': total_distance / 1852.0,
                'main_lines_distance_m': main_lines_total_distance,
                'main_lines_distance_km': main_lines_total_distance / 1000.0,
                'main_lines_distance_nm': main_lines_total_distance / 1852.0,
                'travel_between_lines_distance_m': travel_between_lines_total_distance,
                'travel_between_lines_distance_km': travel_between_lines_total_distance / 1000.0,
                'travel_between_lines_distance_nm': travel_between_lines_total_distance / 1852.0,
                'travel_to_crossline_distance_m': travel_to_crossline_distance,
                'travel_to_crossline_distance_km': travel_to_crossline_distance / 1000.0,
                'travel_to_crossline_distance_nm': travel_to_crossline_distance / 1852.0,
                'crossline_single_pass_distance_m': crossline_single_pass_distance,
                'crossline_single_pass_distance_km': crossline_single_pass_distance / 1000.0,
                'crossline_single_pass_distance_nm': crossline_single_pass_distance / 1852.0,
                'crossline_total_distance_m': crossline_single_pass_distance * num_passes,
                'crossline_total_distance_km': (crossline_single_pass_distance * num_passes) / 1000.0,
                'crossline_total_distance_nm': (crossline_single_pass_distance * num_passes) / 1852.0,
                'num_crossline_passes': num_passes
            }
        except Exception as e:
            print(f"Error calculating total survey time: {e}")
            return {
                'total_minutes': 0,
                'total_hours': 0,
                'main_lines_minutes': 0,
                'crossline_minutes': 0,
                'travel_minutes': 0,
                'travel_to_crossline_minutes': 0,
                'total_distance_m': 0,
                'total_distance_km': 0,
                'total_distance_nm': 0,
                'main_lines_distance_m': 0,
                'main_lines_distance_km': 0,
                'main_lines_distance_nm': 0,
                'travel_between_lines_distance_m': 0,
                'travel_between_lines_distance_km': 0,
                'travel_between_lines_distance_nm': 0,
                'travel_to_crossline_distance_m': 0,
                'travel_to_crossline_distance_km': 0,
                'travel_to_crossline_distance_nm': 0,
                'crossline_single_pass_distance_m': 0,
                'crossline_single_pass_distance_km': 0,
                'crossline_single_pass_distance_nm': 0,
                'crossline_total_distance_m': 0,
                'crossline_total_distance_km': 0,
                'crossline_total_distance_nm': 0,
                'num_crossline_passes': 1
            }

    def _calculate_reference_survey_statistics(self):
        """Calculate comprehensive reference planning survey statistics using the existing calculation function."""
        try:
            stats = self._calculate_total_survey_time()
            if not stats:
                return None

            stats['num_main_lines'] = len(self.survey_lines_data)

            return stats

        except Exception as e:
            print(f"Error calculating reference survey statistics: {e}")
            return None

    def _show_reference_planning_info(self):
        """Display comprehensive reference planning survey statistics in a custom dialog with copy functionality."""
        stats = self._calculate_reference_survey_statistics()
        if not stats:
            self._show_message("warning","Statistics Error", "Unable to calculate reference planning statistics. Please ensure survey lines are generated.")
            return

        # Format the statistics for display
        stats_text = "COMPREHENSIVE REFERENCE PLANNING STATISTICS\n"
        stats_text += "=" * 50 + "\n\n"

        # Survey parameters
        try:
            speed_knots = float(self.survey_speed_entry.text()) if self.survey_speed_entry.text() else 8.0
            stats_text += f"Survey Speed: {speed_knots} knots\n\n"
        except Exception:
            stats_text += "Survey Speed: 8.0 knots (default)\n\n"

        # Main survey lines
        stats_text += "MAIN SURVEY LINES\n"
        stats_text += "-" * 20 + "\n"
        stats_text += f"Number of Main Lines: {stats['num_main_lines']}\n"

        # Add heading information for the first main line
        if self.survey_lines_data and len(self.survey_lines_data) > 0:
            try:
                if pyproj is not None:
                    geod = pyproj.Geod(ellps="WGS84")
                    first_line = self.survey_lines_data[0]
                    lat1, lon1 = first_line[0]
                    lat2, lon2 = first_line[1]
                    fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)

                    heading = fwd_az % 360
                    reciprocal_heading = back_az % 360

                    stats_text += f"Heading: {heading:.1f}°\n"
                    stats_text += f"Reciprocal Heading: {reciprocal_heading:.1f}°\n"
                else:
                    stats_text += f"Heading: pyproj not available\n"
            except Exception:
                stats_text += f"Heading: Unable to calculate\n"

        if stats['num_main_lines'] > 0:
            single_line_length_m = stats['main_lines_distance_m'] / stats['num_main_lines']
            single_line_length_km = single_line_length_m / 1000
            single_line_length_nm = single_line_length_m / 1852
            stats_text += f"Length of Main Line: {single_line_length_m:.1f} m ({single_line_length_km:.3f} km, {single_line_length_nm:.3f} nm)\n"
        stats_text += f"Total Distance: {stats['main_lines_distance_m']:.1f} m ({stats['main_lines_distance_km']:.3f} km, {stats['main_lines_distance_nm']:.3f} nm)\n"
        stats_text += f"Survey Time: {stats['main_lines_minutes']:.1f} min\n\n"

        # Travel between lines
        if stats['travel_between_lines_distance_m'] > 0:
            stats_text += "TRAVEL BETWEEN LINES\n"
            stats_text += "-" * 25 + "\n"
            stats_text += f"Total Travel Distance: {stats['travel_between_lines_distance_m']:.1f} m ({stats['travel_between_lines_distance_km']:.3f} km, {stats['travel_between_lines_distance_nm']:.3f} nm)\n"
            stats_text += f"Travel Time: {stats['travel_minutes']:.1f} min\n\n"

        # Crossline information
        if stats['crossline_single_pass_distance_m'] > 0:
            stats_text += "CROSSLINE\n"
            stats_text += "-" * 10 + "\n"
            stats_text += f"Crossline Passes: {stats['num_crossline_passes']}\n"

            if self.cross_line_data and len(self.cross_line_data) >= 2:
                try:
                    if pyproj is not None:
                        geod = pyproj.Geod(ellps="WGS84")
                        lat1, lon1 = self.cross_line_data[0]
                        lat2, lon2 = self.cross_line_data[1]
                        fwd_az, back_az, _ = geod.inv(lon1, lat1, lon2, lat2)

                        crossline_heading = fwd_az % 360
                        crossline_reciprocal_heading = back_az % 360

                        stats_text += f"Crossline Heading: {crossline_heading:.1f}°\n"
                        stats_text += f"Crossline Reciprocal Heading: {crossline_reciprocal_heading:.1f}°\n"
                    else:
                        stats_text += f"Crossline Heading: pyproj not available\n"
                except Exception:
                    stats_text += f"Crossline Heading: Unable to calculate\n"

            stats_text += f"Crossline (per pass): {stats['crossline_single_pass_distance_m']:.1f} m ({stats['crossline_single_pass_distance_km']:.3f} km, {stats['crossline_single_pass_distance_nm']:.3f} nm)\n"
            stats_text += f"Crossline Total Distance: {stats['crossline_total_distance_m']:.1f} m ({stats['crossline_total_distance_km']:.3f} km, {stats['crossline_total_distance_nm']:.3f} nm)\n"
            stats_text += f"Crossline Time: {stats['crossline_minutes']:.1f} min\n"
            if stats['travel_to_crossline_distance_m'] > 0:
                stats_text += f"Travel to Crossline: {stats['travel_to_crossline_distance_m']:.1f} m ({stats['travel_to_crossline_distance_km']:.3f} km, {stats['travel_to_crossline_distance_nm']:.3f} nm)\n"
            stats_text += f"Travel to Crossline Time: {stats['travel_to_crossline_minutes']:.1f} min\n"
            stats_text += "\n"

        # Total summary
        stats_text += "TOTAL SURVEY SUMMARY\n"
        stats_text += "-" * 25 + "\n"
        stats_text += f"Total Distance: {stats['total_distance_m']:.1f} m ({stats['total_distance_km']:.3f} km, {stats['total_distance_nm']:.3f} nm)\n"
        stats_text += f"Total Time: {stats['total_minutes']:.1f} min ({stats['total_hours']:.2f} hr)\n\n"

        # Survey pattern explanation
        stats_text += "SURVEY PATTERN\n"
        stats_text += "-" * 15 + "\n"
        stats_text += "1. Run main survey lines (zigzag pattern)\n"
        if stats['crossline_single_pass_distance_m'] > 0:
            stats_text += "2. Travel from last main line to crossline start\n"
            stats_text += f"3. Run crossline ({stats['num_crossline_passes']} passes)\n"
        stats_text += "\n"

        # Waypoints with labels
        stats_text += "WAYPOINTS\n"
        stats_text += "-" * 10 + "\n"
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
                stats_text += f"{start_label}: {start_lat_ddm}, {start_lon_ddm} ({start[0]:.6f}, {start[1]:.6f})\n"
                stats_text += f"{end_label}: {end_lat_ddm}, {end_lon_ddm} ({end[0]:.6f}, {end[1]:.6f})\n"
        if self.cross_line_data:
            cls_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][0], is_latitude=True)
            cls_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[0][1], is_latitude=False)
            cle_lat_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][0], is_latitude=True)
            cle_lon_ddm = decimal_degrees_to_ddm(self.cross_line_data[1][1], is_latitude=False)
            stats_text += f"CLS: {cls_lat_ddm}, {cls_lon_ddm} ({self.cross_line_data[0][0]:.6f}, {self.cross_line_data[0][1]:.6f})\n"
            stats_text += f"CLE: {cle_lat_ddm}, {cle_lon_ddm} ({self.cross_line_data[1][0]:.6f}, {self.cross_line_data[1][1]:.6f})\n"

        # Create custom dialog window with copy functionality
        self._show_statistics_dialog("Reference Planning Statistics", stats_text)

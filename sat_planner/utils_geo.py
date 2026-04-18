"""
Pure geometry/coordinate utilities for SAT Planner.
No Qt or app state; safe to unit test.
"""


def _signed_deg_min_sec_to_decimal(degrees, minutes, seconds=0.0):
    """Rebuild decimal degrees from signed degree part and magnitude of min/sec (export_utils convention)."""
    sign = -1 if float(degrees) < 0 else 1
    dec = abs(float(degrees)) + float(minutes) / 60.0 + float(seconds) / 3600.0
    return sign * dec


def lat_lon_decimal_from_survey_csv_row(row):
    """
    Parse one csv.DictReader row from SAT Planner exports: *_DDD.csv, *_DMM.csv, or *_DMS.csv.

    DDD uses 'Latitude' / 'Longitude' (decimal). DMM/DMS use 'Latitude (Deg)', 'Latitude (Min)',
    [ 'Latitude (Sec)' ], and the same for Longitude — matching export_utils.write_dmm_csv / write_dms_csv.

    Returns:
        (lat, lon) as floats.

    Raises:
        ValueError: if the row does not contain a recognizable coordinate set.
    """
    if not isinstance(row, dict):
        raise ValueError("row must be a dict")
    lat_ddd = row.get("Latitude")
    lon_ddd = row.get("Longitude")
    if lat_ddd not in (None, "") and lon_ddd not in (None, ""):
        try:
            return float(lat_ddd), float(lon_ddd)
        except (TypeError, ValueError) as e:
            raise ValueError("invalid decimal Latitude/Longitude") from e

    try:
        lat_deg = float(row["Latitude (Deg)"])
        lat_min = float(row["Latitude (Min)"])
        lon_deg = float(row["Longitude (Deg)"])
        lon_min = float(row["Longitude (Min)"])
    except (KeyError, TypeError, ValueError) as e:
        raise ValueError("missing or invalid DMM/DMS CSV columns") from e

    lat_sec_raw = row.get("Latitude (Sec)")
    lon_sec_raw = row.get("Longitude (Sec)")
    if lat_sec_raw not in (None, "") and str(lat_sec_raw).strip() != "":
        try:
            lat_sec = float(lat_sec_raw)
            lon_sec = float(lon_sec_raw)
        except (TypeError, ValueError) as e:
            raise ValueError("invalid DMS seconds") from e
        return (
            _signed_deg_min_sec_to_decimal(lat_deg, lat_min, lat_sec),
            _signed_deg_min_sec_to_decimal(lon_deg, lon_min, lon_sec),
        )

    return (
        _signed_deg_min_sec_to_decimal(lat_deg, lat_min, 0.0),
        _signed_deg_min_sec_to_decimal(lon_deg, lon_min, 0.0),
    )


def decimal_degrees_to_ddm(decimal_deg, is_latitude=True, decimal_places=3):
    """Convert decimal degrees to degrees and decimal minutes format.

    Args:
        decimal_deg: Decimal degrees value
        is_latitude: True for latitude, False for longitude
        decimal_places: Number of decimal places for minutes (default 3)

    Returns:
        Formatted string like "DD°MM.mmm'N" or "DDD°MM.mmm'E" (minutes as MM.mmm)
    """
    try:
        if is_latitude:
            direction = 'N' if decimal_deg >= 0 else 'S'
            degrees_width = 2
        else:
            direction = 'E' if decimal_deg >= 0 else 'W'
            degrees_width = 3

        abs_deg = abs(decimal_deg)
        degrees = int(abs_deg)
        decimal_minutes = (abs_deg - degrees) * 60.0
        # Format minutes as MM.mmm (2-digit minutes, then decimal places)
        minutes_int = int(decimal_minutes)
        frac = decimal_minutes - minutes_int
        frac_units = round(frac * (10 ** decimal_places))
        if frac_units >= 10 ** decimal_places:
            frac_units = 0
            minutes_int += 1
        if minutes_int >= 60:
            minutes_int = 59
            frac_units = (10 ** decimal_places) - 1
        minutes_str = f"{minutes_int:02d}.{frac_units:0{decimal_places}d}"
        return f"{degrees:0{degrees_width}d}°{minutes_str}'{direction}"
    except Exception:
        return f"{decimal_deg:.6f}"


def decimal_degrees_to_dms_string(decimal_deg, is_latitude=True, decimals_seconds=1):
    """Convert decimal degrees to degrees-minutes-seconds (DMS) format.

    Args:
        decimal_deg: Decimal degrees value
        is_latitude: True for latitude (N/S), False for longitude (E/W)
        decimals_seconds: Number of decimal places for seconds (default 1)

    Returns:
        Formatted string like "42°30'45.2\"N" or "071°15'30.0\"W"
    """
    try:
        if is_latitude:
            direction = 'N' if decimal_deg >= 0 else 'S'
            degrees_width = 2
        else:
            direction = 'E' if decimal_deg >= 0 else 'W'
            degrees_width = 3

        abs_deg = abs(decimal_deg)
        degrees = int(abs_deg)
        minutes_float = (abs_deg - degrees) * 60.0
        minutes = int(minutes_float)
        seconds = (minutes_float - minutes) * 60.0

        if decimals_seconds <= 0:
            seconds_str = f"{int(round(seconds))}"
        else:
            seconds_str = f"{seconds:.{decimals_seconds}f}"
        return f"{degrees:0{degrees_width}d}°{minutes:02d}'{seconds_str}\"{direction}"
    except Exception:
        return f"{decimal_deg:.6f}"

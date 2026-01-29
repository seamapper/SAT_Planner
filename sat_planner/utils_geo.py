"""
Pure geometry/coordinate utilities for SAT Planner.
No Qt or app state; safe to unit test.
"""


def decimal_degrees_to_ddm(decimal_deg, is_latitude=True):
    """Convert decimal degrees to degrees and decimal minutes format.

    Args:
        decimal_deg: Decimal degrees value
        is_latitude: True for latitude, False for longitude

    Returns:
        Formatted string like "DD°MM.mmm'N" or "DDD°MM.mmm'E"
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
        return f"{degrees:0{degrees_width}d}°{decimal_minutes:06.3f}'{direction}"
    except Exception:
        return f"{decimal_deg:.6f}"

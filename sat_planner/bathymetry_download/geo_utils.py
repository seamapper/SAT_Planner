"""Geographic coordinate and extent helpers."""

import pyproj


def transform_bbox(bbox, from_crs, to_crs):
    """Transform (xmin, ymin, xmax, ymax) between CRS definitions."""
    transformer = pyproj.Transformer.from_crs(from_crs, to_crs, always_xy=True)
    xmin, ymin, xmax, ymax = bbox
    return transformer.transform_bounds(xmin, ymin, xmax, ymax)


def bbox_to_meters(bbox_4326):
    """Convert a GCS selection bbox to Web Mercator meters."""
    return transform_bbox(bbox_4326, "EPSG:4326", "EPSG:3857")


def bbox_to_degrees(bbox_3857):
    """Convert a Web Mercator selection bbox to GCS degrees."""
    return transform_bbox(bbox_3857, "EPSG:3857", "EPSG:4326")


def extents_equal(extent_a, extent_b, tol=1e-5):
    """Return True when two extents match within tolerance."""
    return all(abs(a - b) < tol for a, b in zip(extent_a, extent_b))


def clamp_extent_to_bounds(extent, bounds):
    """Clamp a view extent to valid geographic/service bounds."""
    if bounds is None:
        return extent
    xmin, ymin, xmax, ymax = extent
    bxmin, bymin, bxmax, bymax = bounds
    return (
        max(xmin, bxmin),
        max(ymin, bymin),
        min(xmax, bxmax),
        min(ymax, bymax),
    )


def extent_compare_tolerance(web_mercator_map):
    """Tolerance for comparing extents in the active map CRS."""
    return 1.0 if web_mercator_map else 1e-5


def bboxes_overlap(bbox1, bbox2):
    """Return True when two axis-aligned bboxes overlap."""
    x1min, y1min, x1max, y1max = bbox1
    x2min, y2min, x2max, y2max = bbox2
    return not (x1max <= x2min or x2max <= x1min or y1max <= y2min or y2max <= y1min)

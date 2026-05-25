"""
Shared GMRT GeoTIFF topo/bathy splitter.

Splits a downloaded GMRT GeoTIFF into two single-band float32 GeoTIFFs at the
mean-sea-level boundary:

  * topo  file: values >= 0 (0 lands in topo)
  * bathy file: values <  0

Used by both the standalone Download GMRT Grid dialog and the per-tab
"Download GMRT on import" workflow in SAT Planner.

The helper writes the side files, removes any side file whose values are all
NaN (or whose size is zero), and reports back which paths exist and whether
either side was empty. The caller decides whether to keep or delete the
original combined GeoTIFF after inspecting the result.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np

try:
    import rasterio
except ImportError:
    rasterio = None


@dataclass
class SplitResult:
    """Outcome of split_topo_bathy().

    Attributes:
        topo_path: Path to the topography (>=0) GeoTIFF, or None if it was
            all-NaN / not written.
        bathy_path: Path to the bathymetry (<0) GeoTIFF, or None if it was
            all-NaN / not written.
        topo_all_nan: True if the >=0 partition contained no finite values.
        bathy_all_nan: True if the <0 partition contained no finite values.
        error: Error message if the split could not be performed, else None.
    """

    topo_path: Optional[str] = None
    bathy_path: Optional[str] = None
    topo_all_nan: bool = False
    bathy_all_nan: bool = False
    error: Optional[str] = None


def bathy_path_after_split(src_path: str) -> str:
    """Return the bathy filename that split_topo_bathy() would produce for src_path."""
    base, ext = os.path.splitext(src_path)
    output_ext = ext if ext.lower() in (".tif", ".tiff") else ".tif"
    return base + "_bathy" + output_ext


def topo_path_after_split(src_path: str) -> str:
    """Return the topo filename that split_topo_bathy() would produce for src_path."""
    base, ext = os.path.splitext(src_path)
    output_ext = ext if ext.lower() in (".tif", ".tiff") else ".tif"
    return base + "_topo" + output_ext


def split_topo_bathy(src_path: str, *, log: Optional[Callable[[str], None]] = None) -> SplitResult:
    """Split src_path into <base>_topo and <base>_bathy GeoTIFFs.

    The original combined GeoTIFF is left in place; the caller decides whether
    to delete it (typical usage: delete on success, keep on failure).

    Side-file behavior:
      * A side file whose values are entirely NaN is removed from disk and the
        corresponding result path is set to None.

    Args:
        src_path: Path to the combined topo+bathy GeoTIFF.
        log: Optional logger callable invoked with informational strings.

    Returns:
        SplitResult describing the produced paths and per-side empty flags.
    """

    def _log(msg: str) -> None:
        if log is not None:
            try:
                log(msg)
            except Exception:
                pass

    if rasterio is None:
        return SplitResult(error="rasterio is required to split GMRT grids; not installed.")

    if not os.path.exists(src_path):
        return SplitResult(error=f"Source GeoTIFF not found: {src_path}")

    topo_file = topo_path_after_split(src_path)
    bathy_file = bathy_path_after_split(src_path)

    try:
        with rasterio.open(src_path) as src:
            data = src.read(1)
            profile = src.profile

            # 0 lands in topo (>= 0); strictly negative depths in bathy (< 0).
            topo_data = np.where(data >= 0, data, np.nan)
            bathy_data = np.where(data < 0, data, np.nan)

            profile.update(dtype=rasterio.float32, nodata=np.nan)

            with rasterio.open(topo_file, "w", **profile) as dst:
                dst.write(topo_data.astype(np.float32), 1)
            with rasterio.open(bathy_file, "w", **profile) as dst:
                dst.write(bathy_data.astype(np.float32), 1)
    except Exception as exc:
        return SplitResult(error=f"Error splitting grid: {exc}")

    result = SplitResult(topo_path=topo_file, bathy_path=bathy_file)

    for path, arr, label, is_topo in (
        (topo_file, topo_data, "topo", True),
        (bathy_file, bathy_data, "bathy", False),
    ):
        try:
            all_nan = bool(np.isnan(arr).all()) or os.path.getsize(path) == 0
        except OSError:
            all_nan = True
        if all_nan:
            if is_topo:
                result.topo_all_nan = True
                result.topo_path = None
            else:
                result.bathy_all_nan = True
                result.bathy_path = None
            try:
                os.remove(path)
                _log(f"Deleted empty {label} GeoTIFF: {os.path.basename(path)}")
            except OSError as exc:
                _log(f"Failed to delete empty {label} GeoTIFF: {os.path.basename(path)}: {exc}")

    _log(
        "Split GeoTIFF: "
        + ", ".join(
            os.path.basename(p)
            for p in (result.topo_path, result.bathy_path)
            if p is not None
        )
    )
    return result

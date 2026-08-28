# Runtime hook for PyInstaller one-file exe on Windows.
# Configure DLL search path and GDAL/PROJ data dirs before geospatial imports.
import os
import sys

if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
    root = sys._MEIPASS

    # Bundled DLLs must win over any OSGeo/QGIS/ArcGIS install on PATH.
    os.environ["PATH"] = root + os.pathsep + os.environ.get("PATH", "")

    try:
        os.add_dll_directory(root)
        for subdir in ("rasterio", "fiona"):
            pkg_dir = os.path.join(root, subdir)
            if os.path.isdir(pkg_dir):
                os.add_dll_directory(pkg_dir)
    except (AttributeError, OSError):
        pass

    # Preload GDAL dependency chain so rasterio._base gets matching symbols.
    if sys.platform.startswith("win"):
        import ctypes

        preload_names = (
            "zlib.dll",
            "liblzma.dll",
            "sqlite3.dll",
            "proj_9.dll",
            "geos_c.dll",
            "gdal.dll",
        )
        search_dirs = [root]
        for subdir in ("rasterio", "fiona"):
            pkg_dir = os.path.join(root, subdir)
            if os.path.isdir(pkg_dir):
                search_dirs.append(pkg_dir)
        for name in preload_names:
            for directory in search_dirs:
                dll_path = os.path.join(directory, name)
                if os.path.isfile(dll_path):
                    try:
                        ctypes.WinDLL(dll_path)
                    except OSError:
                        pass
                    break

    # GDAL_DATA: match PyInstaller osgeo hook layout (Library/share/gdal).
    if sys.platform.startswith("win"):
        gdal_candidates = (
            os.path.join(root, "Library", "share", "gdal"),
            os.path.join(root, "Library", "data"),
            os.path.join(root, "gdal_data"),  # legacy bundle layout
        )
    else:
        gdal_candidates = (os.path.join(root, "share", "gdal"),)
    for gdal_data in gdal_candidates:
        if os.path.isdir(gdal_data):
            os.environ["GDAL_DATA"] = gdal_data
            break

    # PROJ_LIB: match PyInstaller pyproj hook layout (Library/share/proj).
    if sys.platform.startswith("win"):
        proj_candidates = (
            os.path.join(root, "Library", "share", "proj"),
            os.path.join(root, "proj_data"),  # legacy bundle layout
        )
    else:
        proj_candidates = (os.path.join(root, "share", "proj"),)
    for proj_data in proj_candidates:
        if os.path.isdir(proj_data):
            os.environ["PROJ_LIB"] = proj_data
            os.environ["PROJ_DATA"] = proj_data
            break

"""
Helpers for PyInstaller: bundle GDAL/PROJ data and native DLLs used by rasterio/fiona/pyproj.

Conda/mamba installs place GDAL and PROJ in Library/share and their DLLs in Library/bin.
PyInstaller's collect_data_files() only gathers files inside each Python package, so the
share/proj and share/gdal trees (and many GDAL DLLs) are missed unless added explicitly.
"""
from __future__ import annotations

import os
import sys

# Windows Universal CRT forwarders must come from the OS, not from conda.
_SKIP_DLL_PREFIXES = ("api-ms-win-",)
_SKIP_DLL_NAMES = frozenset(
    {
        "ucrtbase.dll",
        "kernel32.dll",
        "advapi32.dll",
        "user32.dll",
        "ole32.dll",
        "oleaut32.dll",
        "ws2_32.dll",
        "crypt32.dll",
        "odbc32.dll",
    }
)


def conda_library_dirs():
    """Return (bin_dir, share_dir) for a conda/mamba environment, or (None, None)."""
    candidates = [
        os.environ.get("CONDA_PREFIX", ""),
        sys.prefix,
        getattr(sys, "base_prefix", sys.prefix),
    ]
    seen = set()
    for prefix in candidates:
        if not prefix or prefix in seen:
            continue
        seen.add(prefix)
        bin_dir = os.path.join(prefix, "Library", "bin")
        share_dir = os.path.join(prefix, "Library", "share")
        if os.path.isdir(bin_dir) and os.path.isdir(share_dir):
            return bin_dir, share_dir
    return None, None


def collect_gdal_proj_data(share_dir):
    """Return PyInstaller datas entries using conda layout expected by GDAL/PROJ runtime hooks."""
    datas = []
    gdal_data = os.path.join(share_dir, "gdal")
    proj_data = os.path.join(share_dir, "proj")
    if os.path.isdir(gdal_data):
        datas.append((gdal_data, os.path.join("Library", "share", "gdal")))
    if os.path.isdir(proj_data):
        datas.append((proj_data, os.path.join("Library", "share", "proj")))
    return datas


def _should_bundle_dll(dll_name):
    lower = dll_name.lower()
    if lower in _SKIP_DLL_NAMES:
        return False
    return not any(lower.startswith(prefix) for prefix in _SKIP_DLL_PREFIXES)


# DLLs that rasterio/fiona .pyd extensions look for beside the importing module.
_PACKAGE_DLL_NAMES = frozenset(
    {
        "gdal.dll",
        "geos.dll",
        "geos_c.dll",
        "proj_9.dll",
    }
)
_PACKAGE_DESTS = ("rasterio", "fiona")


def conda_macos_prefix():
    """Return the conda prefix when this Mac env has GDAL/PROJ share data."""
    if not sys.platform.startswith("darwin"):
        return None
    # sys.prefix is the interpreter running this build. Prefer it over
    # CONDA_PREFIX so a base-env variable does not hide the env that has rasterio.
    py_tag = f"python{sys.version_info.major}.{sys.version_info.minor}"
    candidates = [
        sys.prefix,
        os.environ.get("CONDA_PREFIX", ""),
        getattr(sys, "base_prefix", ""),
    ]
    seen = set()
    for prefix in candidates:
        if not prefix or prefix in seen:
            continue
        seen.add(prefix)
        gdal_share = os.path.join(prefix, "share", "gdal")
        rasterio_dir = os.path.join(prefix, "lib", py_tag, "site-packages", "rasterio")
        if os.path.isdir(gdal_share) and os.path.isdir(rasterio_dir):
            return prefix
    return None


def collect_macos_gdal_data(prefix):
    """Return PyInstaller datas for the Mac runtime hook (share/gdal, share/proj)."""
    datas = []
    for name in ("gdal", "proj"):
        src = os.path.join(prefix, "share", name)
        if os.path.isdir(src):
            datas.append((src, os.path.join("share", name)))
    plugins = os.path.join(prefix, "lib", "gdalplugins")
    if os.path.isdir(plugins):
        datas.append((plugins, "gdalplugins"))
    return datas


def collect_macos_dylibs(prefix):
    """Bundle conda dylibs linked from rasterio, fiona, shapely, and pyproj.

    On macOS those extensions load libgdal/libproj/libgeos by @rpath from the
    conda env, and libgdal pulls in more conda libraries by name. PyInstaller's
    package scan does not see that tree.
    """
    import glob

    from PyInstaller.depend.bindepend import get_imports

    lib_dir = os.path.join(prefix, "lib")
    py_tag = f"python{sys.version_info.major}.{sys.version_info.minor}"
    site = os.path.join(lib_dir, py_tag, "site-packages")
    seeds = []
    for pattern in ("rasterio/*.so", "fiona/*.so", "shapely/*.so", "pyproj/*.so"):
        for so_path in glob.glob(os.path.join(site, pattern)):
            try:
                deps = get_imports(so_path)
            except Exception:
                continue
            for dep, _ in deps:
                resolved = _resolve_conda_dylib(dep, lib_dir, prefix)
                if resolved:
                    seeds.append(resolved)

    seen = set()
    binaries = []
    queue = list(seeds)
    while queue:
        path = queue.pop()
        key = os.path.normpath(path)
        if key in seen:
            continue
        if not os.path.exists(key):
            continue
        seen.add(key)
        binaries.append((key, "."))
        real = os.path.realpath(key)
        try:
            deps = get_imports(real if os.path.isfile(real) else key)
        except Exception:
            continue
        for dep, _ in deps:
            resolved = _resolve_conda_dylib(dep, lib_dir, prefix)
            if resolved and os.path.normpath(resolved) not in seen:
                queue.append(resolved)
    return binaries


def _resolve_conda_dylib(dep, lib_dir, prefix):
    """Return a conda-prefix dylib path, or None for system libraries."""
    if not dep:
        return None
    name = os.path.basename(dep)
    if name in {"libSystem.B.dylib", "libc++.1.dylib", "libc++abi.dylib"}:
        return None
    if dep.startswith(("/usr/lib", "/System")):
        return None
    if os.path.isabs(dep) and os.path.exists(dep):
        if dep.startswith(prefix + os.sep):
            return dep
        return None
    candidate = os.path.join(lib_dir, name)
    if os.path.exists(candidate):
        return candidate
    return None


def collect_gdal_dlls(bin_dir):
    """Return PyInstaller binaries entries for GDAL/PROJ/GEOS DLLs and dependencies."""
    from PyInstaller.depend.bindepend import get_imports

    seeds = (
        "gdal.dll",
        "geos_c.dll",
        "proj_9.dll",
        "spatialite.dll",
    )
    seen = set()
    queue = []
    for name in seeds:
        path = os.path.join(bin_dir, name)
        if os.path.isfile(path):
            queue.append(os.path.normpath(path))

    binaries = []
    collected_paths = []
    while queue:
        dll_path = queue.pop()
        dll_name = os.path.basename(dll_path)
        key = dll_path.lower()
        if key in seen:
            continue
        if not _should_bundle_dll(dll_name):
            continue
        seen.add(key)
        collected_paths.append(dll_path)
        try:
            for dep, _ in get_imports(dll_path):
                dep_name = os.path.basename(dep)
                if not dep_name.lower().endswith(".dll"):
                    continue
                if not _should_bundle_dll(dep_name):
                    continue
                dep_path = dep if os.path.isabs(dep) else os.path.join(bin_dir, dep_name)
                dep_path = os.path.normpath(dep_path)
                if os.path.isfile(dep_path) and dep_path.lower() not in seen:
                    queue.append(dep_path)
        except Exception:
            pass

    for dll_path in collected_paths:
        binaries.append((dll_path, "."))
        if os.path.basename(dll_path).lower() in _PACKAGE_DLL_NAMES:
            for dest in _PACKAGE_DESTS:
                binaries.append((dll_path, dest))
    return binaries

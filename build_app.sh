#!/bin/bash
# Run from this script's folder so SAT_Planner_macOS.spec and paths resolve from any cwd.
# Build a macOS app for SAT Planner.
# Output: dist/SAT_Planner_v{version}.app (version from sat_planner/constants.py).
# Icon: media/CCOM.icns, generated here from media/CCOM.png.

cd "$(dirname "$0")"

PYTHON_PATH="${PYTHON_PATH:-/Users/pjohnson/miniforge3/envs/python313/bin/python}"

echo "Building SAT Planner macOS app (name and version from code)..."
echo

if [[ ! -x "$PYTHON_PATH" ]]; then
    echo "Python not found: $PYTHON_PATH"
    echo "Set PYTHON_PATH to the interpreter that has this app's dependencies."
    exit 1
fi

if ! "$PYTHON_PATH" -c "import PyQt6.QtWebEngineWidgets" >/dev/null 2>&1; then
    echo "PyQt6 WebEngine is not installed for $PYTHON_PATH"
    echo "Install it, then run this script again:"
    echo "  $PYTHON_PATH -m pip install PyQt6-WebEngine"
    exit 1
fi

if ! "$PYTHON_PATH" -m pip show pyinstaller >/dev/null 2>&1; then
    echo "PyInstaller not found. Installing..."
    "$PYTHON_PATH" -m pip install pyinstaller
fi

if [[ ! -f media/CCOM.png ]]; then
    echo "media/CCOM.png is missing, so the Dock icon cannot be built."
    exit 1
fi

echo "Building media/CCOM.icns from media/CCOM.png..."
ICONSET="$(mktemp -d)/CCOM.iconset"
mkdir -p "$ICONSET"
sips -z 16 16     media/CCOM.png --out "$ICONSET/icon_16x16.png" >/dev/null
sips -z 32 32     media/CCOM.png --out "$ICONSET/icon_16x16@2x.png" >/dev/null
sips -z 32 32     media/CCOM.png --out "$ICONSET/icon_32x32.png" >/dev/null
sips -z 64 64     media/CCOM.png --out "$ICONSET/icon_32x32@2x.png" >/dev/null
sips -z 128 128   media/CCOM.png --out "$ICONSET/icon_128x128.png" >/dev/null
sips -z 256 256   media/CCOM.png --out "$ICONSET/icon_128x128@2x.png" >/dev/null
sips -z 256 256   media/CCOM.png --out "$ICONSET/icon_256x256.png" >/dev/null
sips -z 512 512   media/CCOM.png --out "$ICONSET/icon_256x256@2x.png" >/dev/null
sips -z 512 512   media/CCOM.png --out "$ICONSET/icon_512x512.png" >/dev/null
sips -z 1024 1024 media/CCOM.png --out "$ICONSET/icon_512x512@2x.png" >/dev/null
iconutil -c icns "$ICONSET" -o media/CCOM.icns
rm -rf "$(dirname "$ICONSET")"

"$PYTHON_PATH" -m PyInstaller --noconfirm SAT_Planner_macOS.spec
if [[ $? -ne 0 ]]; then
    echo
    echo "Build failed!"
    exit 1
fi

VERSION="$("$PYTHON_PATH" -c "import re; print(re.search(r'__version__\s*=\s*[\"\\']([^\"\\']+)', open('sat_planner/constants.py', encoding='utf-8').read()).group(1))")"
APP="dist/SAT_Planner_v${VERSION}.app"
if [[ ! -d "$APP" ]]; then
    echo "Build finished, but $APP was not created."
    exit 1
fi

# Ad-hoc sign the bundle, including Qt WebEngine's helper, so the map can start.
codesign --force --deep --sign - --entitlements macos/entitlements.plist "$APP"

echo
echo "Build complete! App is ${APP}"

@echo off
REM Run from this script's folder so ShapefileLinesToGpx.spec resolves from any cwd
cd /d "%~dp0"

REM Build one-file GUI exe: ShapefileLinesToGpx_v{version}.exe (version from __version__ in ShapefileLinesToGpx_PyQt.py)

set PYTHON_PATH=C:\Users\pjohnson\Documents\Python\.venv\Scripts\python.exe

echo Building Shapefile lines to GPX executable...
echo.

"%PYTHON_PATH%" -m pip show pyinstaller >nul 2>&1
if errorlevel 1 (
    echo PyInstaller not found. Installing...
    "%PYTHON_PATH%" -m pip install pyinstaller
)

"%PYTHON_PATH%" -m PyInstaller ShapefileLinesToGpx.spec

if errorlevel 1 (
    echo.
    echo Build failed!
    pause
    exit /b 1
)

echo.
echo Build complete! Executable is in the dist folder.
pause

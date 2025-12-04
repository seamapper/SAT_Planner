@echo off
REM Build executable for SAT Planner
REM This script builds the executable using PyInstaller

set PYTHON_PATH=C:\Users\pjohnson\PycharmProjects\MultibeamToolsMolokai\.venv\Scripts\python.exe

echo Building Sat_Planner_v2025.08.exe...
echo.

REM Check if PyInstaller is installed
"%PYTHON_PATH%" -m pip show pyinstaller >nul 2>&1
if errorlevel 1 (
    echo PyInstaller not found. Installing...
    "%PYTHON_PATH%" -m pip install pyinstaller
)

REM Build the executable using the spec file
"%PYTHON_PATH%" -m PyInstaller Sat_Planner_v2025.08.spec

if errorlevel 1 (
    echo.
    echo Build failed!
    pause
    exit /b 1
)

echo.
echo Build complete! Executable is in the dist folder.
pause


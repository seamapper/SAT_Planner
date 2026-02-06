@echo off
REM Build executable for SAT Planner
REM Output: SAT_Planner_v{version}.exe (version from sat_planner.constants). Icon: media\CCOM.ico

set PYTHON_PATH=C:\Users\pjohnson\PycharmProjects\MultibeamToolsMolokai\.venv\Scripts\python.exe

echo Building SAT Planner executable (name and version from code)...
echo.

REM Check if PyInstaller is installed
"%PYTHON_PATH%" -m pip show pyinstaller >nul 2>&1
if errorlevel 1 (
    echo PyInstaller not found. Installing...
    "%PYTHON_PATH%" -m pip install pyinstaller
)

REM Build using SAT_Planner.spec (exe name = SAT_Planner_v + __version__, icon = media\CCOM.ico)
"%PYTHON_PATH%" -m PyInstaller SAT_Planner.spec

if errorlevel 1 (
    echo.
    echo Build failed!
    pause
    exit /b 1
)

echo.
echo Build complete! Executable is in the dist folder.
pause


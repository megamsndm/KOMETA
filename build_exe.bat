@echo off
REM =============================================================================
REM KOMETA — build_exe.bat
REM Builds KOMETA.exe using PyInstaller
REM Compatible: Python 3.10–3.14, PyInstaller 6.x, Windows 10/11
REM =============================================================================

echo.
echo  ============================================
echo   KOMETA — PyInstaller Build
echo  ============================================
echo.

REM ── Check Python ─────────────────────────────────────────────────────────
python --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Python not found. Install Python 3.10+ and add to PATH.
    pause & exit /b 1
)
for /f "tokens=2 delims= " %%v in ('python --version 2^>^&1') do set PYVER=%%v
echo [INFO] Python %PYVER% detected.

REM ── Install / upgrade PyInstaller ────────────────────────────────────────
python -c "import PyInstaller" >nul 2>&1
if errorlevel 1 (
    echo [INFO] PyInstaller not found — installing...
    pip install pyinstaller
)

REM ── Install required packages if missing ─────────────────────────────────
echo [INFO] Checking required packages...
pip install --quiet pyqt5 matplotlib scipy numpy

REM ── Clean previous build artifacts ───────────────────────────────────────
echo [INFO] Cleaning previous build...
if exist build   rmdir /s /q build
if exist dist\KOMETA.exe del /f /q dist\KOMETA.exe

REM ── Run PyInstaller ──────────────────────────────────────────────────────
echo [INFO] Building KOMETA.exe (this takes 1-3 minutes)...
echo.
pyinstaller KOMETA.spec

if errorlevel 1 (
    echo.
    echo [ERROR] Build failed. Check the output above for details.
    echo.
    echo  Common fixes:
    echo  1. Upgrade PyInstaller:  pip install --upgrade pyinstaller
    echo  2. Install all deps:     pip install pyqt5 matplotlib scipy numpy h5py
    echo  3. Disable antivirus temporarily during build
    echo  4. Delete build\ and dist\ folders, then retry
    echo.
    pause & exit /b 1
)

echo.
echo  ============================================
echo   SUCCESS:  dist\KOMETA.exe  is ready
echo  ============================================
echo.
echo  You can copy dist\KOMETA.exe to any Windows PC and run it directly.
echo  No Python installation is needed on the target machine.
echo.
pause

@echo off
setlocal

cd /d "%~dp0"
echo.
echo === Genomics Platform: Install ===
echo.

where py >nul 2>nul
if %errorlevel% neq 0 (
  echo ERROR: Python launcher "py" not found.
  echo Install Python 3.11+ from https://www.python.org/downloads/ and re-run.
  echo.
  pause
  exit /b 1
)

if not exist ".venv\Scripts\python.exe" (
  echo Creating virtual environment in .venv...
  py -3 -m venv .venv
  if %errorlevel% neq 0 (
    echo ERROR: Failed to create venv.
    echo.
    pause
    exit /b %errorlevel%
  )
) else (
  echo Using existing virtual environment in .venv...
)

call ".venv\Scripts\activate.bat"
if %errorlevel% neq 0 (
  echo ERROR: Failed to activate venv.
  echo.
  pause
  exit /b %errorlevel%
)

echo Upgrading pip tooling...
python -m pip install --upgrade pip setuptools wheel
if %errorlevel% neq 0 (
  echo ERROR: Failed to upgrade pip/setuptools/wheel.
  echo.
  pause
  exit /b %errorlevel%
)

if exist "requirements.txt" (
  echo Installing dependencies from requirements.txt...
  python -m pip install -r requirements.txt
  if %errorlevel% neq 0 (
    echo ERROR: Failed to install requirements.
    echo.
    pause
    exit /b %errorlevel%
  )
) else (
  echo WARNING: requirements.txt not found; skipping dependency install.
)

echo Installing this project (editable)...
python -m pip install -e .
if %errorlevel% neq 0 (
  echo ERROR: Failed to install this project.
  echo.
  pause
  exit /b %errorlevel%
)

echo.
echo Install complete.
echo Run start.bat to launch the GUI.
echo.
pause
endlocal

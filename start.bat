@echo off
setlocal

cd /d "%~dp0"
echo.
echo === Genomics Platform: Start GUI ===
echo.

if not exist ".venv\Scripts\activate.bat" (
  echo ERROR: Virtual environment not found.
  echo Run install.bat first.
  echo.
  pause
  exit /b 1
)

call ".venv\Scripts\activate.bat"
if %errorlevel% neq 0 (
  echo ERROR: Failed to activate venv.
  echo.
  pause
  exit /b %errorlevel%
)

if not exist "app.py" (
  echo ERROR: app.py not found in %cd%
  echo.
  pause
  exit /b 1
)

echo Starting Streamlit...
python -m streamlit run app.py

endlocal

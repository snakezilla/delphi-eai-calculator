@echo off
REM Delphi EAI Calculator - Windows Run Script
REM This script starts the Streamlit web application

cd /d "%~dp0"

REM Check if virtual environment exists
if not exist "venv" (
    echo Creating virtual environment...
    python -m venv venv
)

REM Activate virtual environment
call venv\Scripts\activate.bat

REM Install requirements
pip install -q -r requirements.txt

echo.
echo ==========================================
echo   Delphi EAI Calculator
echo ==========================================
echo.
echo Starting web application...
echo The browser should open automatically.
echo If not, go to: http://localhost:8501
echo.
echo Press Ctrl+C to stop the server.
echo.

streamlit run app.py --server.headless false

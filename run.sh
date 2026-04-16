#!/bin/bash
# Delphi EAI Calculator - Run Script
# This script starts the Streamlit web application

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Check if virtual environment exists
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
fi

# Activate virtual environment
source venv/bin/activate

# Install requirements if needed
pip install -q -r requirements.txt

# Run the app
echo ""
echo "=========================================="
echo "  Delphi EAI Calculator"
echo "=========================================="
echo ""
echo "Starting web application..."
echo "The browser should open automatically."
echo "If not, go to: http://localhost:8501"
echo ""
echo "Press Ctrl+C to stop the server."
echo ""

streamlit run app.py --server.headless false

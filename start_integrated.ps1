# Integrated startup script for React + Flask molecular docking system
Write-Host "🧬 Starting AI-Enabled Molecular Docking System" -ForegroundColor Magenta
Write-Host "=============================================" -ForegroundColor Magenta

# Check if we're in the right directory
if (-not (Test-Path "frontend") -or -not (Test-Path "src")) {
    Write-Host "❌ Please run this script from the project root directory" -ForegroundColor Red
    exit 1
}

# Step 1: Build React frontend
Write-Host "`n📦 Building React frontend..." -ForegroundColor Green
Set-Location "frontend"

# Install dependencies if needed
if (-not (Test-Path "node_modules")) {
    Write-Host "Installing frontend dependencies..." -ForegroundColor Yellow
    npm install
    if ($LASTEXITCODE -ne 0) {
        Write-Host "❌ Failed to install frontend dependencies" -ForegroundColor Red
        exit 1
    }
}

# Build the React app
Write-Host "Building React app..." -ForegroundColor Yellow
npm run build
if ($LASTEXITCODE -ne 0) {
    Write-Host "❌ Frontend build failed!" -ForegroundColor Red
    exit 1
}

if (Test-Path "dist") {
    Write-Host "✅ Frontend build completed successfully!" -ForegroundColor Green
} else {
    Write-Host "❌ Frontend build failed - no dist folder found!" -ForegroundColor Red
    exit 1
}

# Return to root directory
Set-Location ".."

# Step 2: Check Python environment
Write-Host "`n🐍 Checking Python environment..." -ForegroundColor Green

# Check if virtual environment exists
if (Test-Path "venv") {
    Write-Host "Activating virtual environment..." -ForegroundColor Yellow
    & "venv\Scripts\Activate.ps1"
} elseif (Test-Path ".venv") {
    Write-Host "Activating virtual environment..." -ForegroundColor Yellow
    & ".venv\Scripts\Activate.ps1"
} else {
    Write-Host "⚠️  No virtual environment found. Using system Python." -ForegroundColor Yellow
}

# Install Python dependencies if requirements.txt exists
if (Test-Path "requirements.txt") {
    Write-Host "Installing Python dependencies..." -ForegroundColor Yellow
    pip install -r requirements.txt
}

# Step 3: Start the integrated Flask server
Write-Host "`n🚀 Starting integrated Flask server..." -ForegroundColor Green
Write-Host "The server will serve both the React frontend and API endpoints" -ForegroundColor Cyan
Write-Host "Frontend: http://localhost:5000" -ForegroundColor Cyan
Write-Host "API: http://localhost:5000/api" -ForegroundColor Cyan
Write-Host "`nPress Ctrl+C to stop the server" -ForegroundColor Yellow

# Start the Flask application
python -m src.api.app

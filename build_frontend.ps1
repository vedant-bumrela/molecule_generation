# Build script for React frontend
Write-Host "Building React frontend..." -ForegroundColor Green

# Navigate to frontend directory
Set-Location "frontend"

# Install dependencies if node_modules doesn't exist
if (-not (Test-Path "node_modules")) {
    Write-Host "Installing dependencies..." -ForegroundColor Yellow
    npm install
}

# Build the React app
Write-Host "Building React app..." -ForegroundColor Yellow
npm run build

# Check if build was successful
if (Test-Path "dist") {
    Write-Host "✅ Frontend build completed successfully!" -ForegroundColor Green
    Write-Host "Built files are in: frontend/dist" -ForegroundColor Cyan
} else {
    Write-Host "❌ Frontend build failed!" -ForegroundColor Red
    exit 1
}

# Return to root directory
Set-Location ".."

Write-Host "🚀 Ready to start the Flask server!" -ForegroundColor Green
Write-Host "Run: python -m src.api.app" -ForegroundColor Cyan

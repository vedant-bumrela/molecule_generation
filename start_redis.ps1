# Start Redis using Docker
Write-Host "Starting Redis container..." -ForegroundColor Green

# Check if Redis container already exists
$existingContainer = docker ps -a --filter "name=docking-redis" --format "{{.Names}}"

if ($existingContainer -eq "docking-redis") {
    Write-Host "Redis container already exists. Starting it..." -ForegroundColor Yellow
    docker start docking-redis
} else {
    Write-Host "Creating new Redis container..." -ForegroundColor Yellow
    docker run -d `
        --name docking-redis `
        -p 6379:6379 `
        redis:latest
}

# Wait a moment for Redis to start
Start-Sleep -Seconds 2

# Check if Redis is running
$status = docker ps --filter "name=docking-redis" --format "{{.Status}}"
if ($status -like "Up*") {
    Write-Host "[SUCCESS] Redis is running successfully on port 6379" -ForegroundColor Green
} else {
    Write-Host "[ERROR] Failed to start Redis" -ForegroundColor Red
    exit 1
}

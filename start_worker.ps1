# Start the docking worker
Write-Host "Starting docking worker..." -ForegroundColor Green

# Set environment variables
$env:REDIS_HOST = "localhost"
$env:REDIS_PORT = "6379"
$env:KMP_DUPLICATE_LIB_OK = "TRUE"

# Check if Redis is running
Write-Host "Checking Redis connection..." -ForegroundColor Yellow
try {
    $redisCheck = docker ps --filter "name=docking-redis" --format "{{.Status}}"
    if ($redisCheck -like "Up*") {
        Write-Host "[SUCCESS] Redis is running" -ForegroundColor Green
    } else {
        Write-Host "[ERROR] Redis is not running. Please start Redis first using .\start_redis.ps1" -ForegroundColor Red
        exit 1
    }
} catch {
    Write-Host "[ERROR] Could not check Redis status" -ForegroundColor Red
    exit 1
}

# Start the worker
Write-Host "Starting worker process..." -ForegroundColor Yellow
Write-Host "Press Ctrl+C to stop the worker" -ForegroundColor Cyan
Write-Host ""

python -m src.api.worker

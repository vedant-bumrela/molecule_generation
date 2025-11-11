# Start all services for the docking system
Write-Host "=== AI-Enabled Molecular Docking System ===" -ForegroundColor Cyan
Write-Host ""

# Start Redis
Write-Host "[1/3] Starting Redis..." -ForegroundColor Green
& .\start_redis.ps1
if ($LASTEXITCODE -ne 0) {
    Write-Host "Failed to start Redis. Exiting..." -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "[2/3] Starting Flask API..." -ForegroundColor Green
$apiJob = Start-Job -ScriptBlock {
    Set-Location $using:PWD
    $env:REDIS_HOST = "localhost"
    $env:REDIS_PORT = "6379"
    python src/api/app.py
}

# Wait for API to start
Start-Sleep -Seconds 3

Write-Host ""
Write-Host "[3/3] Starting Worker..." -ForegroundColor Green
$workerJob = Start-Job -ScriptBlock {
    Set-Location $using:PWD
    $env:REDIS_HOST = "localhost"
    $env:REDIS_PORT = "6379"
    python -m src.api.worker
}

# Wait for worker to start
Start-Sleep -Seconds 2

Write-Host ""
Write-Host "=== All Services Started ===" -ForegroundColor Green
Write-Host "[OK] Redis:   Running on port 6379" -ForegroundColor Green
Write-Host "[OK] API:     Running on http://localhost:5000" -ForegroundColor Green
Write-Host "[OK] Worker:  Processing jobs from queue" -ForegroundColor Green
Write-Host ""
Write-Host "View logs:" -ForegroundColor Yellow
Write-Host "  API:    Receive-Job -Id $($apiJob.Id) -Keep" -ForegroundColor Gray
Write-Host "  Worker: Receive-Job -Id $($workerJob.Id) -Keep" -ForegroundColor Gray
Write-Host ""
Write-Host "Stop all services:" -ForegroundColor Yellow
Write-Host "  Stop-Job -Id $($apiJob.Id),$($workerJob.Id)" -ForegroundColor Gray
Write-Host "  docker stop docking-redis" -ForegroundColor Gray
Write-Host ""
Write-Host "Press Ctrl+C to stop monitoring..." -ForegroundColor Cyan

# Monitor the jobs
try {
    while ($true) {
        Start-Sleep -Seconds 5
        
        # Check if jobs are still running
        $apiState = (Get-Job -Id $apiJob.Id).State
        $workerState = (Get-Job -Id $workerJob.Id).State
        
        if ($apiState -ne "Running") {
            Write-Host "[WARNING] API job stopped unexpectedly" -ForegroundColor Red
            Receive-Job -Id $apiJob.Id
            break
        }
        
        if ($workerState -ne "Running") {
            Write-Host "[WARNING] Worker job stopped unexpectedly" -ForegroundColor Red
            Receive-Job -Id $workerJob.Id
            break
        }
    }
} finally {
    Write-Host ""
    Write-Host "Cleaning up..." -ForegroundColor Yellow
    Stop-Job -Id $apiJob.Id,$workerJob.Id -ErrorAction SilentlyContinue
    Remove-Job -Id $apiJob.Id,$workerJob.Id -ErrorAction SilentlyContinue
    Write-Host "Jobs stopped. Redis is still running." -ForegroundColor Green
    Write-Host "To stop Redis: docker stop docking-redis" -ForegroundColor Gray
}

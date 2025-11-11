# Startup Guide - AI-Enabled Molecular Docking System

## Quick Start (Recommended)

### Option 1: Start Everything at Once
```powershell
.\start_all.ps1
```
This will start Redis, the Flask API, and the worker in one command.

### Option 2: Start Services Individually

#### Step 1: Start Redis
```powershell
.\start_redis.ps1
```

#### Step 2: Start the Flask API (in a new terminal)
```powershell
$env:REDIS_HOST = "localhost"
$env:REDIS_PORT = "6379"
python src/api/app.py
```

#### Step 3: Start the Worker (in another new terminal)
```powershell
.\start_worker.ps1
```

## Accessing the Application

Once all services are running, open your browser and go to:
- **Web Interface**: http://localhost:5000
- **API Health Check**: http://localhost:5000/api/health

## Stopping Services

### Stop All Services
```powershell
# Stop the Flask API and Worker (if running in background)
# Press Ctrl+C in their respective terminals

# Stop Redis
docker stop docking-redis
```

### Remove Redis Container (Optional)
```powershell
docker rm docking-redis
```

## Troubleshooting

### Redis Connection Issues
If you see "Redis not available" errors:

1. Check if Redis is running:
   ```powershell
   docker ps --filter "name=docking-redis"
   ```

2. Check Redis logs:
   ```powershell
   docker logs docking-redis
   ```

3. Restart Redis:
   ```powershell
   docker restart docking-redis
   ```

### Worker Not Processing Jobs
If jobs stay in "Pending" status:

1. Check if the worker is running and look for errors in the console

2. Check worker logs (if running in background):
   ```powershell
   # Look for worker.log file
   Get-Content worker.log -Tail 50
   ```

3. Verify Redis connection:
   ```powershell
   docker exec -it docking-redis redis-cli ping
   # Should return: PONG
   ```

### Port Already in Use
If port 5000 or 6379 is already in use:

**For Flask API (port 5000):**
- Find and stop the process using port 5000
- Or modify the port in `src/api/app.py` (line 565)

**For Redis (port 6379):**
```powershell
# Find what's using the port
netstat -ano | findstr :6379

# Stop the existing Redis container
docker stop docking-redis
docker rm docking-redis

# Start fresh
.\start_redis.ps1
```

## Service Status Check

### Check All Services
```powershell
# Check Redis
docker ps --filter "name=docking-redis"

# Check if Flask API is responding
curl http://localhost:5000/api/health

# Check if worker is running (look for python process)
Get-Process python | Where-Object {$_.CommandLine -like "*worker*"}
```

## Environment Variables

The following environment variables are used:

- `REDIS_HOST`: Redis server hostname (default: localhost)
- `REDIS_PORT`: Redis server port (default: 6379)
- `UPLOAD_FOLDER`: Directory for uploaded files
- `RESULTS_FOLDER`: Directory for docking results
- `BATCH_FOLDER`: Directory for batch processing

## Development Mode

For development with auto-reload:

```powershell
# Terminal 1: Redis
.\start_redis.ps1

# Terminal 2: Flask API with debug mode
$env:FLASK_ENV = "development"
$env:REDIS_HOST = "localhost"
python src/api/app.py

# Terminal 3: Worker
.\start_worker.ps1
```

## Production Deployment

For production, use Docker Compose:

```powershell
docker-compose up -d
```

This will start all services (Redis, API, Worker) in containers with proper configuration.

## Next Steps

1. Submit a test job through the web interface
2. Monitor the worker console for processing logs
3. Check the Job Status tab to see your job progress
4. View results in the Results tab once completed

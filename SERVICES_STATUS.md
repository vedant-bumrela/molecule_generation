# Services Status - AI-Enabled Molecular Docking System

## ✅ Currently Running Services

### 1. Redis (Docker Container)
- **Status**: ✅ RUNNING
- **Port**: 6379
- **Container Name**: docking-redis
- **Check Status**: `docker ps --filter "name=docking-redis"`

### 2. Worker Process
- **Status**: ✅ RUNNING  
- **Process ID**: Background Job 68
- **Function**: Processes docking jobs from Redis queue
- **Logs**: Check terminal output or `worker.log`

### 3. Flask API
- **Status**: ⚠️ NEEDS TO BE STARTED
- **Port**: 5000
- **URL**: http://localhost:5000

## 🚀 Next Steps

### Start the Flask API

Open a **NEW terminal** and run:

```powershell
cd c:\docking
$env:REDIS_HOST = "localhost"
$env:REDIS_PORT = "6379"
python src/api/app.py
```

Or simply run:
```powershell
python src/api/app.py
```

The API will automatically connect to Redis on localhost:6379.

### Access the Application

Once the Flask API is running, open your browser:
- **Web Interface**: http://localhost:5000
- **API Health Check**: http://localhost:5000/api/health

## 📊 Monitor Services

### Check Worker Status
The worker is running in the background. To see its output:
```powershell
# View worker logs in real-time
Get-Content worker.log -Wait -Tail 50
```

### Check Redis Status
```powershell
docker ps --filter "name=docking-redis"
docker logs docking-redis --tail 50
```

### Test Redis Connection
```powershell
docker exec -it docking-redis redis-cli ping
# Should return: PONG
```

## 🔄 Restart Services

### Restart Worker
```powershell
# Stop the current worker (Ctrl+C in its terminal)
# Then restart:
.\start_worker.ps1
```

### Restart Redis
```powershell
docker restart docking-redis
```

## 🛑 Stop Services

### Stop Worker
Press `Ctrl+C` in the worker terminal

### Stop Flask API
Press `Ctrl+C` in the API terminal

### Stop Redis
```powershell
docker stop docking-redis
```

### Stop All Services
```powershell
# Stop Redis
docker stop docking-redis

# Worker and API will stop when you close their terminals or press Ctrl+C
```

## 🐛 Troubleshooting

### Jobs Still Stuck on Pending?

1. **Verify all services are running**:
   - Redis: `docker ps --filter "name=docking-redis"`
   - Worker: Check if process is running
   - API: Check if http://localhost:5000 responds

2. **Check worker logs for errors**:
   ```powershell
   Get-Content worker.log -Tail 50
   ```

3. **Test Redis connection from worker**:
   ```powershell
   docker exec -it docking-redis redis-cli
   # In Redis CLI:
   LLEN docking_tasks  # Should show number of pending tasks
   KEYS result:*       # Should show completed job results
   ```

4. **Submit a new test job** after all services are running

### Port Conflicts

If you see "port already in use" errors:

**Port 5000 (Flask API)**:
```powershell
# Find process using port 5000
netstat -ano | findstr :5000
# Kill the process if needed
```

**Port 6379 (Redis)**:
```powershell
# Check if Redis is already running
docker ps -a --filter "name=docking-redis"
# Remove old container if needed
docker rm -f docking-redis
```

## 📝 Notes

- The worker is configured to process jobs from the `docking_tasks` queue
- Job results are stored in Redis with key prefix `result:`
- Results expire after 7 days
- The worker uses 8 threads by default (based on your CPU cores)
- OpenMP library conflict is handled with `KMP_DUPLICATE_LIB_OK=TRUE`

## ✨ Ready to Test!

Once you start the Flask API, you can:
1. Open http://localhost:5000 in your browser
2. Submit a docking job
3. Watch it process in real-time
4. The status should change from "Pending" → "Running" → "Completed"

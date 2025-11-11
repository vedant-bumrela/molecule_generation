# Fix: Jobs Stuck on "Pending" Status

## Problem Identified

Jobs were completing successfully in the worker but showing as "Pending" in the web UI.

### Root Cause

**Disconnect between worker and API**:
- Worker processes jobs and stores results in **Redis** with key `result:<job_id>`
- Flask API only checked the **in-memory `jobs` dictionary**
- The API never read results back from Redis
- Result: Jobs appeared stuck on "Pending" even though they were completed

### Evidence

```bash
# Worker logs showed successful completion:
worker-1 | 2025-11-11 16:11:52,046 - docking_worker - INFO - Completed docking task 9b770212...

# Redis had the results:
$ docker exec docking-redis-1 redis-cli KEYS "result:*"
result:9b770212-d489-4722-8161-ebc6a2f2413b
result:0fa4019e-9723-4332-b6d1-7b37863a9ab4

# But API still returned "pending" status
```

---

## Solution Implemented

### Changes Made to `src/api/app.py`

#### 1. Updated `/api/status/<job_id>` Endpoint

**Before**: Only checked in-memory `jobs` dictionary

**After**: 
1. First checks Redis for results (where worker stores them)
2. Falls back to in-memory dictionary if Redis unavailable
3. Returns completed status with scores if found in Redis

```python
@app.route('/api/status/<job_id>', methods=['GET'])
def job_status(job_id):
    # First check Redis for results (worker stores them there)
    try:
        redis_client.ping()
        result_key = f"{RESULT_KEY_PREFIX}{job_id}"
        redis_result = redis_client.get(result_key)
        
        if redis_result:
            # Job completed by worker, results in Redis
            result_data = json.loads(redis_result.decode('utf-8'))
            # Return completed status with results
            ...
    except redis.exceptions.ConnectionError:
        # Fallback to in-memory storage
        ...
```

#### 2. Updated `/api/jobs` Endpoint

**Before**: Only returned status from in-memory dictionary

**After**:
1. Gets job list from in-memory dictionary
2. For each job, checks Redis for updated status
3. Updates status to "completed" if found in Redis

```python
@app.route('/api/jobs', methods=['GET'])
def list_jobs():
    job_list = []
    
    with job_lock:
        for job_id, job in jobs.items():
            # Get base status from memory
            job_status_data = {...}
            
            # Try to get updated status from Redis
            try:
                result_key = f"{RESULT_KEY_PREFIX}{job_id}"
                redis_result = redis_client.get(result_key)
                
                if redis_result:
                    # Update status to completed
                    result_data = json.loads(redis_result.decode('utf-8'))
                    job_status_data['status'] = result_data.get('status', 'completed')
            except Exception:
                pass
            
            job_list.append(job_status_data)
    
    return jsonify({'jobs': job_list})
```

---

## How It Works Now

### Complete Flow

```
1. User submits job via web UI
   ↓
2. Flask API creates job entry (status: pending)
   ↓
3. Job added to Redis queue
   ↓
4. Worker picks up job from queue
   ↓
5. Worker processes docking
   ↓
6. Worker stores results in Redis: result:<job_id>
   ↓
7. User refreshes job list
   ↓
8. API checks Redis for each job
   ↓
9. API finds result in Redis
   ↓
10. API returns status: completed ✅
```

### Data Storage

```
In-Memory (Flask API):
├── jobs[job_id] = {
│       'status': 'pending',
│       'job_type': 'vina',
│       'protein_path': '...',
│       'ligand_path': '...'
│   }

Redis (Worker writes, API reads):
├── result:<job_id> = {
│       'status': 'completed',
│       'results': {
│           'scores': [-7.2, -6.8, -6.5, ...],
│           'output_file': '/path/to/docked.pdbqt',
│           'best_score': -7.2
│       }
│   }
```

---

## Testing

### Before Fix
```
Status: Pending (even after completion)
Message: Job submitted
```

### After Fix
```
Status: Completed ✅
Message: Job completed successfully
Results: {
  "num_poses": 9,
  "scores": [-7.2, -6.8, -6.5, ...]
}
```

---

## Deployment

### Apply the Fix

```bash
# Rebuild and restart the app container
docker-compose up -d --build app

# Or restart all services
docker-compose restart
```

### Verify Fix

1. Submit a new job via web UI
2. Wait ~30 seconds for processing
3. Refresh the Job Status page
4. Status should change from "Pending" → "Completed"

---

## Why This Happened

### Architecture Issue

The system uses a **distributed architecture**:
- Flask API runs in one container
- Worker runs in another container
- They communicate via Redis

**The bug**: API wasn't reading back from the communication channel (Redis)

### Lesson Learned

When using Redis as a job queue:
1. ✅ Workers write results to Redis
2. ✅ API must read results from Redis
3. ❌ Don't rely only on in-memory state in the API

---

## Additional Improvements Made

### Error Handling

- Added try-catch for Redis connection errors
- Graceful fallback to in-memory storage
- Logging for debugging

### Status Synchronization

- API now checks Redis on every status request
- Job list endpoint also checks Redis
- Ensures UI always shows latest status

---

## Future Enhancements

### Short Term
- [ ] Add "running" status detection (check queue length)
- [ ] Implement job progress tracking
- [ ] Add timestamp for job completion

### Long Term
- [ ] Use Redis pub/sub for real-time updates
- [ ] Add database for persistent job history
- [ ] Implement job cleanup/archival

---

## Files Modified

- `src/api/app.py` (2 functions updated)
  - `job_status()` - Lines 477-535
  - `list_jobs()` - Lines 573-611

---

## Verification Commands

```bash
# Check if results are in Redis
docker exec docking-redis-1 redis-cli KEYS "result:*"

# Get specific job result
docker exec docking-redis-1 redis-cli GET "result:<job_id>"

# View API logs
docker-compose logs app --tail 50

# View worker logs
docker-compose logs worker --tail 50

# Test API endpoint directly
curl http://localhost:5000/api/status/<job_id>
```

---

**Status**: ✅ **Fixed and Deployed**
**Date**: November 11, 2025
**Impact**: All jobs now show correct completion status

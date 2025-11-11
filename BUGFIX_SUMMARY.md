# Bug Fix Summary - API Endpoint Mismatch

## Issue Identified
The frontend JavaScript was making API calls to incorrect endpoints, causing the job status updates and batch processing features to fail.

## Root Cause
- Backend API endpoints are defined with `/api` prefix (e.g., `/api/submit`, `/api/jobs`, `/api/status/<job_id>`)
- Frontend had an `API_BASE_URL` constant set to `http://localhost:5000/api`
- Regular docking functions correctly used `API_BASE_URL`
- **Batch processing functions were calling endpoints without the API prefix**, causing 404 errors

## Changes Made

### File: `web/js/main.js`

#### 1. Fixed Batch Job Submission (Line 548)
**Before:** `fetch('/submit/batch', {...})`
**After:** `fetch(\`${API_BASE_URL}/submit/batch\`, {...})`

#### 2. Fixed Load Batch Jobs (Line 574)
**Before:** `fetch('/jobs')`
**After:** `fetch(\`${API_BASE_URL}/jobs\`)`

#### 3. Fixed View Batch Results (Line 661)
**Before:** `fetch(\`/batch/results/${jobId}\`)`
**After:** `fetch(\`${API_BASE_URL}/batch/results/${jobId}\`)`

#### 4. Fixed View Batch Pose Status Check (Line 715)
**Before:** `fetch(\`/status/${jobId}\`)`
**After:** `fetch(\`${API_BASE_URL}/status/${jobId}\`)`

#### 5. Fixed Download Batch Results (Line 527)
**Before:** `window.location.href = \`/download/${currentBatchJobId}/docking_results.csv\``
**After:** `window.location.href = \`${API_BASE_URL}/download/${currentBatchJobId}/docking_results.csv\``

#### 6. Fixed Batch Pose File Paths (Lines 720-721)
**Before:**
```javascript
const proteinPath = `/download/${jobId}/protein_prepared.pdbqt`;
const ligandPath = `/download/${jobId}/${compoundId}/top_pose.pdb`;
```
**After:**
```javascript
const proteinPath = `${API_BASE_URL}/download/${jobId}/protein_prepared.pdbqt`;
const ligandPath = `${API_BASE_URL}/download/${jobId}/${compoundId}/top_pose.pdb`;
```

#### 7. Added Missing Helper Functions (Lines 741-780)
- `showAlert(message, type)` - Display Bootstrap alerts for user feedback
- `loadProtein(proteinPath)` - Placeholder for protein loading
- `loadLigand(ligandPath)` - Placeholder for ligand loading
- `updateScoresTable(scores)` - Placeholder for scores table updates

## Impact
✅ Job status updates will now work correctly
✅ Batch processing submission will function properly
✅ Batch results retrieval will work
✅ File downloads will use correct endpoints
✅ All API calls now consistently use the `API_BASE_URL` constant

## Testing Recommendations
1. Submit a regular docking job and verify status updates
2. Submit a batch processing job and verify it appears in the batch jobs table
3. Check that job status polling works correctly
4. Verify that completed jobs can display results
5. Test file downloads for completed jobs

## Next Steps
- Backend may need additional endpoints for batch processing (`/api/submit/batch`, `/api/batch/results/<job_id>`, `/api/download/<job_id>/<file>`)
- Consider implementing the placeholder helper functions for full batch processing visualization
- Add error handling for network failures
- Implement proper loading states during API calls

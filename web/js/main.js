/**
 * AI-Enabled Molecular Docking System
 * Frontend JavaScript
 */

// Global variables
let stage = null;
let currentJobId = null;
let proteinComponent = null;
let ligandComponent = null;

// API endpoint (change this to match your deployment)
const API_BASE_URL = 'http://localhost:5000/api';

// Initialize the 3D viewer
function initViewer() {
    stage = new NGL.Stage('viewer', { backgroundColor: 'white' });
    stage.setParameters({
        cameraType: 'perspective',
        cameraFov: 40
    });
    
    // Handle window resize
    window.addEventListener('resize', function() {
        stage.handleResize();
    }, false);
}

// Initialize event listeners
function initEventListeners() {
    // Form submission
    document.getElementById('dockingForm').addEventListener('submit', submitJob);
    
    // Ligand input type toggle
    document.querySelectorAll('input[name="ligandInputType"]').forEach(radio => {
        radio.addEventListener('change', toggleLigandInput);
    });
    
    // Advanced options toggle
    document.getElementById('advancedOptions').addEventListener('change', toggleAdvancedOptions);
    
    // Job type change
    document.getElementById('jobType').addEventListener('change', updateMethodInfo);
    
    // Refresh jobs button
    document.getElementById('refreshJobsBtn').addEventListener('click', loadJobs);
    
    // Viewer controls
    document.getElementById('showProtein').addEventListener('click', () => toggleVisibility('protein'));
    document.getElementById('showLigand').addEventListener('click', () => toggleVisibility('ligand'));
    document.getElementById('showBoth').addEventListener('click', () => toggleVisibility('both'));
    document.getElementById('viewCartoon').addEventListener('click', () => changeProteinRepresentation('cartoon'));
    document.getElementById('viewSurface').addEventListener('click', () => changeProteinRepresentation('surface'));
    
    // Tab change events
    document.getElementById('jobs-tab').addEventListener('click', loadJobs);
}

// Toggle ligand input type (file or SMILES)
function toggleLigandInput() {
    const inputType = document.querySelector('input[name="ligandInputType"]:checked').value;
    
    if (inputType === 'file') {
        document.getElementById('ligandFileInput').classList.remove('d-none');
        document.getElementById('ligandSmilesInput').classList.add('d-none');
    } else {
        document.getElementById('ligandFileInput').classList.add('d-none');
        document.getElementById('ligandSmilesInput').classList.remove('d-none');
    }
}

// Toggle advanced options panel
function toggleAdvancedOptions() {
    const showAdvanced = document.getElementById('advancedOptions').checked;
    
    if (showAdvanced) {
        document.getElementById('advancedOptionsPanel').classList.remove('d-none');
    } else {
        document.getElementById('advancedOptionsPanel').classList.add('d-none');
    }
}

// Update method information based on selected job type
function updateMethodInfo() {
    const jobType = document.getElementById('jobType').value;
    const methodInfoDiv = document.getElementById('methodInfo');
    
    // Show/hide method-specific options
    if (jobType === 'vina') {
        document.getElementById('vinaOptions').classList.remove('d-none');
        document.getElementById('diffdockOptions').classList.add('d-none');
    } else if (jobType === 'diffdock') {
        document.getElementById('vinaOptions').classList.add('d-none');
        document.getElementById('diffdockOptions').classList.remove('d-none');
    } else {
        document.getElementById('vinaOptions').classList.add('d-none');
        document.getElementById('diffdockOptions').classList.add('d-none');
    }
    
    // Update method information
    let methodHtml = '';
    
    if (jobType === 'vina') {
        methodHtml = `
            <h5>AutoDock Vina</h5>
            <p>Classical molecular docking engine that uses a scoring function and search algorithm to predict binding poses and affinities.</p>
            <ul>
                <li><strong>Speed:</strong> Fast</li>
                <li><strong>Accuracy:</strong> Good baseline</li>
                <li><strong>Best for:</strong> Standard docking, virtual screening</li>
            </ul>
        `;
    } else if (jobType === 'equibind') {
        methodHtml = `
            <h5>EquiBind</h5>
            <p>SE(3)-equivariant graph neural network for direct pose prediction without sampling.</p>
            <ul>
                <li><strong>Speed:</strong> Very fast</li>
                <li><strong>Accuracy:</strong> Good for initial poses</li>
                <li><strong>Best for:</strong> Blind docking, fast pose generation</li>
            </ul>
        `;
    } else if (jobType === 'diffdock') {
        methodHtml = `
            <h5>DiffDock</h5>
            <p>Diffusion-based generative model for pose prediction with confidence estimates.</p>
            <ul>
                <li><strong>Speed:</strong> Moderate</li>
                <li><strong>Accuracy:</strong> State-of-the-art</li>
                <li><strong>Best for:</strong> High-accuracy pose prediction, confidence estimation</li>
            </ul>
        `;
    }
    
    methodInfoDiv.innerHTML = methodHtml;
}

// Submit docking job
async function submitJob(event) {
    event.preventDefault();
    
    const form = document.getElementById('dockingForm');
    const formData = new FormData(form);
    
    // Check if required fields are filled
    const jobType = formData.get('job_type');
    const proteinFile = document.getElementById('proteinFile').files[0];
    
    if (!proteinFile) {
        alert('Please select a protein file');
        return;
    }
    
    // Check ligand input
    const ligandInputType = document.querySelector('input[name="ligandInputType"]:checked').value;
    
    if (ligandInputType === 'file') {
        const ligandFile = document.getElementById('ligandFileUpload').files[0];
        if (!ligandFile) {
            alert('Please select a ligand file');
            return;
        }
    } else {
        const smiles = document.getElementById('smilesInput').value.trim();
        if (!smiles) {
            alert('Please enter a SMILES string');
            return;
        }
    }
    
    try {
        // Show loading indicator
        form.querySelector('button[type="submit"]').disabled = true;
        form.querySelector('button[type="submit"]').innerHTML = '<span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span> Submitting...';
        
        // Submit job
        const response = await fetch(`${API_BASE_URL}/submit`, {
            method: 'POST',
            body: formData
        });
        
        if (!response.ok) {
            throw new Error(`HTTP error ${response.status}`);
        }
        
        const data = await response.json();
        
        // Reset form
        form.querySelector('button[type="submit"]').disabled = false;
        form.querySelector('button[type="submit"]').textContent = 'Submit Job';
        
        // Show success message
        alert(`Job submitted successfully! Job ID: ${data.job_id}`);
        
        // Switch to jobs tab
        document.getElementById('jobs-tab').click();
        
    } catch (error) {
        console.error('Error submitting job:', error);
        alert('Error submitting job: ' + error.message);
        
        // Reset button
        form.querySelector('button[type="submit"]').disabled = false;
        form.querySelector('button[type="submit"]').textContent = 'Submit Job';
    }
}

// Load jobs list
async function loadJobs() {
    const jobsListDiv = document.getElementById('jobsList');
    
    // Show loading indicator
    jobsListDiv.innerHTML = `
        <div class="loading">
            <div class="spinner-border text-primary" role="status">
                <span class="visually-hidden">Loading...</span>
            </div>
        </div>
    `;
    
    try {
        const response = await fetch(`${API_BASE_URL}/jobs`);
        
        if (!response.ok) {
            throw new Error(`HTTP error ${response.status}`);
        }
        
        const data = await response.json();
        
        // Display jobs
        if (data.jobs && data.jobs.length > 0) {
            let jobsHtml = '';
            
            data.jobs.forEach(job => {
                let statusBadge = '';
                
                switch (job.status) {
                    case 'pending':
                        statusBadge = '<span class="badge bg-warning">Pending</span>';
                        break;
                    case 'running':
                        statusBadge = '<span class="badge bg-primary">Running</span>';
                        break;
                    case 'completed':
                        statusBadge = '<span class="badge bg-success">Completed</span>';
                        break;
                    case 'failed':
                        statusBadge = '<span class="badge bg-danger">Failed</span>';
                        break;
                    default:
                        statusBadge = '<span class="badge bg-secondary">Unknown</span>';
                }
                
                jobsHtml += `
                    <div class="col-md-6">
                        <div class="card job-card" data-job-id="${job.job_id}">
                            <div class="card-body">
                                <div class="d-flex justify-content-between align-items-center">
                                    <h6 class="card-title mb-0">${job.job_type.toUpperCase()} Job</h6>
                                    ${statusBadge}
                                </div>
                                <p class="card-text text-muted small mt-2">Job ID: ${job.job_id}</p>
                                <p class="card-text">${job.message || 'No status message'}</p>
                                <button class="btn btn-sm btn-outline-primary check-status-btn">Check Status</button>
                                ${job.status === 'completed' ? '<button class="btn btn-sm btn-success ms-2 view-results-btn">View Results</button>' : ''}
                            </div>
                        </div>
                    </div>
                `;
            });
            
            jobsListDiv.innerHTML = jobsHtml;
            
            // Add event listeners to job cards
            document.querySelectorAll('.check-status-btn').forEach(button => {
                button.addEventListener('click', function() {
                    const jobId = this.closest('.job-card').dataset.jobId;
                    checkJobStatus(jobId);
                });
            });
            
            document.querySelectorAll('.view-results-btn').forEach(button => {
                button.addEventListener('click', function() {
                    const jobId = this.closest('.job-card').dataset.jobId;
                    loadJobResults(jobId);
                });
            });
            
        } else {
            jobsListDiv.innerHTML = '<p class="text-muted">No jobs found</p>';
        }
        
    } catch (error) {
        console.error('Error loading jobs:', error);
        jobsListDiv.innerHTML = `<p class="text-danger">Error loading jobs: ${error.message}</p>`;
    }
}

// Check job status
async function checkJobStatus(jobId) {
    try {
        const response = await fetch(`${API_BASE_URL}/status/${jobId}`);
        
        if (!response.ok) {
            throw new Error(`HTTP error ${response.status}`);
        }
        
        const data = await response.json();
        
        // Show status in alert
        alert(`Job Status: ${data.status}\nMessage: ${data.message}`);
        
        // Refresh jobs list
        loadJobs();
        
    } catch (error) {
        console.error('Error checking job status:', error);
        alert('Error checking job status: ' + error.message);
    }
}

// Load job results
async function loadJobResults(jobId) {
    currentJobId = jobId;
    
    try {
        const response = await fetch(`${API_BASE_URL}/status/${jobId}`);
        
        if (!response.ok) {
            throw new Error(`HTTP error ${response.status}`);
        }
        
        const data = await response.json();
        
        if (data.status !== 'completed') {
            alert('Job is not completed yet');
            return;
        }
        
        // Switch to results tab
        document.getElementById('results-tab').click();
        
        // Update result info
        const resultInfoDiv = document.getElementById('resultInfo');
        resultInfoDiv.innerHTML = `
            <h6>Job Details</h6>
            <p><strong>Job ID:</strong> ${jobId}</p>
            <p><strong>Status:</strong> ${data.status}</p>
            <p><strong>Message:</strong> ${data.message}</p>
        `;
        
        // Update download links
        const downloadPdbqt = document.getElementById('downloadPdbqt');
        const downloadPdb = document.getElementById('downloadPdb');
        
        downloadPdbqt.href = `${API_BASE_URL}/result/${jobId}/pdbqt`;
        downloadPdbqt.classList.remove('d-none');
        
        downloadPdb.href = `${API_BASE_URL}/result/${jobId}/pdb`;
        downloadPdb.classList.remove('d-none');
        
        // Update score table if available
        const scoreTable = document.getElementById('scoreTable');
        const scoreTableBody = document.getElementById('scoreTableBody');
        
        if (data.results && data.results.scores) {
            scoreTable.classList.remove('d-none');
            
            let scoresHtml = '';
            data.results.scores.forEach((score, index) => {
                scoresHtml += `
                    <tr>
                        <td>${index + 1}</td>
                        <td>${score.toFixed(2)}</td>
                        <td><button class="btn btn-sm btn-outline-primary view-pose-btn" data-pose="${index + 1}">View</button></td>
                    </tr>
                `;
            });
            
            scoreTableBody.innerHTML = scoresHtml;
            
            // Add event listeners to view pose buttons
            document.querySelectorAll('.view-pose-btn').forEach(button => {
                button.addEventListener('click', function() {
                    const poseNum = this.dataset.pose;
                    viewPose(jobId, poseNum);
                });
            });
        } else {
            scoreTable.classList.add('d-none');
        }
        
        // Load 3D visualization
        loadVisualization(jobId);
        
    } catch (error) {
        console.error('Error loading job results:', error);
        alert('Error loading job results: ' + error.message);
    }
}

// Load 3D visualization
async function loadVisualization(jobId) {
    // Clear existing components
    if (stage) {
        stage.removeAllComponents();
    } else {
        initViewer();
    }
    
    try {
        // Load protein
        const proteinUrl = `${API_BASE_URL}/result/${jobId}/protein`;
        
        try {
            proteinComponent = await stage.loadFile(proteinUrl, { ext: 'pdb' });
            proteinComponent.addRepresentation('cartoon', { color: 'chainname' });
            proteinComponent.addRepresentation('licorice', { 
                sele: 'hetero and not water', 
                multipleBond: true 
            });
        } catch (error) {
            console.warn('Error loading protein:', error);
        }
        
        // Load ligand
        const ligandUrl = `${API_BASE_URL}/result/${jobId}/pdb`;
        
        try {
            ligandComponent = await stage.loadFile(ligandUrl, { ext: 'pdb' });
            ligandComponent.addRepresentation('licorice', { 
                multipleBond: true,
                colorScheme: 'element'
            });
        } catch (error) {
            console.warn('Error loading ligand:', error);
        }
        
        // Adjust view
        stage.autoView();
        
    } catch (error) {
        console.error('Error loading visualization:', error);
    }
}

// View specific pose
async function viewPose(jobId, poseNum) {
    // TODO: Implement pose selection
    alert(`Viewing pose ${poseNum} is not implemented yet`);
}

// Toggle visibility of protein/ligand
function toggleVisibility(what) {
    if (!stage) return;
    
    if (what === 'protein' || what === 'both') {
        if (proteinComponent) {
            proteinComponent.setVisibility(true);
        }
    } else {
        if (proteinComponent) {
            proteinComponent.setVisibility(false);
        }
    }
    
    if (what === 'ligand' || what === 'both') {
        if (ligandComponent) {
            ligandComponent.setVisibility(true);
        }
    } else {
        if (ligandComponent) {
            ligandComponent.setVisibility(false);
        }
    }
}

// Change protein representation
function changeProteinRepresentation(type) {
    if (!proteinComponent) return;
    
    proteinComponent.removeAllRepresentations();
    
    if (type === 'cartoon') {
        proteinComponent.addRepresentation('cartoon', { color: 'chainname' });
        proteinComponent.addRepresentation('licorice', { 
            sele: 'hetero and not water', 
            multipleBond: true 
        });
    } else if (type === 'surface') {
        proteinComponent.addRepresentation('surface', { 
            opacity: 0.7,
            colorScheme: 'bfactor'
        });
    }
}

// Initialize when DOM is loaded
document.addEventListener('DOMContentLoaded', function() {
    initViewer();
    initEventListeners();
    updateMethodInfo();
    
    // Batch form submission
    const batchForm = document.getElementById('batch-form');
    if (batchForm) {
        batchForm.addEventListener('submit', function(e) {
            e.preventDefault();
            submitBatchJob();
        });
    }

    // Load batch jobs on tab activation
    const batchTab = document.getElementById('batch-tab');
    if (batchTab) {
        batchTab.addEventListener('shown.bs.tab', function() {
            loadBatchJobs();
        });
    }

    // Download batch results
    const downloadBatchResults = document.getElementById('download-batch-results');
    if (downloadBatchResults) {
        downloadBatchResults.addEventListener('click', function(e) {
            e.preventDefault();
            if (currentBatchJobId) {
                window.location.href = `${API_BASE_URL}/download/${currentBatchJobId}/docking_results.csv`;
            }
        });
    }
});

// Global variables for batch processing
let batchJobsPollingInterval = null;
let currentBatchJobId = null;

// Submit batch job
function submitBatchJob() {
    const batchForm = document.getElementById('batch-form');
    const formData = new FormData(batchForm);
    
    // Disable submit button
    const submitBtn = document.getElementById('batch-submit-btn');
    submitBtn.disabled = true;
    submitBtn.innerHTML = '<span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span> Submitting...';
    
    // Submit job
    fetch(`${API_BASE_URL}/submit/batch`, {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        if (data.job_id) {
            showAlert('Batch job submitted successfully!', 'success');
            loadBatchJobs();
            startBatchJobsPolling();
        } else {
            showAlert('Error submitting batch job: ' + (data.error || 'Unknown error'), 'danger');
        }
    })
    .catch(error => {
        showAlert('Error submitting batch job: ' + error, 'danger');
    })
    .finally(() => {
        // Re-enable submit button
        submitBtn.disabled = false;
        submitBtn.innerHTML = 'Submit Batch Job';
    });
}

// Load batch jobs
function loadBatchJobs() {
    fetch(`${API_BASE_URL}/jobs`)
    .then(response => response.json())
    .then(data => {
        const batchJobsTable = document.getElementById('batch-jobs-table').getElementsByTagName('tbody')[0];
        batchJobsTable.innerHTML = '';
        
        // Filter for batch jobs
        const batchJobs = data.jobs.filter(job => job.job_type === 'batch');
        
        if (batchJobs.length === 0) {
            const row = batchJobsTable.insertRow();
            const cell = row.insertCell();
            cell.colSpan = 5;
            cell.textContent = 'No batch jobs found.';
            cell.className = 'text-center';
        } else {
            batchJobs.forEach(job => {
                const row = batchJobsTable.insertRow();
                
                // Job ID
                const idCell = row.insertCell();
                idCell.textContent = job.job_id.substring(0, 8) + '...';
                idCell.title = job.job_id;
                
                // Method
                const methodCell = row.insertCell();
                methodCell.textContent = job.method;
                
                // Status
                const statusCell = row.insertCell();
                let statusBadgeClass = 'bg-secondary';
                if (job.status === 'completed') statusBadgeClass = 'bg-success';
                if (job.status === 'failed') statusBadgeClass = 'bg-danger';
                if (job.status === 'processing') statusBadgeClass = 'bg-primary';
                statusCell.innerHTML = `<span class="badge ${statusBadgeClass}">${job.status}</span>`;
                
                // Progress
                const progressCell = row.insertCell();
                if (job.progress !== undefined) {
                    progressCell.innerHTML = `
                        <div class="progress">
                            <div class="progress-bar" role="progressbar" style="width: ${job.progress}%;" 
                                aria-valuenow="${job.progress}" aria-valuemin="0" aria-valuemax="100">
                                ${job.progress}%
                            </div>
                        </div>
                        <small>${job.processed_compounds || 0} / ${job.total_compounds || 0} compounds</small>
                    `;
                } else {
                    progressCell.textContent = 'N/A';
                }
                
                // Actions
                const actionsCell = row.insertCell();
                let actionsHtml = '';
                
                if (job.status === 'completed') {
                    actionsHtml += `<button class="btn btn-sm btn-primary me-1" onclick="viewBatchResults('${job.job_id}')">View Results</button>`;
                }
                
                actionsCell.innerHTML = actionsHtml || 'No actions available';
            });
            
            // Start polling for updates
            startBatchJobsPolling();
        }
    })
    .catch(error => {
        showAlert('Error loading batch jobs: ' + error, 'danger');
    });
}

// Start polling for batch job updates
function startBatchJobsPolling() {
    // Clear existing interval
    if (batchJobsPollingInterval) {
        clearInterval(batchJobsPollingInterval);
    }
    
    // Start new polling interval
    batchJobsPollingInterval = setInterval(loadBatchJobs, 10000); // Poll every 10 seconds
}

// View batch results
function viewBatchResults(jobId) {
    currentBatchJobId = jobId;
    
    fetch(`${API_BASE_URL}/batch/results/${jobId}`)
    .then(response => response.json())
    .then(data => {
        if (data.results) {
            // Show results container
            document.getElementById('batch-results-container').style.display = 'block';
            
            // Populate results table
            const resultsTable = document.getElementById('batch-results-table').getElementsByTagName('tbody')[0];
            resultsTable.innerHTML = '';
            
            data.results.forEach(result => {
                const row = resultsTable.insertRow();
                
                // Compound ID
                const idCell = row.insertCell();
                idCell.textContent = result.id;
                
                // Score
                const scoreCell = row.insertCell();
                if (result.top_score !== undefined) {
                    scoreCell.textContent = result.top_score.toFixed(2);
                } else if (result.top_confidence !== undefined) {
                    scoreCell.textContent = `Confidence: ${result.top_confidence.toFixed(2)}`;
                } else {
                    scoreCell.textContent = 'N/A';
                }
                
                // Actions
                const actionsCell = row.insertCell();
                if (result.top_pose_path) {
                    actionsCell.innerHTML = `<button class="btn btn-sm btn-primary" onclick="viewBatchPose('${jobId}', '${result.id}')">View Pose</button>`;
                } else {
                    actionsCell.textContent = 'No pose available';
                }
            });
            
            // Scroll to results
            document.getElementById('batch-results-container').scrollIntoView({ behavior: 'smooth' });
        } else {
            showAlert('Error loading batch results: ' + (data.error || 'Unknown error'), 'danger');
        }
    })
    .catch(error => {
        showAlert('Error loading batch results: ' + error, 'danger');
    });
}

// View batch pose
function viewBatchPose(jobId, compoundId) {
    // Switch to results tab
    document.getElementById('results-tab').click();
    
    // Load protein and ligand
    fetch(`${API_BASE_URL}/status/${jobId}`)
    .then(response => response.json())
    .then(data => {
        if (data.status === 'completed') {
            // Get paths to files
            const proteinPath = `${API_BASE_URL}/download/${jobId}/protein_prepared.pdbqt`;
            const ligandPath = `${API_BASE_URL}/download/${jobId}/${compoundId}/top_pose.pdb`;
            
            // Load structures
            loadProtein(proteinPath);
            loadLigand(ligandPath);
            
            // Update scores table
            updateScoresTable([{
                pose: 1,
                score: data.results?.scores?.[0] || 'N/A'
            }]);
        } else {
            showAlert('Job is not completed yet.', 'warning');
        }
    })
    .catch(error => {
        showAlert('Error loading pose: ' + error, 'danger');
    });
}

// Helper function to show alerts
function showAlert(message, type) {
    // Create alert element
    const alertDiv = document.createElement('div');
    alertDiv.className = `alert alert-${type} alert-dismissible fade show position-fixed top-0 start-50 translate-middle-x mt-3`;
    alertDiv.style.zIndex = '9999';
    alertDiv.innerHTML = `
        ${message}
        <button type="button" class="btn-close" data-bs-dismiss="alert" aria-label="Close"></button>
    `;
    
    // Add to body
    document.body.appendChild(alertDiv);
    
    // Auto-dismiss after 5 seconds
    setTimeout(() => {
        alertDiv.remove();
    }, 5000);
}

// Helper function to load protein (placeholder - can be enhanced)
function loadProtein(proteinPath) {
    console.log('Loading protein from:', proteinPath);
    // This would integrate with the existing loadVisualization function
    // For now, just log the path
}

// Helper function to load ligand (placeholder - can be enhanced)
function loadLigand(ligandPath) {
    console.log('Loading ligand from:', ligandPath);
    // This would integrate with the existing loadVisualization function
    // For now, just log the path
}

// Helper function to update scores table (placeholder - can be enhanced)
function updateScoresTable(scores) {
    console.log('Updating scores table:', scores);
    // This would update the scores table in the UI
    // For now, just log the scores
}
#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Flask API for molecular docking system.
Provides endpoints for docking, ML pose prediction, and result retrieval.
"""

import os
import uuid
import json
import tempfile
import shutil
from pathlib import Path
from typing import Dict, Any, List, Optional, Union

from flask import Flask, request, jsonify, send_file, send_from_directory
from flask_cors import CORS
from werkzeug.utils import secure_filename
import threading
import logging
import redis

# Import docking components
from ..preprocessing.protein_prep import ProteinPreparation
from ..preprocessing.ligand_prep import LigandPreparation
from ..docking.vina_docking import VinaDocking
from ..ml.gnina_rescoring import GNINAScoring
from ..ml.equibind_pose import EquiBindPose
from ..ml.diffdock_pose import DiffDockPose
from .batch import BatchProcessor
from .hdi_routes import hdi_bp

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Initialize Flask app with static file serving from React build directory
frontend_build_dir = Path(__file__).resolve().parents[2] / 'frontend' / 'dist'
app = Flask(__name__, static_folder=str(frontend_build_dir), static_url_path='')

# Enable CORS
CORS(app)

# Register blueprints
app.register_blueprint(hdi_bp)

# Configuration
app.config['UPLOAD_FOLDER'] = os.environ.get('UPLOAD_FOLDER', os.path.join(tempfile.gettempdir(), 'docking_uploads'))
app.config['RESULTS_FOLDER'] = os.environ.get('RESULTS_FOLDER', os.path.join(tempfile.gettempdir(), 'docking_results'))
app.config['BATCH_FOLDER'] = os.environ.get('BATCH_FOLDER', os.path.join(tempfile.gettempdir(), 'docking_batch'))
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024  # 50 MB max upload size

# Create shortcuts for easier access
UPLOAD_FOLDER = app.config['UPLOAD_FOLDER']
RESULTS_FOLDER = app.config['RESULTS_FOLDER']
BATCH_FOLDER = app.config['BATCH_FOLDER']

# Configure Redis
REDIS_HOST = os.environ.get('REDIS_HOST', 'localhost')
REDIS_PORT = int(os.environ.get('REDIS_PORT', 6379))
redis_client = redis.Redis(host=REDIS_HOST, port=REDIS_PORT)
DOCKING_QUEUE = 'docking_tasks'
BATCH_QUEUE = 'batch_tasks'
RESULT_KEY_PREFIX = 'result:'

# Create necessary directories
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['RESULTS_FOLDER'], exist_ok=True)

# Job storage
jobs = {}
job_lock = threading.Lock()

# Allowed file extensions
ALLOWED_EXTENSIONS = {
    'protein': {'pdb'},
    'ligand': {'sdf', 'mol2', 'pdb', 'pdbqt', 'smi'}
}

def allowed_file(filename: str, file_type: str) -> bool:
    """Check if file has an allowed extension.
    
    Args:
        filename: Name of the file to check
        file_type: Type of file ('protein' or 'ligand')
        
    Returns:
        True if file extension is allowed, False otherwise
    """
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS.get(file_type, set())

def generate_job_id() -> str:
    """Generate a unique job ID.
    
    Returns:
        Unique job ID string
    """
    return str(uuid.uuid4())

def update_job_status(job_id: str, status: str, message: Optional[str] = None, 
                     results: Optional[Dict[str, Any]] = None) -> None:
    """Update job status in the jobs dictionary.
    
    Args:
        job_id: ID of the job to update
        status: New status ('pending', 'running', 'completed', 'failed')
        message: Optional status message
        results: Optional results dictionary
    """
    with job_lock:
        if job_id in jobs:
            jobs[job_id]['status'] = status
            if message:
                jobs[job_id]['message'] = message
            if results:
                jobs[job_id]['results'] = results

def run_docking_job(job_id: str, job_data: Dict[str, Any]) -> None:
    """Run docking job in a separate thread.
    
    Args:
        job_id: ID of the job to run
        job_data: Job data dictionary
    """
    try:
        update_job_status(job_id, 'running', 'Processing files')
        
        # Create job directory
        job_dir = os.path.join(app.config['RESULTS_FOLDER'], job_id)
        os.makedirs(job_dir, exist_ok=True)
        
        # Get job parameters
        protein_path = job_data['protein_path']
        ligand_path = job_data['ligand_path']
        job_type = job_data.get('job_type', 'vina')
        params = job_data.get('params', {})
        
        # Prepare output paths
        prepared_protein_path = os.path.join(job_dir, 'protein_prepared.pdbqt')
        prepared_ligand_path = os.path.join(job_dir, 'ligand_prepared.pdbqt')
        output_path = os.path.join(job_dir, 'docking_results.pdbqt')
        
        # Step 1: Prepare protein
        update_job_status(job_id, 'running', 'Preparing protein')
        protein_prep = ProteinPreparation()
        prepared_protein_path = protein_prep.prepare_protein(
            protein_path, 
            output_path=prepared_protein_path,
            add_hydrogens=params.get('add_hydrogens', True),
            assign_charges=params.get('assign_charges', True)
        )
        
        # Step 2: Prepare ligand
        update_job_status(job_id, 'running', 'Preparing ligand')
        ligand_prep = LigandPreparation()
        
        # Check if input is SMILES
        if ligand_path.endswith('.smi'):
            with open(ligand_path, 'r') as f:
                smiles = f.read().strip()
            prepared_ligand_path = ligand_prep.prepare_ligand_from_smiles(
                smiles, 
                output_path=prepared_ligand_path
            )
        else:
            prepared_ligand_path = ligand_prep.prepare_ligand_from_file(
                ligand_path, 
                output_path=prepared_ligand_path
            )
        
        results = {}
        
        # Step 3: Run appropriate docking/prediction method
        if job_type == 'vina':
            update_job_status(job_id, 'running', 'Running Vina docking')
            vina = VinaDocking()
            
            # Get box parameters
            box_params = params.get('box', {})
            if not box_params:
                # Auto-detect binding site
                binding_site = vina.detect_binding_site(prepared_protein_path)
                box_center = [
                    binding_site['center_x'],
                    binding_site['center_y'],
                    binding_site['center_z']
                ]
                box_size = [
                    binding_site['size_x'],
                    binding_site['size_y'],
                    binding_site['size_z']
                ]
            else:
                box_center = [
                    box_params.get('center_x', 0),
                    box_params.get('center_y', 0),
                    box_params.get('center_z', 0)
                ]
                box_size = [
                    box_params.get('size_x', 20),
                    box_params.get('size_y', 20),
                    box_params.get('size_z', 20)
                ]
            
            # Run docking
            binding_site_params = {
                'center_x': box_center[0],
                'center_y': box_center[1],
                'center_z': box_center[2],
                'size_x': box_size[0],
                'size_y': box_size[1],
                'size_z': box_size[2]
            }
            output_path = vina.run_docking(
                prepared_protein_path,
                prepared_ligand_path,
                output_path=output_path,
                binding_site=binding_site_params,
                exhaustiveness=params.get('exhaustiveness', 8),
                num_modes=params.get('num_modes', 9)
            )
            
            # Extract poses and convert top pose to PDB for visualization
            pose_files = vina.extract_poses(output_path, os.path.join(job_dir, 'poses'))
            top_pose_path = os.path.join(job_dir, 'top_pose.pdb')
            if pose_files:
                vina.convert_pose_to_pdb(pose_files[0], top_pose_path)
            
            # Load scores from JSON file created by run_docking
            scores_file = output_path.replace('.pdbqt', '.json')
            scores = {}
            if os.path.exists(scores_file):
                import json
                with open(scores_file, 'r') as f:
                    scores = json.load(f)
            
            results = {
                'output_path': output_path,
                'top_pose_path': top_pose_path,
                'scores': scores,
                'num_poses': len(pose_files)
            }
            
            # Optional: Run GNINA rescoring
            if params.get('run_rescoring', False):
                update_job_status(job_id, 'running', 'Running GNINA rescoring')
                gnina = GNINAScoring()
                rescoring_results = gnina.rescore_poses(
                    prepared_protein_path,
                    output_path
                )
                results['rescoring'] = rescoring_results
        
        elif job_type == 'equibind':
            update_job_status(job_id, 'running', 'Running EquiBind pose prediction')
            equibind = EquiBindPose()
            
            # Convert PDBQT to PDB for EquiBind
            protein_pdb = os.path.join(job_dir, 'protein.pdb')
            ligand_sdf = os.path.join(job_dir, 'ligand.sdf')
            
            # Convert files if needed
            if protein_path.endswith('.pdbqt'):
                cmd = ["obabel", protein_path, "-O", protein_pdb]
                os.system(" ".join(cmd))
            else:
                shutil.copy(protein_path, protein_pdb)
                
            if ligand_path.endswith('.pdbqt'):
                cmd = ["obabel", ligand_path, "-O", ligand_sdf]
                os.system(" ".join(cmd))
            else:
                shutil.copy(ligand_path, ligand_sdf)
            
            # Run EquiBind
            output_path = equibind.predict_pose(
                protein_pdb,
                ligand_sdf,
                os.path.join(job_dir, 'equibind_pose.pdbqt')
            )
            
            # Convert to PDB for visualization
            top_pose_path = os.path.join(job_dir, 'equibind_pose.pdb')
            cmd = ["obabel", output_path, "-O", top_pose_path]
            os.system(" ".join(cmd))
            
            results = {
                'output_path': output_path,
                'top_pose_path': top_pose_path
            }
            
        elif job_type == 'diffdock':
            update_job_status(job_id, 'running', 'Running DiffDock pose prediction')
            diffdock = DiffDockPose({
                'num_samples': params.get('num_samples', 10),
                'confidence_model': params.get('confidence_model', True)
            })
            
            # Convert PDBQT to PDB for DiffDock
            protein_pdb = os.path.join(job_dir, 'protein.pdb')
            ligand_sdf = os.path.join(job_dir, 'ligand.sdf')
            
            # Convert files if needed
            if protein_path.endswith('.pdbqt'):
                cmd = ["obabel", protein_path, "-O", protein_pdb]
                os.system(" ".join(cmd))
            else:
                shutil.copy(protein_path, protein_pdb)
                
            if ligand_path.endswith('.pdbqt'):
                cmd = ["obabel", ligand_path, "-O", ligand_sdf]
                os.system(" ".join(cmd))
            else:
                shutil.copy(ligand_path, ligand_sdf)
            
            # Run DiffDock
            output_paths, confidence_scores = diffdock.predict_pose(
                protein_pdb,
                ligand_sdf,
                os.path.join(job_dir, 'diffdock_output')
            )
            
            # Convert top pose to PDB for visualization
            top_pose_path = os.path.join(job_dir, 'diffdock_top_pose.pdb')
            if output_paths:
                cmd = ["obabel", output_paths[0], "-O", top_pose_path]
                os.system(" ".join(cmd))
            
            results = {
                'output_paths': output_paths,
                'top_pose_path': top_pose_path if output_paths else None,
                'confidence_scores': confidence_scores
            }
        
        # Update job status to completed
        update_job_status(job_id, 'completed', 'Job completed successfully', results)
        
    except Exception as e:
        logger.error(f"Error in job {job_id}: {str(e)}")
        update_job_status(job_id, 'failed', f"Error: {str(e)}")

@app.route('/')
def index():
    """Serve the main UI page.
    
    Returns:
        HTML page
    """
    return send_from_directory(app.static_folder, 'index.html')

@app.route('/api/health', methods=['GET'])
def health_check():
    """Health check endpoint.
    
    Returns:
        JSON response with status
    """
    return jsonify({'status': 'ok'})

@app.route('/api/submit', methods=['POST'])
def submit_job():
    """Submit a new docking job.
    
    Returns:
        JSON response with job ID and status
    """
    try:
        # Check if protein file is provided
        if 'protein' not in request.files:
            logger.error('No protein file in request')
            return jsonify({'error': 'No protein file provided'}), 400
        
        protein_file = request.files['protein']
        if protein_file.filename == '':
            return jsonify({'error': 'No protein file selected'}), 400
        
        if not allowed_file(protein_file.filename, 'protein'):
            return jsonify({'error': f'Invalid protein file format. Allowed formats: {", ".join(ALLOWED_EXTENSIONS["protein"])}'}), 400
        
        # Check if ligand file or SMILES is provided
        ligand_file = None
        smiles = request.form.get('smiles')
        
        if 'ligand' in request.files:
            ligand_file = request.files['ligand']
            if ligand_file.filename != '' and not allowed_file(ligand_file.filename, 'ligand'):
                return jsonify({'error': f'Invalid ligand file format. Allowed formats: {", ".join(ALLOWED_EXTENSIONS["ligand"])}'}), 400
        
        if not ligand_file and not smiles:
            return jsonify({'error': 'No ligand file or SMILES provided'}), 400
        
        # Generate job ID
        job_id = generate_job_id()
        
        # Create job directory
        job_dir = os.path.join(app.config['UPLOAD_FOLDER'], job_id)
        os.makedirs(job_dir, exist_ok=True)
        
        # Save protein file
        protein_filename = secure_filename(protein_file.filename)
        protein_path = os.path.join(job_dir, protein_filename)
        protein_file.save(protein_path)
        
        # Save ligand file or SMILES
        if ligand_file and ligand_file.filename != '':
            ligand_filename = secure_filename(ligand_file.filename)
            ligand_path = os.path.join(job_dir, ligand_filename)
            ligand_file.save(ligand_path)
        elif smiles:
            ligand_path = os.path.join(job_dir, 'ligand.smi')
            with open(ligand_path, 'w') as f:
                f.write(smiles)
        
        # Get job parameters
        job_type = request.form.get('job_type', 'vina')
        params = {}
        
        # Parse box parameters if provided (skip if empty strings)
        box_fields = ['center_x', 'center_y', 'center_z', 'size_x', 'size_y', 'size_z']
        if all(f'box_{p}' in request.form and request.form.get(f'box_{p}', '').strip() for p in box_fields):
            params['box'] = {
                'center_x': float(request.form.get('box_center_x')),
                'center_y': float(request.form.get('box_center_y')),
                'center_z': float(request.form.get('box_center_z')),
                'size_x': float(request.form.get('box_size_x', 20)),
                'size_y': float(request.form.get('box_size_y', 20)),
                'size_z': float(request.form.get('box_size_z', 20))
            }
        
        # Parse other parameters
        for param in ['exhaustiveness', 'num_modes', 'run_rescoring', 'num_samples', 'confidence_model']:
            if param in request.form:
                if param in ['exhaustiveness', 'num_modes', 'num_samples']:
                    params[param] = int(request.form.get(param))
                elif param in ['run_rescoring', 'confidence_model']:
                    params[param] = request.form.get(param).lower() == 'true'
        
        # Create job data
        job_data = {
            'protein_path': protein_path,
            'ligand_path': ligand_path,
            'job_type': job_type,
            'params': params,
            'status': 'pending',
            'message': 'Job submitted'
        }
        
        # Store job data
        with job_lock:
            jobs[job_id] = job_data
        
        # Create task data for Redis queue (align with worker schema)
        output_dir = os.path.join(app.config['RESULTS_FOLDER'], job_id)
        os.makedirs(output_dir, exist_ok=True)

        task_data = {
            'job_id': job_id,
            'method': job_type,
            'protein_file': protein_path,
            'output_dir': output_dir,
            'exhaustiveness': params.get('exhaustiveness', 8),
            'num_modes': params.get('num_modes', 9),
        }

        # Provide ligand as file or SMILES
        if 'ligand' in request.files and request.files['ligand'].filename != '':
            task_data['ligand_file'] = ligand_path
        elif smiles:
            task_data['ligand_smiles'] = smiles

        # Flatten box parameters if provided
        box_params = params.get('box', {})
        if box_params:
            task_data.update({
                'center_x': box_params.get('center_x'),
                'center_y': box_params.get('center_y'),
                'center_z': box_params.get('center_z'),
                'size_x': box_params.get('size_x', 20),
                'size_y': box_params.get('size_y', 20),
                'size_z': box_params.get('size_z', 20)
            })
        
        # Check if Redis is available
        try:
            redis_client.ping()
            # Add task to Redis queue
            redis_client.rpush(DOCKING_QUEUE, json.dumps(task_data))
            logger.info(f"Job {job_id} added to Redis queue")
        except redis.exceptions.ConnectionError:
            logger.warning("Redis not available, falling back to local processing")
            # Start job in a separate thread (fallback)
            threading.Thread(target=run_docking_job, args=(job_id, job_data)).start()
        
        # Return job ID
        return jsonify({
            'job_id': job_id,
            'status': 'pending',
            'message': 'Job submitted successfully'
        })
    
    except Exception as e:
        logger.error(f'Error in submit_job: {str(e)}', exc_info=True)
        return jsonify({'error': f'Server error: {str(e)}'}), 500

@app.route('/api/status/<job_id>', methods=['GET'])
def job_status(job_id):
    """Get job status.
    
    Args:
        job_id: ID of the job to check
        
    Returns:
        JSON response with job status
    """
    # First check Redis for results (worker stores them there)
    try:
        redis_client.ping()
        result_key = f"{RESULT_KEY_PREFIX}{job_id}"
        redis_result = redis_client.get(result_key)
        
        if redis_result:
            # Job completed by worker, results in Redis
            result_data = json.loads(redis_result.decode('utf-8'))
            
            response = {
                'job_id': job_id,
                'status': result_data.get('status', 'completed'),
                'message': 'Job completed successfully' if result_data.get('status') == 'completed' else result_data.get('error', '')
            }
            
            # Include results if available
            if result_data.get('status') == 'completed' and 'results' in result_data:
                results = result_data['results']
                response['results'] = {
                    'num_poses': len(results.get('scores', [])),
                    'scores': results.get('scores', [])
                }
            
            return jsonify(response)
    except redis.exceptions.ConnectionError:
        logger.warning("Redis not available, checking local job storage")
    
    # Fallback to in-memory jobs dictionary
    with job_lock:
        if job_id not in jobs:
            return jsonify({'error': 'Job not found'}), 404
        
        job = jobs[job_id]
        
        response = {
            'job_id': job_id,
            'status': job.get('status', 'unknown'),
            'message': job.get('message', '')
        }
        
        # Include results if job is completed
        if job.get('status') == 'completed' and 'results' in job:
            response['results'] = {
                'num_poses': job['results'].get('num_poses', 1),
                'scores': job['results'].get('scores', [])
            }
        
        return jsonify(response)

@app.route('/api/result/<job_id>/<file_type>', methods=['GET'])
def get_result_file(job_id, file_type):
    """Get result file.
    
    Args:
        job_id: ID of the job
        file_type: Type of file to retrieve ('pdbqt', 'pdb', 'sdf')
        
    Returns:
        File download
    """
    with job_lock:
        if job_id not in jobs:
            return jsonify({'error': 'Job not found'}), 404
        
        job = jobs[job_id]
        
        if job.get('status') != 'completed':
            return jsonify({'error': 'Job not completed'}), 400
        
        if 'results' not in job:
            return jsonify({'error': 'No results available'}), 400
        
        results = job['results']
        
        if file_type == 'pdbqt':
            if 'output_path' in results:
                return send_file(results['output_path'], as_attachment=True)
            elif 'output_paths' in results and results['output_paths']:
                return send_file(results['output_paths'][0], as_attachment=True)
        elif file_type == 'pdb':
            if 'top_pose_path' in results:
                return send_file(results['top_pose_path'], as_attachment=True)
        
        return jsonify({'error': f'No {file_type} file available'}), 404

@app.route('/api/download/<job_id>/<path:filename>', methods=['GET'])
def download_file(job_id, filename):
    """Download a specific file from job results.
    
    Args:
        job_id: ID of the job
        filename: Name of the file to download
        
    Returns:
        File download
    """
    import os
    
    # Construct the file path
    results_dir = os.path.join(RESULTS_FOLDER, job_id)
    file_path = os.path.join(results_dir, filename)
    
    # Security check: ensure the file is within the results directory
    if not os.path.abspath(file_path).startswith(os.path.abspath(results_dir)):
        return jsonify({'error': 'Invalid file path'}), 403
    
    # Check if file exists
    if not os.path.exists(file_path):
        logger.warning(f"File not found: {file_path}")
        return jsonify({'error': f'File not found: {filename}'}), 404
    
    try:
        return send_file(file_path, as_attachment=False, mimetype='chemical/x-pdb')
    except Exception as e:
        logger.error(f"Error sending file {file_path}: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/jobs', methods=['GET'])
def list_jobs():
    """List all jobs.
    
    Returns:
        JSON response with list of jobs
    """
    job_list = []
    
    # Get jobs from in-memory storage
    with job_lock:
        for job_id, job in jobs.items():
            # Check Redis for updated status
            job_status_data = {
                'job_id': job_id,
                'status': job.get('status', 'pending'),
                'message': job.get('message', ''),
                'job_type': job.get('job_type', 'vina')
            }
            
            # Try to get updated status from Redis
            try:
                redis_client.ping()
                result_key = f"{RESULT_KEY_PREFIX}{job_id}"
                redis_result = redis_client.get(result_key)
                
                if redis_result:
                    result_data = json.loads(redis_result.decode('utf-8'))
                    job_status_data['status'] = result_data.get('status', 'completed')
                    if result_data.get('status') == 'completed':
                        job_status_data['message'] = 'Job completed successfully'
                    elif result_data.get('status') == 'failed':
                        job_status_data['message'] = result_data.get('error', 'Job failed')
            except (redis.exceptions.ConnectionError, Exception) as e:
                logger.debug(f"Could not check Redis for job {job_id}: {e}")
            
            job_list.append(job_status_data)
    
    return jsonify({'jobs': job_list})

@app.route('/', defaults={'path': ''})
@app.route('/<path:path>')
def serve_react_app(path):
    """Serve React app for all non-API routes."""
    # If it's an API route, let it be handled by the API endpoints
    if path.startswith('api/'):
        return jsonify({'error': 'API endpoint not found'}), 404
    
    # For all other routes, serve the React app
    try:
        if path != "" and os.path.exists(os.path.join(app.static_folder, path)):
            return send_from_directory(app.static_folder, path)
        else:
            return send_from_directory(app.static_folder, 'index.html')
    except Exception as e:
        logger.error(f"Error serving React app: {e}")
        return jsonify({'error': 'Frontend not available'}), 500

def main():
    """Run the Flask application."""
    app.run(host='0.0.0.0', port=5000, debug=True)

if __name__ == '__main__':
    main()
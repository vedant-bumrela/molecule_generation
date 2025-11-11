#!/usr/bin/env python3
"""
Worker script for processing molecular docking tasks from Redis queue.
This enables distributed high-throughput screening (HTS) capabilities.
"""

import os
import sys
import json
import time
import redis
import logging
from pathlib import Path

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.preprocessing.protein_prep import ProteinPreparation
from src.preprocessing.ligand_prep import LigandPreparation
from src.docking.vina_docking import VinaDocking
from src.ml.gnina_rescoring import GNINAScoring
from src.ml.equibind_pose import EquiBindPose
from src.ml.diffdock_pose import DiffDockPose
from src.api.batch import BatchProcessor

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('worker.log')
    ]
)
logger = logging.getLogger('docking_worker')

# Redis connection
REDIS_HOST = os.environ.get('REDIS_HOST', 'localhost')
REDIS_PORT = int(os.environ.get('REDIS_PORT', 6379))
redis_client = redis.Redis(host=REDIS_HOST, port=REDIS_PORT)

# Queue names
DOCKING_QUEUE = 'docking_tasks'
BATCH_QUEUE = 'batch_tasks'
RESULT_KEY_PREFIX = 'result:'

def process_docking_task(task_data):
    """Process a single docking task from the queue"""
    try:
        job_id = task_data.get('job_id')
        logger.info(f"Processing docking task {job_id}")
        
        # Extract parameters
        method = task_data.get('method', 'vina')
        protein_file = task_data.get('protein_file')
        ligand_file = task_data.get('ligand_file')
        ligand_smiles = task_data.get('ligand_smiles')
        output_dir = task_data.get('output_dir')
        
        # Create output directory if it doesn't exist
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Prepare protein
        protein_prep = ProteinPreparation()
        protein_output = os.path.join(output_dir, 'protein_prepared.pdbqt')
        prepared_protein = protein_prep.prepare_protein(
            protein_file, 
            output_path=protein_output
        )
        
        # Prepare ligand
        ligand_prep = LigandPreparation()
        ligand_output = os.path.join(output_dir, 'ligand_prepared.pdbqt')
        if ligand_file:
            prepared_ligand = ligand_prep.prepare_ligand_from_file(
                ligand_file, 
                output_path=ligand_output
            )
        elif ligand_smiles:
            prepared_ligand = ligand_prep.prepare_ligand_from_smiles(
                ligand_smiles, 
                output_path=ligand_output
            )
        else:
            raise ValueError("Either ligand_file or ligand_smiles must be provided")
        
        # Perform docking based on method
        results = {}
        
        if method == 'vina':
            # AutoDock Vina docking
            vina = VinaDocking()
            
            # Prepare binding site parameters
            binding_site = None
            if all(task_data.get(k) is not None for k in ['center_x', 'center_y', 'center_z']):
                binding_site = {
                    'center_x': task_data.get('center_x'),
                    'center_y': task_data.get('center_y'),
                    'center_z': task_data.get('center_z'),
                    'size_x': task_data.get('size_x', 20),
                    'size_y': task_data.get('size_y', 20),
                    'size_z': task_data.get('size_z', 20)
                }
            
            # Run docking
            output_pdbqt = os.path.join(output_dir, 'docked.pdbqt')
            result_path = vina.run_docking(
                prepared_protein,
                prepared_ligand,
                output_path=output_pdbqt,
                binding_site=binding_site,
                exhaustiveness=task_data.get('exhaustiveness', 8),
                num_modes=task_data.get('num_modes', 9)
            )
            
            # Load scores from JSON file
            scores_path = Path(result_path).with_suffix('.json')
            with open(scores_path, 'r') as f:
                scores_data = json.load(f)
            
            # Extract affinity scores from modes
            scores = []
            if 'modes' in scores_data:
                scores = [mode['affinity'] for mode in scores_data['modes']]
            
            results = {
                'output_file': result_path,
                'scores': scores,
                'best_score': scores_data.get('best_affinity', scores[0] if scores else None)
            }
            
        elif method == 'gnina':
            # GNINA ML-based rescoring
            gnina = GNINAScoring()
            rescoring_results = gnina.rescore(
                prepared_protein,
                prepared_ligand,
                output_dir=output_dir
            )
            results = {
                'cnn_scores': rescoring_results['cnn_scores'],
                'affinities': rescoring_results['affinities'],
                'poses': rescoring_results['poses']
            }
            
        elif method == 'equibind':
            # EquiBind ML-based pose prediction
            equibind = EquiBindPose()
            pose_results = equibind.predict_pose(
                prepared_protein,
                prepared_ligand,
                output_dir=output_dir
            )
            results = {
                'poses': pose_results['poses']
            }
            
        elif method == 'diffdock':
            # DiffDock ML-based pose prediction
            diffdock = DiffDockPose()
            pose_results = diffdock.predict_pose(
                prepared_protein,
                prepared_ligand,
                output_dir=output_dir
            )
            results = {
                'poses': pose_results['poses'],
                'confidence': pose_results.get('confidence', [])
            }
        
        # Store results in Redis
        result_key = f"{RESULT_KEY_PREFIX}{job_id}"
        redis_client.set(result_key, json.dumps({
            'status': 'completed',
            'results': results
        }))
        redis_client.expire(result_key, 86400 * 7)  # Expire after 7 days
        
        logger.info(f"Completed docking task {job_id}")
        return True
        
    except Exception as e:
        logger.error(f"Error processing docking task: {str(e)}", exc_info=True)
        
        # Store error in Redis
        if job_id:
            result_key = f"{RESULT_KEY_PREFIX}{job_id}"
            redis_client.set(result_key, json.dumps({
                'status': 'failed',
                'error': str(e)
            }))
            redis_client.expire(result_key, 86400 * 7)  # Expire after 7 days
        
        return False

def process_batch_task(task_data):
    """Process a batch docking task from the queue"""
    try:
        job_id = task_data.get('job_id')
        logger.info(f"Processing batch task {job_id}")
        
        # Extract parameters
        method = task_data.get('method', 'vina')
        protein_file = task_data.get('protein_file')
        compound_library = task_data.get('compound_library')
        is_smiles = task_data.get('is_smiles', True)
        output_dir = task_data.get('output_dir')
        num_workers = task_data.get('num_workers', 4)
        
        # Create output directory if it doesn't exist
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Initialize batch processor
        batch_processor = BatchProcessor(
            method=method,
            output_dir=output_dir,
            num_workers=num_workers,
            center_x=task_data.get('center_x'),
            center_y=task_data.get('center_y'),
            center_z=task_data.get('center_z'),
            size_x=task_data.get('size_x', 20),
            size_y=task_data.get('size_y', 20),
            size_z=task_data.get('size_z', 20),
            exhaustiveness=task_data.get('exhaustiveness', 8),
            energy_range=task_data.get('energy_range', 3),
            num_modes=task_data.get('num_modes', 9)
        )
        
        # Set up progress callback to update Redis
        def progress_callback(processed, total, current_result=None):
            progress = int((processed / total) * 100) if total > 0 else 0
            redis_client.set(f"progress:{job_id}", json.dumps({
                'processed': processed,
                'total': total,
                'progress': progress,
                'current_result': current_result
            }))
        
        # Process batch
        results = None
        if is_smiles:
            results = batch_processor.process_smiles_library(
                protein_file,
                compound_library,
                progress_callback=progress_callback
            )
        else:
            results = batch_processor.process_sdf_library(
                protein_file,
                compound_library,
                progress_callback=progress_callback
            )
        
        # Store results in Redis
        result_key = f"{RESULT_KEY_PREFIX}{job_id}"
        redis_client.set(result_key, json.dumps({
            'status': 'completed',
            'results': results
        }))
        redis_client.expire(result_key, 86400 * 7)  # Expire after 7 days
        
        logger.info(f"Completed batch task {job_id}")
        return True
        
    except Exception as e:
        logger.error(f"Error processing batch task: {str(e)}", exc_info=True)
        
        # Store error in Redis
        if job_id:
            result_key = f"{RESULT_KEY_PREFIX}{job_id}"
            redis_client.set(result_key, json.dumps({
                'status': 'failed',
                'error': str(e)
            }))
            redis_client.expire(result_key, 86400 * 7)  # Expire after 7 days
        
        return False

def main():
    """Main worker loop to process tasks from Redis queues"""
    logger.info("Starting docking worker")
    
    while True:
        try:
            # Try to get a task from the docking queue first (higher priority)
            task = redis_client.blpop([DOCKING_QUEUE, BATCH_QUEUE], timeout=1)
            
            if task:
                queue_name, task_data = task
                queue_name = queue_name.decode('utf-8')
                task_data = json.loads(task_data.decode('utf-8'))
                
                logger.info(f"Received task from {queue_name}")
                
                if queue_name == DOCKING_QUEUE:
                    process_docking_task(task_data)
                elif queue_name == BATCH_QUEUE:
                    process_batch_task(task_data)
            
            # Small delay to prevent CPU hogging
            time.sleep(0.1)
            
        except KeyboardInterrupt:
            logger.info("Worker shutting down")
            break
        except Exception as e:
            logger.error(f"Error in worker main loop: {str(e)}", exc_info=True)
            time.sleep(5)  # Wait before retrying

if __name__ == "__main__":
    main()
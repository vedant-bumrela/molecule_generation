#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Batch processing module for high-throughput screening (HTS) capabilities.
Provides functionality for processing multiple ligands against a target protein.
"""

import os
import csv
import json
import time
import logging
import threading
import queue
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple

import pandas as pd
from rdkit import Chem

from ..preprocessing.protein_prep import ProteinPreparation
from ..preprocessing.ligand_prep import LigandPreparation
from ..docking.vina_docking import VinaDocking
from ..ml.gnina_rescoring import GNINAScoring
from ..ml.equibind_pose import EquiBindPose
from ..ml.diffdock_pose import DiffDockPose

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class BatchProcessor:
    """Class for batch processing of molecular docking jobs."""
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """Initialize batch processor with configuration.
        
        Args:
            config: Configuration dictionary with parameters for batch processing
        """
        self.config = config or {}
        self.output_dir = self.config.get("output_dir", "./batch_results")
        self.num_workers = self.config.get("num_workers", 4)
        self.job_queue = queue.Queue()
        self.results = {}
        self.workers = []
        self.stop_event = threading.Event()
        
        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)
    
    def prepare_protein(self, protein_path: Union[str, Path], 
                       output_path: Optional[Union[str, Path]] = None) -> str:
        """Prepare protein for docking.
        
        Args:
            protein_path: Path to protein PDB file
            output_path: Path to output PDBQT file (optional)
            
        Returns:
            Path to prepared protein PDBQT file
        """
        logger.info(f"Preparing protein: {protein_path}")
        
        protein_prep = ProteinPreparation()
        prepared_path = protein_prep.prepare_protein(
            protein_path,
            output_path=output_path,
            add_hydrogens=self.config.get("add_hydrogens", True),
            assign_charges=self.config.get("assign_charges", True)
        )
        
        logger.info(f"Protein preparation completed: {prepared_path}")
        return prepared_path
    
    def process_smiles_library(self, smiles_file: Union[str, Path], 
                              protein_path: Union[str, Path],
                              output_dir: Optional[Union[str, Path]] = None,
                              method: str = "vina") -> str:
        """Process a library of SMILES strings for docking.
        
        Args:
            smiles_file: Path to CSV file with SMILES strings (must have 'SMILES' column)
            protein_path: Path to protein PDB file
            output_dir: Directory to save results (optional)
            method: Docking method to use ('vina', 'equibind', 'diffdock')
            
        Returns:
            Path to results CSV file
        """
        if output_dir is None:
            output_dir = os.path.join(self.output_dir, f"smiles_library_{int(time.time())}")
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Prepare protein once
        prepared_protein = self.prepare_protein(
            protein_path,
            os.path.join(output_dir, "protein_prepared.pdbqt")
        )
        
        # Read SMILES file
        try:
            df = pd.read_csv(smiles_file)
            if 'SMILES' not in df.columns:
                raise ValueError("SMILES file must contain a 'SMILES' column")
            
            # Add ID column if not present
            if 'ID' not in df.columns:
                df['ID'] = [f"compound_{i}" for i in range(len(df))]
            
            # Create jobs
            for _, row in df.iterrows():
                job = {
                    "id": row['ID'],
                    "smiles": row['SMILES'],
                    "protein_path": prepared_protein,
                    "output_dir": os.path.join(output_dir, row['ID']),
                    "method": method
                }
                
                # Add any additional properties from the dataframe
                for col in df.columns:
                    if col not in ['SMILES', 'ID']:
                        job[col] = row[col]
                
                self.job_queue.put(job)
            
            # Start workers
            self._start_workers()
            
            # Wait for all jobs to complete
            self.job_queue.join()
            
            # Stop workers
            self._stop_workers()
            
            # Compile results
            results_file = os.path.join(output_dir, "docking_results.csv")
            self._compile_results(results_file)
            
            logger.info(f"Batch processing completed: {results_file}")
            return results_file
            
        except Exception as e:
            logger.error(f"Error processing SMILES library: {e}")
            raise
    
    def process_sdf_library(self, sdf_file: Union[str, Path], 
                           protein_path: Union[str, Path],
                           output_dir: Optional[Union[str, Path]] = None,
                           method: str = "vina") -> str:
        """Process a library of compounds from an SDF file for docking.
        
        Args:
            sdf_file: Path to SDF file with compounds
            protein_path: Path to protein PDB file
            output_dir: Directory to save results (optional)
            method: Docking method to use ('vina', 'equibind', 'diffdock')
            
        Returns:
            Path to results CSV file
        """
        if output_dir is None:
            output_dir = os.path.join(self.output_dir, f"sdf_library_{int(time.time())}")
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Prepare protein once
        prepared_protein = self.prepare_protein(
            protein_path,
            os.path.join(output_dir, "protein_prepared.pdbqt")
        )
        
        # Read SDF file
        try:
            suppl = Chem.SDMolSupplier(str(sdf_file))
            
            # Create jobs
            for i, mol in enumerate(suppl):
                if mol is None:
                    continue
                
                # Get molecule name or create one
                name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"compound_{i}"
                
                # Create job
                job = {
                    "id": name,
                    "mol": mol,
                    "protein_path": prepared_protein,
                    "output_dir": os.path.join(output_dir, name),
                    "method": method
                }
                
                # Add any additional properties from the molecule
                for prop_name in mol.GetPropNames():
                    if prop_name != "_Name":
                        job[prop_name] = mol.GetProp(prop_name)
                
                self.job_queue.put(job)
            
            # Start workers
            self._start_workers()
            
            # Wait for all jobs to complete
            self.job_queue.join()
            
            # Stop workers
            self._stop_workers()
            
            # Compile results
            results_file = os.path.join(output_dir, "docking_results.csv")
            self._compile_results(results_file)
            
            logger.info(f"Batch processing completed: {results_file}")
            return results_file
            
        except Exception as e:
            logger.error(f"Error processing SDF library: {e}")
            raise
    
    def _worker(self):
        """Worker thread for processing jobs."""
        while not self.stop_event.is_set():
            try:
                # Get job with timeout to check stop_event periodically
                job = self.job_queue.get(timeout=1)
                
                try:
                    # Process job
                    result = self._process_job(job)
                    
                    # Store result
                    self.results[job["id"]] = result
                    
                except Exception as e:
                    logger.error(f"Error processing job {job['id']}: {e}")
                    self.results[job["id"]] = {"error": str(e)}
                
                finally:
                    # Mark job as done
                    self.job_queue.task_done()
                    
            except queue.Empty:
                continue
    
    def _start_workers(self):
        """Start worker threads."""
        self.stop_event.clear()
        self.workers = []
        
        for _ in range(self.num_workers):
            worker = threading.Thread(target=self._worker)
            worker.daemon = True
            worker.start()
            self.workers.append(worker)
    
    def _stop_workers(self):
        """Stop worker threads."""
        self.stop_event.set()
        
        for worker in self.workers:
            worker.join()
    
    def _process_job(self, job: Dict[str, Any]) -> Dict[str, Any]:
        """Process a single docking job.
        
        Args:
            job: Job dictionary with parameters
            
        Returns:
            Dictionary with job results
        """
        job_id = job["id"]
        output_dir = job["output_dir"]
        method = job["method"]
        
        os.makedirs(output_dir, exist_ok=True)
        
        logger.info(f"Processing job {job_id} with method {method}")
        
        # Prepare ligand
        ligand_prep = LigandPreparation()
        
        if "smiles" in job:
            # Prepare from SMILES
            prepared_ligand = ligand_prep.prepare_from_smiles(
                job["smiles"],
                output_path=os.path.join(output_dir, "ligand_prepared.pdbqt")
            )
        elif "mol" in job:
            # Prepare from RDKit molecule
            # First save to SDF
            mol_file = os.path.join(output_dir, "ligand.sdf")
            writer = Chem.SDWriter(mol_file)
            writer.write(job["mol"])
            writer.close()
            
            prepared_ligand = ligand_prep.prepare_from_file(
                mol_file,
                output_path=os.path.join(output_dir, "ligand_prepared.pdbqt")
            )
        else:
            raise ValueError("Job must contain either 'smiles' or 'mol'")
        
        result = {
            "id": job_id,
            "method": method,
            "output_dir": output_dir
        }
        
        # Run appropriate docking method
        if method == "vina":
            result.update(self._run_vina(job["protein_path"], prepared_ligand, output_dir))
        elif method == "equibind":
            result.update(self._run_equibind(job["protein_path"], prepared_ligand, output_dir))
        elif method == "diffdock":
            result.update(self._run_diffdock(job["protein_path"], prepared_ligand, output_dir))
        else:
            raise ValueError(f"Unknown method: {method}")
        
        return result
    
    def _run_vina(self, protein_path: str, ligand_path: str, output_dir: str) -> Dict[str, Any]:
        """Run AutoDock Vina docking.
        
        Args:
            protein_path: Path to prepared protein PDBQT file
            ligand_path: Path to prepared ligand PDBQT file
            output_dir: Directory to save results
            
        Returns:
            Dictionary with docking results
        """
        vina = VinaDocking()
        
        # Auto-detect binding site if not specified
        box_params = self.config.get("box", {})
        if not box_params:
            binding_site = vina.detect_binding_site(protein_path)
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
                box_params.get("center_x", 0),
                box_params.get("center_y", 0),
                box_params.get("center_z", 0)
            ]
            box_size = [
                box_params.get("size_x", 20),
                box_params.get("size_y", 20),
                box_params.get("size_z", 20)
            ]
        
        # Run docking
        output_path = os.path.join(output_dir, "docking_results.pdbqt")
        binding_site_params = {
            'center_x': box_center[0],
            'center_y': box_center[1],
            'center_z': box_center[2],
            'size_x': box_size[0],
            'size_y': box_size[1],
            'size_z': box_size[2]
        }
        output_path = vina.run_docking(
            protein_path,
            ligand_path,
            output_path=output_path,
            binding_site=binding_site_params,
            exhaustiveness=self.config.get("exhaustiveness", 8),
            num_modes=self.config.get("num_modes", 9)
        )
        
        # Extract poses and convert top pose to PDB
        pose_files = vina.extract_poses(output_path, os.path.join(output_dir, "poses"))
        top_pose_path = os.path.join(output_dir, "top_pose.pdb")
        if pose_files:
            vina.convert_pose_to_pdb(pose_files[0], top_pose_path)
        
        # Load scores from JSON file created by run_docking
        scores_file = output_path.replace('.pdbqt', '.json')
        scores = []
        if os.path.exists(scores_file):
            import json
            with open(scores_file, 'r') as f:
                score_data = json.load(f)
                scores = [mode['affinity'] for mode in score_data.get('modes', [])]
        
        # Run GNINA rescoring if requested
        rescoring_results = None
        if self.config.get("run_rescoring", False):
            gnina = GNINAScoring()
            rescoring_results = gnina.rescore_poses(protein_path, output_path)
        
        return {
            "output_path": output_path,
            "top_pose_path": top_pose_path,
            "scores": scores,
            "top_score": scores[0] if scores else None,
            "num_poses": len(pose_files),
            "rescoring": rescoring_results
        }
    
    def _run_equibind(self, protein_path: str, ligand_path: str, output_dir: str) -> Dict[str, Any]:
        """Run EquiBind pose prediction.
        
        Args:
            protein_path: Path to prepared protein PDBQT file
            ligand_path: Path to prepared ligand PDBQT file
            output_dir: Directory to save results
            
        Returns:
            Dictionary with pose prediction results
        """
        equibind = EquiBindPose(self.config.get("equibind_config", {}))
        
        # Convert PDBQT to PDB/SDF for EquiBind
        protein_pdb = os.path.join(output_dir, "protein.pdb")
        ligand_sdf = os.path.join(output_dir, "ligand.sdf")
        
        # Convert files using Open Babel
        import subprocess
        
        subprocess.run(["obabel", protein_path, "-O", protein_pdb], check=True)
        subprocess.run(["obabel", ligand_path, "-O", ligand_sdf], check=True)
        
        # Run EquiBind
        output_path = equibind.predict_pose(
            protein_pdb,
            ligand_sdf,
            os.path.join(output_dir, "equibind_pose.pdbqt")
        )
        
        # Convert to PDB for visualization
        top_pose_path = os.path.join(output_dir, "equibind_pose.pdb")
        subprocess.run(["obabel", output_path, "-O", top_pose_path], check=True)
        
        return {
            "output_path": output_path,
            "top_pose_path": top_pose_path
        }
    
    def _run_diffdock(self, protein_path: str, ligand_path: str, output_dir: str) -> Dict[str, Any]:
        """Run DiffDock pose prediction.
        
        Args:
            protein_path: Path to prepared protein PDBQT file
            ligand_path: Path to prepared ligand PDBQT file
            output_dir: Directory to save results
            
        Returns:
            Dictionary with pose prediction results
        """
        diffdock_config = self.config.get("diffdock_config", {})
        diffdock_config.update({
            "num_samples": self.config.get("num_samples", 10),
            "confidence_model": self.config.get("confidence_model", True)
        })
        
        diffdock = DiffDockPose(diffdock_config)
        
        # Convert PDBQT to PDB/SDF for DiffDock
        protein_pdb = os.path.join(output_dir, "protein.pdb")
        ligand_sdf = os.path.join(output_dir, "ligand.sdf")
        
        # Convert files using Open Babel
        import subprocess
        
        subprocess.run(["obabel", protein_path, "-O", protein_pdb], check=True)
        subprocess.run(["obabel", ligand_path, "-O", ligand_sdf], check=True)
        
        # Run DiffDock
        output_paths, confidence_scores = diffdock.predict_pose(
            protein_pdb,
            ligand_sdf,
            os.path.join(output_dir, "diffdock_output")
        )
        
        # Convert top pose to PDB for visualization
        top_pose_path = None
        if output_paths:
            top_pose_path = os.path.join(output_dir, "diffdock_top_pose.pdb")
            subprocess.run(["obabel", output_paths[0], "-O", top_pose_path], check=True)
        
        return {
            "output_paths": output_paths,
            "top_pose_path": top_pose_path,
            "confidence_scores": confidence_scores,
            "top_confidence": list(confidence_scores.values())[0] if confidence_scores else None
        }
    
    def _compile_results(self, output_file: Union[str, Path]) -> None:
        """Compile results from all jobs into a CSV file.
        
        Args:
            output_file: Path to output CSV file
        """
        # Collect all possible fields
        fields = set()
        for result in self.results.values():
            fields.update(result.keys())
        
        # Ensure ID is first
        fields = ["id"] + [f for f in fields if f != "id"]
        
        # Write results to CSV
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            
            for job_id, result in self.results.items():
                # Handle nested dictionaries
                flat_result = {}
                for key, value in result.items():
                    if isinstance(value, dict):
                        for subkey, subvalue in value.items():
                            flat_result[f"{key}_{subkey}"] = subvalue
                    else:
                        flat_result[key] = value
                
                writer.writerow(flat_result)


def main():
    """Command line interface for batch processing."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Batch processing for molecular docking")
    parser.add_argument("--protein", required=True, help="Path to protein PDB file")
    parser.add_argument("--library", required=True, help="Path to compound library (CSV with SMILES column or SDF)")
    parser.add_argument("--output-dir", help="Directory to save results")
    parser.add_argument("--method", choices=["vina", "equibind", "diffdock"], default="vina", help="Docking method")
    parser.add_argument("--num-workers", type=int, default=4, help="Number of worker threads")
    parser.add_argument("--exhaustiveness", type=int, default=8, help="Vina exhaustiveness")
    parser.add_argument("--num-modes", type=int, default=9, help="Number of binding modes to generate")
    parser.add_argument("--run-rescoring", action="store_true", help="Run GNINA rescoring")
    parser.add_argument("--num-samples", type=int, default=10, help="Number of DiffDock samples")
    
    args = parser.parse_args()
    
    # Configure batch processor
    config = {
        "output_dir": args.output_dir,
        "num_workers": args.num_workers,
        "exhaustiveness": args.exhaustiveness,
        "num_modes": args.num_modes,
        "run_rescoring": args.run_rescoring,
        "num_samples": args.num_samples
    }
    
    batch_processor = BatchProcessor(config)
    
    # Process library
    if args.library.lower().endswith('.csv'):
        results_file = batch_processor.process_smiles_library(
            args.library,
            args.protein,
            args.output_dir,
            args.method
        )
    elif args.library.lower().endswith('.sdf'):
        results_file = batch_processor.process_sdf_library(
            args.library,
            args.protein,
            args.output_dir,
            args.method
        )
    else:
        raise ValueError("Library file must be CSV or SDF")
    
    print(f"Batch processing completed. Results saved to: {results_file}")


if __name__ == "__main__":
    main()
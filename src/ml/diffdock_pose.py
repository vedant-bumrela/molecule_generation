#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DiffDock ML pose prediction module for molecular docking.
Handles generative pose prediction using DiffDock's diffusion-based model.
"""

import os
import logging
import subprocess
import json
import tempfile
from pathlib import Path
import shutil
from typing import Optional, Union, Dict, Any, List, Tuple

import numpy as np
import torch

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class DiffDockPose:
    """Class for generative pose prediction using DiffDock."""
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """Initialize DiffDock pose prediction with configuration.
        
        Args:
            config: Configuration dictionary with parameters for pose prediction
        """
        self.config = config or {}
        self.diffdock_path = self.config.get("diffdock_path", "./DiffDock")
        self.model_path = self.config.get("model_path", os.path.join(self.diffdock_path, "workdir", "pretrained"))
        self.device = self.config.get("device", "cuda" if torch.cuda.is_available() else "cpu")
        self.num_samples = self.config.get("num_samples", 10)
        self.confidence_model = self.config.get("confidence_model", True)
        
    def setup_diffdock(self):
        """Set up DiffDock environment if not already set up."""
        if not os.path.exists(self.diffdock_path):
            logger.info("Setting up DiffDock environment")
            
            # Clone DiffDock repository
            cmd = [
                "git", "clone", "https://github.com/gcorso/DiffDock.git", self.diffdock_path
            ]
            
            try:
                subprocess.run(cmd, check=True, capture_output=True)
                logger.info(f"DiffDock repository cloned to {self.diffdock_path}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to clone DiffDock repository: {e.stderr.decode()}")
                raise
            
            # Create model directory
            os.makedirs(os.path.join(self.diffdock_path, "workdir", "pretrained"), exist_ok=True)
            
            # Note: In a real implementation, you would download the pretrained models here
            # For this example, we'll assume they're already available or will be downloaded manually
            logger.info("Please download pretrained models manually from the DiffDock repository")
        
        logger.info("DiffDock environment is set up")
    
    def prepare_input(self, protein_path: Union[str, Path], ligand_path: Union[str, Path],
                     output_dir: Optional[Union[str, Path]] = None) -> str:
        """Prepare input files for DiffDock.
        
        Args:
            protein_path: Path to protein PDB file
            ligand_path: Path to ligand file (SDF, MOL2, etc.)
            output_dir: Directory to save prepared input files (optional)
            
        Returns:
            Path to directory containing prepared input files
        """
        protein_path = Path(protein_path)
        ligand_path = Path(ligand_path)
        
        if output_dir is None:
            output_dir = Path(tempfile.mkdtemp(prefix="diffdock_input_"))
        else:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Preparing input files for DiffDock in {output_dir}")
        
        # Create data directory structure
        data_dir = output_dir / "data"
        data_dir.mkdir(exist_ok=True)
        
        # Copy protein and ligand files
        protein_dest = data_dir / f"{protein_path.stem}.pdb"
        ligand_dest = data_dir / f"{ligand_path.stem}.sdf"
        
        # Convert ligand to SDF if needed
        if ligand_path.suffix.lower() != ".sdf":
            cmd = [
                "obabel", 
                str(ligand_path), 
                "-O", str(ligand_dest)
            ]
            
            try:
                subprocess.run(cmd, check=True, capture_output=True)
                logger.info(f"Converted ligand to SDF format: {ligand_dest}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to convert ligand to SDF: {e.stderr.decode()}")
                raise
        else:
            shutil.copy(ligand_path, ligand_dest)
        
        # Copy protein file
        shutil.copy(protein_path, protein_dest)
        
        # Create dataset file for DiffDock
        dataset_file = output_dir / "input_data.csv"
        with open(dataset_file, 'w') as f:
            f.write("protein_path,ligand_path,complex_name\n")
            f.write(f"{protein_dest},{ligand_dest},{protein_path.stem}_{ligand_path.stem}\n")
        
        logger.info(f"Input files prepared for DiffDock: {output_dir}")
        return str(output_dir)
    
    def run_prediction(self, input_dir: Union[str, Path], output_dir: Optional[Union[str, Path]] = None) -> str:
        """Run DiffDock pose prediction.
        
        Args:
            input_dir: Directory containing prepared input files
            output_dir: Directory to save prediction results (optional)
            
        Returns:
            Path to directory containing prediction results
        """
        input_dir = Path(input_dir)
        
        if output_dir is None:
            output_dir = input_dir / "results"
        else:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Running DiffDock pose prediction with input from {input_dir}")
        
        # Change to DiffDock directory
        original_dir = os.getcwd()
        os.chdir(self.diffdock_path)
        
        try:
            # Run DiffDock inference
            cmd = [
                "python", "inference.py", 
                "--protein_path", str(input_dir / "data" / f"{Path(input_dir).stem}.pdb"),
                "--ligand_path", str(input_dir / "data" / f"{Path(input_dir).stem}.sdf"),
                "--out_dir", str(output_dir),
                "--samples_per_complex", str(self.num_samples),
                "--model_dir", str(self.model_path),
                "--device", self.device
            ]
            
            if self.confidence_model:
                cmd.append("--confidence_model")
            
            subprocess.run(cmd, check=True, capture_output=True)
            logger.info(f"DiffDock pose prediction completed successfully: {output_dir}")
            
            # Return to original directory
            os.chdir(original_dir)
            return str(output_dir)
        
        except subprocess.CalledProcessError as e:
            logger.error(f"DiffDock pose prediction failed: {e.stderr.decode()}")
            # Return to original directory
            os.chdir(original_dir)
            raise
        
        except Exception as e:
            logger.error(f"Error during DiffDock pose prediction: {e}")
            # Return to original directory
            os.chdir(original_dir)
            raise
    
    def convert_output_to_pdbqt(self, output_dir: Union[str, Path], 
                               output_path: Optional[Union[str, Path]] = None) -> List[str]:
        """Convert DiffDock output to PDBQT format for docking.
        
        Args:
            output_dir: Directory containing DiffDock prediction results
            output_path: Base path for output PDBQT files (optional)
            
        Returns:
            List of paths to PDBQT files
        """
        output_dir = Path(output_dir)
        
        # Find predicted ligand files
        ligand_files = list(output_dir.glob("*_ligand_*.sdf"))
        
        if not ligand_files:
            raise ValueError(f"No predicted ligand files found in {output_dir}")
        
        # Sort by confidence score if available
        ligand_files.sort(key=lambda x: int(x.stem.split('_')[-1]))
        
        pdbqt_files = []
        for i, ligand_file in enumerate(ligand_files):
            if output_path is None:
                pdbqt_file = ligand_file.with_suffix(".pdbqt")
            else:
                pdbqt_file = Path(f"{output_path}_{i}.pdbqt")
            
            logger.info(f"Converting DiffDock output {ligand_file} to PDBQT format")
            
            # Use Open Babel to convert to PDBQT
            cmd = [
                "obabel", 
                str(ligand_file), 
                "-O", str(pdbqt_file), 
                "-xh"  # Add AutoDock atom types and optimize hydrogens
            ]
            
            try:
                subprocess.run(cmd, check=True, capture_output=True)
                logger.info(f"Converted DiffDock output to PDBQT format: {pdbqt_file}")
                pdbqt_files.append(str(pdbqt_file))
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to convert to PDBQT: {e.stderr.decode()}")
                raise
        
        return pdbqt_files
    
    def extract_confidence_scores(self, output_dir: Union[str, Path]) -> Dict[str, float]:
        """Extract confidence scores from DiffDock output.
        
        Args:
            output_dir: Directory containing DiffDock prediction results
            
        Returns:
            Dictionary mapping pose file paths to confidence scores
        """
        output_dir = Path(output_dir)
        confidence_file = output_dir / "confidence_scores.txt"
        
        if not confidence_file.exists():
            logger.warning(f"Confidence scores file not found: {confidence_file}")
            return {}
        
        confidence_scores = {}
        try:
            with open(confidence_file, 'r') as f:
                for line in f:
                    if line.strip():
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            file_path = parts[0]
                            score = float(parts[1])
                            confidence_scores[file_path] = score
            
            logger.info(f"Extracted {len(confidence_scores)} confidence scores")
            return confidence_scores
        
        except Exception as e:
            logger.error(f"Failed to extract confidence scores: {e}")
            return {}
    
    def predict_pose(self, protein_path: Union[str, Path], ligand_path: Union[str, Path],
                    output_dir: Optional[Union[str, Path]] = None) -> Tuple[List[str], Dict[str, float]]:
        """Full DiffDock pose prediction pipeline.
        
        Args:
            protein_path: Path to protein PDB file
            ligand_path: Path to ligand file
            output_dir: Directory to save output files (optional)
            
        Returns:
            Tuple of (list of paths to predicted pose PDBQT files, confidence scores)
        """
        protein_path = Path(protein_path)
        ligand_path = Path(ligand_path)
        
        if output_dir is None:
            output_dir = ligand_path.parent / f"{ligand_path.stem}_diffdock"
            output_dir.mkdir(parents=True, exist_ok=True)
        else:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Starting full DiffDock pose prediction pipeline")
        
        # Step 1: Set up DiffDock
        self.setup_diffdock()
        
        # Step 2: Prepare input files
        input_dir = self.prepare_input(protein_path, ligand_path)
        
        # Step 3: Run prediction
        results_dir = self.run_prediction(input_dir)
        
        # Step 4: Convert output to PDBQT
        output_base = output_dir / f"{ligand_path.stem}_diffdock"
        pdbqt_paths = self.convert_output_to_pdbqt(results_dir, output_base)
        
        # Step 5: Extract confidence scores
        confidence_scores = self.extract_confidence_scores(results_dir)
        
        # Clean up temporary files if needed
        if self.config.get("cleanup_intermediates", True):
            if Path(input_dir).exists() and Path(input_dir).is_dir() and "diffdock_input_" in input_dir:
                shutil.rmtree(input_dir)
        
        logger.info(f"Completed DiffDock pose prediction: {len(pdbqt_paths)} poses")
        return pdbqt_paths, confidence_scores


def main():
    """Command line interface for DiffDock pose prediction."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Predict ligand poses with DiffDock")
    parser.add_argument("--protein", required=True, help="Path to protein PDB file")
    parser.add_argument("--ligand", required=True, help="Path to ligand file")
    parser.add_argument("--out-dir", dest="output_dir", help="Directory to save output files")
    parser.add_argument("--diffdock-path", help="Path to DiffDock repository")
    parser.add_argument("--model-path", help="Path to pretrained DiffDock model")
    parser.add_argument("--device", choices=["cpu", "cuda"], help="Device to use for inference")
    parser.add_argument("--num-samples", type=int, default=10, help="Number of poses to generate")
    parser.add_argument("--no-confidence", action="store_true", help="Disable confidence model")
    parser.add_argument("--keep-intermediates", action="store_true", help="Keep intermediate files")
    
    args = parser.parse_args()
    
    # Configure pose prediction
    config = {
        "cleanup_intermediates": not args.keep_intermediates,
        "num_samples": args.num_samples,
        "confidence_model": not args.no_confidence
    }
    if args.diffdock_path:
        config["diffdock_path"] = args.diffdock_path
    if args.model_path:
        config["model_path"] = args.model_path
    if args.device:
        config["device"] = args.device
    
    # Run pose prediction
    diffdock = DiffDockPose(config)
    pdbqt_paths, confidence_scores = diffdock.predict_pose(
        args.protein,
        args.ligand,
        args.output_dir
    )
    
    print(f"Generated {len(pdbqt_paths)} poses:")
    for i, path in enumerate(pdbqt_paths):
        confidence = confidence_scores.get(Path(path).stem, "N/A")
        print(f"  Pose {i+1}: {path} (Confidence: {confidence})")


if __name__ == "__main__":
    main()
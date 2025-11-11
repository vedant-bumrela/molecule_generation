#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
EquiBind ML pose prediction module for molecular docking.
Handles direct pose prediction using EquiBind's SE(3)-equivariant GNN.
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

class EquiBindPose:
    """Class for direct pose prediction using EquiBind."""
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """Initialize EquiBind pose prediction with configuration.
        
        Args:
            config: Configuration dictionary with parameters for pose prediction
        """
        self.config = config or {}
        self.equibind_path = self.config.get("equibind_path", "./EquiBind")
        self.model_path = self.config.get("model_path", os.path.join(self.equibind_path, "pretrained_models", "default"))
        self.device = self.config.get("device", "cuda" if torch.cuda.is_available() else "cpu")
        
    def setup_equibind(self):
        """Set up EquiBind environment if not already set up."""
        if not os.path.exists(self.equibind_path):
            logger.info("Setting up EquiBind environment")
            
            # Clone EquiBind repository
            cmd = [
                "git", "clone", "https://github.com/HannesStark/EquiBind.git", self.equibind_path
            ]
            
            try:
                subprocess.run(cmd, check=True, capture_output=True)
                logger.info(f"EquiBind repository cloned to {self.equibind_path}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to clone EquiBind repository: {e.stderr.decode()}")
                raise
            
            # Download pretrained models
            os.makedirs(os.path.join(self.equibind_path, "pretrained_models"), exist_ok=True)
            
            # Note: In a real implementation, you would download the pretrained models here
            # For this example, we'll assume they're already available or will be downloaded manually
            logger.info("Please download pretrained models manually from the EquiBind repository")
        
        logger.info("EquiBind environment is set up")
    
    def prepare_input(self, protein_path: Union[str, Path], ligand_path: Union[str, Path],
                     output_dir: Optional[Union[str, Path]] = None) -> str:
        """Prepare input files for EquiBind.
        
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
            output_dir = Path(tempfile.mkdtemp(prefix="equibind_input_"))
        else:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Preparing input files for EquiBind in {output_dir}")
        
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
        
        # Create dataset file
        dataset_file = output_dir / "data_to_predict.csv"
        with open(dataset_file, 'w') as f:
            f.write("protein_path,ligand_path\n")
            f.write(f"{protein_dest.name},{ligand_dest.name}\n")
        
        logger.info(f"Input files prepared for EquiBind: {output_dir}")
        return str(output_dir)
    
    def run_prediction(self, input_dir: Union[str, Path], output_dir: Optional[Union[str, Path]] = None) -> str:
        """Run EquiBind pose prediction.
        
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
        
        logger.info(f"Running EquiBind pose prediction with input from {input_dir}")
        
        # Change to EquiBind directory
        original_dir = os.getcwd()
        os.chdir(self.equibind_path)
        
        try:
            # Run EquiBind inference
            cmd = [
                "python", "inference.py", 
                "--config", "configs/inference.yml",
                "--input_path", str(input_dir / "data_to_predict.csv"),
                "--output_directory", str(output_dir),
                "--run_corrections", "True",
                "--model_path", str(self.model_path),
                "--device", self.device
            ]
            
            subprocess.run(cmd, check=True, capture_output=True)
            logger.info(f"EquiBind pose prediction completed successfully: {output_dir}")
            
            # Return to original directory
            os.chdir(original_dir)
            return str(output_dir)
        
        except subprocess.CalledProcessError as e:
            logger.error(f"EquiBind pose prediction failed: {e.stderr.decode()}")
            # Return to original directory
            os.chdir(original_dir)
            raise
        
        except Exception as e:
            logger.error(f"Error during EquiBind pose prediction: {e}")
            # Return to original directory
            os.chdir(original_dir)
            raise
    
    def convert_output_to_pdbqt(self, output_dir: Union[str, Path], 
                               output_path: Optional[Union[str, Path]] = None) -> str:
        """Convert EquiBind output to PDBQT format for docking.
        
        Args:
            output_dir: Directory containing EquiBind prediction results
            output_path: Path to output PDBQT file (optional)
            
        Returns:
            Path to PDBQT file
        """
        output_dir = Path(output_dir)
        
        # Find predicted ligand file
        ligand_files = list(output_dir.glob("*_lig_equibind_corrected.sdf"))
        if not ligand_files:
            ligand_files = list(output_dir.glob("*_lig_equibind.sdf"))
        
        if not ligand_files:
            raise ValueError(f"No predicted ligand files found in {output_dir}")
        
        ligand_file = ligand_files[0]
        
        if output_path is None:
            output_path = ligand_file.with_suffix(".pdbqt")
        else:
            output_path = Path(output_path)
        
        logger.info(f"Converting EquiBind output {ligand_file} to PDBQT format")
        
        # Use Open Babel to convert to PDBQT
        cmd = [
            "obabel", 
            str(ligand_file), 
            "-O", str(output_path), 
            "-xh"  # Add AutoDock atom types and optimize hydrogens
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            logger.info(f"Converted EquiBind output to PDBQT format: {output_path}")
            return str(output_path)
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to convert to PDBQT: {e.stderr.decode()}")
            raise
    
    def predict_pose(self, protein_path: Union[str, Path], ligand_path: Union[str, Path],
                    output_path: Optional[Union[str, Path]] = None) -> str:
        """Full EquiBind pose prediction pipeline.
        
        Args:
            protein_path: Path to protein PDB file
            ligand_path: Path to ligand file
            output_path: Path to output PDBQT file (optional)
            
        Returns:
            Path to predicted pose PDBQT file
        """
        protein_path = Path(protein_path)
        ligand_path = Path(ligand_path)
        
        if output_path is None:
            output_path = ligand_path.parent / f"{ligand_path.stem}_equibind.pdbqt"
        else:
            output_path = Path(output_path)
        
        logger.info(f"Starting full EquiBind pose prediction pipeline")
        
        # Step 1: Set up EquiBind
        self.setup_equibind()
        
        # Step 2: Prepare input files
        input_dir = self.prepare_input(protein_path, ligand_path)
        
        # Step 3: Run prediction
        output_dir = self.run_prediction(input_dir)
        
        # Step 4: Convert output to PDBQT
        pdbqt_path = self.convert_output_to_pdbqt(output_dir, output_path)
        
        # Clean up temporary files if needed
        if self.config.get("cleanup_intermediates", True):
            if Path(input_dir).exists() and Path(input_dir).is_dir() and "equibind_input_" in input_dir:
                shutil.rmtree(input_dir)
        
        logger.info(f"Completed EquiBind pose prediction: {pdbqt_path}")
        return pdbqt_path


def main():
    """Command line interface for EquiBind pose prediction."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Predict ligand poses with EquiBind")
    parser.add_argument("--protein", required=True, help="Path to protein PDB file")
    parser.add_argument("--ligand", required=True, help="Path to ligand file")
    parser.add_argument("--out", dest="output_path", help="Path to output PDBQT file")
    parser.add_argument("--equibind-path", help="Path to EquiBind repository")
    parser.add_argument("--model-path", help="Path to pretrained EquiBind model")
    parser.add_argument("--device", choices=["cpu", "cuda"], help="Device to use for inference")
    parser.add_argument("--keep-intermediates", action="store_true", help="Keep intermediate files")
    
    args = parser.parse_args()
    
    # Configure pose prediction
    config = {"cleanup_intermediates": not args.keep_intermediates}
    if args.equibind_path:
        config["equibind_path"] = args.equibind_path
    if args.model_path:
        config["model_path"] = args.model_path
    if args.device:
        config["device"] = args.device
    
    # Run pose prediction
    equibind = EquiBindPose(config)
    output_path = equibind.predict_pose(
        args.protein,
        args.ligand,
        args.output_path
    )
    
    print(f"Predicted pose saved to: {output_path}")


if __name__ == "__main__":
    main()
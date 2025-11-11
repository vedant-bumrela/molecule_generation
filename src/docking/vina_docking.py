#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
AutoDock Vina docking module for molecular docking.
Handles docking of ligands to proteins using AutoDock Vina.
"""

import os
import logging
import subprocess
import json
from pathlib import Path
from typing import Optional, Union, Dict, Any, List, Tuple

import numpy as np

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class VinaDocking:
    """Class for molecular docking using AutoDock Vina."""
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """Initialize Vina docking with configuration.
        
        Args:
            config: Configuration dictionary with parameters for docking
        """
        self.config = config or {}
        self.vina_executable = self.config.get("vina_executable", "vina")
        
    def detect_binding_site(self, protein_path: Union[str, Path], 
                           ligand_path: Optional[Union[str, Path]] = None) -> Dict[str, float]:
        """Detect binding site in protein, optionally using reference ligand.
        
        Args:
            protein_path: Path to protein PDBQT file
            ligand_path: Path to reference ligand PDBQT file (optional)
            
        Returns:
            Dictionary with binding site center and size
        """
        logger.info(f"Detecting binding site for protein: {protein_path}")
        
        # If reference ligand is provided, use its center as binding site
        if ligand_path:
            # Parse ligand PDBQT to get coordinates
            coords = []
            with open(ligand_path, 'r') as f:
                for line in f:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coords.append((x, y, z))
            
            # Calculate center and size
            coords = np.array(coords)
            center = coords.mean(axis=0)
            size = coords.max(axis=0) - coords.min(axis=0) + 10.0  # Add 10 Å padding
            
            binding_site = {
                "center_x": float(center[0]),
                "center_y": float(center[1]),
                "center_z": float(center[2]),
                "size_x": float(size[0]),
                "size_y": float(size[1]),
                "size_z": float(size[2])
            }
        else:
            # Without reference ligand, use protein center and a default size
            # Parse protein PDBQT to get coordinates
            coords = []
            with open(protein_path, 'r') as f:
                for line in f:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coords.append((x, y, z))
            
            # Calculate center and use default size
            coords = np.array(coords)
            center = coords.mean(axis=0)
            
            binding_site = {
                "center_x": float(center[0]),
                "center_y": float(center[1]),
                "center_z": float(center[2]),
                "size_x": 20.0,  # Default size
                "size_y": 20.0,
                "size_z": 20.0
            }
        
        logger.info(f"Binding site detected: {binding_site}")
        return binding_site
    
    def run_docking(self, receptor_path: Union[str, Path], ligand_path: Union[str, Path],
                   output_path: Optional[Union[str, Path]] = None,
                   binding_site: Optional[Dict[str, float]] = None,
                   exhaustiveness: int = 8, num_modes: int = 9) -> str:
        """Run docking with AutoDock Vina.
        
        Args:
            receptor_path: Path to receptor PDBQT file
            ligand_path: Path to ligand PDBQT file
            output_path: Path to output PDBQT file (optional)
            binding_site: Dictionary with binding site center and size
            exhaustiveness: Exhaustiveness of the global search
            num_modes: Maximum number of binding modes to generate
            
        Returns:
            Path to docking results PDBQT file
        """
        receptor_path = Path(receptor_path)
        ligand_path = Path(ligand_path)
        
        if output_path is None:
            output_path = ligand_path.parent / f"{ligand_path.stem}_docked.pdbqt"
        else:
            output_path = Path(output_path)
        
        # Ensure output directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Detect binding site if not provided
        if binding_site is None:
            binding_site = self.detect_binding_site(receptor_path)
        
        logger.info(f"Running docking for receptor {receptor_path} and ligand {ligand_path}")
        
        # Build Vina command
        cmd = [
            self.vina_executable,
            "--receptor", str(receptor_path),
            "--ligand", str(ligand_path),
            "--center_x", str(binding_site["center_x"]),
            "--center_y", str(binding_site["center_y"]),
            "--center_z", str(binding_site["center_z"]),
            "--size_x", str(binding_site["size_x"]),
            "--size_y", str(binding_site["size_y"]),
            "--size_z", str(binding_site["size_z"]),
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(num_modes),
            "--out", str(output_path)
        ]
        
        # Add optional parameters
        if "cpu" in self.config:
            cmd.extend(["--cpu", str(self.config["cpu"])])
        
        # Run Vina
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"Docking completed successfully: {output_path}")
            
            # Parse Vina output to get scores
            scores = self._parse_vina_output(result.stdout)
            
            # Save scores to JSON file
            scores_path = output_path.with_suffix(".json")
            with open(scores_path, 'w') as f:
                json.dump(scores, f, indent=2)
            
            logger.info(f"Docking scores saved to: {scores_path}")
            return str(output_path)
        
        except subprocess.CalledProcessError as e:
            logger.error(f"Docking failed: {e.stderr}")
            raise
    
    def _parse_vina_output(self, output: str) -> Dict[str, Any]:
        """Parse Vina output to extract scores.
        
        Args:
            output: Vina output text
            
        Returns:
            Dictionary with docking scores
        """
        scores = {"modes": []}
        
        # Extract binding modes and scores
        mode_section = False
        for line in output.split('\n'):
            if "-----+------------+----------+----------" in line:
                mode_section = True
                continue
            
            if mode_section and line.strip() and "Writing output" not in line:
                parts = line.split()
                if len(parts) >= 4 and parts[0].isdigit():
                    mode = {
                        "mode": int(parts[0]),
                        "affinity": float(parts[1]),
                        "rmsd_lb": float(parts[2]),
                        "rmsd_ub": float(parts[3])
                    }
                    scores["modes"].append(mode)
        
        # Add best score
        if scores["modes"]:
            scores["best_affinity"] = scores["modes"][0]["affinity"]
        
        return scores
    
    def extract_poses(self, docked_pdbqt: Union[str, Path], 
                     output_dir: Optional[Union[str, Path]] = None) -> List[str]:
        """Extract individual poses from docked PDBQT file.
        
        Args:
            docked_pdbqt: Path to docked PDBQT file
            output_dir: Directory to save individual poses (optional)
            
        Returns:
            List of paths to individual pose files
        """
        docked_pdbqt = Path(docked_pdbqt)
        
        if output_dir is None:
            output_dir = docked_pdbqt.parent / f"{docked_pdbqt.stem}_poses"
        else:
            output_dir = Path(output_dir)
        
        # Ensure output directory exists
        output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Extracting poses from {docked_pdbqt}")
        
        # Read docked PDBQT file
        with open(docked_pdbqt, 'r') as f:
            content = f.read()
        
        # Split by MODEL/ENDMDL tags
        pose_blocks = content.split("MODEL")[1:]  # Skip first empty part
        
        pose_files = []
        for i, block in enumerate(pose_blocks):
            # Add MODEL tag back
            pose_content = f"MODEL{block}"
            
            # Save to file
            pose_file = output_dir / f"pose_{i+1}.pdbqt"
            with open(pose_file, 'w') as f:
                f.write(pose_content)
            
            pose_files.append(str(pose_file))
        
        logger.info(f"Extracted {len(pose_files)} poses to {output_dir}")
        return pose_files
    
    def convert_pose_to_pdb(self, pdbqt_path: Union[str, Path], 
                           output_path: Optional[Union[str, Path]] = None) -> str:
        """Convert PDBQT pose to PDB format.
        
        Args:
            pdbqt_path: Path to PDBQT pose file
            output_path: Path to output PDB file (optional)
            
        Returns:
            Path to PDB file
        """
        pdbqt_path = Path(pdbqt_path)
        
        if output_path is None:
            output_path = pdbqt_path.with_suffix(".pdb")
        else:
            output_path = Path(output_path)
        
        logger.info(f"Converting {pdbqt_path} to PDB format")
        
        # Use Open Babel to convert PDBQT to PDB
        cmd = [
            "obabel", 
            str(pdbqt_path), 
            "-O", str(output_path)
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            logger.info(f"Converted to PDB format: {output_path}")
            return str(output_path)
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to convert to PDB: {e.stderr.decode()}")
            raise


def main():
    """Command line interface for Vina docking."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Run molecular docking with AutoDock Vina")
    parser.add_argument("--receptor", required=True, help="Path to receptor PDBQT file")
    parser.add_argument("--ligand", required=True, help="Path to ligand PDBQT file")
    parser.add_argument("--out", dest="output_path", help="Path to output PDBQT file")
    parser.add_argument("--center_x", type=float, help="X coordinate of binding site center")
    parser.add_argument("--center_y", type=float, help="Y coordinate of binding site center")
    parser.add_argument("--center_z", type=float, help="Z coordinate of binding site center")
    parser.add_argument("--size_x", type=float, help="X dimension of binding site")
    parser.add_argument("--size_y", type=float, help="Y dimension of binding site")
    parser.add_argument("--size_z", type=float, help="Z dimension of binding site")
    parser.add_argument("--exhaustiveness", type=int, default=8, help="Exhaustiveness of the global search")
    parser.add_argument("--num_modes", type=int, default=9, help="Maximum number of binding modes to generate")
    parser.add_argument("--cpu", type=int, help="Number of CPUs to use")
    parser.add_argument("--extract-poses", action="store_true", help="Extract individual poses")
    
    args = parser.parse_args()
    
    # Check if binding site is provided
    binding_site = None
    if all(getattr(args, f"center_{axis}") is not None for axis in "xyz") and \
       all(getattr(args, f"size_{axis}") is not None for axis in "xyz"):
        binding_site = {
            "center_x": args.center_x,
            "center_y": args.center_y,
            "center_z": args.center_z,
            "size_x": args.size_x,
            "size_y": args.size_y,
            "size_z": args.size_z
        }
    
    # Configure docking
    config = {}
    if args.cpu:
        config["cpu"] = args.cpu
    
    # Run docking
    docking = VinaDocking(config)
    output_path = docking.run_docking(
        args.receptor,
        args.ligand,
        args.output_path,
        binding_site,
        args.exhaustiveness,
        args.num_modes
    )
    
    print(f"Docking results saved to: {output_path}")
    
    # Extract poses if requested
    if args.extract_poses:
        pose_files = docking.extract_poses(output_path)
        print(f"Extracted {len(pose_files)} poses")


if __name__ == "__main__":
    main()
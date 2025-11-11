#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
GNINA ML rescoring module for molecular docking.
Handles rescoring of docking poses using GNINA's CNN scoring function.
"""

import os
import logging
import subprocess
import json
from pathlib import Path
from typing import Optional, Union, Dict, Any, List, Tuple

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class GNINAScoring:
    """Class for rescoring docking poses using GNINA's CNN scoring function."""
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """Initialize GNINA scoring with configuration.
        
        Args:
            config: Configuration dictionary with parameters for scoring
        """
        self.config = config or {}
        self.gnina_executable = self.config.get("gnina_executable", "gnina")
        self.cnn_model = self.config.get("cnn_model", "default")  # Use default GNINA model
        
    def rescore_pose(self, receptor_path: Union[str, Path], ligand_path: Union[str, Path],
                    output_path: Optional[Union[str, Path]] = None) -> Dict[str, Any]:
        """Rescore a single docking pose using GNINA's CNN scoring function.
        
        Args:
            receptor_path: Path to receptor file (PDBQT or PDB)
            ligand_path: Path to ligand pose file (PDBQT or PDB)
            output_path: Path to output JSON file with scores (optional)
            
        Returns:
            Dictionary with CNN scores
        """
        receptor_path = Path(receptor_path)
        ligand_path = Path(ligand_path)
        
        if output_path is None:
            output_path = ligand_path.parent / f"{ligand_path.stem}_gnina_scores.json"
        else:
            output_path = Path(output_path)
        
        logger.info(f"Rescoring pose {ligand_path} with GNINA CNN scoring")
        
        # Build GNINA command for rescoring
        cmd = [
            self.gnina_executable,
            "--receptor", str(receptor_path),
            "--ligand", str(ligand_path),
            "--score_only",  # Only compute scores, don't perform docking
            "--cnn_scoring", "only",  # Use only CNN scoring
            "--cnn_model", self.cnn_model
        ]
        
        # Add optional parameters
        if "cpu" in self.config:
            cmd.extend(["--cpu", str(self.config["cpu"])])
        
        # Run GNINA
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info("GNINA rescoring completed successfully")
            
            # Parse GNINA output to get scores
            scores = self._parse_gnina_output(result.stdout)
            
            # Save scores to JSON file
            with open(output_path, 'w') as f:
                json.dump(scores, f, indent=2)
            
            logger.info(f"GNINA scores saved to: {output_path}")
            return scores
        
        except subprocess.CalledProcessError as e:
            logger.error(f"GNINA rescoring failed: {e.stderr}")
            raise
    
    def rescore_poses(self, receptor_path: Union[str, Path], pose_dir: Union[str, Path],
                     output_path: Optional[Union[str, Path]] = None) -> Dict[str, Any]:
        """Rescore multiple docking poses using GNINA's CNN scoring function.
        
        Args:
            receptor_path: Path to receptor file (PDBQT or PDB)
            pose_dir: Directory containing pose files (PDBQT or PDB)
            output_path: Path to output JSON file with scores (optional)
            
        Returns:
            Dictionary with CNN scores for all poses
        """
        receptor_path = Path(receptor_path)
        pose_dir = Path(pose_dir)
        
        if output_path is None:
            output_path = pose_dir.parent / f"{pose_dir.name}_gnina_scores.json"
        else:
            output_path = Path(output_path)
        
        logger.info(f"Rescoring poses in {pose_dir} with GNINA CNN scoring")
        
        # Find all pose files
        pose_files = list(pose_dir.glob("*.pdbqt")) + list(pose_dir.glob("*.pdb"))
        if not pose_files:
            raise ValueError(f"No pose files found in {pose_dir}")
        
        # Rescore each pose
        all_scores = {"poses": []}
        for pose_file in pose_files:
            try:
                scores = self.rescore_pose(receptor_path, pose_file)
                scores["pose_file"] = str(pose_file)
                all_scores["poses"].append(scores)
            except Exception as e:
                logger.error(f"Failed to rescore {pose_file}: {e}")
        
        # Sort poses by CNN score
        all_scores["poses"].sort(key=lambda x: x.get("cnn_score", float("inf")))
        
        # Add best score
        if all_scores["poses"]:
            all_scores["best_pose"] = all_scores["poses"][0]["pose_file"]
            all_scores["best_cnn_score"] = all_scores["poses"][0].get("cnn_score")
        
        # Save all scores to JSON file
        with open(output_path, 'w') as f:
            json.dump(all_scores, f, indent=2)
        
        logger.info(f"All GNINA scores saved to: {output_path}")
        return all_scores
    
    def _parse_gnina_output(self, output: str) -> Dict[str, Any]:
        """Parse GNINA output to extract CNN scores.
        
        Args:
            output: GNINA output text
            
        Returns:
            Dictionary with CNN scores
        """
        scores = {}
        
        # Extract CNN score
        for line in output.split('\n'):
            if "CNNscore" in line:
                parts = line.split()
                if len(parts) >= 2:
                    scores["cnn_score"] = float(parts[-1])
            elif "CNNaffinity" in line:
                parts = line.split()
                if len(parts) >= 2:
                    scores["cnn_affinity"] = float(parts[-1])
        
        return scores
    
    def run_gnina_docking(self, receptor_path: Union[str, Path], ligand_path: Union[str, Path],
                         output_path: Optional[Union[str, Path]] = None,
                         binding_site: Optional[Dict[str, float]] = None,
                         exhaustiveness: int = 8, num_modes: int = 9) -> str:
        """Run docking with GNINA (includes both sampling and CNN scoring).
        
        Args:
            receptor_path: Path to receptor file (PDBQT or PDB)
            ligand_path: Path to ligand file (PDBQT or PDB)
            output_path: Path to output file (optional)
            binding_site: Dictionary with binding site center and size
            exhaustiveness: Exhaustiveness of the global search
            num_modes: Maximum number of binding modes to generate
            
        Returns:
            Path to docking results file
        """
        receptor_path = Path(receptor_path)
        ligand_path = Path(ligand_path)
        
        if output_path is None:
            output_path = ligand_path.parent / f"{ligand_path.stem}_gnina_docked.sdf"
        else:
            output_path = Path(output_path)
        
        logger.info(f"Running GNINA docking for receptor {receptor_path} and ligand {ligand_path}")
        
        # Build GNINA command
        cmd = [
            self.gnina_executable,
            "--receptor", str(receptor_path),
            "--ligand", str(ligand_path),
            "--out", str(output_path),
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(num_modes),
            "--cnn_scoring", "all"  # Use CNN scoring for both ranking and refinement
        ]
        
        # Add binding site if provided
        if binding_site:
            cmd.extend([
                "--center_x", str(binding_site["center_x"]),
                "--center_y", str(binding_site["center_y"]),
                "--center_z", str(binding_site["center_z"]),
                "--size_x", str(binding_site["size_x"]),
                "--size_y", str(binding_site["size_y"]),
                "--size_z", str(binding_site["size_z"])
            ])
        
        # Add optional parameters
        if "cpu" in self.config:
            cmd.extend(["--cpu", str(self.config["cpu"])])
        
        # Run GNINA
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"GNINA docking completed successfully: {output_path}")
            
            # Parse GNINA output to get scores
            scores = self._parse_gnina_docking_output(result.stdout)
            
            # Save scores to JSON file
            scores_path = output_path.with_suffix(".json")
            with open(scores_path, 'w') as f:
                json.dump(scores, f, indent=2)
            
            logger.info(f"GNINA docking scores saved to: {scores_path}")
            return str(output_path)
        
        except subprocess.CalledProcessError as e:
            logger.error(f"GNINA docking failed: {e.stderr}")
            raise
    
    def _parse_gnina_docking_output(self, output: str) -> Dict[str, Any]:
        """Parse GNINA docking output to extract scores.
        
        Args:
            output: GNINA output text
            
        Returns:
            Dictionary with docking scores
        """
        scores = {"modes": []}
        
        # Extract binding modes and scores
        mode_section = False
        current_mode = None
        
        for line in output.split('\n'):
            if "mode" in line and "affinity" in line:
                mode_section = True
                continue
            
            if mode_section and line.strip():
                if line.startswith("mode"):
                    if current_mode:
                        scores["modes"].append(current_mode)
                    
                    parts = line.split()
                    current_mode = {
                        "mode": int(parts[1]),
                        "affinity": float(parts[3])
                    }
                elif "CNNscore" in line and current_mode:
                    parts = line.split()
                    current_mode["cnn_score"] = float(parts[-1])
                elif "CNNaffinity" in line and current_mode:
                    parts = line.split()
                    current_mode["cnn_affinity"] = float(parts[-1])
        
        # Add last mode
        if current_mode:
            scores["modes"].append(current_mode)
        
        # Add best scores
        if scores["modes"]:
            scores["best_affinity"] = scores["modes"][0]["affinity"]
            scores["best_cnn_score"] = scores["modes"][0].get("cnn_score")
            scores["best_cnn_affinity"] = scores["modes"][0].get("cnn_affinity")
        
        return scores


def main():
    """Command line interface for GNINA scoring."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Rescore docking poses with GNINA CNN scoring")
    parser.add_argument("--receptor", required=True, help="Path to receptor file")
    parser.add_argument("--ligand", help="Path to ligand pose file")
    parser.add_argument("--poses-dir", help="Directory containing pose files")
    parser.add_argument("--out", dest="output_path", help="Path to output file")
    parser.add_argument("--cnn-model", help="CNN model to use for scoring")
    parser.add_argument("--cpu", type=int, help="Number of CPUs to use")
    parser.add_argument("--run-docking", action="store_true", help="Run full GNINA docking instead of just rescoring")
    
    # Docking options (only used with --run-docking)
    parser.add_argument("--center_x", type=float, help="X coordinate of binding site center")
    parser.add_argument("--center_y", type=float, help="Y coordinate of binding site center")
    parser.add_argument("--center_z", type=float, help="Z coordinate of binding site center")
    parser.add_argument("--size_x", type=float, default=20.0, help="X dimension of binding site")
    parser.add_argument("--size_y", type=float, default=20.0, help="Y dimension of binding site")
    parser.add_argument("--size_z", type=float, default=20.0, help="Z dimension of binding site")
    parser.add_argument("--exhaustiveness", type=int, default=8, help="Exhaustiveness of the global search")
    parser.add_argument("--num-modes", type=int, default=9, help="Maximum number of binding modes")
    
    args = parser.parse_args()
    
    # Configure scoring
    config = {}
    if args.cnn_model:
        config["cnn_model"] = args.cnn_model
    if args.cpu:
        config["cpu"] = args.cpu
    
    gnina = GNINAScoring(config)
    
    if args.run_docking:
        # Check if binding site is provided
        binding_site = None
        if args.center_x is not None and args.center_y is not None and args.center_z is not None:
            binding_site = {
                "center_x": args.center_x,
                "center_y": args.center_y,
                "center_z": args.center_z,
                "size_x": args.size_x,
                "size_y": args.size_y,
                "size_z": args.size_z
            }
        
        # Run GNINA docking
        output_path = gnina.run_gnina_docking(
            args.receptor,
            args.ligand,
            args.output_path,
            binding_site,
            args.exhaustiveness,
            args.num_modes
        )
        print(f"GNINA docking results saved to: {output_path}")
    
    elif args.poses_dir:
        # Rescore multiple poses
        scores = gnina.rescore_poses(
            args.receptor,
            args.poses_dir,
            args.output_path
        )
        print(f"Best pose: {scores.get('best_pose')}")
        print(f"Best CNN score: {scores.get('best_cnn_score')}")
    
    elif args.ligand:
        # Rescore single pose
        scores = gnina.rescore_pose(
            args.receptor,
            args.ligand,
            args.output_path
        )
        print(f"CNN score: {scores.get('cnn_score')}")
        print(f"CNN affinity: {scores.get('cnn_affinity')}")
    
    else:
        parser.error("Either --ligand or --poses-dir must be provided")


if __name__ == "__main__":
    main()
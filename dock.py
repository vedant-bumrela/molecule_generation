#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Main CLI script for molecular docking pipeline.
Integrates protein preparation, ligand preparation, and docking.
"""

import os
import sys
import logging
import argparse
from pathlib import Path

from src.preprocessing import ProteinPreparation, LigandPreparation
from src.docking import VinaDocking

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="AI-enabled molecular docking pipeline")
    
    # Input options
    input_group = parser.add_argument_group("Input Options")
    input_group.add_argument("--protein", required=True, help="Path to protein PDB file")
    
    # Mutually exclusive ligand input options
    ligand_group = input_group.add_mutually_exclusive_group(required=True)
    ligand_group.add_argument("--ligand", help="Path to ligand file (SDF, MOL2, etc.)")
    ligand_group.add_argument("--smiles", help="SMILES string of ligand")
    
    # Output options
    output_group = parser.add_argument_group("Output Options")
    output_group.add_argument("--out", dest="output_dir", default="./output", help="Output directory")
    output_group.add_argument("--name", help="Base name for output files")
    
    # Protein preparation options
    protein_group = parser.add_argument_group("Protein Preparation Options")
    protein_group.add_argument("--ph", type=float, default=7.4, help="pH for hydrogen addition (default: 7.4)")
    protein_group.add_argument("--keep-water", action="store_true", help="Keep water molecules in protein")
    
    # Ligand preparation options
    ligand_group = parser.add_argument_group("Ligand Preparation Options")
    ligand_group.add_argument("--num-confs", type=int, default=10, help="Number of conformers to generate (default: 10)")
    
    # Docking options
    docking_group = parser.add_argument_group("Docking Options")
    docking_group.add_argument("--center_x", type=float, help="X coordinate of binding site center")
    docking_group.add_argument("--center_y", type=float, help="Y coordinate of binding site center")
    docking_group.add_argument("--center_z", type=float, help="Z coordinate of binding site center")
    docking_group.add_argument("--size_x", type=float, default=20.0, help="X dimension of binding site (default: 20.0)")
    docking_group.add_argument("--size_y", type=float, default=20.0, help="Y dimension of binding site (default: 20.0)")
    docking_group.add_argument("--size_z", type=float, default=20.0, help="Z dimension of binding site (default: 20.0)")
    docking_group.add_argument("--exhaustiveness", type=int, default=8, help="Exhaustiveness of the global search (default: 8)")
    docking_group.add_argument("--num-modes", type=int, default=9, help="Maximum number of binding modes (default: 9)")
    docking_group.add_argument("--cpu", type=int, help="Number of CPUs to use")
    
    # Other options
    parser.add_argument("--keep-intermediates", action="store_true", help="Keep intermediate files")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    
    return parser.parse_args()

def main():
    """Main function for molecular docking pipeline."""
    # Parse arguments
    args = parse_args()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Determine base name for output files
    if args.name:
        base_name = args.name
    elif args.ligand:
        base_name = Path(args.ligand).stem
    else:
        base_name = "ligand"
    
    # Configuration
    config = {"cleanup_intermediates": not args.keep_intermediates}
    
    logger.info("Starting molecular docking pipeline")
    
    # Step 1: Prepare protein
    logger.info("Step 1: Preparing protein")
    protein_prep = ProteinPreparation(config)
    protein_pdbqt = output_dir / f"{Path(args.protein).stem}_prepared.pdbqt"
    protein_pdbqt = protein_prep.prepare_protein(
        args.protein,
        protein_pdbqt,
        ph=args.ph
    )
    logger.info(f"Prepared protein: {protein_pdbqt}")
    
    # Step 2: Prepare ligand
    logger.info("Step 2: Preparing ligand")
    ligand_prep = LigandPreparation(config)
    ligand_pdbqt = output_dir / f"{base_name}_prepared.pdbqt"
    
    if args.ligand:
        ligand_pdbqt = ligand_prep.prepare_ligand_from_file(
            args.ligand,
            ligand_pdbqt
        )
    else:
        ligand_pdbqt = ligand_prep.prepare_ligand_from_smiles(
            args.smiles,
            ligand_pdbqt,
            num_confs=args.num_confs
        )
    logger.info(f"Prepared ligand: {ligand_pdbqt}")
    
    # Step 3: Run docking
    logger.info("Step 3: Running docking")
    docking = VinaDocking(config)
    
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
    
    # Run docking
    docked_pdbqt = output_dir / f"{base_name}_docked.pdbqt"
    docked_pdbqt = docking.run_docking(
        protein_pdbqt,
        ligand_pdbqt,
        docked_pdbqt,
        binding_site=binding_site,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes
    )
    logger.info(f"Docking completed: {docked_pdbqt}")
    
    # Extract poses
    logger.info("Extracting poses")
    poses_dir = output_dir / f"{base_name}_poses"
    pose_files = docking.extract_poses(docked_pdbqt, poses_dir)
    logger.info(f"Extracted {len(pose_files)} poses to {poses_dir}")
    
    # Convert best pose to PDB
    best_pose_pdb = docking.convert_pose_to_pdb(pose_files[0])
    logger.info(f"Best pose saved as PDB: {best_pose_pdb}")
    
    logger.info("Molecular docking pipeline completed successfully")
    print(f"\nResults summary:")
    print(f"  Prepared protein: {protein_pdbqt}")
    print(f"  Prepared ligand: {ligand_pdbqt}")
    print(f"  Docking results: {docked_pdbqt}")
    print(f"  Best pose: {best_pose_pdb}")
    print(f"  All poses: {poses_dir}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.error(f"Error: {e}", exc_info=True)
        sys.exit(1)
#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Protein preparation module for molecular docking.
Handles tasks like fixing missing atoms, adding hydrogens, and assigning protonation states.
"""

import os
import logging
from pathlib import Path
import subprocess
from typing import Optional, Union, Dict, Any

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ProteinPreparation:
    """Class for preparing protein structures for molecular docking."""
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """Initialize protein preparation with configuration.
        
        Args:
            config: Configuration dictionary with parameters for preparation
        """
        self.config = config or {}
        self.parser = PDBParser(QUIET=True)
        self.io = PDBIO()
        
    def fix_protein(self, pdb_path: Union[str, Path], output_path: Optional[Union[str, Path]] = None) -> str:
        """Fix missing atoms and other issues in protein structure.
        
        Args:
            pdb_path: Path to input PDB file
            output_path: Path to output fixed PDB file (optional)
            
        Returns:
            Path to fixed PDB file
        """
        pdb_path = Path(pdb_path)
        if output_path is None:
            output_path = pdb_path.parent / f"{pdb_path.stem}_fixed.pdb"
        else:
            output_path = Path(output_path)
            
        logger.info(f"Fixing protein structure: {pdb_path}")
        
        # Parse PDB file
        structure = self.parser.get_structure("protein", pdb_path)
        
        # Basic cleanup - remove hetero atoms, waters, etc.
        class ProteinSelect(Select):
            def accept_residue(self, residue):
                return is_aa(residue)
        
        # Write fixed structure
        self.io.set_structure(structure)
        self.io.save(str(output_path), ProteinSelect())
        
        logger.info(f"Fixed protein structure saved to: {output_path}")
        return str(output_path)
    
    def add_hydrogens(self, pdb_path: Union[str, Path], output_path: Optional[Union[str, Path]] = None,
                     ph: float = 7.4) -> str:
        """Add hydrogens to protein structure using Open Babel.
        
        Args:
            pdb_path: Path to input PDB file
            output_path: Path to output PDB file with hydrogens (optional)
            ph: pH value for hydrogen addition
            
        Returns:
            Path to PDB file with hydrogens
        """
        pdb_path = Path(pdb_path)
        if output_path is None:
            output_path = pdb_path.parent / f"{pdb_path.stem}_h.pdb"
        else:
            output_path = Path(output_path)
            
        logger.info(f"Adding hydrogens to protein structure: {pdb_path} at pH {ph}")
        
        # Use Open Babel to add hydrogens
        cmd = [
            "obabel", 
            str(pdb_path), 
            "-O", str(output_path), 
            "-h", # Add hydrogens
            f"-p {ph}"  # Set pH
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            logger.info(f"Added hydrogens to protein structure: {output_path}")
            return str(output_path)
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to add hydrogens: {e.stderr.decode()}")
            raise
    
    def assign_charges(self, pdb_path: Union[str, Path], output_path: Optional[Union[str, Path]] = None) -> str:
        """Assign charges to protein atoms.
        
        Args:
            pdb_path: Path to input PDB file
            output_path: Path to output PDB file with charges (optional)
            
        Returns:
            Path to PDB file with charges
        """
        pdb_path = Path(pdb_path)
        if output_path is None:
            output_path = pdb_path.parent / f"{pdb_path.stem}_charged.pdb"
        else:
            output_path = Path(output_path)
            
        logger.info(f"Assigning charges to protein structure: {pdb_path}")
        
        # Use Open Babel to assign charges (Gasteiger charges)
        cmd = [
            "obabel", 
            str(pdb_path), 
            "-O", str(output_path), 
            "-xg"  # Add Gasteiger charges
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            logger.info(f"Assigned charges to protein structure: {output_path}")
            return str(output_path)
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to assign charges: {e.stderr.decode()}")
            raise
    
    def convert_to_pdbqt(self, pdb_path: Union[str, Path], output_path: Optional[Union[str, Path]] = None) -> str:
        """Convert PDB file to PDBQT format for docking with AutoDock Vina.
        
        Args:
            pdb_path: Path to input PDB file
            output_path: Path to output PDBQT file (optional)
            
        Returns:
            Path to PDBQT file
        """
        pdb_path = Path(pdb_path)
        if output_path is None:
            output_path = pdb_path.parent / f"{pdb_path.stem}.pdbqt"
        else:
            output_path = Path(output_path)
            
        logger.info(f"Converting protein to PDBQT format: {pdb_path}")
        
        # Use Open Babel to convert to PDBQT
        cmd = [
            "obabel", 
            str(pdb_path), 
            "-O", str(output_path), 
            "-xr"  # Add AutoDock atom types
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            logger.info(f"Converted protein to PDBQT format: {output_path}")
            return str(output_path)
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to convert to PDBQT: {e.stderr.decode()}")
            raise
    
    def prepare_protein(self, pdb_path: Union[str, Path], output_path: Optional[Union[str, Path]] = None,
                       ph: float = 7.4) -> str:
        """Full protein preparation pipeline.
        
        Args:
            pdb_path: Path to input PDB file
            output_path: Path to output PDBQT file (optional)
            ph: pH value for hydrogen addition
            
        Returns:
            Path to prepared PDBQT file
        """
        pdb_path = Path(pdb_path)
        if output_path is None:
            output_path = pdb_path.parent / f"{pdb_path.stem}_prepared.pdbqt"
        else:
            output_path = Path(output_path)
            
        logger.info(f"Starting full protein preparation pipeline for: {pdb_path}")
        
        # Step 1: Fix protein structure
        fixed_pdb = self.fix_protein(pdb_path)
        
        # Step 2: Add hydrogens
        h_pdb = self.add_hydrogens(fixed_pdb, ph=ph)
        
        # Step 3: Assign charges
        charged_pdb = self.assign_charges(h_pdb)
        
        # Step 4: Convert to PDBQT
        pdbqt_path = self.convert_to_pdbqt(charged_pdb, output_path)
        
        # Clean up intermediate files if needed
        if self.config.get("cleanup_intermediates", True):
            for path in [fixed_pdb, h_pdb, charged_pdb]:
                if Path(path).exists():
                    os.remove(path)
        
        logger.info(f"Completed protein preparation: {pdbqt_path}")
        return pdbqt_path


def main():
    """Command line interface for protein preparation."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Prepare protein structure for molecular docking")
    parser.add_argument("--in", dest="input_path", required=True, help="Input PDB file path")
    parser.add_argument("--out", dest="output_path", help="Output PDBQT file path")
    parser.add_argument("--ph", type=float, default=7.4, help="pH for hydrogen addition (default: 7.4)")
    parser.add_argument("--keep-intermediates", action="store_true", help="Keep intermediate files")
    
    args = parser.parse_args()
    
    config = {"cleanup_intermediates": not args.keep_intermediates}
    prep = ProteinPreparation(config)
    
    output_path = prep.prepare_protein(
        args.input_path,
        args.output_path,
        ph=args.ph
    )
    
    print(f"Prepared protein saved to: {output_path}")


if __name__ == "__main__":
    main()
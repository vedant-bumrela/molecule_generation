#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Ligand preparation module for molecular docking.
Handles tasks like SMILES sanitization, 3D conformer generation, and charge assignment.
"""

import os
import logging
from pathlib import Path
import subprocess
from typing import Optional, Union, Dict, Any, List, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class LigandPreparation:
    """Class for preparing ligand structures for molecular docking."""
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """Initialize ligand preparation with configuration.
        
        Args:
            config: Configuration dictionary with parameters for preparation
        """
        self.config = config or {}
        
    def sanitize_smiles(self, smiles: str) -> Chem.Mol:
        """Sanitize SMILES string and return RDKit molecule.
        
        Args:
            smiles: SMILES string of ligand
            
        Returns:
            RDKit molecule object
        """
        logger.info(f"Sanitizing SMILES: {smiles}")
        
        # Parse SMILES and sanitize
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        
        # Add explicit hydrogens
        mol = Chem.AddHs(mol)
        
        logger.info("SMILES sanitization completed")
        return mol
    
    def generate_3d_conformers(self, mol: Chem.Mol, num_confs: int = 10) -> Chem.Mol:
        """Generate 3D conformers for molecule.
        
        Args:
            mol: RDKit molecule object
            num_confs: Number of conformers to generate
            
        Returns:
            RDKit molecule with 3D conformers
        """
        logger.info(f"Generating {num_confs} 3D conformers")
        
        # Generate 3D coordinates
        AllChem.EmbedMultipleConfs(
            mol, 
            numConfs=num_confs,
            useExpTorsionAnglePrefs=True,
            useBasicKnowledge=True,
            randomSeed=42,
            numThreads=0  # Use all available CPUs
        )
        
        # Optimize conformers with MMFF
        for conf_id in range(mol.GetNumConformers()):
            AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
        
        logger.info(f"Generated {mol.GetNumConformers()} conformers")
        return mol
    
    def assign_charges(self, mol: Chem.Mol) -> Chem.Mol:
        """Assign partial charges to molecule atoms.
        
        Args:
            mol: RDKit molecule object
            
        Returns:
            RDKit molecule with charges
        """
        logger.info("Assigning partial charges to ligand")
        
        # Use MMFF94 charges
        props = AllChem.MMFFGetMoleculeProperties(mol)
        if props is None:
            logger.warning("Could not compute MMFF charges, using Gasteiger charges instead")
            # Fallback to Gasteiger charges
            Chem.ComputeGasteigerCharges(mol)
        
        logger.info("Charge assignment completed")
        return mol
    
    def save_mol_to_file(self, mol: Chem.Mol, output_path: Union[str, Path], 
                         conf_id: int = 0, file_format: str = "sdf") -> str:
        """Save molecule to file.
        
        Args:
            mol: RDKit molecule object
            output_path: Path to output file
            conf_id: Conformer ID to save
            file_format: Output file format (sdf, mol2, pdb)
            
        Returns:
            Path to saved file
        """
        output_path = Path(output_path)
        
        # Ensure directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Saving molecule to {output_path}")
        
        # Save molecule based on format
        if file_format.lower() == "sdf":
            writer = Chem.SDWriter(str(output_path))
            writer.write(mol, conf_id)
            writer.close()
        elif file_format.lower() == "pdb":
            Chem.MolToPDBFile(mol, str(output_path), conf_id)
        else:
            # For other formats, use Open Babel conversion
            temp_sdf = output_path.with_suffix(".sdf")
            writer = Chem.SDWriter(str(temp_sdf))
            writer.write(mol, conf_id)
            writer.close()
            
            # Convert using Open Babel
            cmd = [
                "obabel", 
                str(temp_sdf), 
                "-O", str(output_path)
            ]
            
            try:
                subprocess.run(cmd, check=True, capture_output=True)
                # Remove temporary file
                os.remove(temp_sdf)
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to convert file format: {e.stderr.decode()}")
                raise
        
        logger.info(f"Molecule saved to {output_path}")
        return str(output_path)
    
    def convert_to_pdbqt(self, input_path: Union[str, Path], 
                         output_path: Optional[Union[str, Path]] = None) -> str:
        """Convert molecule file to PDBQT format for docking with AutoDock Vina.
        
        Args:
            input_path: Path to input molecule file
            output_path: Path to output PDBQT file (optional)
            
        Returns:
            Path to PDBQT file
        """
        input_path = Path(input_path)
        if output_path is None:
            output_path = input_path.with_suffix(".pdbqt")
        else:
            output_path = Path(output_path)
            
        logger.info(f"Converting ligand to PDBQT format: {input_path}")
        
        # Use Open Babel to convert to PDBQT
        cmd = [
            "obabel", 
            str(input_path), 
            "-O", str(output_path), 
            "-xh"  # Add AutoDock atom types and optimize hydrogens
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            logger.info(f"Converted ligand to PDBQT format: {output_path}")
            return str(output_path)
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to convert to PDBQT: {e.stderr.decode()}")
            raise
    
    def prepare_ligand_from_smiles(self, smiles: str, output_path: Optional[Union[str, Path]] = None,
                                  num_confs: int = 10) -> str:
        """Full ligand preparation pipeline from SMILES.
        
        Args:
            smiles: SMILES string of ligand
            output_path: Path to output PDBQT file (optional)
            num_confs: Number of conformers to generate
            
        Returns:
            Path to prepared PDBQT file
        """
        if output_path is None:
            output_path = Path("ligand_prepared.pdbqt")
        else:
            output_path = Path(output_path)
            
        logger.info(f"Starting full ligand preparation pipeline for SMILES: {smiles}")
        
        # Step 1: Sanitize SMILES
        mol = self.sanitize_smiles(smiles)
        
        # Step 2: Generate 3D conformers
        mol = self.generate_3d_conformers(mol, num_confs=num_confs)
        
        # Step 3: Assign charges
        mol = self.assign_charges(mol)
        
        # Step 4: Save to SDF file
        sdf_path = output_path.with_suffix(".sdf")
        self.save_mol_to_file(mol, sdf_path, conf_id=0, file_format="sdf")
        
        # Step 5: Convert to PDBQT
        pdbqt_path = self.convert_to_pdbqt(sdf_path, output_path)
        
        # Clean up intermediate files if needed
        if self.config.get("cleanup_intermediates", True):
            if Path(sdf_path).exists():
                os.remove(sdf_path)
        
        logger.info(f"Completed ligand preparation: {pdbqt_path}")
        return pdbqt_path
    
    def prepare_ligand_from_file(self, input_path: Union[str, Path], 
                                output_path: Optional[Union[str, Path]] = None) -> str:
        """Full ligand preparation pipeline from file.
        
        Args:
            input_path: Path to input molecule file (SDF, MOL2, etc.)
            output_path: Path to output PDBQT file (optional)
            
        Returns:
            Path to prepared PDBQT file
        """
        input_path = Path(input_path)
        if output_path is None:
            output_path = input_path.with_suffix(".pdbqt")
        else:
            output_path = Path(output_path)
            
        logger.info(f"Starting full ligand preparation pipeline for file: {input_path}")
        
        # Load molecule from file
        if input_path.suffix.lower() == ".sdf":
            mol = Chem.SDMolSupplier(str(input_path))[0]
        elif input_path.suffix.lower() == ".mol":
            mol = Chem.MolFromMolFile(str(input_path))
        else:
            # For other formats, use Open Babel to convert to SDF first
            temp_sdf = input_path.with_suffix(".sdf")
            cmd = [
                "obabel", 
                str(input_path), 
                "-O", str(temp_sdf)
            ]
            
            try:
                subprocess.run(cmd, check=True, capture_output=True)
                mol = Chem.SDMolSupplier(str(temp_sdf))[0]
                # Remove temporary file
                os.remove(temp_sdf)
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to convert file format: {e.stderr.decode()}")
                raise
        
        if mol is None:
            raise ValueError(f"Could not load molecule from file: {input_path}")
        
        # Add hydrogens if not present
        if Chem.AddHs(mol).GetNumAtoms() > mol.GetNumAtoms():
            mol = Chem.AddHs(mol)
        
        # Assign charges
        mol = self.assign_charges(mol)
        
        # Save to SDF file
        sdf_path = output_path.with_suffix(".sdf")
        self.save_mol_to_file(mol, sdf_path, conf_id=0, file_format="sdf")
        
        # Convert to PDBQT
        pdbqt_path = self.convert_to_pdbqt(sdf_path, output_path)
        
        # Clean up intermediate files if needed
        if self.config.get("cleanup_intermediates", True):
            if Path(sdf_path).exists():
                os.remove(sdf_path)
        
        logger.info(f"Completed ligand preparation: {pdbqt_path}")
        return pdbqt_path


def main():
    """Command line interface for ligand preparation."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Prepare ligand for molecular docking")
    parser.add_argument("--smiles", help="SMILES string of ligand")
    parser.add_argument("--in", dest="input_path", help="Input molecule file path")
    parser.add_argument("--out", dest="output_path", help="Output PDBQT file path")
    parser.add_argument("--num-confs", type=int, default=10, help="Number of conformers to generate (default: 10)")
    parser.add_argument("--keep-intermediates", action="store_true", help="Keep intermediate files")
    
    args = parser.parse_args()
    
    if not args.smiles and not args.input_path:
        parser.error("Either --smiles or --in must be provided")
    
    config = {"cleanup_intermediates": not args.keep_intermediates}
    prep = LigandPreparation(config)
    
    if args.smiles:
        output_path = prep.prepare_ligand_from_smiles(
            args.smiles,
            args.output_path,
            num_confs=args.num_confs
        )
    else:
        output_path = prep.prepare_ligand_from_file(
            args.input_path,
            args.output_path
        )
    
    print(f"Prepared ligand saved to: {output_path}")


if __name__ == "__main__":
    main()
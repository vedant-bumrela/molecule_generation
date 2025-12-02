#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Molecule utility functions for ChemVerse.
Handles SMILES conversion, validation, and property calculation.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from typing import Tuple, Dict


class MoleculeConverter:
    """Convert between different molecular formats."""
    
    def smiles_to_pdb(self, smiles: str) -> str:
        """
        Convert SMILES to PDB format with 3D coordinates.
        
        Args:
            smiles: SMILES string
            
        Returns:
            PDB format string
        """
        try:
            # Create molecule from SMILES
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return None
            
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
            
            # Convert to PDB
            pdb_string = Chem.MolToPDBBlock(mol)
            
            return pdb_string
            
        except Exception as e:
            print(f"Error converting SMILES to PDB: {e}")
            return None
    
    def smiles_to_mol(self, smiles: str):
        """Convert SMILES to RDKit molecule object."""
        try:
            return Chem.MolFromSmiles(smiles)
        except:
            return None


class MoleculeValidator:
    """Validate molecules and calculate properties."""
    
    def validate_smiles(self, smiles: str) -> Tuple[bool, str]:
        """
        Validate SMILES string.
        
        Returns:
            (is_valid, message)
        """
        if not smiles or not isinstance(smiles, str):
            return False, "SMILES string is required"
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False, "Invalid SMILES string"
            
            # Check for reasonable size
            num_atoms = mol.GetNumAtoms()
            if num_atoms == 0:
                return False, "Molecule has no atoms"
            if num_atoms > 200:
                return False, "Molecule too large (max 200 atoms)"
            
            return True, "Valid"
            
        except Exception as e:
            return False, f"Validation error: {str(e)}"
    
    def get_properties(self, smiles: str) -> Dict:
        """
        Calculate molecular properties.
        
        Returns:
            Dictionary with MW, formula, logP, etc.
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return {}
            
            properties = {
                'molecular_weight': round(Descriptors.MolWt(mol), 2),
                'formula': Chem.rdMolDescriptors.CalcMolFormula(mol),
                'logp': round(Descriptors.MolLogP(mol), 2),
                'hbd': Descriptors.NumHDonors(mol),
                'hba': Descriptors.NumHAcceptors(mol),
                'tpsa': round(Descriptors.TPSA(mol), 2),
                'num_atoms': mol.GetNumAtoms(),
                'num_bonds': mol.GetNumBonds(),
                'num_rings': Descriptors.RingCount(mol),
                'num_aromatic_rings': Descriptors.NumAromaticRings(mol)
            }
            
            return properties
            
        except Exception as e:
            print(f"Error calculating properties: {e}")
            return {}

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Enzyme structure management for HDI analysis.
Handles downloading and preparing CYP450 enzyme structures.
"""

import os
import requests
from pathlib import Path
from typing import Optional, Dict

# CYP450 enzyme database
CYP450_ENZYMES = {
    'CYP3A4': {
        'pdb_id': '1TQN',
        'name': 'Cytochrome P450 3A4',
        'percentage': '~50% of drugs'
    },
    'CYP2D6': {
        'pdb_id': '2F9Q',
        'name': 'Cytochrome P450 2D6',
        'percentage': '~25% of drugs'
    },
    'CYP2C9': {
        'pdb_id': '1R9O',
        'name': 'Cytochrome P450 2C9',
        'percentage': 'Warfarin, NSAIDs'
    },
    'CYP2C19': {
        'pdb_id': '4GQS',
        'name': 'Cytochrome P450 2C19',
        'percentage': 'PPIs, Clopidogrel'
    },
    'CYP1A2': {
        'pdb_id': '2HI4',
        'name': 'Cytochrome P450 1A2',
        'percentage': 'Caffeine metabolism'
    }
}


class EnzymeStructureManager:
    """Manage CYP450 enzyme structures for HDI analysis."""
    
    def __init__(self, base_dir: Optional[str] = None):
        """
        Initialize enzyme structure manager.
        
        Args:
            base_dir: Base directory for enzyme structures
        """
        if base_dir is None:
            # Default to data/hdi/enzymes directory
            base_dir = Path(__file__).parents[2] / 'data' / 'hdi' / 'enzymes'
        
        self.base_dir = Path(base_dir)
        self.prepared_dir = self.base_dir / 'prepared'
        
        # Create directories if they don't exist
        self.base_dir.mkdir(parents=True, exist_ok=True)
        self.prepared_dir.mkdir(parents=True, exist_ok=True)
    
    def get_enzyme_structure(self, enzyme_name: str) -> Optional[str]:
        """
        Get path to prepared enzyme structure.
        
        Args:
            enzyme_name: Name of enzyme (e.g., 'CYP3A4')
            
        Returns:
            Path to prepared PDBQT file, or None if not available
        """
        prepared_path = self.prepared_dir / f"{enzyme_name.lower()}.pdbqt"
        
        if prepared_path.exists():
            return str(prepared_path)
        
        # If prepared file doesn't exist, check for raw PDB
        raw_path = self.base_dir / f"{enzyme_name.lower()}.pdb"
        if raw_path.exists():
            # Prepare the structure
            return self._prepare_enzyme(enzyme_name, str(raw_path))
        
        return None
    
    def _prepare_enzyme(self, enzyme_name: str, pdb_path: str) -> str:
        """
        Prepare enzyme structure for docking (convert to PDBQT).
        
        Args:
            enzyme_name: Name of enzyme
            pdb_path: Path to raw PDB file
            
        Returns:
            Path to prepared PDBQT file
        """
        from ..preprocessing.protein_prep import ProteinPreparation
        
        output_path = self.prepared_dir / f"{enzyme_name.lower()}.pdbqt"
        
        protein_prep = ProteinPreparation()
        prepared = protein_prep.prepare_protein(
            pdb_path,
            output_path=str(output_path),
            add_hydrogens=True,
            assign_charges=True
        )
        
        return prepared
    
    def validate_structures(self) -> Dict[str, bool]:
        """
        Check which enzyme structures are available.
        
        Returns:
            Dictionary mapping enzyme names to availability
        """
        availability = {}
        
        for enzyme_name in CYP450_ENZYMES.keys():
            structure_path = self.get_enzyme_structure(enzyme_name)
            availability[enzyme_name] = structure_path is not None
        
        return availability
    
    def download_structure(self, enzyme_name: str) -> Optional[str]:
        """
        Download enzyme structure from RCSB PDB.
        
        Args:
            enzyme_name: Name of enzyme
            
        Returns:
            Path to downloaded PDB file, or None if failed
        """
        if enzyme_name not in CYP450_ENZYMES:
            return None
        
        pdb_id = CYP450_ENZYMES[enzyme_name]['pdb_id']
        output_path = self.base_dir / f"{enzyme_name.lower()}.pdb"
        
        try:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            response = requests.get(url)
            response.raise_for_status()
            
            with open(output_path, 'w') as f:
                f.write(response.text)
            
            return str(output_path)
        except Exception as e:
            print(f"Failed to download {pdb_id}: {e}")
            return None
    
    def get_available_enzymes(self) -> list:
        """
        Get list of available enzyme names.
        
        Returns:
            List of enzyme names with available structures
        """
        availability = self.validate_structures()
        return [name for name, available in availability.items() if available]
    
    def setup_all_enzymes(self) -> Dict[str, str]:
        """
        Download and prepare all CYP450 enzyme structures.
        
        Returns:
            Dictionary mapping enzyme names to prepared PDBQT paths
        """
        results = {}
        
        for enzyme_name in CYP450_ENZYMES.keys():
            # Check if already available
            structure_path = self.get_enzyme_structure(enzyme_name)
            
            if structure_path is None:
                # Try to download
                print(f"Downloading {enzyme_name}...")
                pdb_path = self.download_structure(enzyme_name)
                
                if pdb_path:
                    # Prepare the structure
                    structure_path = self._prepare_enzyme(enzyme_name, pdb_path)
            
            if structure_path:
                results[enzyme_name] = structure_path
                print(f"✓ {enzyme_name} ready: {structure_path}")
            else:
                print(f"✗ {enzyme_name} failed")
        
        return results

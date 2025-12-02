#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Reaction prediction engine for ChemVerse.
Uses RDKit reaction templates to predict products.
"""

from typing import List, Dict
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import re


class ReactionEngine:
    """Predict chemical reactions using SMARTS templates."""
    
    # Reaction templates (SMARTS based)
    TEMPLATES = {
        'acid_base': {
            'name': 'Acid-Base Neutralization',
            'description': 'Proton transfer from acid to base',
            'smarts': '[CX3](=O)[OX2H1].[OH-]>>[CX3](=O)[O-].[OH2]',
            'example': 'CH3COOH + OH- → CH3COO- + H2O'
        },
        'esterification': {
            'name': 'Fischer Esterification',
            'description': 'Formation of ester from carboxylic acid and alcohol',
            'smarts': '[CX3](=O)[OX2H1].[OX2H1]>>[CX3](=O)[OX2].[OH2]',
            'example': 'RCOOH + R\'OH → RCOOR\' + H2O'
        },
        'addition': {
            'name': 'Addition Reaction',
            'description': 'Addition across double bond',
            'smarts': '[CX3]=[CX3].[H][H]>>[CX4][CX4]',
            'example': 'C=C + H2 → C-C'
        },
        'substitution': {
            'name': 'Nucleophilic Substitution',
            'description': 'Substitution of leaving group',
            'smarts': '[CX4][Cl].[ O-]>>[CX4][O].[Cl-]',
            'example': 'R-Cl + OH- → R-OH + Cl-'
        }
    }
    
    def predict_products(self, reactant_smiles: List[str], reaction_type: str) -> Dict:
        """
        Predict reaction products using template matching.
        
        Args:
            reactant_smiles: List of SMILES strings for reactants
            reaction_type: Type of reaction (matches TEMPLATES keys)
            
        Returns:
            {
                'products': [SMILES],
                'feasible': bool,
                'energy_change': float,
                'mechanism': str,
                'safety_warnings': [str]
            }
        """
        if reaction_type not in self.TEMPLATES:
            return {'error': f'Unknown reaction type: {reaction_type}'}
        
        template = self.TEMPLATES[reaction_type]
        
        # Simple product generation for demonstration
        # In production, use RDKit AllChem.ReactionFromSmarts
        products = self._apply_template(reactant_smiles, template)
        
        # Check feasibility
        feasibility = self._check_feasibility(reactant_smiles, products)
        
        # Safety check
        safety_warnings = self._check_safety(reactant_smiles + products)
        
        return {
            'products': products,
            'feasible': feasibility['feasible'],
            'energy_change': feasibility.get('delta_g', 0),
            'mechanism': template['name'],
            'description': template['description'],
            'safety_warnings': safety_warnings
        }
    
    def _apply_template(self, reactants: List[str], template: Dict) -> List[str]:
        """Apply reaction template to get products."""
        try:
            # Parse reaction SMARTS
            rxn = AllChem.ReactionFromSmarts(template['smarts'])
            
            # Convert SMILES to molecules
            reactant_mols = [Chem.MolFromSmiles(s) for s in reactants]
            reactant_mols = [mol for mol in reactant_mols if mol is not None]
            
            if len(reactant_mols) < len(reactants):
                return []  # Invalid reactants
            
            # Run reaction
            products = []
            try:
                product_sets = rxn.RunReactants(tuple(reactant_mols))
                
                for product_set in product_sets:
                    for product_mol in product_set:
                        smiles = Chem.MolToSmiles(product_mol)
                        if smiles and smiles not in products:
                            products.append(smiles)
            except:
                # If reaction fails, return simple concatenation as placeholder
                products = [reactants[0] + reactants[1][:10]]  # Simplified
            
            return products if products else ['CC(=O)O']  # Placeholder
            
        except Exception as e:
            print(f"Reaction error: {e}")
            return ['CC(=O)O']  # Return placeholder product
    
    def _check_feasibility(self, reactants: List[str], products: List[str]) -> Dict:
        """
        Simple thermodynamic feasibility check.
        In production, use quantum chemistry or ML models.
        """
        try:
            # Estimate energy based on molecular weight (very simplified!)
            reactant_energy = sum([self._estimate_energy(r) for r in reactants])
            product_energy = sum([self._estimate_energy(p) for p in products])
            
            delta_g = product_energy - reactant_energy
            
            # Simple rule: ΔG < 0 = feasible
            feasible = delta_g < 0
            
            return {
                'feasible': feasible,
                'delta_g': round(delta_g, 2),
                'note': 'Simplified estimation'
            }
        except:
            return {'feasible': True, 'delta_g': -5.0}  # Default
    
    def _estimate_energy(self, smiles: str) -> float:
        """
        Rough energy estimation (placeholder).
        In production, use: DFT calculations, ML models, or thermodynamic tables.
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return 0
            
            # Use molecular weight as very rough proxy
            mw = Descriptors.MolWt(mol)
            
            # Simplified: larger molecules = lower energy (more stable)
            return -0.5 * mw
        except:
            return 0
    
    def _check_safety(self, smiles_list: List[str]) -> List[str]:
        """
        Check for hazardous functional groups.
        """
        warnings = []
        
        hazard_patterns = {
            'Explosive': ['N(=O)(=O)', '[N+](=O)[O-]'],  # Nitro groups
            'Toxic': ['C#N', '[As]', '[Hg]'],  # Cyanides, heavy metals
            'Corrosive': ['S(=O)(=O)[OH]', 'OS(=O)(=O)O'],  # Strong acids
            'Flammable': ['C#C', 'CC#C'],  # Alkynes
        }
        
        for smiles in smiles_list:
            try:
                mol = Chem.MolFromSmiles(smiles)
                if not mol:
                    continue
                
                for hazard_type, patterns in hazard_patterns.items():
                    for pattern in patterns:
                        pattern_mol = Chem.MolFromSmarts(pattern)
                        if pattern_mol and mol.HasSubstructMatch(pattern_mol):
                            warnings.append(f'{hazard_type} - handles with care')
                            break
            except:
                continue
        
        return list(set(warnings))  # Remove duplicates

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Feature extraction for CYP450 ML models.
Combines molecular fingerprints, descriptors, and pharmacophore features.
"""

import numpy as np
from typing import Dict, List, Optional
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski, Crippen, MACCSkeys
from rdkit.Chem import rdMolDescriptors
import warnings
warnings.filterwarnings('ignore')


class MolecularFeatureExtractor:
    """Extract comprehensive molecular features for ML models."""
    
    def __init__(self, 
                 morgan_radius: int = 2,
                 morgan_bits: int = 2048,
                 include_counts: bool = False):
        """
        Initialize feature extractor.
        
        Args:
            morgan_radius: Radius for Morgan fingerprint
            morgan_bits: Number of bits for Morgan fingerprint
            include_counts: Whether to include count-based fingerprints
        """
        self.morgan_radius = morgan_radius
        self.morgan_bits = morgan_bits
        self.include_counts = include_counts
        
        # Feature names for interpretability
        self.descriptor_names = []
        self.fingerprint_names = []
        
    def extract_all_features(self, smiles: str) -> Optional[np.ndarray]:
        """
        Extract all features from SMILES string.
        
        Args:
            smiles: SMILES string
            
        Returns:
            Feature vector as numpy array, or None if invalid
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        features = []
        
        # 1. Molecular descriptors (20 features)
        descriptors = self.extract_descriptors(mol)
        features.extend(descriptors.values())
        
        # 2. Morgan fingerprint (2048 bits)
        morgan_fp = self.extract_morgan_fingerprint(mol)
        features.extend(morgan_fp)
        
        # 3. MACCS keys (166 bits)
        maccs = self.extract_maccs_keys(mol)
        features.extend(maccs)
        
        # 4. Pharmacophore features (10 features)
        pharmacophore = self.extract_pharmacophore_features(mol)
        features.extend(pharmacophore.values())
        
        return np.array(features, dtype=np.float32)
    
    def extract_descriptors(self, mol) -> Dict[str, float]:
        """
        Extract molecular descriptors.
        
        Returns:
            Dictionary of descriptor name -> value
        """
        try:
            descriptors = {
                # Size and shape
                'mol_weight': Descriptors.MolWt(mol),
                'heavy_atom_count': Lipinski.HeavyAtomCount(mol),
                'num_atoms': mol.GetNumAtoms(),
                'num_rotatable_bonds': Lipinski.NumRotatableBonds(mol),
                
                # Lipophilicity
                'logp': Crippen.MolLogP(mol),
                'molar_refractivity': Crippen.MolMR(mol),
                
                # Hydrogen bonding
                'hbd': Lipinski.NumHDonors(mol),
                'hba': Lipinski.NumHAcceptors(mol),
                
                # Aromaticity
                'aromatic_rings': Lipinski.NumAromaticRings(mol),
                'aliphatic_rings': Lipinski.NumAliphaticRings(mol),
                'num_aromatic_atoms': sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic()),
                
                # Polarity
                'tpsa': Descriptors.TPSA(mol),
                
                # Complexity
                'num_heteroatoms': Lipinski.NumHeteroatoms(mol),
                'fraction_csp3': Lipinski.FractionCSP3(mol),
                'num_rings': Lipinski.RingCount(mol),
                
                # Additional descriptors
                'num_saturated_rings': Lipinski.NumSaturatedRings(mol),
                'num_valence_electrons': Descriptors.NumValenceElectrons(mol),
                'bertz_complexity': Descriptors.BertzCT(mol),
                
                # Charge-related
                'formal_charge': Chem.GetFormalCharge(mol),
                'num_radical_electrons': Descriptors.NumRadicalElectrons(mol),
            }
            
            self.descriptor_names = list(descriptors.keys())
            return descriptors
            
        except Exception as e:
            print(f"Error extracting descriptors: {e}")
            return {name: 0.0 for name in self.descriptor_names or range(20)}
    
    def extract_morgan_fingerprint(self, mol) -> List[int]:
        """
        Extract Morgan (ECFP) fingerprint.
        
        Returns:
            List of bits (0 or 1)
        """
        try:
            fp = AllChem.GetMorganFingerprintAsBitVect(
                mol,
                radius=self.morgan_radius,
                nBits=self.morgan_bits
            )
            return list(fp)
        except Exception as e:
            print(f"Error extracting Morgan fingerprint: {e}")
            return [0] * self.morgan_bits
    
    def extract_maccs_keys(self, mol) -> List[int]:
        """
        Extract MACCS keys (166 structural keys).
        
        Returns:
            List of bits (0 or 1)
        """
        try:
            maccs = MACCSkeys.GenMACCSKeys(mol)
            return list(maccs)
        except Exception as e:
            print(f"Error extracting MACCS keys: {e}")
            return [0] * 167  # MACCS has 167 bits (0-indexed)
    
    def extract_pharmacophore_features(self, mol) -> Dict[str, float]:
        """
        Extract custom pharmacophore features relevant to CYP450 binding.
        
        Returns:
            Dictionary of feature name -> value
        """
        try:
            features = {
                'has_basic_nitrogen': int(self._has_basic_nitrogen(mol)),
                'has_aromatic_n': int(self._has_aromatic_nitrogen(mol)),
                'num_benzene_rings': self._count_benzene_rings(mol),
                'has_halogen': int(self._has_halogen(mol)),
                'num_hydroxyl': self._count_functional_group(mol, '[OH]'),
                'num_carbonyl': self._count_functional_group(mol, '[C](=O)'),
                'num_amine': self._count_functional_group(mol, '[N;!H0]'),
                'has_sulfur': int(any(atom.GetSymbol() == 'S' for atom in mol.GetAtoms())),
                'aromatic_ratio': self._calculate_aromatic_ratio(mol),
                'lipinski_violations': self._count_lipinski_violations(mol),
            }
            return features
            
        except Exception as e:
            print(f"Error extracting pharmacophore features: {e}")
            return {f'pharm_feat_{i}': 0.0 for i in range(10)}
    
    def _has_basic_nitrogen(self, mol) -> bool:
        """Check if molecule has basic nitrogen (potential H-bond donor)."""
        pattern = Chem.MolFromSmarts('[N;!$(N-C=O);!$(N-S(=O)=O)]')
        return mol.HasSubstructMatch(pattern) if pattern else False
    
    def _has_aromatic_nitrogen(self, mol) -> bool:
        """Check if molecule has aromatic nitrogen."""
        pattern = Chem.MolFromSmarts('[n]')
        return mol.HasSubstructMatch(pattern) if pattern else False
    
    def _count_benzene_rings(self, mol) -> int:
        """Count number of benzene rings."""
        pattern = Chem.MolFromSmarts('c1ccccc1')
        if pattern:
            return len(mol.GetSubstructMatches(pattern))
        return 0
    
    def _has_halogen(self, mol) -> bool:
        """Check if molecule contains halogen atoms."""
        return any(atom.GetSymbol() in ['F', 'Cl', 'Br', 'I'] for atom in mol.GetAtoms())
    
    def _count_functional_group(self, mol, smarts: str) -> int:
        """Count occurrences of a functional group defined by SMARTS."""
        pattern = Chem.MolFromSmarts(smarts)
        if pattern:
            return len(mol.GetSubstructMatches(pattern))
        return 0
    
    def _calculate_aromatic_ratio(self, mol) -> float:
        """Calculate ratio of aromatic atoms to total atoms."""
        total_atoms = mol.GetNumAtoms()
        if total_atoms == 0:
            return 0.0
        aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
        return aromatic_atoms / total_atoms
    
    def _count_lipinski_violations(self, mol) -> int:
        """Count number of Lipinski rule violations."""
        violations = 0
        
        # Rule of 5
        if Descriptors.MolWt(mol) > 500:
            violations += 1
        if Crippen.MolLogP(mol) > 5:
            violations += 1
        if Lipinski.NumHDonors(mol) > 5:
            violations += 1
        if Lipinski.NumHAcceptors(mol) > 10:
            violations += 1
            
        return violations
    
    def get_feature_names(self) -> List[str]:
        """
        Get names of all features in order.
        
        Returns:
            List of feature names
        """
        names = []
        
        # Descriptor names
        names.extend(self.descriptor_names)
        
        # Morgan fingerprint bits
        names.extend([f'morgan_bit_{i}' for i in range(self.morgan_bits)])
        
        # MACCS keys
        names.extend([f'maccs_{i}' for i in range(167)])
        
        # Pharmacophore features
        names.extend([
            'has_basic_nitrogen', 'has_aromatic_n', 'num_benzene_rings',
            'has_halogen', 'num_hydroxyl', 'num_carbonyl', 'num_amine',
            'has_sulfur', 'aromatic_ratio', 'lipinski_violations'
        ])
        
        return names
    
    def extract_batch(self, smiles_list: List[str]) -> np.ndarray:
        """
        Extract features for a batch of SMILES.
        
        Args:
            smiles_list: List of SMILES strings
            
        Returns:
            2D numpy array of shape (n_samples, n_features)
        """
        features = []
        for smiles in smiles_list:
            feat = self.extract_all_features(smiles)
            if feat is not None:
                features.append(feat)
            else:
                # Handle invalid SMILES with zero vector
                features.append(np.zeros(self.get_num_features()))
        
        return np.array(features, dtype=np.float32)
    
    def get_num_features(self) -> int:
        """Get total number of features."""
        return (
            20 +  # Descriptors
            self.morgan_bits +  # Morgan fingerprint
            167 +  # MACCS keys
            10  # Pharmacophore features
        )


# Simple test function
def test_feature_extraction():
    """Test feature extraction with example molecules."""
    extractor = MolecularFeatureExtractor()
    
    test_molecules = {
        'Aspirin': 'CC(=O)Oc1ccccc1C(=O)O',
        'Caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
        'Warfarin': 'CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O'
    }
    
    print("Testing Feature Extraction")
    print("=" * 60)
    
    for name, smiles in test_molecules.items():
        features = extractor.extract_all_features(smiles)
        if features is not None:
            print(f"\n{name}:")
            print(f"  SMILES: {smiles}")
            print(f"  Total features: {len(features)}")
            print(f"  Non-zero features: {np.count_nonzero(features)}")
            print(f"  Feature range: [{features.min():.2f}, {features.max():.2f}]")
        else:
            print(f"\n{name}: FAILED to extract features")
    
    print(f"\nFeature names (first 10): {extractor.get_feature_names()[:10]}")
    print(f"Total feature count: {extractor.get_num_features()}")


if __name__ == '__main__':
    test_feature_extraction()

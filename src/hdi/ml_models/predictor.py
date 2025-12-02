#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Simple CYP450 predictor using rule-based and ML-hybrid approach.
This is a prototype that combines structural alerts with lightweight ML.
"""

import numpy as np
from typing import Dict, Optional
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen
import pickle
import os


class SimpleCYP450Predictor:
    """
    Simple CYP450 inhibition/induction predictor.
    Uses rule-based heuristics until ML models are trained.
    """
    
    def __init__(self, model_dir: Optional[str] = None):
        """
        Initialize predictor.
        
        Args:
            model_dir: Directory containing trained models (optional)
        """
        self.model_dir = model_dir or 'src/hdi/ml_models/trained'
        self.models = {}
        
        # Try to load trained models if they exist
        self._load_models()
    
    def _load_models(self):
        """Load trained models from disk if available."""
        if not os.path.exists(self.model_dir):
            print(f"Model directory {self.model_dir} not found. Using rule-based predictions.")
            return
        
        # Try to load CYP3A4 models
        try:
            inh_path = f"{self.model_dir}/CYP3A4_inhibition.pkl"
            if os.path.exists(inh_path):
                with open(inh_path, 'rb') as f:
                    data = pickle.load(f)
                    self.models['CYP3A4_inhibition'] = data
                print(f"✓ Loaded CYP3A4 inhibition model")
            
            ind_path = f"{self.model_dir}/CYP3A4_induction.pkl"
            if os.path.exists(ind_path):
                with open(ind_path, 'rb') as f:
                    data = pickle.load(f)
                    self.models['CYP3A4_induction'] = data
                print(f"✓ Loaded CYP3A4 induction model")
                
        except Exception as e:
            print(f"Error loading models: {e}")
    
    def _predict_with_model(self, mol, enzyme: str, task: str = 'inhibition') -> Dict:
        """
        Predict using trained ML model.
        
        Args:
            mol: RDKit molecule
            enzyme: Enzyme name
            task: 'inhibition' or 'induction'
            
        Returns:
            Dictionary with prediction results
        """
        model_key = f"{enzyme}_{task}"
        
        if model_key not in self.models:
            return {'error': 'Model not loaded', 'probability': 0, 'risk': 'UNKNOWN'}
        
        try:
            from .features import MolecularFeatureExtractor
            
            # Extract features
            extractor = MolecularFeatureExtractor()
            smiles = Chem.MolToSmiles(mol)
            features = extractor.extract_all_features(smiles)
            
            if features is None:
                return {'error': 'Feature extraction failed', 'probability': 0, 'risk': 'UNKNOWN'}
            
            # Load model and scaler
            model_data = self.models[model_key]
            model = model_data['model']
            scaler = model_data['scaler']
            
            # Scale features
            features_scaled = scaler.transform(features.reshape(1, -1))
            
            # Predict
            prob = model.predict_proba(features_scaled)[0, 1]
            
            # Determine risk
            if prob > 0.7:
                risk = 'HIGH'
            elif prob > 0.4:
                risk = 'MEDIUM'
            else:
                risk = 'LOW'
            
            return {
                'probability': float(prob),
                'risk': risk,
                'confidence': 0.9,  # High confidence for ML models
                'method': 'ML-trained'
            }
            
        except Exception as e:
            print(f"Error in ML prediction: {e}")
            return {'error': str(e), 'probability': 0, 'risk': 'UNKNOWN'}
    
    def predict_inhibition(self, smiles: str, enzyme: str = 'CYP3A4') -> Dict:
        """
        Predict CYP450 inhibition probability.
        
        Args:
            smiles: SMILES string
            enzyme: Enzyme name
            
        Returns:
            Dictionary with prediction results
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES', 'probability': 0, 'risk': 'UNKNOWN'}
        
        # If trained model exists, use it
        model_key = f"{enzyme}_inhibition"
        if model_key in self.models:
            return self._predict_with_model(mol, enzyme, 'inhibition')
        
        # Otherwise, use rule-based heuristics
        return self._predict_with_rules_inhibition(mol, enzyme)
    
    def predict_induction(self, smiles: str, enzyme: str = 'CYP3A4') -> Dict:
        """
        Predict CYP450 induction probability.
        
        Args:
            smiles: SMILES string
            enzyme: Enzyme name
            
        Returns:
            Dictionary with prediction results
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES', 'probability': 0, 'risk': 'UNKNOWN'}
        
        # If trained model exists, use it
        model_key = f"{enzyme}_induction"
        if model_key in self.models:
            result = self._predict_with_model(mol, enzyme, 'induction')
            # Add mechanism for induction
            if result['probability'] > 0.5:
                result['mechanism'] = 'PXR/CAR activation (ML predicted)'
            return result
        
        # Use rule-based heuristics
        return self._predict_with_rules_induction(mol, enzyme)
    
    def predict_induction(self, smiles: str, enzyme: str = 'CYP3A4') -> Dict:
        """
        Predict CYP450 induction probability.
        
        Args:
            smiles: SMILES string
            enzyme: Enzyme name
            
        Returns:
            Dictionary with prediction results
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES', 'probability': 0, 'risk': 'UNKNOWN'}
        
        # Use rule-based heuristics for induction
        return self._predict_with_rules_induction(mol, enzyme)
    
    def _predict_with_rules_inhibition(self, mol, enzyme: str) -> Dict:
        """
        Rule-based CYP450 inhibition prediction.
        Based on known pharmacophores and molecular properties.
        """
        score = 0
        features = []
        
        # Calculate molecular properties
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        aromatic_rings = Lipinski.NumAromaticRings(mol)
        
        # Rule 1: Size and lipophilicity (common in CYP inhibitors)
        if mw > 300 and logp > 2:
            score += 0.2
            features.append('Large lipophilic molecule')
        
        # Rule 2: Multiple aromatic rings
        if aromatic_rings >= 2:
            score += 0.15 * aromatic_rings
            features.append(f'{aromatic_rings} aromatic rings')
        
        # Rule 3: Nitrogen-containing heterocycles (common in CYP inhibitors)
        if self._has_nitrogen_heterocycle(mol):
            score += 0.25
            features.append('Nitrogen heterocycle')
        
        # Rule 4: Azole groups (strong CYP3A4 inhibitors)
        if self._has_azole(mol):
            score += 0.3
            features.append('Azole group (strong inhibitor)')
        
        # Enzyme-specific rules
        if enzyme == 'CYP3A4':
            # CYP3A4 prefers larger, more lipophilic molecules
            if mw > 400 and logp > 3:
                score += 0.15
        elif enzyme == 'CYP2D6':
            # CYP2D6 often inhibited by basic amines
            if self._has_basic_amine(mol):
                score += 0.2
                features.append('Basic amine')
        
        # Cap score at 1.0
        probability = min(score, 1.0)
        
        # Determine risk level
        if probability > 0.6:
            risk = 'HIGH'
        elif probability > 0.3:
            risk = 'MEDIUM'
        else:
            risk = 'LOW'
        
        return {
            'probability': probability,
            'risk': risk,
            'confidence': 0.6,  # Lower confidence for rule-based
            'features': features,
            'method': 'rule-based'
        }
    
    def _predict_with_rules_induction(self, mol, enzyme: str) -> Dict:
        """
        Rule-based CYP450 induction prediction.
        Focuses on PXR/CAR activators.
        """
        score = 0
        features = []
        
        # Calculate properties
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        aromatic_rings = Lipinski.NumAromaticRings(mol)
        
        # Rule 1: Large, lipophilic molecules (PXR activators)
        if mw > 400 and logp > 4:
            score += 0.35
            features.append('Large lipophilic (PXR activator profile)')
        
        # Rule 2: Multiple aromatic systems
        if aromatic_rings >= 3:
            score += 0.25
            features.append('Polyaromatic system')
        
        # Rule 3: Steroid-like scaffold
        if self._has_steroid_scaffold(mol):
            score += 0.3
            features.append('Steroid-like scaffold')
        
        # Rule 4: Phenobarbital-like (CAR activator)
        if self._has_barbiturate(mol):
            score += 0.35
            features.append('Barbiturate (CAR activator)')
        
        # Enzyme-specific
        if enzyme == 'CYP3A4':
            # CYP3A4 is the primary target of PXR
            score *= 1.2  # Boost score for CYP3A4
        
        probability = min(score, 1.0)
        
        # Determine risk
        if probability > 0.5:
            risk = 'HIGH'
        elif probability > 0.25:
            risk = 'MEDIUM'
        else:
            risk = 'LOW'
        
        return {
            'probability': probability,
            'risk': risk,
            'confidence': 0.5,  # Lower confidence for rule-based
            'features': features,
            'method': 'rule-based',
            'mechanism': 'PXR/CAR activation' if probability > 0.3 else 'Unknown'
        }
    
    # Helper methods for structural pattern detection
    
    def _has_nitrogen_heterocycle(self, mol) -> bool:
        """Check for nitrogen-containing heterocycles."""
        pattern = Chem.MolFromSmarts('[n;r]')  # Aromatic nitrogen in ring
        return mol.HasSubstructMatch(pattern) if pattern else False
    
    def _has_azole(self, mol) -> bool:
        """Check for azole groups (imidazole, triazole, etc.)."""
        patterns = [
            'c1ncnc1',  # Imidazole
            'c1nncn1',  # Triazole
        ]
        for smarts in patterns:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                return True
        return False
    
    def _has_basic_amine(self, mol) -> bool:
        """Check for basic amine groups."""
        pattern = Chem.MolFromSmarts('[N;!$(N-C=O);!$(N-S(=O)=O)]')
        return mol.HasSubstructMatch(pattern) if pattern else False
    
    def _has_steroid_scaffold(self, mol) -> bool:
        """Check for steroid-like four-ring scaffold."""
        # Simplified steroid pattern
        pattern = Chem.MolFromSmarts('C1CC2CCC3CCCC3C2C1')
        return mol.HasSubstructMatch(pattern) if pattern else False
    
    def _has_barbiturate(self, mol) -> bool:
        """Check for barbiturate core."""
        pattern = Chem.MolFromSmarts('C(=O)NC(=O)NC(=O)')
        return mol.HasSubstructMatch(pattern) if pattern else False


# Test function
def test_predictor():
    """Test the predictor with known molecules."""
    predictor = SimpleCYP450Predictor()
    
    test_cases = {
        'Ketoconazole (strong CYP3A4 inhibitor)': 'CC1COC(O1)(CN2C=NC=N2)C3=C(C=C(C=C3)Cl)Cl',
        'Rifampin (CYP3A4 inducer)': 'CC1C=CC=C(C(=O)NC2C(C(C3C(C2O)C(C(C4C3CC(O4)C5C(C(C(C(C5O)C)O)C)OC)C)O)O)O)C',
        'Aspirin (weak interaction)': 'CC(=O)Oc1ccccc1C(=O)O',
    }
    
    print("Testing Simple CYP450 Predictor")
    print("=" * 70)
    
    for name, smiles in test_cases.items():
        print(f"\n{name}:")
        print(f"  SMILES: {smiles[:60]}...")
        
        # Test inhibition
        inh_result = predictor.predict_inhibition(smiles, 'CYP3A4')
        print(f"  Inhibition: {inh_result['risk']} (prob: {inh_result['probability']:.2f})")
        if inh_result.get('features'):
            print(f"    Features: {', '.join(inh_result['features'][:3])}")
        
        # Test induction
        ind_result = predictor.predict_induction(smiles, 'CYP3A4')
        print(f"  Induction: {ind_result['risk']} (prob: {ind_result['probability']:.2f})")
        if ind_result.get('features'):
            print(f"    Features: {', '.join(ind_result['features'][:3])}")


if __name__ == '__main__':
    test_predictor()

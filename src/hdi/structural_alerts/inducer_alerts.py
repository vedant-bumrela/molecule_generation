#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Structural alerts for CYP450 inducers.
Detects molecular patterns known to induce CYP450 enzyme expression.
"""

from typing import List, Dict
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

# SMARTS patterns for known CYP450 inducers
# Based on literature and known pharmacophores

INDUCER_SMARTS = {
    # Multiple aromatic rings (common in CYP inducers)
    'multiple_aromatic_rings': 'c1ccccc1-c2ccccc2',  # Biphenyl
    'polycyclic_aromatic': 'c1ccc2ccccc2c1',  # Naphthalene-like
    
    # Large hydrophobic scaffold
    'large_aromatic_system': 'c1cc2ccc3ccccc3cc2cc1',  # Anthracene
    
    # Steroid-like (simplified)
    'steroid_scaffold': 'C1CC2CCC3CCCC3C2C1',  # Simplified steroid core
    
    # Phenobarbital-like
    'barbiturate': 'C(=O)NC(=O)NC(=O)',  # Barbiturate core
}

# Known inducer scaffolds (as SMILES for exact matching)
KNOWN_INDUCER_SMILES = {
    'rifampin': 'CC1C=CC=C(C(=O)NC2=C(C(=C3C(=C2O)C(=C(C4=C3C(=O)C(O4)(OC=CC(C(C(C(C(C(C1O)C)O)C)OC(=O)C)C)OC)C)C)O)O)C)O',
    'hyperforin': 'CC1=C2C(=C(C(=C1O)C(C)CC(C)(C)C)C)C(=O)C3=C(C2=O)C(=C(C(=C3)O)C(C)CC(C)(C)C)O',
    'phenobarbital': 'CCC1(C(=O)NC(=O)NC1=O)C2=CC=CC=C2',
    'carbamazepine': 'C1=CC=C(C=C1)C2=CC=CC=C2C(=O)N',
}


class CYP450InducerAlerts:
    """Detect structural features that indicate CYP450 induction potential."""
    
    def __init__(self):
        """Initialize with compiled SMARTS patterns."""
        self.compiled_patterns = {}
        for name, smarts in INDUCER_SMARTS.items():
            try:
                self.compiled_patterns[name] = Chem.MolFromSmarts(smarts)
            except:
                print(f"Warning: Could not compile SMARTS pattern {name}: {smarts}")
    
    def scan_molecule(self, smiles: str) -> Dict:
        """
        Scan molecule for CYP450 inducer structural alerts.
        
        Args:
            smiles: SMILES string of compound
            
        Returns:
            Dictionary with alerts and risk assessment
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES', 'alerts': [], 'risk': 'UNKNOWN'}
        
        alerts = []
        
        # Check SMARTS patterns
        for pattern_name, pattern_mol in self.compiled_patterns.items():
            if mol.HasSubstructMatch(pattern_mol):
                alerts.append({
                    'type': 'structural_pattern',
                    'pattern': pattern_name,
                    'description': self._get_pattern_description(pattern_name),
                    'severity': self._get_pattern_severity(pattern_name)
                })
        
        # Check molecular properties indicative of inducers
        property_alerts = self._check_physicochemical_properties(mol)
        alerts.extend(property_alerts)
        
        # Check against known inducer library
        similarity_alerts = self._check_known_inducers(smiles)
        alerts.extend(similarity_alerts)
        
        # Calculate overall induction risk
        risk_level = self._calculate_induction_risk(alerts)
        
        return {
            'alerts': alerts,
            'num_alerts': len(alerts),
            'risk': risk_level,
            'likely_mechanisms': self._infer_mechanisms(alerts),
            'affected_enzymes': self._predict_affected_enzymes(alerts)
        }
    
    def _check_physicochemical_properties(self, mol) -> List[Dict]:
        """Check if physicochemical properties match inducer profile."""
        alerts = []
        
        # Calculate properties
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        aromatic_rings = Lipinski.NumAromaticRings(mol)
        
        # PXR activators tend to be large and lipophilic
        if mw > 400 and logp > 4:
            alerts.append({
                'type': 'physicochemical',
                'pattern': 'large_lipophilic',
                'description': 'Large molecular weight and high lipophilicity (typical of PXR activators)',
                'severity': 'MEDIUM',
                'values': {'MW': mw, 'LogP': logp}
            })
        
        # Multiple aromatic rings
        if aromatic_rings >= 3:
            alerts.append({
                'type': 'physicochemical',
                'pattern': 'polyaromatic',
                'description': f'Multiple aromatic rings ({aromatic_rings}) - common in CYP inducers',
                'severity': 'LOW'
            })
        
        return alerts
    
    def _check_known_inducers(self, smiles: str) -> List[Dict]:
        """Check structural similarity to known inducers."""
        from rdkit import DataStructs
        from rdkit.Chem import AllChem
        
        alerts = []
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return alerts
        
        fp_query = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        
        for name, known_smiles in KNOWN_INDUCER_SMILES.items():
            known_mol = Chem.MolFromSmiles(known_smiles)
            if known_mol is None:
                continue
            
            fp_known = AllChem.GetMorganFingerprintAsBitVect(known_mol, 2, nBits=2048)
            similarity = DataStructs.TanimotoSimilarity(fp_query, fp_known)
            
            if similarity > 0.7:  # High similarity
                alerts.append({
                    'type': 'similarity',
                    'pattern': f'similar_to_{name}',
                    'description': f'Structurally similar to {name} (similarity: {similarity:.2f})',
                    'severity': 'HIGH',
                    'similarity': similarity
                })
            elif similarity > 0.5:  # Moderate similarity
                alerts.append({
                    'type': 'similarity',
                    'pattern': f'moderately_similar_to_{name}',
                    'description': f'Moderately similar to known inducer {name} (similarity: {similarity:.2f})',
                    'severity': 'MEDIUM',
                    'similarity': similarity
                })
        
        return alerts
    
    def _calculate_induction_risk(self, alerts: List[Dict]) -> str:
        """Calculate overall CYP450 induction risk from alerts."""
        if not alerts:
            return 'LOW'
        
        # Count severity levels
        high_severity = sum(1 for a in alerts if a.get('severity') == 'HIGH')
        medium_severity = sum(1 for a in alerts if a.get('severity') == 'MEDIUM')
        
        if high_severity >= 2 or (high_severity >= 1 and medium_severity >= 2):
            return 'HIGH'
        elif high_severity >= 1 or medium_severity >= 2:
            return 'MEDIUM'
        else:
            return 'LOW'
    
    def _infer_mechanisms(self, alerts: List[Dict]) -> List[str]:
        """Infer likely induction mechanisms from alerts."""
        mechanisms = set()
        
        for alert in alerts:
            pattern = alert.get('pattern', '')
            if 'PXR' in pattern or 'steroid' in pattern or 'large_lipophilic' in pattern:
                mechanisms.add('PXR activation (CYP3A4 induction)')
            if 'CAR' in pattern or 'phenobarbital' in pattern:
                mechanisms.add('CAR activation (CYP2B6, CYP3A4 induction)')
            if 'hyperforin' in pattern or 'similar_to_hyperforin' in pattern:
                mechanisms.add('PXR/CAR activation (CYP3A4, CYP2C9 induction)')
            if 'rifampin' in pattern or 'similar_to_rifampin' in pattern:
                mechanisms.add('Strong PXR activation (multiple CYP450s)')
        
        return list(mechanisms) if mechanisms else ['Unknown mechanism']
    
    def _predict_affected_enzymes(self, alerts: List[Dict]) -> List[str]:
        """Predict which CYP450 enzymes are likely to be induced."""
        enzymes = set()
        
        for alert in alerts:
            pattern = alert.get('pattern', '')
            # Most inducers affect CYP3A4
            if any(x in pattern for x in ['PXR', 'steroid', 'large_lipophilic', 'hyperforin', 'rifampin']):
                enzymes.add('CYP3A4')
            # Some also affect CYP2C9
            if any(x in pattern for x in ['hyperforin', 'rifampin', 'CAR']):
                enzymes.add('CYP2C9')
            # CAR activators induce CYP2B6
            if 'CAR' in pattern or 'phenobarbital' in pattern:
                enzymes.add('CYP2B6')
        
        return list(enzymes) if enzymes else ['CYP3A4']  # Default to CYP3A4
    
    def _get_pattern_description(self, pattern_name: str) -> str:
        """Get human-readable description of pattern."""
        descriptions = {
            'PXR_activator_steroid': 'Steroid-like scaffold (PXR activator)',
            'PXR_activator_large_hydrophobic': 'Large hydrophobic core structure',
            'hyperforin_scaffold': 'Hyperforin-like core (St. John\'s Wort active compound)',
            'CAR_activator_phenobarbital': 'Phenobarbital-like structure',
            'rifampin_like': 'Rifampin-like pharmacophore',
            'multiple_aromatic_rings': 'Multiple aromatic ring system',
            'polycyclic_aromatic': 'Polycyclic aromatic hydrocarbon'
        }
        return descriptions.get(pattern_name, pattern_name)
    
    def _get_pattern_severity(self, pattern_name: str) -> str:
        """Get severity level for pattern."""
        high_severity = ['hyperforin_scaffold', 'rifampin_like', 'PXR_activator_steroid']
        medium_severity = ['CAR_activator_phenobarbital', 'PXR_activator_large_hydrophobic']
        
        if pattern_name in high_severity:
            return 'HIGH'
        elif pattern_name in medium_severity:
            return 'MEDIUM'
        else:
            return 'LOW'

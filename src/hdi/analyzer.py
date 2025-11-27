#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
HDI (Herb-Drug Interaction) Analyzer.
Core analysis logic for predicting herb-drug interactions through molecular docking.
"""

import os
import json
from typing import Dict, List, Optional, Tuple
from pathlib import Path

from ..preprocessing.ligand_prep import LigandPreparation
from ..docking.vina_docking import VinaDocking
from .enzyme_structures import EnzymeStructureManager, CYP450_ENZYMES


class HDIAnalyzer:
    """Analyze herb-drug interactions using molecular docking against CYP450 enzymes."""
    
    # Risk thresholds based on binding affinity (kcal/mol)
    HIGH_RISK_THRESHOLD = -8.0    # Strong binding
    MEDIUM_RISK_THRESHOLD = -6.0  # Moderate binding
    
    def __init__(self):
        """Initialize HDI analyzer."""
        self.enzyme_manager = EnzymeStructureManager()
        self.ligand_prep = LigandPreparation()
        self.vina = VinaDocking()
    
    def analyze_hdi(
        self,
        herb_smiles: str,
        drug_smiles: str,
        herb_name: str = "Unknown Herb",
        drug_name: str = "Unknown Drug",
        output_dir: str = None
    ) -> Dict:
        """
        Analyze herb-drug interaction.
        
        Args:
            herb_smiles: SMILES string of herbal compound
            drug_smiles: SMILES string of pharmaceutical drug
            herb_name: Name of herb (for reporting)
            drug_name: Name of drug (for reporting)
            output_dir: Directory to save analysis results
            
        Returns:
            Dictionary containing analysis results
        """
        if output_dir is None:
            output_dir = Path('uploads') / 'hdi_analysis' / 'temp'
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Prepare ligands from SMILES
        herb_ligand_path = output_dir / 'herb_ligand.pdbqt'
        drug_ligand_path = output_dir / 'drug_ligand.pdbqt'
        
        herb_ligand = self.ligand_prep.prepare_ligand_from_smiles(
            herb_smiles,
            output_path=str(herb_ligand_path)
        )
        
        drug_ligand = self.ligand_prep.prepare_ligand_from_smiles(
            drug_smiles,
            output_path=str(drug_ligand_path)
        )
        
        # Dock against CYP450 enzyme panel
        herb_results = self._dock_against_cyp450_panel(herb_ligand, output_dir / 'herb_docking')
        drug_results = self._dock_against_cyp450_panel(drug_ligand, output_dir / 'drug_docking')
        
        # Calculate interaction risk
        interaction_analysis = self._analyze_competition(herb_results, drug_results)
        
        # Assess clinical significance
        clinical_assessment = self._assess_clinical_significance(
            interaction_analysis,
            herb_name,
            drug_name
        )
        
        # Compile results
        results = {
            'herb_info': {
                'name': herb_name,
                'smiles': herb_smiles
            },
            'drug_info': {
                'name': drug_name,
                'smiles': drug_smiles
            },
            'enzyme_binding': {
                'herb': herb_results,
                'drug': drug_results
            },
            'interaction_risk': interaction_analysis['overall_risk'],
            'affected_enzymes': interaction_analysis['affected_enzymes'],
            'clinical_significance': clinical_assessment['significance'],
            'mechanism': clinical_assessment['mechanism'],
            'recommendation': clinical_assessment['recommendation']
        }
        
        # Save results
        results_file = output_dir / 'results.json'
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        return results
    
    def _dock_against_cyp450_panel(
        self,
        ligand_path: str,
        output_dir: Path
    ) -> Dict[str, Dict]:
        """
        Dock ligand against all available CYP450 enzymes.
        
        Args:
            ligand_path: Path to prepared ligand (PDBQT)
            output_dir: Directory to save docking results
            
        Returns:
            Dictionary mapping enzyme names to docking results
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        results = {}
        
        # Get available enzymes
        available_enzymes = self.enzyme_manager.get_available_enzymes()
        
        if not available_enzymes:
            # Try to set up enzymes
            self.enzyme_manager.setup_all_enzymes()
            available_enzymes = self.enzyme_manager.get_available_enzymes()
        
        for enzyme_name in available_enzymes:
            enzyme_path = self.enzyme_manager.get_enzyme_structure(enzyme_name)
            
            if enzyme_path is None:
                continue
            
            # Create enzyme-specific output directory
            enzyme_output_dir = output_dir / enzyme_name.lower()
            enzyme_output_dir.mkdir(parents=True, exist_ok=True)
            
            output_file = enzyme_output_dir / 'docking_result.pdbqt'
            
            try:
                # Auto-detect binding site
                binding_site = self.vina.detect_binding_site(enzyme_path)
                
                binding_site_params = {
                    'center_x': binding_site['center_x'],
                    'center_y': binding_site['center_y'],
                    'center_z': binding_site['center_z'],
                    'size_x': binding_site['size_x'],
                    'size_y': binding_site['size_y'],
                    'size_z': binding_site['size_z']
                }
                
                # Run docking
                result_path = self.vina.run_docking(
                    enzyme_path,
                    ligand_path,
                    output_path=str(output_file),
                    binding_site=binding_site_params,
                    exhaustiveness=8,
                    num_modes=3  # Fewer modes for speed
                )
                
                # Load scores
                scores_file = str(output_file).replace('.pdbqt', '.json')
                if os.path.exists(scores_file):
                    with open(scores_file, 'r') as f:
                        scores_data = json.load(f)
                        scores = scores_data.get('scores', [])
                else:
                    scores = []
                
                results[enzyme_name] = {
                    'success': True,
                    'output_path': str(result_path),
                    'scores': scores,
                    'best_affinity': scores[0]['affinity'] if scores else None,
                    'enzyme_info': CYP450_ENZYMES.get(enzyme_name, {})
                }
                
            except Exception as e:
                results[enzyme_name] = {
                    'success': False,
                    'error': str(e)
                }
        
        return results
    
    def _analyze_competition(
        self,
        herb_results: Dict[str, Dict],
        drug_results: Dict[str, Dict]
    ) -> Dict:
        """
        Analyze competitive binding between herb and drug.
        
        Args:
            herb_results: Herb docking results
            drug_results: Drug docking results
            
        Returns:
            Dictionary with interaction analysis
        """
        affected_enzymes = []
        max_risk_level = "LOW"
        
        for enzyme_name in herb_results.keys():
            if enzyme_name not in drug_results:
                continue
            
            herb_data = herb_results[enzyme_name]
            drug_data = drug_results[enzyme_name]
            
            if not herb_data.get('success') or not drug_data.get('success'):
                continue
            
            herb_affinity = herb_data.get('best_affinity')
            drug_affinity = drug_data.get('best_affinity')
            
            if herb_affinity is None or drug_affinity is None:
                continue
            
            # Determine if there's competitive binding
            # Strong herb binding indicates potential interference
            risk_level = self._calculate_risk_level(herb_affinity)
            
            # Check if drug also binds (would be affected)
            drug_binds = drug_affinity < self.MEDIUM_RISK_THRESHOLD
            
            if drug_binds:
                affected_enzymes.append({
                    'name': enzyme_name,
                    'herb_affinity': herb_affinity,
                    'drug_affinity': drug_affinity,
                    'risk_level': risk_level,
                    'competition': True,
                    'enzyme_info': CYP450_ENZYMES.get(enzyme_name, {})
                })
                
                # Update overall risk
                if risk_level == "HIGH":
                    max_risk_level = "HIGH"
                elif risk_level == "MEDIUM" and max_risk_level != "HIGH":
                    max_risk_level = "MEDIUM"
        
        return {
            'overall_risk': max_risk_level,
            'affected_enzymes': affected_enzymes,
            'num_affected': len(affected_enzymes)
        }
    
    def _calculate_risk_level(self, affinity: float) -> str:
        """
        Calculate risk level from binding affinity.
        
        Args:
            affinity: Binding affinity in kcal/mol
            
        Returns:
            Risk level string ('HIGH', 'MEDIUM', 'LOW')
        """
        if affinity < self.HIGH_RISK_THRESHOLD:
            return "HIGH"
        elif affinity < self.MEDIUM_RISK_THRESHOLD:
            return "MEDIUM"
        else:
            return "LOW"
    
    def _assess_clinical_significance(
        self,
        interaction_analysis: Dict,
        herb_name: str,
        drug_name: str
    ) -> Dict:
        """
        Assess clinical significance of interaction.
        
        Args:
            interaction_analysis: Results from competition analysis
            herb_name: Name of herb
            drug_name: Name of drug
            
        Returns:
            Dictionary with clinical assessment
        """
        risk_level = interaction_analysis['overall_risk']
        affected_enzymes = interaction_analysis['affected_enzymes']
        
        if risk_level == "HIGH":
            significance = f"HIGH RISK: {herb_name} shows strong binding to CYP450 enzymes that metabolize {drug_name}. "
            significance += f"This may significantly alter {drug_name} metabolism, potentially leading to "
            significance += "either reduced efficacy (if enzyme is induced) or increased toxicity (if enzyme is inhibited)."
            
            mechanism = "Competitive inhibition and/or enzyme induction at CYP450 active sites"
            recommendation = f"AVOID concurrent use of {herb_name} and {drug_name}. "
            recommendation += "If unavoidable, monitor drug levels closely and adjust dosage as needed under medical supervision."
            
        elif risk_level == "MEDIUM":
            significance = f"MODERATE RISK: {herb_name} may interfere with {drug_name} metabolism through CYP450 enzymes. "
            significance += "Clinical monitoring is recommended."
            
            mechanism = "Potential competitive binding at CYP450 active sites"
            recommendation = f"Use caution when combining {herb_name} with {drug_name}. "
            recommendation += "Monitor for changes in drug efficacy or adverse effects. Consult healthcare provider."
            
        else:
            significance = f"LOW RISK: Minimal predicted interaction between {herb_name} and {drug_name} "
            significance += "based on CYP450 enzyme binding analysis."
            
            mechanism = "No significant competitive binding detected"
            recommendation = "No specific precautions required based on current analysis. However, always inform healthcare providers of all supplements and medications."
        
        # Add specific enzyme information
        if affected_enzymes:
            enzyme_names = [e['name'] for e in affected_enzymes]
            significance += f"\n\nAffected enzymes: {', '.join(enzyme_names)}"
        
        return {
            'significance': significance,
            'mechanism': mechanism,
            'recommendation': recommendation
        }

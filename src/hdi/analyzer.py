#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
HDI (Herb-Drug Interaction) Analyzer.
Enhanced multi-layer analysis for predicting herb-drug interactions.
"""

import os
import json
from typing import Dict, List, Optional, Tuple
from pathlib import Path

from ..preprocessing.ligand_prep import LigandPreparation
from ..docking.vina_docking import VinaDocking
from .enzyme_structures import EnzymeStructureManager, CYP450_ENZYMES
from .structural_alerts import CYP450InducerAlerts
from .risk_aggregator import RiskAggregator
from .ml_models.predictor import SimpleCYP450Predictor


class HDIAnalyzer:
    """Analyze herb-drug interactions using multi-layer prediction system."""
    
    # Risk thresholds based on binding affinity (kcal/mol)
    HIGH_RISK_THRESHOLD = -8.0    # Strong binding
    MEDIUM_RISK_THRESHOLD = -6.0  # Moderate binding
    
    def __init__(self):
        """Initialize HDI analyzer with all prediction layers."""
        self.enzyme_manager = EnzymeStructureManager()
        self.ligand_prep = LigandPreparation()
        self.vina = VinaDocking()
        self.inducer_alerts = CYP450InducerAlerts()
        self.risk_aggregator = RiskAggregator()
        self.ml_predictor = SimpleCYP450Predictor()  # 🆕 ML predictor
    
    def analyze_hdi(
        self,
        herb_smiles: str,
        drug_smiles: str,
        herb_name: str = "Unknown Herb",
        drug_name: str = "Unknown Drug",
        output_dir: str = None
    ) -> Dict:
        """
        Analyze herb-drug interaction using multi-layer prediction.
        
        Args:
            herb_smiles: SMILES string of herbal compound
            drug_smiles: SMILES string of pharmaceutical drug
            herb_name: Name of herb (for reporting)
            drug_name: Name of drug (for reporting)
            output_dir: Directory to save analysis results
            
        Returns:
            Dictionary containing comprehensive analysis results
        """
        if output_dir is None:
            output_dir = Path('uploads') / 'hdi_analysis' / 'temp'
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # ==== LAYER 1: Structural Alerts ====
        print("🔍 Layer 1: Scanning structural alerts...")
        herb_alerts = self.inducer_alerts.scan_molecule(herb_smiles)
        drug_alerts = self.inducer_alerts.scan_molecule(drug_smiles)
        
        # ==== LAYER 2: ML Predictions ====
        print("🤖 Layer 2: ML-based predictions...")
        herb_ml_inhibition = self.ml_predictor.predict_inhibition(herb_smiles, 'CYP3A4')
        herb_ml_induction = self.ml_predictor.predict_induction(herb_smiles, 'CYP3A4')
        
        # ==== LAYER 3: Molecular Docking ====
        print("🔬 Layer 3: Running molecular docking...")
        
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
        herb_docking = self._dock_against_cyp450_panel(herb_ligand, output_dir / 'herb_docking')
        drug_docking = self._dock_against_cyp450_panel(drug_ligand, output_dir / 'drug_docking')
        
        # Analyze competitive binding
        docking_analysis = self._analyze_competition(herb_docking, drug_docking)
        
        # ==== LAYER 4: Risk Aggregation ====
        print("⚖️ Layer 4: Aggregating risk from all layers...")
        
        layer_results = {
            'docking': docking_analysis,
            'structural_alerts': herb_alerts,  # Herb alerts are primary concern
            'ml_inhibition': herb_ml_inhibition,  # 🆕 ML inhibition prediction
            'ml_induction': herb_ml_induction,     # 🆕 ML induction prediction
        }
        
        aggregated_risk = self.risk_aggregator.aggregate(layer_results)
        
        # ==== Enhanced Clinical Assessment ====
        clinical_assessment = self._assess_clinical_significance_enhanced(
            aggregated_risk,
            herb_alerts,
            drug_alerts,
            docking_analysis,
            herb_name,
            drug_name
        )
        
        # ==== Compile Enhanced Results ====
        results = {
            'herb_info': {
                'name': herb_name,
                'smiles': herb_smiles
            },
            'drug_info': {
                'name': drug_name,
                'smiles': drug_smiles
            },
            # Aggregated risk (from all layers)
            'interaction_risk': aggregated_risk['interaction_risk'],
            'confidence': aggregated_risk.get('confidence', 0.5),
            
            # Affected enzymes (combined from all sources)
            'affected_enzymes': self._combine_affected_enzymes(
                docking_analysis.get('affected_enzymes', []),
                herb_alerts.get('affected_enzymes', [])
            ),
            
            # Mechanisms (from all layers)
            'mechanisms': aggregated_risk.get('mechanisms', []),
            
            # Clinical assessment
            'clinical_significance': clinical_assessment['significance'],
            'mechanism': clinical_assessment['mechanism'],
            'recommendation': clinical_assessment['recommendation'],
            
            # Detailed layer results
            'analysis_layers': {
                'structural_alerts': {
                    'herb': herb_alerts,
                    'drug': drug_alerts
                },
                'docking': {
                    'herb': herb_docking,
                    'drug': drug_docking,
                    'competition_analysis': docking_analysis
                },
                'risk_breakdown': aggregated_risk.get('breakdown', {})
            },
            
            # Metadata
            'evidence_layers': aggregated_risk.get('evidence_layers', []),
            'layer_count': aggregated_risk.get('layer_count', 1)
        }
        
        # Save results
        results_file = output_dir / 'results.json'
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"✅ Analysis complete! Risk: {results['interaction_risk']} (Confidence: {results['confidence']:.0%})")
        
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
    
    def _combine_affected_enzymes(
        self, 
        docking_enzymes: List[Dict], 
        alert_enzymes: List[str]
    ) -> List[Dict]:
        """
        Combine affected enzymes from docking and structural alerts.
        
        Args:
            docking_enzymes: List of enzyme dicts from docking analysis
            alert_enzymes: List of enzyme names from structural alerts
            
        Returns:
            Combined list of affected enzymes with enhanced information
        """
        # Create dict for easy lookup
        combined = {}
        
        # Add docking results
        for enzyme in docking_enzymes:
            name = enzyme.get('name')
            if name:
                combined[name] = enzyme.copy()
        
        # Add alert-predicted enzymes
        for enzyme_name in alert_enzymes:
            if enzyme_name in combined:
                # Mark that it was also predicted by alerts
                combined[enzyme_name]['predicted_by_alerts'] = True
                combined[enzyme_name]['risk_level'] = 'HIGH'  # Upgrade risk if predicted by both
            else:
                # Add as new enzyme (from induction prediction)
                combined[enzyme_name] = {
                    'name': enzyme_name,
                    'risk_level': 'HIGH',
                    'mechanism': 'Predicted enzyme induction', 
                    'predicted_by_alerts': True,
                    'herb_affinity': None,
                    'drug_affinity': None
                }
        
        return list(combined.values())
    
    def _assess_clinical_significance_enhanced(
        self,
        aggregated_risk: Dict,
        herb_alerts: Dict,
        drug_alerts: Dict,
        docking_analysis: Dict,
        herb_name: str,
        drug_name: str
    ) -> Dict:
        """
        Enhanced clinical significance assessment using all analysis layers.
        
        Args:
            aggregated_risk: Aggregated risk from all layers
            herb_alerts: Structural alerts for herb
            drug_alerts: Structural alerts for drug  
            docking_analysis: Docking competition analysis
            herb_name: Name of herb
            drug_name: Name of drug
            
        Returns:
            Dictionary with clinical assessment
        """
        risk_level = aggregated_risk['interaction_risk']
        mechanisms = aggregated_risk.get('mechanisms', [])
        confidence = aggregated_risk.get('confidence', 0.5)
        
        # Build significance message
        if risk_level == "HIGH":
            significance = f"⚠️ HIGH RISK: {herb_name} shows strong potential for interaction with {drug_name}. "
            
            # Add specific mechanisms
            if any('induction' in m.lower() for m in mechanisms):
                significance += f"The herbal compound may induce CYP450 enzymes, potentially reducing {drug_name} effectiveness. "
            if any('inhibition' in m.lower() or 'competitive' in m.lower() for m in mechanisms):
                significance += f"Competitive enzyme inhibition may increase {drug_name} levels, risking toxicity. "
            
            # Add structural alert details
            if herb_alerts.get('num_alerts', 0) > 0:
                significance += f"\n\n🔬 Structural Analysis: Detected {herb_alerts['num_alerts']} inducer/inhibitor patterns. "
                top_alerts = herb_alerts.get('alerts', [])[:3]
                for alert in top_alerts:
                    significance += f"\n  • {alert.get('description', '')}"
            
            mechanism = " | ".join(mechanisms[:3]) if mechanisms else "Multiple metabolic pathways"
            
            recommendation = f"⛔ AVOID concurrent use of {herb_name} and {drug_name}. "
            recommendation += f"If unavoidable, close monitoring of {drug_name} levels and clinical response is essential. "
            recommendation += "Dose adjustments may be required under medical supervision."
            
        elif risk_level == "MEDIUM":
            significance = f"⚡ MODERATE RISK: {herb_name} may interfere with {drug_name} metabolism. "
            significance += "Clinical monitoring recommended."
            
            if mechanisms:
                significance += f"\n\nPotential mechanisms: {', '.join(mechanisms[:2])}"
            
            mechanism = " | ".join(mechanisms[:2]) if mechanisms else "CYP450 interference"
            
            recommendation = f"⚠️ Use caution when combining {herb_name} with {drug_name}. "
            recommendation += "Monitor for changes in drug efficacy or adverse effects. Consult healthcare provider before use."
            
        else:
            significance = f"✓ LOW RISK: Minimal predicted interaction between {herb_name} and {drug_name} "
            significance += "based on current analysis."
            
            if confidence < 0.5:
                significance += "\n\nNote: Limited data available. Exercise caution and inform healthcare providers of all supplements."
            
            mechanism = "No significant metabolic interference detected"
            recommendation = "No specific precautions required based on current analysis. However, always inform healthcare providers of all supplements and medications."
        
        # Add confidence level disclaimer
        confidence_text = "high" if confidence > 0.8 else "moderate" if confidence > 0.5 else "limited"
        significance += f"\n\n📊 Confidence: {confidence_text} ({confidence:.0%}) | Evidence from {aggregated_risk.get('layer_count', 1)} analysis layer(s)"
        
        return {
            'significance': significance,
            'mechanism': mechanism,
            'recommendation': recommendation,
            'confidence': confidence
        }

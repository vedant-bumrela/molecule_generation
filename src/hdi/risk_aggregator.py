#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Risk Aggregation Engine for HDI Analysis.
Combines multiple prediction layers into a final risk assessment.
"""

from typing import Dict, List, Tuple
from enum import Enum


class RiskLevel(Enum):
    """Risk level enum for consistency."""
    LOW = "LOW"
    MEDIUM = "MEDIUM"
    HIGH = "HIGH"
    UNKNOWN = "UNKNOWN"


class RiskAggregator:
    """Aggregate risk scores from multiple analysis layers."""
    
    # Weight for each prediction layer
    LAYER_WEIGHTS = {
        'clinical_evidence': 1.0,      # Highest priority - known interactions
        'structural_alerts': 0.85,     # Very reliable - pattern matching
        'ml_induction': 0.80,          # ML model for induction (when available)
        'ml_inhibition': 0.70,         # ML model for inhibition (when available)
        'docking': 0.60,               # Molecular docking
        'pk_prediction': 0.40,         # PK parameter prediction
    }
    
    def __init__(self):
        """Initialize risk aggregator."""
        pass
    
    def aggregate(self, layer_results: Dict) -> Dict:
        """
        Aggregate risk from multiple analysis layers.
        
        Args:
            layer_results: Dictionary containing results from each layer:
                - docking: {overall_risk: str, affected_enzymes: List}
                - structural_alerts: {risk: str, alerts: List}
                - clinical_evidence: {is_known: bool, severity: str} (optional)
                - ml_induction: {risk: str, confidence: float} (optional)
                - ml_inhibition: {risk: str, confidence: float} (optional)
        
        Returns:
            Dictionary with aggregated risk assessment
        """
        
        # If we have clinical evidence, prioritize it
        clinical = layer_results.get('clinical_evidence', {})
        if clinical.get('is_known') and clinical.get('severity'):
            return self._create_clinical_evidence_result(clinical)
        
        # Collect risk scores from all layers
        risk_scores = []
        mechanisms = []
        affected_enzymes = set()
        confidence_factors = []
        
        # Layer 1: Structural Alerts (HIGH PRIORITY)
        alerts_result = layer_results.get('structural_alerts', {})
        if alerts_result and alerts_result.get('risk') != 'UNKNOWN':
            risk_level = alerts_result['risk']
            weight = self.LAYER_WEIGHTS['structural_alerts']
            
            # Induction alerts are critical!
            if alerts_result.get('alerts'):
                risk_scores.append((risk_level, weight, 'Structural alerts'))
                mechanisms.extend(alerts_result.get('likely_mechanisms', []))
                affected_enzymes.update(alerts_result.get('affected_enzymes', []))
                confidence_factors.append(('Structural alerts', weight))
        
        # Layer 2: Molecular Docking
        docking_result = layer_results.get('docking', {})
        if docking_result and docking_result.get('overall_risk'):
            risk_level = docking_result['overall_risk']
            weight = self.LAYER_WEIGHTS['docking']
            
            if risk_level != 'LOW' or not risk_scores:  # Include if significant or only option
                risk_scores.append((risk_level, weight, 'Competitive binding'))
                
                # Add affected enzymes from docking
                for enzyme in docking_result.get('affected_enzymes', []):
                    affected_enzymes.add(enzyme.get('name', ''))
                    if enzyme.get('competition'):
                        mechanisms.append(f"Competitive inhibition at {enzyme.get('name')}")
                
                confidence_factors.append(('Molecular docking', weight))
        
        # Layer 3: ML Induction Prediction (when available)
        ml_induction = layer_results.get('ml_induction', {})
        if ml_induction and ml_induction.get('risk'):
            risk_level = ml_induction['risk']
            confidence = ml_induction.get('confidence', 0.5)
            weight = self.LAYER_WEIGHTS['ml_induction'] * confidence
            
            risk_scores.append((risk_level, weight, 'ML induction prediction'))
            confidence_factors.append(('ML induction model', weight))
        
        # Layer 4: ML Inhibition Prediction (when available)
        ml_inhibition = layer_results.get('ml_inhibition', {})
        if ml_inhibition and ml_inhibition.get('risk'):
            risk_level = ml_inhibition['risk']
            confidence = ml_inhibition.get('confidence', 0.5)
            weight = self.LAYER_WEIGHTS['ml_inhibition'] * confidence
            
            risk_scores.append((risk_level, weight, 'ML inhibition prediction'))
            confidence_factors.append(('ML inhibition model', weight))
        
        # Calculate aggregated risk
        if not risk_scores:
            return self._create_unknown_result()
        
        final_risk = self._calculate_weighted_risk(risk_scores)
        overall_confidence = self._calculate_confidence(confidence_factors)
        
        return {
            'interaction_risk': final_risk,
            'confidence': overall_confidence,
            'mechanisms': list(set(mechanisms)),
            'affected_enzymes': list(affected_enzymes),
            'evidence_layers': [source for _, _, source in risk_scores],
            'layer_count': len(risk_scores),
            'breakdown': self._create_breakdown(risk_scores)
        }
    
    def _calculate_weighted_risk(self, risk_scores: List[Tuple[str, float, str]]) -> str:
        """
        Calculate weighted risk level from multiple scores.
        
        Args:
            risk_scores: List of (risk_level, weight, source) tuples
        
        Returns:
            Final risk level string
        """
        # Convert risk levels to numeric scores
        risk_values = {
            'HIGH': 3,
            'MEDIUM': 2,
            'LOW': 1,
            'UNKNOWN': 0
        }
        
        # Calculate weighted average
        total_weight = 0
        weighted_sum = 0
        
        for risk_level, weight, _ in risk_scores:
            risk_value = risk_values.get(risk_level, 0)
            weighted_sum += risk_value * weight
            total_weight += weight
        
        if total_weight == 0:
            return 'UNKNOWN'
        
        weighted_avg = weighted_sum / total_weight
        
        # Convert back to risk level
        # If any HIGH with significant weight, or weighted average > 2.5, return HIGH
        high_count = sum(1 for r, w, _ in risk_scores if r == 'HIGH' and w > 0.5)
        if high_count >= 1 or weighted_avg >= 2.5:
            return 'HIGH'
        elif weighted_avg >= 1.5:
            return 'MEDIUM'
        else:
            return 'LOW'
    
    def _calculate_confidence(self, confidence_factors: List[Tuple[str, float]]) -> float:
        """Calculate overall confidence score."""
        if not confidence_factors:
            return 0.3  # Low confidence if no evidence
        
        # Confidence increases with more evidence layers
        layer_count = len(confidence_factors)
        weights = [w for _, w in confidence_factors]
        
        # Base confidence from strongest evidence
        max_weight = max(weights) if weights else 0.3
        
        # Bonus for multiple layers
        layer_bonus = min(0.2 * (layer_count - 1), 0.4)
        
        confidence = min(max_weight + layer_bonus, 1.0)
        return round(confidence, 2)
    
    def _create_breakdown(self, risk_scores: List[Tuple[str, float, str]]) -> Dict:
        """Create detailed breakdown of risk assessment."""
        breakdown = {}
        for risk_level, weight, source in risk_scores:
            breakdown[source] = {
                'risk': risk_level,
                'weight': round(weight, 2)
            }
        return breakdown
    
    def _create_clinical_evidence_result(self, clinical: Dict) -> Dict:
        """Create result based on clinical evidence."""
        return {
            'interaction_risk': clinical['severity'].upper(),
            'confidence': 1.0,
            'mechanisms': [clinical.get('mechanism', 'Clinically documented')],
            'affected_enzymes': clinical.get('affected_enzymes', []),
            'evidence_layers': ['Clinical database'],
            'layer_count': 1,
            'clinical_cases': clinical.get('num_cases', 0),
            'references': clinical.get('references', []),
            'breakdown': {
                'Clinical evidence': {
                    'risk': clinical['severity'].upper(),
                    'weight': 1.0
                }
            }
        }
    
    def _create_unknown_result(self) -> Dict:
        """Create result when no evidence is available."""
        return {
            'interaction_risk': 'UNKNOWN',
            'confidence': 0.1,
            'mechanisms': ['Insufficient data'],
            'affected_enzymes': [],
            'evidence_layers': [],
            'layer_count': 0,
            'breakdown': {}
        }

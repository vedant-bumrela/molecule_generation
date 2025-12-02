#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ChemVerse backend API routes.
Handles molecule visualization, reaction prediction, and validation.
"""

from flask import Blueprint, request, jsonify
from pathlib import Path
import os

from .reaction_engine import ReactionEngine
from .molecule_utils import MoleculeConverter, MoleculeValidator

chemverse_bp = Blueprint('chemverse', __name__, url_prefix='/api/chemverse')

# Initialize engines
reaction_engine = ReactionEngine()
molecule_converter = MoleculeConverter()
molecule_validator = MoleculeValidator()


@chemverse_bp.route('/smiles-to-3d', methods=['POST'])
def smiles_to_3d():
    """
    Convert SMILES to 3D PDB structure.
    
    Input: { "smiles": "CCO" }
    Output: { "pdb": "..." }
    """
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES required'}), 400
        
        # Validate SMILES
        is_valid, message = molecule_validator.validate_smiles(smiles)
        if not is_valid:
            return jsonify({'error': message}), 400
        
        # Convert to 3D
        pdb_data = molecule_converter.smiles_to_pdb(smiles)
        
        if not pdb_data:
            return jsonify({'error': 'Failed to generate 3D structure'}), 500
        
        return jsonify({
            'pdb': pdb_data,
            'smiles': smiles
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@chemverse_bp.route('/react', methods=['POST'])
def predict_reaction():
    """
    Predict chemical reaction products.
    
    Input: {
        "reactants": ["CCO", "CC(=O)O"],
        "reaction_type": "esterification"
    }
    Output: {
        "products": ["..."],
        "feasible": true,
        "energy_change": -15.2,
        "mechanism": "..."
    }
    """
    try:
        data = request.get_json()
        reactants = data.get('reactants', [])
        reaction_type = data.get('reaction_type')
        
        if not reactants or len(reactants) < 2:
            return jsonify({'error': 'At least 2 reactants required'}), 400
        
        if not reaction_type:
            return jsonify({'error': 'Reaction type required'}), 400
        
        # Predict products
        result = reaction_engine.predict_products(reactants, reaction_type)
        
        if 'error' in result:
            return jsonify(result), 400
        
        return jsonify(result)
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@chemverse_bp.route('/validate', methods=['POST'])
def validate_molecule():
    """
    Validate molecule and get properties.
    
    Input: { "smiles": "CCO" }
    Output: {
        "valid": true,
        "molecular_weight": 46.07,
        "formula": "C2H6O",
        "warnings": []
    }
    """
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES required'}), 400
        
        # Validate and get properties
        is_valid, message = molecule_validator.validate_smiles(smiles)
        
        if not is_valid:
            return jsonify({
                'valid': False,
                'error': message
            }), 200
        
        # Get molecular properties
        properties = molecule_validator.get_properties(smiles)
        
        return jsonify({
            'valid': True,
            **properties
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@chemverse_bp.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return jsonify({'status': 'ok', 'service': 'ChemVerse API'})

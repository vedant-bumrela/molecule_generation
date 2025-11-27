#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
HDI (Herb-Drug Interaction) API routes.
Provides endpoints for analyzing interactions between herbal compounds and pharmaceutical drugs.
"""

import os
import uuid
import json
from datetime import datetime
from flask import Blueprint, request, jsonify
from werkzeug.utils import secure_filename

# Create blueprint
hdi_bp = Blueprint('hdi', __name__, url_prefix='/api/hdi')

# Allowed file extensions
ALLOWED_EXTENSIONS = {'pdb', 'mol', 'sdf', 'mol2'}

def allowed_file(filename):
    """Check if file extension is allowed."""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def generate_hdi_job_id():
    """Generate unique job ID for HDI analysis."""
    return f"hdi_{uuid.uuid4().hex[:12]}"

@hdi_bp.route('/analyze', methods=['POST'])
def analyze_hdi():
    """
    Analyze herb-drug interaction.
    
    Accepts herbal and drug compound information and performs docking
    against CYP450 enzyme panel to assess interaction risk.
    """
    try:
        # Get form data
        herb_name = request.form.get('herb_name', 'Unknown Herb')
        herb_compound = request.form.get('herb_compound', '')
        herb_smiles = request.form.get('herb_smiles', '')
        
        drug_name = request.form.get('drug_name', 'Unknown Drug')
        drug_compound = request.form.get('drug_compound', '')
        drug_smiles = request.form.get('drug_smiles', '')
        
        # Validate input - need at least SMILES for both
        if not herb_smiles:
            return jsonify({
                'error': 'Herb SMILES is required for analysis'
            }), 400
            
        if not drug_smiles:
            return jsonify({
                'error': 'Drug SMILES is required for analysis'
            }), 400
        
        # Generate job ID
        job_id = generate_hdi_job_id()
        
        # Create job directory
        job_dir = os.path.join('uploads', 'hdi_jobs', job_id)
        os.makedirs(job_dir, exist_ok=True)
        
        # Handle file uploads if provided
        herb_file = request.files.get('herb_file')
        drug_file = request.files.get('drug_file')
        
        herb_structure_path = None
        drug_structure_path = None
        
        if herb_file and allowed_file(herb_file.filename):
            filename = secure_filename(herb_file.filename)
            herb_structure_path = os.path.join(job_dir, f'herb_{filename}')
            herb_file.save(herb_structure_path)
        
        if drug_file and allowed_file(drug_file.filename):
            filename = secure_filename(drug_file.filename)
            drug_structure_path = os.path.join(job_dir, f'drug_{filename}')
            drug_file.save(drug_structure_path)
        
        # Store job metadata
        job_data = {
            'job_id': job_id,
            'status': 'pending',
            'herb_info': {
                'name': herb_name,
                'compound': herb_compound,
                'smiles': herb_smiles,
                'structure_file': herb_structure_path
            },
            'drug_info': {
                'name': drug_name,
                'compound': drug_compound,
                'smiles': drug_smiles,
                'structure_file': drug_structure_path
            },
            'timestamp': datetime.utcnow().isoformat(),
            'job_dir': job_dir
        }
        
        # Save job metadata
        metadata_path = os.path.join(job_dir, 'metadata.json')
        with open(metadata_path, 'w') as f:
            json.dump(job_data, f, indent=2)
        
        # Run HDI analysis in background thread
        import threading
        from ..hdi import HDIAnalyzer
        
        def run_hdi_analysis():
            """Run HDI analysis in background."""
            try:
                # Update status to running
                job_data['status'] = 'running'
                with open(metadata_path, 'w') as f:
                    json.dump(job_data, f, indent=2)
                
                # Initialize analyzer
                analyzer = HDIAnalyzer()
                
                # Run analysis
                results = analyzer.analyze_hdi(
                    herb_smiles=herb_smiles,
                    drug_smiles=drug_smiles,
                    herb_name=herb_name,
                    drug_name=drug_name,
                    output_dir=job_dir
                )
                
                # Update status to completed
                job_data['status'] = 'completed'
                job_data['results'] = results
                with open(metadata_path, 'w') as f:
                    json.dump(job_data, f, indent=2)
                
                # Also save standalone results file
                results_path = os.path.join(job_dir, 'results.json')
                with open(results_path, 'w') as f:
                    json.dump(results, f, indent=2)
                    
            except Exception as e:
                # Update status to failed
                job_data['status'] = 'failed'
                job_data['error'] = str(e)
                with open(metadata_path, 'w') as f:
                    json.dump(job_data, f, indent=2)
        
        # Start analysis in background
        thread = threading.Thread(target=run_hdi_analysis)
        thread.daemon = True
        thread.start()
        
        return jsonify({
            'job_id': job_id,
            'status': 'pending',
            'message': 'HDI analysis job submitted successfully',
            'herb_name': herb_name,
            'drug_name': drug_name
        }), 202
        
    except Exception as e:
        return jsonify({
            'error': f'Failed to submit HDI analysis: {str(e)}'
        }), 500


@hdi_bp.route('/status/<job_id>', methods=['GET'])
def get_hdi_status(job_id):
    """Get status of HDI analysis job."""
    try:
        job_dir = os.path.join('uploads', 'hdi_jobs', job_id)
        metadata_path = os.path.join(job_dir, 'metadata.json')
        
        if not os.path.exists(metadata_path):
            return jsonify({
                'error': 'Job not found'
            }), 404
        
        with open(metadata_path, 'r') as f:
            job_data = json.load(f)
        
        return jsonify({
            'job_id': job_id,
            'status': job_data.get('status', 'unknown'),
            'herb_name': job_data['herb_info']['name'],
            'drug_name': job_data['drug_info']['name'],
            'timestamp': job_data.get('timestamp')
        }), 200
        
    except Exception as e:
        return jsonify({
            'error': f'Failed to get job status: {str(e)}'
        }), 500


@hdi_bp.route('/results/<job_id>', methods=['GET'])
def get_hdi_results(job_id):
    """Get results of completed HDI analysis."""
    try:
        job_dir = os.path.join('uploads', 'hdi_jobs', job_id)
        results_path = os.path.join(job_dir, 'results.json')
        
        if not os.path.exists(results_path):
            # Check if job exists
            metadata_path = os.path.join(job_dir, 'metadata.json')
            if not os.path.exists(metadata_path):
                return jsonify({'error': 'Job not found'}), 404
            
            # Job exists but no results yet
            with open(metadata_path, 'r') as f:
                job_data = json.load(f)
            
            return jsonify({
                'job_id': job_id,
                'status': job_data.get('status', 'pending'),
                'message': 'Analysis in progress'
            }), 202
        
        # Read results
        with open(results_path, 'r') as f:
            results = json.load(f)
        
        return jsonify(results), 200
        
    except Exception as e:
        return jsonify({
            'error': f'Failed to retrieve results: {str(e)}'
        }), 500

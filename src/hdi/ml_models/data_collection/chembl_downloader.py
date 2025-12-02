#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Data collection from ChEMBL for CYP450 bioactivity.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Optional
import requests
import time


class ChEMBLDataCollector:
    """Collect CYP450 bioactivity data from ChEMBL database."""
    
    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
    
    CYP450_TARGETS = {
        'CYP3A4': 'CHEMBL340',
        'CYP2D6': 'CHEMBL299',
        'CYP2C9': 'CHEMBL259',
        'CYP2C19': 'CHEMBL4432',
        'CYP1A2': 'CHEMBL3356'
    }
    
    def __init__(self, output_dir: str = 'data/chembl'):
        """
        Initialize ChEMBL data collector.
        
        Args:
            output_dir: Directory to save downloaded data
        """
        self.output_dir = output_dir
        self.session = requests.Session()
        
    def download_cyp450_inhibition(self, enzyme: str = 'CYP3A4') -> pd.DataFrame:
        """
        Download CYP450 inhibition data from ChEMBL.
        
        Args:
            enzyme: Enzyme name (CYP3A4, CYP2D6, etc.)
            
        Returns:
            DataFrame with SMILES, activity values, and metadata
        """
        if enzyme not in self.CYP450_TARGETS:
            raise ValueError(f"Unknown enzyme: {enzyme}. Choose from {list(self.CYP450_TARGETS.keys())}")
        
        target_id = self.CYP450_TARGETS[enzyme]
        print(f"Downloading {enzyme} data from ChEMBL (target: {target_id})...")
        
        # Query ChEMBL API for bioactivity data
        url = f"{self.BASE_URL}/activity.json"
        params = {
            'target_chembl_id': target_id,
            'standard_type__in': 'IC50,Ki,Kd',  # Inhibition assays
            'standard_relation': '=',
            'limit': 10000  # Max results per page
        }
        
        all_activities = []
        offset = 0
        
        while True:
            params['offset'] = offset
            response = self.session.get(url, params=params)
            
            if response.status_code != 200:
                print(f"Error {response.status_code}: {response.text}")
                break
            
            data = response.json()
            activities = data.get('activities', [])
            
            if not activities:
                break
            
            all_activities.extend(activities)
            print(f"  Downloaded {len(all_activities)} activities...")
            
            # Check if there are more pages
            if len(activities) < params['limit']:
                break
            
            offset += params['limit']
            time.sleep(0.5)  # Be nice to the API
        
        print(f"Total activities downloaded: {len(all_activities)}")
        
        # Convert to DataFrame
        df = self._parse_activities(all_activities)
        df['enzyme'] = enzyme
        
        return df
    
    def _parse_activities(self, activities: List[Dict]) -> pd.DataFrame:
        """
        Parse ChEMBL activities into clean DataFrame.
        
        Args:
            activities: List of activity dictionaries from ChEMBL
            
        Returns:
            Cleaned DataFrame
        """
        records = []
        
        for act in activities:
            try:
                # Extract key information
                record = {
                    'chembl_id': act.get('molecule_chembl_id'),
                    'smiles': act.get('canonical_smiles'),
                    'activity_type': act.get('standard_type'),
                    'activity_value': act.get('standard_value'),
                    'activity_units': act.get('standard_units'),
                    'activity_relation': act.get('standard_relation'),
                    'assay_chembl_id': act.get('assay_chembl_id'),
                    'confidence_score': act.get('confidence_score'),
                    'pchembl_value': act.get('pchembl_value')
                }
                
                records.append(record)
            except Exception as e:
                print(f"Warning: Skipping activity due to error: {e}")
                continue
        
        df = pd.DataFrame(records)
        return df
    
    def clean_and_label_data(self, df: pd.DataFrame, threshold_nm: float = 10000) -> pd.DataFrame:
        """
        Clean data and create binary labels for inhibition.
        
        Args:
            df: Raw activity DataFrame
            threshold_nm: Threshold for inhibitor classification (nM)
            
        Returns:
            Cleaned DataFrame with 'is_inhibitor' label
        """
        # Filter for quality
        df_clean = df.copy()
        
        # Remove rows without SMILES
        df_clean = df_clean[df_clean['smiles'].notna()]
        
        # Remove rows without activity value
        df_clean = df_clean[df_clean['activity_value'].notna()]
        
        # Convert activity to numeric
        df_clean['activity_value'] = pd.to_numeric(df_clean['activity_value'], errors='coerce')
        df_clean = df_clean[df_clean['activity_value'].notna()]
        
        # Convert to nM if needed
        df_clean['activity_value_nm'] = df_clean.apply(
            lambda row: self._convert_to_nm(row['activity_value'], row['activity_units']),
            axis=1
        )
        
        # Create binary label (1 = inhibitor, 0 = non-inhibitor)
        df_clean['is_inhibitor'] = (df_clean['activity_value_nm'] < threshold_nm).astype(int)
        
        # Remove duplicates (keep strongest activity)
        df_clean = df_clean.sort_values('activity_value_nm')
        df_clean = df_clean.drop_duplicates(subset='smiles', keep='first')
        
        print(f"Cleaned dataset: {len(df_clean)} compounds")
        print(f"  Inhibitors: {df_clean['is_inhibitor'].sum()}")
        print(f"  Non-inhibitors: {(~df_clean['is_inhibitor'].astype(bool)).sum()}")
        
        return df_clean
    
    def _convert_to_nm(self, value: float, units: str) -> float:
        """Convert activity value to nM."""
        if pd.isna(units):
            return value  # Assume nM if no units
        
        units = str(units).lower()
        
        if 'nm' in units or 'nanomolar' in units:
            return value
        elif 'um' in units or 'micromolar' in units or 'μm' in units:
            return value * 1000
        elif 'mm' in units or 'millimolar' in units:
            return value * 1000000
        elif 'm' in units and 'molar' in units:
            return value * 1000000000
        else:
            return value  # Unknown units, assume nM
    
    def save_dataset(self, df: pd.DataFrame, enzyme: str):
        """
        Save dataset to CSV.
        
        Args:
            df: DataFrame to save
            enzyme: Enzyme name for filename
        """
        import os
        os.makedirs(self.output_dir, exist_ok=True)
        
        filename = f"{self.output_dir}/{enzyme}_inhibition_data.csv"
        df.to_csv(filename, index=False)
        print(f"Saved to {filename}")
        
        return filename


# Example usage and testing
def example_download():
    """Example: Download CYP3A4 data."""
    collector = ChEMBLDataCollector()
    
    # Download CYP3A4 inhibition data
    df_raw = collector.download_cyp450_inhibition('CYP3A4')
    
    # Clean and label
    df_clean = collector.clean_and_label_data(df_raw, threshold_nm=10000)
    
    # Save
    collector.save_dataset(df_clean, 'CYP3A4')
    
    print("\nDataset Summary:")
    print(df_clean[['smiles', 'activity_value_nm', 'is_inhibitor']].head(10))
    
    return df_clean


if __name__ == '__main__':
    # Note: This will make real API calls to ChEMBL
    print("ChEMBL Data Collector")
    print("=" * 60)
    print("\nThis script downloads real data from ChEMBL.")
    print("To run full download, uncomment the line below.\n")
    
    # Uncomment to run:
    # df = example_download()
    
    print("Setup complete. Ready to download CYP450 data when needed.")

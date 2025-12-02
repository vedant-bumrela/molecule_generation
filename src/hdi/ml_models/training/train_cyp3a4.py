#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Comprehensive training script for CYP3A4 models.
Trains both inhibition and induction classifiers.
"""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.metrics import (
    roc_auc_score, precision_score, recall_score, f1_score,
    confusion_matrix, classification_report
)
from sklearn.preprocessing import StandardScaler
import pickle
import os

from ..features import MolecularFeatureExtractor


class CYP3A4ModelTrainer:
    """Train CYP3A4 inhibition and induction models."""
    
    def __init__(self, output_dir='src/hdi/ml_models/trained'):
        """
        Initialize trainer.
        
        Args:
            output_dir: Directory to save trained models
        """
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        self.feature_extractor = MolecularFeatureExtractor()
        self.scaler = StandardScaler()
        
        self.inhibition_model = None
        self.induction_model = None
    
    def create_synthetic_dataset(self, n_inhibitors=200, n_non_inhibitors=200, n_inducers=50):
        """
        Create synthetic training data based on known CYP3A4 patterns.
        This is a placeholder until real ChEMBL data is downloaded.
        
        Args:
            n_inhibitors: Number of inhibitor samples
            n_non_inhibitors: Number of non-inhibitor samples
            n_inducers: Number of inducer samples
            
        Returns:
            DataFrame with SMILES and labels
        """
        print("Creating synthetic dataset...")
        
        # Known CYP3A4 inhibitors (high affinity)
        known_inhibitors = [
            'CC1COC(O1)(CN2C=NC=N2)C3=C(C=C(C=C3)Cl)Cl',  # Ketoconazole
            'CN1C=NC2=C1C=C(C=C2)NC(=O)C3=CC=CC=C3',  # Similar
            'CC(C)C1=NC2=CC=CC=C2C(=C1)C3=CC=CC=C3',  # Nifedipine-like
            'CC1=CC2=C(C=C1)NC(=O)C3=CC=CC=C32',  # Aromatic
            'C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N',  # Tryptophan-like
        ]
        
        # Known non-inhibitors (weak/no interaction)
        known_non_inhibitors = [
            'CC',
            'CC#N',  # Acetonitrile
            'CC(=O)Oc1ccccc1C(=O)O',  # Aspirin
            'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',  # Caffeine
            'CC(C)NCC(COC1=CC=C(C=C1)COCCOC)O',  # Metoprolol
            'CC(C)O',  # Isopropanol
            'CCO',  # Ethanol
            'CCCCO',  # Butanol
        ]
        
        # Known inducers
        known_inducers = [
            'CC1=C2C(=C(C(=C1O)C(C)CC(C)(C)C)C)C(=O)C3=C(C2=O)C(=C(C(=C3)O)C(C)CC(C)(C)C)O',  # Hyperforin
            'C1CC2CCC3CCCC3C2C1',  # Steroid-like
        ]
        
        data = []
        
        # Add inhibitors
        for _ in range(n_inhibitors):
            smiles = known_inhibitors[_ % len(known_inhibitors)]
            data.append({'smiles': smiles, 'is_inhibitor': 1, 'is_inducer': 0})
        
        # Add non-inhibitors
        for _ in range(n_non_inhibitors):
            smiles = known_non_inhibitors[_ % len(known_non_inhibitors)]
            data.append({'smiles': smiles, 'is_inhibitor': 0, 'is_inducer': 0})
        
        # Add inducers
        for _ in range(n_inducers):
            smiles = known_inducers[_ % len(known_inducers)]
            data.append({'smiles': smiles, 'is_inhibitor': 0, 'is_inducer': 1})
       
        # Also add some non-inducers
        for _ in range(n_inducers):
            smiles = known_non_inhibitors[_ % len(known_non_inhibitors)]
            data.append({'smiles': smiles, 'is_inhibitor': 0, 'is_inducer': 0})
        
        df = pd.DataFrame(data)
        
        print(f"Created dataset with {len(df)} samples")
        print(f"  Inhibitors: {df['is_inhibitor'].sum()}")
        print(f"  Inducers: {df['is_inducer'].sum()}")
        
        return df
    
    def extract_features(self, smiles_list):
        """Extract features from SMILES list."""
        print("Extracting features...")
        features = []
        valid_indices = []
        
        for idx, smiles in enumerate(smiles_list):
            feat = self.feature_extractor.extract_all_features(smiles)
            if feat is not None:
                features.append(feat)
                valid_indices.append(idx)
        
        print(f"Extracted features for {len(features)}/{len(smiles_list)} molecules")
        return np.array(features), valid_indices
    
    def train_inhibition_model(self, X_train, y_train, X_test, y_test):
        """
        Train CYP3A4 inhibition classifier.
        
        Returns:
            Trained model and performance metrics
        """
        print("\n" + "="*60)
        print("Training CYP3A4 Inhibition Model")
        print("="*60)
        
        # Scale features
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)
        
        # Train Random Forest model
        model = RandomForestClassifier(
            n_estimators=100,
            max_depth=20,
            min_samples_split=5,
            class_weight='balanced',
            random_state=42,
            n_jobs=-1
        )
        
        print("Training model...")
        model.fit(X_train_scaled, y_train)
        
        # Cross-validation
        print("\nCross-validation scores:")
        cv_scores = cross_val_score(
            model, X_train_scaled, y_train,
            cv=5, scoring='roc_auc', n_jobs=-1
        )
        print(f"  ROC-AUC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
        
        # Test set evaluation
        y_pred = model.predict(X_test_scaled)
        y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
        
        print("\nTest Set Performance:")
        print(f"  ROC-AUC: {roc_auc_score(y_test, y_pred_proba):.3f}")
        print(f"  Precision: {precision_score(y_test, y_pred):.3f}")
        print(f"  Recall: {recall_score(y_test, y_pred):.3f}")
        print(f"  F1-Score: {f1_score(y_test, y_pred):.3f}")
        
        print("\nConfusion Matrix:")
        print(confusion_matrix(y_test, y_pred))
        
        self.inhibition_model = model
        return model, {
            'cv_roc_auc': cv_scores.mean(),
            'test_roc_auc': roc_auc_score(y_test, y_pred_proba),
            'test_precision': precision_score(y_test, y_pred),
            'test_recall': recall_score(y_test, y_pred),
            'test_f1': f1_score(y_test, y_pred)
        }
    
    def train_induction_model(self, X_train, y_train, X_test, y_test):
        """
        Train CYP3A4 induction classifier.
        
        Returns:
            Trained model and performance metrics
        """
        print("\n" + "="*60)
        print("Training CYP3A4 Induction Model")
        print("="*60)
        
        # Use scaled features (already fitted from inhibition training)
        X_train_scaled = self.scaler.transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)
        
        # Train model
        model = RandomForestClassifier(
            n_estimators=100,
            max_depth=15,
            class_weight='balanced',
            random_state=42,
            n_jobs=-1
        )
        
        print("Training model...")
        model.fit(X_train_scaled, y_train)
        
        # Evaluation
        y_pred = model.predict(X_test_scaled)
        y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
        
        print("\nTest Set Performance:")
        print(f"  ROC-AUC: {roc_auc_score(y_test, y_pred_proba):.3f}")
        print(f"  Precision: {precision_score(y_test, y_pred, zero_division=0):.3f}")
        print(f"  Recall: {recall_score(y_test, y_pred, zero_division=0):.3f}")
        print(f"  F1-Score: {f1_score(y_test, y_pred, zero_division=0):.3f}")
        
        self.induction_model = model
        return model, {
            'test_roc_auc': roc_auc_score(y_test, y_pred_proba),
            'test_precision': precision_score(y_test, y_pred, zero_division=0),
            'test_recall': recall_score(y_test, y_pred, zero_division=0),
            'test_f1': f1_score(y_test, y_pred, zero_division=0)
        }
    
    def save_models(self):
        """Save trained models to disk."""
        print("\nSaving models...")
        
        # Save inhibition model
        if self.inhibition_model:
            model_path = f"{self.output_dir}/CYP3A4_inhibition.pkl"
            with open(model_path, 'wb') as f:
                pickle.dump({
                    'model': self.inhibition_model,
                    'scaler': self.scaler,
                    'feature_extractor_params': {
                        'morgan_radius': self.feature_extractor.morgan_radius,
                        'morgan_bits': self.feature_extractor.morgan_bits
                    }
                }, f)
            print(f"  Saved inhibition model to {model_path}")
        
        # Save induction model
        if self.induction_model:
            model_path = f"{self.output_dir}/CYP3A4_induction.pkl"
            with open(model_path, 'wb') as f:
                pickle.dump({
                    'model': self.induction_model,
                    'scaler': self.scaler
                }, f)
            print(f"  Saved induction model to {model_path}")
    
    def train_all(self):
        """Complete training pipeline."""
        print("="*70)
        print("CYP3A4 ML Model Training Pipeline")
        print("="*70)
        
        # Step 1: Create/load dataset
        df = self.create_synthetic_dataset(n_inhibitors=100, n_non_inhibitors=100, n_inducers=30)
        
        # Step 2: Extract features
        X, valid_indices = self.extract_features(df['smiles'].tolist())
        df_valid = df.iloc[valid_indices].reset_index(drop=True)
        
        # Step 3: Train inhibition model
        y_inhibition = df_valid['is_inhibitor'].values
        X_train, X_test, y_train_inh, y_test_inh = train_test_split(
            X, y_inhibition, test_size=0.2, random_state=42
        )
        
        inh_model, inh_metrics = self.train_inhibition_model(
            X_train, y_train_inh, X_test, y_test_inh
        )
        
        # Step 4: Train induction model
        y_induction = df_valid['is_inducer'].values
        _, _, y_train_ind, y_test_ind = train_test_split(
            X, y_induction, test_size=0.2, random_state=42
        )
        
        ind_model, ind_metrics = self.train_induction_model(
            X_train, y_train_ind, X_test, y_test_ind
        )
        
        # Step 5: Save models
        self.save_models()
        
        print("\n" + "="*70)
        print("Training Complete!")
        print("="*70)
        print("\nModel Performance Summary:")
        print(f"  Inhibition ROC-AUC: {inh_metrics['test_roc_auc']:.3f}")
        print(f"  Induction ROC-AUC: {ind_metrics['test_roc_auc']:.3f}")
        
        return inh_metrics, ind_metrics


if __name__ == '__main__':
    trainer = CYP3A4ModelTrainer()
    inh_metrics, ind_metrics = trainer.train_all()

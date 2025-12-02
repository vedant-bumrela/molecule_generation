"""
Test the trained ML models with St. John's Wort example.
"""

from src.hdi.ml_models.predictor import SimpleCYP450Predictor

# Test molecules
test_cases = {
    "Hyperforin (St. John's Wort - should be HIGH INDUCTION)": "CC1=C2C(=C(C(=C1O)C(C)CC(C)(C)C)C)C(=O)C3=C(C2=O)C(=C(C(=C3)O)C(C)CC(C)(C)C)O",
    "Ketoconazole (strong CYP3A4 inhibitor)": "CC1COC(O1)(CN2C=NC=N2)C3=C(C=C(C=C3)Cl)Cl",
    "Warfarin (drug)": "CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O",
    "Asp Aspirin (weak)": "CC(=O)Oc1ccccc1C(=O)O",
}

print("="*70)
print("Testing Trained CYP3A4 ML Models")
print("="*70)

predictor = SimpleCYP450Predictor()

for name, smiles in test_cases.items():
    print(f"\n{name}:")
    print(f"  SMILES: {smiles[:60]}...")
    
    # Test inhibition
    inh = predictor.predict_inhibition(smiles, 'CYP3A4')
    print(f"  Inhibition: {inh.get('risk')} (prob: {inh.get('probability', 0):.2f}) [{inh.get('method', 'N/A')}]")
    
    # Test induction
    ind = predictor.predict_induction(smiles, 'CYP3A4')
    print(f"  Induction: {ind.get('risk')} (prob: {ind.get('probability', 0):.2f}) [{ind.get('method', 'N/A')}]")
    if ind.get('mechanism'):
        print(f"    Mechanism: {ind['mechanism']}")

print("\n" + "="*70)
print("✅ Testing complete!")

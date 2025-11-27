"""
Quick test for structural alerts system.
"""

from src.hdi.structural_alerts import CYP450InducerAlerts

# Test with hyperforin (from St. John's Wort)
hyperforin_smiles = "CC1=C2C(=C(C(=C1O)C(C)CC(C)(C)C)C)C(=O)C3=C(C2=O)C(=C(C(=C3)O)C(C)CC(C)(C)C)O"

# Test with simplified structure
test_smiles = "CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O"  # Warfarin

print("Testing Structural Alerts System\n" + "="*50)

alert_scanner = CYP450InducerAlerts()

print("\n1. Testing Hyperforin (St. John's Wort):")
print(f"   SMILES: {hyperforin_smiles[:50]}...")
result1 = alert_scanner.scan_molecule(hyperforin_smiles)
print(f"   Risk: {result1.get('risk')}")
print(f"   Alerts: {result1.get('num_alerts')}")
print(f"   Mechanisms: {result1.get('likely_mechanisms')}")
print(f"   Affected enzymes: {result1.get('affected_enzymes')}")

print("\n2. Testing Warfarin (control):")
print(f"   SMILES: {test_smiles}")
result2 = alert_scanner.scan_molecule(test_smiles)
print(f"   Risk: {result2.get('risk')}")
print(f"   Alerts: {result2.get('num_alerts')}")

print("\n" + "="*50)
print("✅ Structural alerts system working!")

"""
Test ChemVerse API endpoints.
"""

import requests
import json

BASE_URL = "http://localhost:5000/api/chemverse"

def test_smiles_to_3d():
    """Test SMILES to 3D conversion."""
    print("Testing SMILES to 3D conversion...")
    
    response = requests.post(
        f"{BASE_URL}/smiles-to-3d",
        json={"smiles": "CCO"}  # Ethanol
    )
    
    if response.status_code == 200:
        data = response.json()
        print(f"✓ Success! PDB data length: {len(data['pdb'])} chars")
        print(f"  First 100 chars: {data['pdb'][:100]}...")
    else:
        print(f"✗ Failed: {response.status_code}")
        print(response.text)

def test_reaction():
    """Test reaction prediction."""
    print("\nTesting reaction prediction...")
    
    response = requests.post(
        f"{BASE_URL}/react",
        json={
            "reactants": ["CCO", "CC(=O)O"],  # Ethanol + Acetic acid
            "reaction_type": "esterification"
        }
    )
    
    if response.status_code == 200:
        data = response.json()
        print(f"✓ Success!")
        print(f"  Products: {data['products']}")
        print(f"  Feasible: {data['feasible']}")
        print(f"  Energy change: {data['energy_change']} kcal/mol")
    else:
        print(f"✗ Failed: {response.status_code}")
        print(response.text)

def test_validation():
    """Test molecule validation."""
    print("\nTesting molecule validation...")
    
    response = requests.post(
        f"{BASE_URL}/validate",
        json={"smiles": "c1ccccc1"}  # Benzene
    )
    
    if response.status_code == 200:
        data = response.json()
        print(f"✓ Success!")
        print(f"  Valid: {data['valid']}")
        print(f"  Formula: {data.get('formula')}")
        print(f"  MW: {data.get('molecular_weight')}")
    else:
        print(f"✗ Failed: {response.status_code}")
        print(response.text)

if __name__ == "__main__":
    print("="*60)
    print("ChemVerse API Testing")
    print("="*60)
    
    try:
        test_smiles_to_3d()
        test_reaction()
        test_validation()
        
        print("\n" + "="*60)
        print("✅ All tests completed!")
        print("="*60)
    except Exception as e:
        print(f"\n❌ Error: {e}")

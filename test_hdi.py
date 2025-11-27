import requests

# Test HDI API
response = requests.post(
    'http://localhost:5000/api/hdi/analyze',
    data={
        'herb_smiles': 'CCO',
        'drug_smiles': 'CC(C)O',
        'herb_name': 'Test Herb',
        'drug_name': 'Test Drug'
    }
)

print(f"Status Code: {response.status_code}")
print(f"Response: {response.text}")

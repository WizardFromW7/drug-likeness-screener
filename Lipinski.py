from rdkit import Chem, name
from rdkit.Chem import Descriptors
import pubchempy as pcp
import csv
def get_smiles(drug_name):
    results = pcp.get_compounds(drug_name,'name')
    if results:
        return results [0].smiles
    else:
        return None
    
def check_lipinski(smiles, name):
    mol = Chem.MolFromSmiles(smiles)
    
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)

    print(f"\n--- {name} ---")
    print(f"Molecular Weight: {mw:.2f}")
    print(f"LogP: {logp:.2f}")
    print(f"H-Bond Donors: {hbd}")
    print(f"H-Bond Acceptors: {hba}")
    passes = mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10
    print(f"Passes Lipinski's Rule of Five: {passes}")
    return {"Drug": name, "MW": round(mw,2), "LogP": round(logp,2), "HBD": hbd, "HBA": hba, "Passes": passes}

results = []
drugs = ["aspirin", "ibuprofen", "paracetamol", "viagra", "penicillin"]
for drug in drugs:
    smiles = get_smiles(drug)
    if smiles:
        result = check_lipinski(smiles, drug)
        results.append(result)
    else:
        print(f"{drug} not found!")

# Results will save to the same folder as this script
with open("drug_results.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["Drug","MW","LogP","HBD","HBA","Passes"])
    writer.writeheader()
    writer.writerows(results)

print("\nResults saved to drug_results.csv!")
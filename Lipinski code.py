from rdkit import Chem
from rdkit.Chem import Descriptors
import pubchempy as pcp  # PubChem API wrapper for fetching molecular data
import csv

def get_smiles(drug_name):
    """Fetch the SMILES string for a drug from PubChem by name."""
    results = pcp.get_compounds(drug_name, 'name')
    if results:
        return results[0].smiles  # Return the first match
    else:
        return None  # Drug not found in PubChem

def check_lipinski(smiles, name):
    """
    Calculate molecular descriptors and check Lipinski's Rule of Five.
    Lipinski's rules predict whether a compound is likely to be orally bioavailable:
    - Molecular weight <= 500 Da
    - LogP <= 5 (lipophilicity)
    - H-bond donors <= 5
    - H-bond acceptors <= 10
    """
    mol = Chem.MolFromSmiles(smiles)  # Convert SMILES string to RDKit molecule object

    # Calculate molecular descriptors
    mw = Descriptors.MolWt(mol)        # Molecular weight
    logp = Descriptors.MolLogP(mol)    # Lipophilicity
    hbd = Descriptors.NumHDonors(mol)  # Hydrogen bond donors
    hba = Descriptors.NumHAcceptors(mol)  # Hydrogen bond acceptors

    # Print results for this compound
    print(f"\n--- {name} ---")
    print(f"Molecular Weight: {mw:.2f}")
    print(f"LogP: {logp:.2f}")
    print(f"H-Bond Donors: {hbd}")
    print(f"H-Bond Acceptors: {hba}")

    # Check all four Lipinski rules
    passes = mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10
    print(f"Passes Lipinski's Rule of Five: {passes}")

    # Return results as a dictionary for CSV export
    return {"Drug": name, "MW": round(mw, 2), "LogP": round(logp, 2), "HBD": hbd, "HBA": hba, "Passes": passes}

# List of drugs to screen
drugs = ["aspirin", "ibuprofen", "paracetamol", "viagra", "penicillin"]

results = []
for drug in drugs:
    smiles = get_smiles(drug)
    if smiles:
        result = check_lipinski(smiles, drug)
        results.append(result)
    else:
        print(f"{drug} not found!")

# Export results to CSV — file saves in the same folder as this script
with open("drug_results.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["Drug", "MW", "LogP", "HBD", "HBA", "Passes"])
    writer.writeheader()
    writer.writerows(results)

print("\nResults saved to drug_results.csv!")
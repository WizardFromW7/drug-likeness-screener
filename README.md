# Drug-Likeness Screening Tool

A Python tool that looks up any drug by name, retrieves its molecular data from PubChem, and checks whether it meets Lipinski's Rule of Five — a key filter used in early drug discovery to predict oral bioavailability.

## What it does
- Takes a drug name as input
- Fetches SMILES data from PubChem via PubChemPy
- Calculates molecular descriptors using RDKit (MW, LogP, H-bond donors/acceptors)
- Evaluates against Lipinski's Rule of Five
- Exports results to CSV

## Requirements
pip install rdkit pubchempy

## Usage
python Lipinski.py

Enter a drug name when prompted. To screen multiple drugs at once, edit the `drugs` list in the script.

## Example output
| Drug | MW | LogP | HBD | HBA | Passes |
|---|---|---|---|---|---|
| Aspirin | 180.16 | 1.31 | 1 | 3 | True |
| Viagra | 666.71 | 0.36 | 5 | 11 | False |

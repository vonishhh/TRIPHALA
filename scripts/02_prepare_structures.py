import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import os

# ─────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────

BASE_DIR = os.path.join(os.path.dirname(__file__), "..")
INPUT_CSV = os.path.join(BASE_DIR, "01_compounds", "raw", "compound_library.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "01_compounds", "prepared")

# ─────────────────────────────────────────────────────
# FUNCTION: Convert one SMILES string into a 3D structure
# ─────────────────────────────────────────────────────

def generate_3d_structure(smiles, compound_id):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  Could not parse SMILES for {compound_id}")
        return None

    # 3D embedding needs explicit hydrogens, not just implicit ones
    mol = Chem.AddHs(mol)

    # Generate an initial 3D shape (randomSeed makes it reproducible)
    result = AllChem.EmbedMolecule(mol, randomSeed=42)
    if result != 0:
        print(f"  3D embedding failed for {compound_id}")
        return None

    # Energy-minimize so bond lengths/angles are physically realistic
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        AllChem.UFFOptimizeMolecule(mol)  # fallback force field

    return mol

# ─────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────

if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df = pd.read_csv(INPUT_CSV)
    print(f"Loaded {len(df)} compounds from {INPUT_CSV}")

    success_count = 0
    fail_count = 0

    for idx, row in df.iterrows():
        cid = row.get("CID", idx)
        smiles = row.get("SMILES")

        if pd.isna(smiles):
            print(f"Skipping CID {cid}: no SMILES available")
            fail_count += 1
            continue

        mol = generate_3d_structure(smiles, cid)

        if mol is None:
            fail_count += 1
            continue

        mol.SetProp("_Name", str(row.get("IUPACName", cid)))

        out_path = os.path.join(OUTPUT_DIR, f"{cid}.sdf")
        writer = Chem.SDWriter(out_path)
        writer.write(mol)
        writer.close()

        success_count += 1
        print(f"Saved 3D structure: {cid}.sdf")

    print(f"\nDone. {success_count} structures generated, {fail_count} failed.")
    print(f"Saved to: {OUTPUT_DIR}")
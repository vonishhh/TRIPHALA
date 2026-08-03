import os
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RAW_DIR = os.path.join(BASE_DIR, "02_proteins", "raw")
PREPARED_DIR = os.path.join(BASE_DIR, "02_proteins", "prepared")

PROTEINS = {
    "keap1.pdb": "keap1",
    "bcl2.pdb": "bcl2",
    "caspase3.pdb": "caspase3"
}

class ProteinOnly(Select):
    """Keep only standard protein residues — drops water and any ligand/ion."""
    def accept_residue(self, residue):
        return residue.id[0] == " "

def get_ligand_center(structure):
    """Finds the co-crystallized ligand (if any) and returns its center (x,y,z)."""
    coords = []
    for residue in structure.get_residues():
        het_flag = residue.id[0]
        if het_flag.startswith("H_") and residue.resname != "HOH":
            for atom in residue:
                coords.append(atom.coord)
    if coords:
        return np.array(coords).mean(axis=0)
    return None

if __name__ == "__main__":
    os.makedirs(PREPARED_DIR, exist_ok=True)
    parser = PDBParser(QUIET=True)
    io = PDBIO()

    for filename, label in PROTEINS.items():
        raw_path = os.path.join(RAW_DIR, filename)
        if not os.path.exists(raw_path):
            print(f"Missing file: {raw_path}, skipping.")
            continue

        print(f"\nProcessing {label}...")
        structure = parser.get_structure(label, raw_path)

        center = get_ligand_center(structure)
        if center is not None:
            print(f"  Ligand-based binding site center (x,y,z): {center}")
        else:
            print(f"  No co-crystallized ligand found in this structure.")

        out_path = os.path.join(PREPARED_DIR, f"{label}_protein_only.pdb")
        io.set_structure(structure)
        io.save(out_path, ProteinOnly())
        print(f"  Saved cleaned protein to: {out_path}")
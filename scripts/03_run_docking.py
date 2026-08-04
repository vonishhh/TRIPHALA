import os
import glob
import csv
from vina import Vina

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
LIGAND_DIR = os.path.join(BASE_DIR, "01_compounds", "pdbqt")
RECEPTOR_DIR = os.path.join(BASE_DIR, "02_proteins", "pdbqt")
DOCKING_DIR = os.path.join(BASE_DIR, "03_docking", "vina")
SCORES_DIR = os.path.join(BASE_DIR, "04_analysis", "scores")

# Each protein: receptor file + binding pocket center (from your ligand-based detection)
RECEPTORS = {
    "keap1": {
        "pdbqt": os.path.join(RECEPTOR_DIR, "keap1.pdbqt"),
        "center": [3.8096485, -0.66292655, -22.942854]
    },
    "bcl2": {
        "pdbqt": os.path.join(RECEPTOR_DIR, "bcl2.pdbqt"),
        "center": [15.966812, 2.049132, 16.19613]
    },
    "caspase3": {
        "pdbqt": os.path.join(RECEPTOR_DIR, "caspase3.pdbqt"),
        "center": [37.208946, 93.54983, 19.316647]
    }
}

BOX_SIZE = [20, 20, 20]  # Angstroms — widen later if PyMOL shows the pocket needs more room

def get_best_affinity(pdbqt_path):
    """Vina writes the top score directly into the output file's first REMARK line."""
    with open(pdbqt_path, "r") as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT"):
                return float(line.split()[3])
    return None

if __name__ == "__main__":
    os.makedirs(SCORES_DIR, exist_ok=True)
    ligand_files = sorted(glob.glob(os.path.join(LIGAND_DIR, "*.pdbqt")))
    print(f"Found {len(ligand_files)} ligands to dock.")

    results = []

    for protein_name, info in RECEPTORS.items():
        print(f"\n{'='*50}")
        print(f"Docking against: {protein_name}")
        print(f"{'='*50}")

        out_dir = os.path.join(DOCKING_DIR, protein_name)
        os.makedirs(out_dir, exist_ok=True)

        # Set receptor + compute grid maps ONCE per protein (reused for all ligands — big time saver)
        v = Vina(sf_name="vina")
        v.set_receptor(info["pdbqt"])
        v.compute_vina_maps(center=info["center"], box_size=BOX_SIZE)

        for i, ligand_path in enumerate(ligand_files):
            compound_id = os.path.splitext(os.path.basename(ligand_path))[0]
            out_path = os.path.join(out_dir, f"{compound_id}_out.pdbqt")

            try:
                v.set_ligand_from_file(ligand_path)
                v.dock(exhaustiveness=8, n_poses=9)  # Vina defaults — fast first-pass screen
                v.write_poses(out_path, n_poses=1, overwrite=True)

                affinity = get_best_affinity(out_path)
                results.append({"Protein": protein_name, "CID": compound_id, "Affinity_kcal_mol": affinity})
                print(f"  [{i+1}/{len(ligand_files)}] {compound_id}: {affinity} kcal/mol")

            except Exception as e:
                print(f"  [{i+1}/{len(ligand_files)}] {compound_id}: FAILED ({e})")
                results.append({"Protein": protein_name, "CID": compound_id, "Affinity_kcal_mol": None})

    out_csv = os.path.join(SCORES_DIR, "vina_scores.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Protein", "CID", "Affinity_kcal_mol"])
        writer.writeheader()
        writer.writerows(results)

    print(f"\nDone. Saved all scores to: {out_csv}")
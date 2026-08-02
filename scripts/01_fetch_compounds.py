import requests
import pandas as pd
import time
import os

# ─────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────

# Maps each plant to its CSV file (must match your saved filenames exactly)
PLANT_FILES = {
    "Phyllanthus emblica": "phyllanthus_emblica_compounds.csv",
    "Terminalia chebula": "terminalia_chebula_compounds.csv",
    "Terminalia bellirica": "terminalia_bellirica_compounds.csv"
}

BASE_DIR = os.path.join(os.path.dirname(__file__), "..")
OUTPUT_DIR = os.path.join(BASE_DIR, "01_compounds", "raw")

# ─────────────────────────────────────────────────────
# FUNCTION 1: Load compound names from a plant's CSV
# ─────────────────────────────────────────────────────
def load_compound_names(csv_path):
    if not os.path.exists(csv_path):
        print(f"File not found: {csv_path}")
        return []
    with open(csv_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    names = [line.strip() for line in lines if line.strip()]
    # Remove a header row if present
    if names and names[0].lower() in ("compound_name", "name", "compound"):
        names = names[1:]
    return names
# ─────────────────────────────────────────────────────
# FUNCTION 2: Look up ONE compound's PubChem CID by name
# ─────────────────────────────────────────────────────

def get_cid_for_compound(compound_name):
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
        f"/name/{compound_name.replace(' ', '%20')}"
        f"/cids/JSON"
    )
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            cids = data.get("IdentifierList", {}).get("CID", [])
            if cids:
                return cids[0]
        return None
    except Exception as e:
        print(f"  Lookup failed for {compound_name}: {e}")
        return None

# ─────────────────────────────────────────────────────
# FUNCTION 3: Get compound properties in bulk (unchanged, already works)
# ─────────────────────────────────────────────────────

def fetch_properties(cid_list, plant_name):
    if not cid_list:
        return pd.DataFrame()

    properties = ",".join([
        "MolecularWeight", "XLogP", "HBondDonorCount", "HBondAcceptorCount",
        "IsomericSMILES", "IUPACName", "RotatableBondCount", "TPSA"
    ])

    all_data = []
    for i in range(0, len(cid_list), 100):
        batch = cid_list[i:i+100]
        cids_str = ",".join(map(str, batch))
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
            f"/cid/{cids_str}/property/{properties}/JSON"
        )
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                batch_data = response.json()["PropertyTable"]["Properties"]
                all_data.extend(batch_data)
                print(f"  Fetched batch {i//100 + 1}: {len(batch_data)} compounds")
        except Exception as e:
            print(f"  Batch {i//100+1} failed: {e}")
        time.sleep(0.5)

    df = pd.DataFrame(all_data)
    df["Source_Plant"] = plant_name
    return df

# ─────────────────────────────────────────────────────
# FUNCTION 4: Lipinski filter (unchanged)
# ─────────────────────────────────────────────────────

def apply_filters(df):
    required = ["MolecularWeight", "XLogP", "HBondDonorCount", "HBondAcceptorCount"]
    df = df.dropna(subset=required)
    for col in required:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    strict = df[
        (df["MolecularWeight"] <= 500) & (df["XLogP"] <= 5) &
        (df["HBondDonorCount"] <= 5) & (df["HBondAcceptorCount"] <= 10)
    ].copy()

    relaxed = df[
        (df["MolecularWeight"] <= 900) & (df["XLogP"] <= 7) &
        (df["HBondDonorCount"] <= 10) & (df["HBondAcceptorCount"] <= 20)
    ].copy()

    strict["filter_passed"] = "strict_lipinski"
    relaxed_only = relaxed[~relaxed["CID"].isin(strict["CID"])].copy()
    relaxed_only["filter_passed"] = "relaxed"

    final = pd.concat([strict, relaxed_only], ignore_index=True)
    print(f"Total compounds: {len(df)}")
    print(f"Passed strict Lipinski: {len(strict)}")
    print(f"Passed relaxed filter: {len(relaxed_only)}")
    print(f"Final compound pool: {len(final)}")
    return final

# ─────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────

if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    all_dfs = []

    for plant_name, filename in PLANT_FILES.items():
        print(f"\n{'='*50}")
        print(f"Processing: {plant_name}")
        print(f"{'='*50}")

        csv_path = os.path.join(OUTPUT_DIR, filename)
        compound_names = load_compound_names(csv_path)

        if not compound_names:
            print(f"No compound names loaded for {plant_name}, skipping.")
            continue

        print(f"Loaded {len(compound_names)} names. Looking up CIDs...")
        cids = []
        for name in compound_names:
            cid = get_cid_for_compound(name)
            if cid:
                cids.append(cid)
            else:
                print(f"  Not found on PubChem: {name}")
            time.sleep(0.3)

        print(f"Resolved {len(cids)} of {len(compound_names)} compounds.")

        if not cids:
            continue

        df = fetch_properties(cids, plant_name)
        if not df.empty:
            all_dfs.append(df)
        time.sleep(1)

    if not all_dfs:
        print("\nNo data collected. Stopping.")
    else:
        master_df = pd.concat(all_dfs, ignore_index=True)
        master_df = master_df.drop_duplicates(subset="CID")
        print(f"\nTotal unique compounds across all plants: {len(master_df)}")

        filtered_df = apply_filters(master_df)
        out_path = os.path.join(OUTPUT_DIR, "compound_library.csv")
        filtered_df.to_csv(out_path, index=False)
        print(f"\nSaved compound library to: {out_path}")

        print("\nTop 10 compounds by molecular weight:")
        cols = ["IUPACName", "MolecularWeight", "XLogP", "Source_Plant", "filter_passed"]
        print(filtered_df[cols].head(10).to_string(index=False))
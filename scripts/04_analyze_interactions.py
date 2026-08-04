import pandas as pd
import os

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SCORES_DIR = os.path.join(BASE_DIR, "04_analysis", "scores")
COMPOUND_LIB = os.path.join(BASE_DIR, "01_compounds", "raw", "compound_library.csv")
VINA_SCORES = os.path.join(SCORES_DIR, "vina_scores.csv")

TOP_N = 10  # how many top hits to keep per protein

if __name__ == "__main__":
    scores_df = pd.read_csv(VINA_SCORES)
    compounds_df = pd.read_csv(COMPOUND_LIB)

    # Attach compound names + source plant to each score
    merged = scores_df.merge(
        compounds_df[["CID", "IUPACName", "Source_Plant"]],
        on="CID",
        how="left"
    )

    # Keep the top N strongest binders per protein (most negative = strongest)
    top_hits = (
        merged.sort_values("Affinity_kcal_mol")
        .groupby("Protein")
        .head(TOP_N)
        .sort_values(["Protein", "Affinity_kcal_mol"])
    )

    out_path = os.path.join(SCORES_DIR, "top_hits_vina.csv")
    top_hits.to_csv(out_path, index=False)
    print(f"Saved top {TOP_N} hits per protein to: {out_path}")
    print(top_hits.to_string(index=False))
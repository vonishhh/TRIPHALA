#!/bin/bash

echo "Creating Triphala project structure..."

mkdir -p 01_compounds/{raw,prepared,pdbqt} 02_proteins/{raw,prepared,pdbqt} 03_docking/vina/{keap1,bcl2,caspase3} 03_docking/gnina/{keap1,bcl2,caspase3} 04_analysis/{interactions,scores,plots} 05_ml/{features,models,synergy} 06_md/{topology,trajectories} scripts notebooks results

echo "# Triphala Dual-Target Docking Project" > README.md
echo "Computational study of Triphala polyphenols against Keap1-Nrf2, Bcl-2, Caspase-3" >> README.md

touch scripts/01_fetch_compounds.py
touch scripts/02_prepare_structures.py
touch scripts/03_run_docking.py
touch scripts/04_analyze_interactions.py
touch scripts/05_ml_synergy.py

echo "Done! Folder structure created inside $(pwd)"

conda create -n triphala python=3.10 -y

# echo "Setup complete. Always run 'conda activate triphala' before working."
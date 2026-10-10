#!/bin/bash
# 06 — Neoantigen prediction (MHCflurry across full diploid profile, all 8–11mers)
set -euo pipefail
mkdir -p results/neoantigens logs
python3 workflow/scripts/vep_to_peptides.py results/vep/SRR30677899.vep.txt results/peptides/tiled_9mers.csv || \
python3 workflow/scripts/vep_to_peptides.py
echo "mhcflurry-predict results/peptides/tiled_9mers.csv --alleles A*02:06 A*30:01 B*08:01 B*53:01 C*07:01 ..."
mhcflurry-predict results/peptides/tiled_9mers.csv --alleles "HLA-A*02:06" \
  --out results/neoantigens/mhcflurry_A0206.csv 2>&1 | tee logs/06_mhcflurry.log || true
python3 workflow/scripts/summarize_binders.py || true

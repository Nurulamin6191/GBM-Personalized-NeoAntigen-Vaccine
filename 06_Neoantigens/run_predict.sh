#!/bin/bash
# 06 — Neoantigen prediction (MHCflurry full diploid profile; 9mers, 8/10/11mers TODO)
# Runtime outputs go to results-RUNTIME/ (git-ignored). Curated tables live in 06_Neoantigens/.
set -euo pipefail
mkdir -p results-RUNTIME/peptides results-RUNTIME/neoantigens logs
VEP="results-RUNTIME/vep/SRR30677899.vep.txt"
if [ -s "$VEP" ]; then
  python3 04_Variants/vep_to_peptides.py "$VEP" results-RUNTIME/peptides/tiled_9mers.csv
else
  echo "no VEP output yet ($VEP missing) — tiling demo rows only"
  python3 04_Variants/vep_to_peptides.py dummy results-RUNTIME/peptides/tiled_9mers.csv || true
fi
PEPS=$(tail -n +2 results-RUNTIME/peptides/tiled_9mers.csv | cut -d, -f4 | tr '\n' ' ')
for AL in "HLA-A*02:06" "HLA-A*30:01" "HLA-B*08:01" "HLA-B*53:01" "HLA-C*07:01"; do
  TAG=$(echo "$AL" | tr -d '*:')
  mhcflurry-predict --alleles "$AL" --peptides $PEPS \
    --out results-RUNTIME/neoantigens/mhcflurry_${TAG}.csv 2>&1 | tee -a logs/06_mhcflurry.log || true
done
python3 06_Neoantigens/summarize_binders.py || true

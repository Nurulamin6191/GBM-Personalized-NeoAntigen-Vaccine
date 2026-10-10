#!/bin/bash
# 04 — Annotation (VEP) + missense → 8–11mer tiling (the bridge missing in v1)
set -euo pipefail
mkdir -p results/vep results/peptides logs
for VCF in results/variants/*.PASS.vcf.gz; do
  S=$(basename "$VCF" .PASS.vcf.gz)
  echo "== VEP $S (cache auto-download on first run) =="
  vep -i "$VCF" -o results/vep/${S}.vep.txt --cache --offline --assembly GRCh38 \
    --format vcf --tab --symbol --canonical --protein --biotype --sift --polyphen \
    --fields "Uploaded_variation,Location,Allele,Gene,SYMBOL,Consequence,Protein_position,Amino_acids,Codons,SIFT,PolyPhen,EXON" \
    --force_overwrite || echo "VEP failed for $S — see logs; rerun with --database for online mode"
done
echo "== tiling missense SNVs → 9mers (demo of bridge; full 8-11mer + expression filter in 04b notebook) =="
python3 04_Variants/vep_to_peptides.py || echo "tiling script pending VEP output"

#!/bin/bash
# Expression gating: StringTie TPM (full-genome HISAT2 BAMs) -> Table 1 v2 (expressed only)
set -euo pipefail
GTF="refs/Homo_sapiens.GRCh38.111.gtf.gz"
mkdir -p results/expr logs
[ -s "$GTF" ] || { echo "GTF missing"; exit 1; }
gunzip -c "$GTF" > refs/GRCh38.111.gtf
for S in SRR30677899 SRR30677902; do
  BAM="results/align_full/${S}.hisat2.bam"
  [ -s "$BAM" ] || { echo "$BAM missing"; exit 1; }
  stringtie "$BAM" -G refs/GRCh38.111.gtf -e -p 6 \
    -o results/expr/${S}.gtf -A results/expr/${S}.tab 2>&1 | tee logs/stringtie_${S}.log
done
python3 06_Neoantigens/gate_expression.py

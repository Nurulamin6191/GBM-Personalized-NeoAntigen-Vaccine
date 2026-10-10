#!/bin/bash
# 05 — HLA typing (OptiType 1.5.0, RNA mode, full diploid 6-allele call)
# Env: conda env gbm-pub provides optitype + razers3 + cbc/glpk (check: optitype check-deps)
set -euo pipefail
mkdir -p 05_HLA/optitype_rerun logs
# v1 reported 5 alleles with 2nd HLA-C missing; prior OptiType tsv showed C*07:01 homozygous.
# This re-derives from trimmed RNA-seq FASTQs:
for S in SRR30677899 SRR30677902; do
  R1=$(ls results-RUNTIME/qc/trimmed/${S}_1_val_1.fq.gz 2>/dev/null || ls data-RUNTIME/fastq/${S}_1.fastq.gz)
  R2=$(ls results-RUNTIME/qc/trimmed/${S}_2_val_2.fq.gz 2>/dev/null || ls data-RUNTIME/fastq/${S}_2.fastq.gz)
  optitype run -i "$R1" -i "$R2" --rna -o 05_HLA/optitype_rerun -p "$S" \
    --solver cbc --threads 4 2>&1 | tee logs/05_optitype_${S}.log
done
cat 05_HLA/optitype_rerun/*_result.tsv 2>/dev/null | tee 05_HLA/rerun_profile.tsv || \
  echo "rerun pending — prior profile: A*02:06/A*30:01 B*08:01/B*53:01 C*07:01/C*07:01 (homozygous)"

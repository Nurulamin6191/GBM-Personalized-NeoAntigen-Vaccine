#!/bin/bash
# Expression gating: StringTie TPM (full-genome HISAT2 BAMs) -> Table 1 v2 (expressed only)
# Prerequisites: 01_Data/fetch.sh, 02_QC/run_qc.sh, HISAT2 index+align (see 03_Alignment/run_align.sh),
#   GTF at refs-RUNTIME/Homo_sapiens.GRCh38.111.gtf.gz (see 01_Data/fetch.sh).
set -euo pipefail
GTF="refs-RUNTIME/Homo_sapiens.GRCh38.111.gtf.gz"
mkdir -p results-RUNTIME/expr logs
[ -s "$GTF" ] || { echo "GTF missing: fetch ENSEMBL GRCh38 r111 GTF to $GTF"; exit 1; }
[ -s refs-RUNTIME/GRCh38.111.gtf ] || gunzip -c "$GTF" > refs-RUNTIME/GRCh38.111.gtf
for S in SRR30677899 SRR30677902; do
  BAM="results-RUNTIME/align_full/${S}.hisat2.bam"
  [ -s "$BAM" ] || { echo "$BAM missing: run HISAT2 full-genome alignment first"; exit 1; }
  stringtie "$BAM" -G refs-RUNTIME/GRCh38.111.gtf -e -p 6 \
    -o results-RUNTIME/expr/${S}.gtf -A results-RUNTIME/expr/${S}.tab 2>&1 | tee logs/stringtie_${S}.log
done
python3 06_Neoantigens/gate_expression.py

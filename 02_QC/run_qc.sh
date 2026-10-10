#!/bin/bash
# 01 — QC: FastQC (raw) → Trim Galore (Q20) → FastQC (trimmed) → MultiQC
set -euo pipefail
FQ="data-RUNTIME/fastq"
QC="results-RUNTIME/qc"
mkdir -p "$QC/raw" "$QC/trimmed" "$QC/multiqc" logs
echo "[QC] FastQC raw..."
fastqc --outdir "$QC/raw" -t 8 "$FQ"/*.fastq.gz
echo "[QC] Trim Galore Q20, min-len 20 (R1 raw=29bp; 50 would delete all R1)..."
# NovaSeq 2-colour poly-G + SMARTer oligo in R2: --2colour 20 trims poly-G by quality
trim_galore --paired --quality 20 --length 20 --trim-n --2colour 20 \
  --cores 4 -o "$QC/trimmed" "$FQ"/*_1.fastq.gz "$FQ"/*_2.fastq.gz || \
trim_galore --quality 20 --length 20 --trim-n --cores 4 -o "$QC/trimmed" "$FQ"/*.fastq.gz
echo "[QC] MultiQC..."
multiqc "$QC" -o "$QC/multiqc" -f

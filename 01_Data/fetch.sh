#!/bin/bash
# 00 — Data acquisition (SRA → FASTQ)
# Runs after: conda activate gbm-pub
set -euo pipefail
ACC_LIST="SRR30677899 SRR30677902"
OUT="data-RUNTIME/fastq"
LOG="logs/00_fetch.log"
mkdir -p "$OUT" logs
: > "$LOG"
for ACC in $ACC_LIST; do
  echo "[$(date -Is)] $ACC prefetch" | tee -a "$LOG"
  prefetch "$ACC" -O data-RUNTIME/sra 2>&1 | tee -a "$LOG"
  echo "[$(date -Is)] $ACC fasterq-dump" | tee -a "$LOG"
  fasterq-dump "data-RUNTIME/sra/$ACC/$ACC.sra" -O "$OUT" -e 8 --split-files 2>&1 | tee -a "$LOG"
  gzip -f "$OUT"/${ACC}*.fastq 2>&1 | tee -a "$LOG"
done
echo "[$(date -Is)] md5 + counts" | tee -a "$LOG"
md5sum "$OUT"/*.gz | tee -a "$LOG"
for f in "$OUT"/*.gz; do echo "$f $(zcat "$f" | wc -l | awk '{print $1/4}') reads" | tee -a "$LOG"; done

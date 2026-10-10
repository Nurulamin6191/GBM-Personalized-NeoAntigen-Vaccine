#!/bin/bash
# 02 — Alignment
# A (pilot, reproduces original): STAR vs GRCh38 chr1
# B (publication): HISAT2 vs GRCh38 full primary assembly (fits 15 GB RAM; STAR full index needs ~30 GB)
set -euo pipefail
mkdir -p results/align logs refs
CHR1_FA="refs/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz"
FULL_FA="refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF="refs/Homo_sapiens.GRCh38.111.gtf.gz"

echo "== download refs if missing =="
for U in \
 "https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz" \
 "https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
 "https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz"; do
  F="refs/$(basename $U)"; [ -s "$F" ] || { echo "GET $U"; curl -L -o "$F" "$U"; }
done

echo "== A: STAR chr1 index (if missing) =="
if [ ! -d refs/star_chr1 ]; then
  mkdir -p refs/star_chr1 && gunzip -c "$CHR1_FA" > refs/chr1.fa
  STAR --runMode genomeGenerate --genomeDir refs/star_chr1 \
    --genomeFastaFiles refs/chr1.fa --runThreadN 8 --genomeSAindexNbases 12
fi
echo "== B: HISAT2 full index (if missing) =="
if [ ! -f refs/hisat2_full/genome.1.ht2 ]; then
  mkdir -p refs/hisat2_full && gunzip -c "$FULL_FA" > refs/full.fa
  hisat2-build -p 8 refs/full.fa refs/hisat2_full/genome
fi
echo "== align each sample, both strategies =="
for R1 in data/fastq/*_1.fastq.gz; do
  S=$(basename "$R1" _1.fastq.gz)
  R2="${R1/_1.fastq.gz/_2.fastq.gz}"
  TRIM1=$(ls results/qc/trimmed/${S}*val_1.fq.gz 2>/dev/null | head -1 || echo "$R1")
  TRIM2=$(ls results/qc/trimmed/${S}*val_2.fq.gz 2>/dev/null | head -1 || echo "$R2")
  echo "-- $S STAR-chr1 --"
  STAR --genomeDir refs/star_chr1 --readFilesIn "$TRIM1" "$TRIM2" \
    --readFilesCommand zcat --runThreadN 8 \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix results/align/${S}.star_chr1.
  echo "-- $S HISAT2-full --"
  hisat2 -p 8 -x refs/hisat2_full/genome -1 "$TRIM1" -2 "$TRIM2" \
    | samtools sort -@8 -o results/align/${S}.hisat2.bam
  samtools index results/align/${S}.hisat2.bam
  samtools flagstat results/align/${S}.hisat2.bam > results/align/${S}.hisat2.flagstat.txt
done

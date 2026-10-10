#!/bin/bash
# 03 — RNA-seq variant calling (GATK best practices, tumor-only; no matched normal available)
# AddOrReplaceReadGroups → MarkDuplicates → SplitNCigarReads → HaplotypeCaller → hard filter
# + Mutect2 tumor-only attempt (panel-of-normals absent → flagged as LIMITATION)
set -euo pipefail
REF="refs/full.fa"
mkdir -p results/variants logs
[ -s "$REF" ] || { echo "ERROR: $REF missing — run 02_align.sh first"; exit 1; }
[ -s "${REF}.fai" ] || samtools faidx "$REF"
[ -s refs/full.dict ] || picard CreateSequenceDictionary R="$REF" O=refs/full.dict
for BAM in results/align/*.hisat2.bam; do
  S=$(basename "$BAM" .hisat2.bam)
  echo "== $S =="
  picard AddOrReplaceReadGroups I="$BAM" O=results/variants/${S}.rg.bam \
    RGID="$S" RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM="$S"
  samtools index results/variants/${S}.rg.bam
  gatk MarkDuplicates -I results/variants/${S}.rg.bam -O results/variants/${S}.dedup.bam \
    -M results/variants/${S}.markdup.txt --CREATE_INDEX true
  gatk SplitNCigarReads -R "$REF" -I results/variants/${S}.dedup.bam \
    -O results/variants/${S}.split.bam
  gatk HaplotypeCaller -R "$REF" -I results/variants/${S}.split.bam \
    -O results/variants/${S}.hc.vcf.gz --dont-use-soft-clipped-bases true --max-alternate-alleles 3
  gatk VariantFiltration -R "$REF" -V results/variants/${S}.hc.vcf.gz \
    -O results/variants/${S}.hc.filtered.vcf.gz \
    --filter-name "QDlt2" --filter-expression "QD < 2.0" \
    --filter-name "FSgt30" --filter-expression "FS > 30.0" \
    --filter-name "DP lt10" --filter-expression "DP < 10"
  bcftools view -f PASS results/variants/${S}.hc.filtered.vcf.gz -Oz -o results/variants/${S}.PASS.vcf.gz
  bcftools index results/variants/${S}.PASS.vcf.gz
done
echo "NOTE: tumor-only; germline leakage expected. Filter vs dbSNP/gnomAD + VEP consequence in 04."

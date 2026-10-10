# Methods — full re-run (publication draft)

## Design
Tumor-only RNA-seq neoantigen discovery on PDX glioblastoma under temozolomide
(PRJNA1161693, SAMN43786660; runs SRR30677899: 15.4M reads, SRR30677902: 18.0M reads;
Illumina NovaSeq 6000, paired). No matched normal — germline leakage is a stated limitation.

## Steps
1. **Acquisition:** `prefetch` + `fasterq-dump --split-files` (script `00_fetch.sh`), md5 + read counts logged.
2. **QC:** FastQC raw → Trim Galore Q20/l≥50 → FastQC trimmed → MultiQC (`01_qc.sh`).
3. **Alignment:** HISAT2 vs GRCh38 primary assembly (publication primary; fits 15 GB RAM)
   + STAR vs chr1 pilot (reproduces original). `samtools flagstat` recorded (`02_align.sh`).
4. **Variants:** Picard AddOrReplaceReadGroups → MarkDuplicates → GATK SplitNCigarReads →
   HaplotypeCaller → VariantFiltration (QD<2, FS>30, DP<10) → PASS set (`03_variants.sh`).
   Mutect2 tumor-only attempted; absent PoN documented.
5. **Annotation:** VEP 111 (canonical, protein, SIFT/PolyPhen) → missense filter (`04_annotate.sh`).
6. **Peptide tiling (bridge added):** missense SNVs → 8–11mers (`vep_to_peptides.py`). v1 skipped this.
7. **HLA:** OptiType on trimmed FASTQs, full diploid 6-allele call (`05_hla.sh`).
   v1 listed 5 alleles (2nd HLA-C missing) — re-derived here.
8. **Binding:** MHCflurry 2.1 across all 6 alleles + lengths 8–11; NetMHCpan-4.1 external
   validation; IEDB cutoffs strong <50 / weak <500 nM (`06_predict.sh`, `summarize_binders.py`).
9. **Prioritisation (new):** expression (StringTie TPM), VAF×DP, agretopicity (mut/WT),
   hydrophobicity, cleavage — absent in v1.

## Corrections to v1
- "~13% uniquely mapped = Excellent" is **expected** for chr1-only alignment (~8% of genome);
  not a quality signal. Full-genome rate reported in full run.
- SYFPEITHI is a textbook HLA-A*02:01 **binder control**, mislabelled "negative control" in v1.
- Single-allele (A*02:06), 2-peptide test → full-allele, proteome-wide screen in rebuild.

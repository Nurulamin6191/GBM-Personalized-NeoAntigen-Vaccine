# Manuscript outline — GBM personalised neoantigen discovery (tumor RNA-seq, chr1 pilot + genome-wide FreeBayes)

Target: brief communication / methods-with-results journal (e.g. BMC Bioinformatics,
Frontiers in Immunology - Cancer Immunity, JITC tech report).

## Title
Evidence-first neoantigen prioritisation in glioblastoma PDX RNA-seq:
reproducibility audit, two-caller variant check and full-diploid MHCflurry screen

## Abstract (5 sentences)
1. Personalised GBM vaccines need ranked neoantigens from tumor RNA-seq.
2. We audited + rebuilt a chr1-pilot pipeline (PRJNA1161693, PDX + temozolomide, 2 runs).
3. STAR mapping reproduced 3x (13.9%); HLA re-typed (A*02:06/A*30:01 B*08:01/B*53:01 C*07:01 hom).
4. FreeBayes 165,818 -> 4,224 missense; GATK HC 43,461 -> 8,260 PASS (19% overlap);
   SDF4 A120V is FreeBayes-only (2 ALT reads, QUAL~0, GATK-absent) and its 9mers do not
   bind A*02:06 (all >6,800 nM) -> excluded as lead; batch top hits ~10-12 nM + IEDB
   NetMHCpan cross-validation reported as the prioritised set.
5. All code, env, figures and validation tables are versioned; limitations + required
   wet-lab confirmations listed.

## Figures
- Fig1 pipeline + added tiling bridge; Fig2 preliminary affinity; Fig3 chr1 mapping math;
  Fig4 top-20 batch binders; Fig5 consequence spectrum; Fig6 adapter burden;
  Fig7 SDF4 full-allele matrix; Fig8 replicate consistency + IEDB validation.

## Tables
- Table 1: top-10 batch candidates (gene symbol, mutation, peptide, A*02:06 affinity,
  dbSNP status, expression status TBD from StringTie).
- Table S1: sample/run metadata; S2: QC metrics; S3: STAR logs; S4: missense counts;
  S5: SDF4 pileup evidence; S6: full-allele SDF4 matrix; S7: per-sample A*02:06 re-score;
  S8: IEDB NetMHCpan values; S9: HLA runs (899 vs 902); S10: tool versions.

## Required wet-lab confirmations (state explicitly)
Sanger/IGV of any lead variant; DNA-based HLA; immunogenicity (ELISpot) for Table 1;
tumour-vs-normal WES to remove germline leakage.

## Data/code availability
SRA PRJNA1161693; reference ENSEMBL GRCh38 r111; env 08_Reproducibility/gbm-pub.yml;
repo GBM-Personalized-NeoAntigen-Vaccine.

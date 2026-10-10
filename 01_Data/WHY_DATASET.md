# 01_Data — why this dataset

## Objective of the pipeline
Identify and rank personalised neoantigen candidates for glioblastoma from
tumor RNA-seq, with every filtering step evidenced and every negative result reported.

## Why PRJNA1161693 (Gliocidin / PDX-161010-1 + temozolomide)
- Disease-relevant: patient-derived xenograft glioblastoma under standard-of-care
  chemotherapy (temozolomide) — the exact clinical context where a personalised
  vaccine would be given.
- Controlled: both runs share BioSample SAMN43786660 (technical replicates), so
  variant calls can be checked for replicate consistency without confounding
  inter-patient HLA differences.
- Practical: ~33M spots total fits academic hardware (this rebuild ran on 8 CPU /
  15 GB RAM / 151 GB free disk); full-genome STAR indexing (~30 GB RAM) does not,
  hence the chr1-pilot + HISAT2-full design documented in 08_Reproducibility/METHODS.md.

## Why chr1 pilot
- The original work aligned to GRCh38 chromosome 1 only. We reproduce it exactly
  (13.9% uniquely mapped, 3 TIMES: Galaxy 13.90 / local 899 13.91 / local 902 13.83)
  and show why 13% is EXPECTED (chr1 ≈ 8% of genome), correcting the original
  "Excellent" claim. Full-genome extension is scripted (03_Alignment/run_align.sh).

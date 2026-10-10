# Workspace organization (2026-10-10)

```
/mnt/nvme0n1p5/AGY/previous/
├── Last/              33G  READ-ONLY provenance (user's previous Galaxy work). NEVER delete/edit.
│   ├── galaxy/             473M  FreeBayes raw VCF (165,818), SnpSift table, final_results.csv
│   ├── Galaxy_Pipelines/   163M  annotated VCF, candidates/batch/high-affinity CSVs, 5 py scripts
│   └── Pipeline_outcome_Data/  ~32G  per-sample FastQC/BAM/VCF/HLA (large binaries stay here)
├── rebuild/           8.7G  active publication build (this is what syncs to repo/)
│   ├── data/          3.3G  SRA + FASTQ (both samples; excluded from git)
│   ├── results/       3.1G  QC/align/variants/sdf4-validation (CSVs + summaries to git; BAMs excluded)
│   ├── refs/          2.4G  chr1 FASTA + STAR index (excluded from git; documented)
│   ├── evidence_prev/ 7.6M  curated previous-work evidence (IN git)
│   ├── figures/       800K  Fig1-7 PNG (IN git)
│   ├── docs/          METHODS, AUDIT, VALIDATION, EVIDENCE_MANIFEST (IN git)
│   ├── config/ workflow/ CITATION.cff LICENSE (IN git)
├── repo/              15M  git checkout (pushed to GitHub main @84b60b8)
├── .conda_envs/       4.5G  gbm-pub env (python3.10 + tools + mhcflurry; on /mnt, root disk was full)
└── .conda_pkgs/       2.3G  conda cache (cleaned 900M tarballs 2026-10-10)
```

## Removed as irrelevant
- `rebuild/**/ _STARtmp/` (STAR temp dirs), `/tmp/opencode` scratch (hc/fb pos lists),
  top-level empty `refs/` + `logs/` (consolidated to `rebuild/logs/`), 900 MB conda tarballs.
- Excluded from git (kept on disk, listed in EVIDENCE_MANIFEST): BAMs, FASTQs, .sra,
  STAR index binaries, VEP tar, .idx/.tbi, 82M raw counts, file1 duplicates, PI_request.txt.

## What syncs to GitHub
`repo/` = curated subset of `rebuild/` (docs, workflow, config, figures, evidence_prev,
results/*.csv + Log.final/flagstat/stats summaries). Large binaries never staged (verified: no files >50M).

# Audit — v1 → publication rebuild

Source: `https://github.com/Nurulamin6191/00_Discovery_Neoantigens_Glioblastoma` (README-only, 8 commits).
Samples: SRR30677899 (15,434,896 reads) + SRR30677902 (17,969,839 reads), PRJNA1161693,
"Gliocidin is a nicotinamide mimetic prodrug that targets glioblastoma", PDX-161010-1 + temozolomide,
NovaSeq 6000 paired RNA-seq. Both runs map to same BioSample SAMN43786660 (replicates, not tumor/normal pair).

## Gaps (each now ticketed in rebuild)

| # | v1 claim / state | Problem | Fix in rebuild |
|---|---|---|---|
| 1 | "~13% uniquely mapped (Excellent for chr1 pilot)" | Expected: chr1 ≈ 8% of genome; 13% reflects off-target mapping, not quality | HISAT2 full-genome primary + STAR-chr1 pilot; flagstat both (`02_align.sh`), Fig 3 |
| 2 | GATK HaplotypeCaller → somatic variants | Germline caller, tumor-only, no normal/PoN, no dbSNP filter | RNA-seq best practices + hard filters + dbSNP flag; limitation stated (`03_variants.sh`, METHODS) |
| 3 | VCF → 2 peptides | Missing VEP-missense → 8–11mer tiling bridge | `vep_to_peptides.py` tiles all missense 8–11mers; Fig 1 yellow box |
| 4 | HLA: A*02:06/A*30:01/B*08:01/B*53:01/C*07:01 | Only 5 alleles; 2nd HLA-C missing; no version/coverage | `05_hla.sh` re-derives diploid 6-allele OptiType call |
| 5 | MHCflurry on A*02:06: LLFGYPVYV 10.31 nM, SYFPEITHI 746.85 nM "negative" | Single allele, 2 peptides; SYFPEITHI is a classic A*02:01 binder control | Full-allele 8–11mer screen, 50/500 nM stratification, NetMHCpan validation (`06_predict.sh`), Fig 2 |
| 6 | No env/versions/logs/figures | Unreproducible | `workflow/envs/gbm-pub.yml`, pinned scripts, MultiQC, Figs 1–3 |

## Status
- [x] Scaffold, scripts 00–06, tiling + ranking scripts, METHODS, Figs 1–3 (preliminary values)
- [ ] Lean conda env (`gbm-pub` on /mnt disk) + chr1 FASTA downloading
- [ ] FASTQ fetch → QC → align → variants → VEP → HLA → MHCflurry → final figures/tables → push

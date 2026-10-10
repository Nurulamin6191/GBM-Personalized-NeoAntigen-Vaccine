# VALIDATION — end-to-end (Galaxy previous work + local re-run)

## Chain of custody (verified counts)
| Step | Source | Value |
|---|---|---|
| Raw reads | ENA / fasterq-dump | SRR30677899: 15,434,896 spots / 30,869,792 reads; SRR30677902: 17,969,839 spots / 35,939,678 reads |
| Trim (Galaxy, single-end Q20) | trimming reports | 899: R1 24.6% adapters; 902: R1 11.9%, R2 40.2%; pairs removed <1% |
| STAR chr1 | Galaxy Log.final vs local re-run | Galaxy 13.90% uniquely mapped / local 899: 13.91% / local 902: 13.83% — MATCH (too-short ~83.9%). Validates reproducibility |
| Splices | Galaxy (with GTF) vs local (no GTF) | 159,301 annotated vs 0 — expected difference; mapping rate unaffected |
| FreeBayes raw | candidate.vcf / SRR_annotated.vcf | 165,818 variants |
| Missense | grep missense_variant | 4,224 (≈4.2k; report said 4,265 — 41 multi-annotation lines) |
| Tiling | batch_peptides.csv | 36,203 9mers from 4,115 candidates |
| Strong binders | high_affinity_candidates.csv | 1,953 (<500 nM, HLA-A*02:06) |
| Focus test | focused_SDF4_vs_control.csv | LLFGYPVYV 10.31 nM; ALAAAVAAV 16.94 nM |

## Resolved issues
- **HLA-C "missing" allele:** OptiType result is C*07:01 / C*07:01 (homozygous, 29 reads). Diploid profile complete.
- **13% mapping:** expected for chr1-only; flagstat 37% properly paired / 63% singletons confirms.

## Fresh MHCflurry 2.1.2 re-run (this session, models_class1, env gbm-pub)
- Galaxy values reproduced to 5 decimals: LLFGYPVYV 10.31053, ALAAAVAAV 16.94121 (HLA-A*02:06).
- SDF4 A120V 9mer tiling × full diploid profile (`results/sdf4_validation_out.csv`,
  `results/sdf4_validation_other_alleles.csv`, Fig 7): all A*02:06 >6,800 nM (non-binders);
  best overall EMQRWVMEK/A*30:01 = 192.05 nM (weak); QRWVMEKTA/C*07:01 = 857.77 nM.
  No strong binder. SDF4-as-9mer-vaccine-candidate is REFUTED for this allele set;
  10/11-mers + expression follow-up remain open, not claimed.
- ALAAAVAAV binds ONLY A*02:06 (elsewhere >3,500 nM); gene origin still unattributed
  (0 hits in batch_peptides/high_affinity, absent from SDF4). Do not label it SDF4.

## Flags (honest limitations, must stay in publication)
1. **ALAAAVAAV provenance UNCONFIRMED:** absent from SDF4 protein (all isoforms checked via Ensembl REST)
   and from batch tilings. The validated SDF4 tiling mutant is SAKEMQRWV-family (I120V window).
   ALAAAVAAV 16.94 nM must be re-derived or relabelled, not presented as SDF4 epitope without proof.
2. **SDF4 call quality LOW:** chr1:1223936 G>A QUAL≈7e-9 (FreeBayes ~0), DP=24. Requires IGV/Sanger validation.
3. **Ranking honesty:** batch top hits (FLIGSTSYV 10.14 nM, …) outrank SDF4 focus (16.94 nM).
   SDF4 is biologically prioritized (PI-validated), not top-affinity. Table 1 must show full top-20.
4. **Tumor-only:** no matched normal; germline leakage expected. Single-allele screen (A*02:06);
   full-diploid + NetMHCpan validation pending.
5. **R1=29bp:** asymmetric SMART-seq-type library; min_len 50 would delete R1 (fixed to 20).

## Hardware notes
- Root `/` filled (137G, 99%) during conda solve → env moved to `/mnt` (pkgs+envs dirs).
- STAR full-genome index needs ~30GB RAM (have 15GB) → chr1 STAR pilot + HISAT2-full deferred.
- VEP via Docker deferred (conda VEP pulls full Perl stack); SnpEff/FreeBayes Galaxy outputs used as evidence.
- Dataset unchanged (both SRRs); no downsampling needed — 151G free on /mnt sufficed.

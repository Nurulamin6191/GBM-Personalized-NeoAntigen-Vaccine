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

## Local GATK validation (this session, SRR30677899, chr1)
- HaplotypeCaller raw 43,461 → PASS 8,260 (7,832 SNPs + 430 indels) after QD<2/FS>30/DP<10 filtration.
- Overlap with FreeBayes raw (165,818, chr-normalized): **1,568 / 8,260 (19%)**.
  Low overlap expected: different callers + ENSEMBL "1" vs UCSC "chr1" reference naming (harmonization TODO).
- **SDF4 chr1:1223936 absent from GATK raw AND PASS; zero GATK calls within ±1 kb.**
  FreeBayes-only call at QUAL≈7e-9. Two-caller non-replication — SDF4 A120V is FreeBayes-only, low-confidence.

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

## Replicate HLA + per-sample re-score (this session)
- OptiType 1.5.0 RNA re-run: 899 EXACT match to prior (A*02:06/A*30:01 B*08:01/B*53:01
  C*07:01 hom, 29 reads, obj 28.478). 902 (20 reads): A*02:01/A*30:01 B*08:01/B*58:01
  C*07:01 hom. Same BioSample, different call -> RNA-based HLA unstable at 20-29 reads.
  Primary profile = 899 consensus (higher reads, independently reproduced 2x);
  DNA-based HLA confirmation listed as required follow-up.
- BUG FIX: per-sample high-affinity tables were run vs HLA-A*02:01 (script default),
  not patient A*02:06. Re-scored all 453+430 peptides vs A*02:06 locally
  (06_Neoantigens/per_sample_A0206/): values reproduce pooled run; 899∩902 share 155 peptides.
- IEDB tables = independent NetMHCpan EL/BA validation (BA IC50 2.5-12 nM, all strong),
  largely non-overlapping peptide set -> orthogonal support (Fig 8).

## Genome-wide VCF provenance flag (this session)
- FreeBayes VCF spans 165 contigs genome-wide (chr1: 14,411; chr2: 13,194; ...),
  yet Galaxy BAMs in that folder are chr1-only (single @SQ contig 1).
  The VCF therefore came from a FULL-genome alignment whose BAM is not archived.
  Consequence: variant evidence is genome-wide (good for discovery) but alignment
  params for the VCF are unrecoverable -> our local chr1 GATK chain is the only
  fully-provenanced variant call. Full-genome HISAT2 re-alignment (refs/full/)
  will close this gap.

## Proteome attribution (this session, UniProt UP000005640 reviewed, 20,652 proteins)
- SAKEMQRWI (wild-type SDF4 window) -> CAB45_HUMAN (SDF4 gene product). Tiling provenance CONFIRMED.
- Mutant SAKEMQRWV -> 0 hits (expected: somatic mutant absent from reference).
- ALAAAVAAV -> 0 hits: non-human/synthetic, gene unattributed. Treat as synthetic binder control.
- LLFGYPVYV -> 0 hits: synthetic positive control (A*02:06-restricted, reproduced 10.31 nM).
- SYFPEITHI -> JAK1_HUMAN exact 9mer: must NOT be used as negative control (v1 error confirmed).
  (06_Neoantigens/proteome_attribution.tsv)

## Germline filter + Table 1 (this session)
- dbSNP rsID filter reconciled: 4,265 = 2,126 COMMON proxy (49.9%) + 1,815 RARE (42.6%) + 324 UNRESOLVED.
- SDF4 chr1:1223936 carries rs1384393162 -> COMMON germline. THIRD strike (QUAL~0, GATK-absent, germline): DEFINITIVELY EXCLUDED.
- IEDB-validated YLIDSAHKI (GSTM1 N85S) also germline-filtered out (GSTM1 famously polymorphic) -> filter demonstrably matters.
- Table 1 (06_Neoantigens/Table1_top10_rare_binders.csv): 1,083 rare strong binders; top-10 ECT2/CUL4A/AP2M1/MACROH2A1/UTP14C/MCUB/SLC16A4/CYB561/VPS13B/CAND1, all <13 nM. Expression (StringTie TPM) + Sanger still required.

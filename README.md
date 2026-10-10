# In Silico Identification of Neoantigens for Personalised Immunotherapy in Glioblastoma

> Repository renamed to `GBM-Personalized-NeoAntigen-Vaccine` (short, searchable;
> old URL redirects automatically).

## Objective
Build an **evidence-first, end-to-end validated** neoantigen discovery pipeline for
glioblastoma from tumor RNA-seq: raw reads → QC → alignment → variant calling →
HLA typing → peptide tiling → binding prediction → ranked candidates — with every
filtering step counted, every negative result reported, and every figure
regenerated from checked-in tables (07_Figures/).

## Why this dataset
- **PRJNA1161693** (Gliocidin PDX glioblastoma + temozolomide, the clinical setting
  for a personalised vaccine): SRR30677899 (15.4M spots) + SRR30677902 (18.0M),
  Illumina NovaSeq 6000, paired. See 01_Data/SAMPLES.csv + WHY_DATASET.md.
- Both runs share BioSample SAMN43786660 (replicates) → variant calls checked for
  consistency without inter-patient HLA confounding.
- ~33M spots fits 8-CPU/15-GB-RAM hardware (full-genome STAR index needs ~30 GB,
  hence chr1-pilot + scripted HISAT2-full design).

## Why this pipeline (each tool earns its place)
| Stage | Tool | Why |
|---|---|---|
| QC | FastQC + Trim Galore Q20/len20/2-colour | R1=29bp/R2=89bp SMART-seq-type; len50 would delete R1; NovaSeq poly-G + SMARTer oligo in R2 |
| Alignment | STAR chr1 pilot (reproduced 3×: 13.90/13.91/13.83%) | reproduces original exactly; proves 13% is expected chr1-only behaviour, not "Excellent" |
| Variants | FreeBayes (Galaxy, 165,818 → 4,224 missense) + GATK HC local (43,461 → 8,260 PASS) | two-caller check; 19% overlap; SDF4 call is FreeBayes-only QUAL~0 |
| Peptides | missense → 8–11mer tiling (36,203 9mers) | the bridge v1 skipped (it tested 2 hand-picked peptides) |
| HLA | OptiType → A*02:06/A*30:01 B*08:01/B*53:01 C*07:01×2 | C is homozygous — resolves v1's "missing" 6th allele |
| Binding | MHCflurry 2.1.2 full-diploid, 50/500 nM IEDB cutoffs | single-allele 2-peptide test → 60 peptide-allele matrix |

## Honest headline results
- Batch screen top hits bind strongly (e.g. FLIGSTSYV 10.14 nM); full top-20 in Fig4.
  They await gene-annotation + expression filtering for Table 1 — not claimed as vaccines.
- **SDF4 A120V 9mers do NOT bind HLA-A*02:06** (all >6,800 nM; best overall
  EMQRWVMEK/A*30:01 192 nM weak). SDF4 call itself is FreeBayes-only, absent from
  GATK ±1 kb. The 9mer vaccine claim is refuted two independent ways (Fig7).
- ALAAAVAAV (16.94 nM, A*02:06-only) reproduced exactly but is **not SDF4-derived**
  (absent from SDF4 protein + all tilings) — gene attribution required.
- SYFPEITHI is a known A*02:01 binder control, not a "negative control" (v1 error).

## Layout
01_Data → 02_QC → 03_Alignment → 04_Variants → 05_HLA → 06_Neoantigens →
07_Figures → 08_Reproducibility (config, env, METHODS, AUDIT, VALIDATION, manifests).
Run stage scripts from repo root (`bash 01_Data/fetch.sh`, …); runtime dirs
(data/, refs/, results-RUNTIME/) are git-ignored and documented.

## Reproduce
1. `conda env create -f 08_Reproducibility/gbm-pub.yml` (or prefix-install on /mnt if root is full)
2. `bash 01_Data/fetch.sh && bash 02_QC/run_qc.sh && bash 03_Alignment/run_align.sh`
3. `bash 04_Variants/run_variants.sh && bash 05_HLA/run_hla.sh && bash 06_Neoantigens/run_predict.sh`
4. `python3 07_Figures/make_figures_prev.py && python3 07_Figures/make_fig07.py`

## Limitations
Tumor-only (no matched normal), chr1 pilot (full-genome scripted, not run here),
FreeBayes/GATK 19% overlap (caller + reference-naming differences), SDF4 needs
IGV/Sanger confirmation, NetMHCpan cross-validation pending.

## Repository address
Canonical: `https://github.com/Nurulamin6191/GBM-Personalized-NeoAntigen-Vaccine`.
Then: `git remote set-url origin https://github.com/Nurulamin6191/GBM-Personalized-NeoAntigen-Vaccine.git`.

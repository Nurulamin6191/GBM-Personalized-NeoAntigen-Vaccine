# Alternative best methodology (what the field would use today)

Reference points: Keskin et al. Nature 2019 (GBM NeoVax), Hilf et al. Nature 2019 (GAPVAC),
Wells et al. Cell 2020 (TESLA), pVACtools (Griffith lab), 2023 RNA-HLA benchmarks
(BMC Genomics; bioRxiv 2023.05.22), 2024 GBM field (Nat Commun 173-patient real-world,
CVGBM mRNA, DNA+PD-1 NCT05743595, O.R.A.N/TOFU mRNA platform).

## Per-stage best practice vs this repo

| Stage | Field standard (2024–2026) | This repo | Gap |
|---|---|---|---|
| Input | Tumor + matched-normal WES + tumor RNA (Keskin; TESLA) | Tumor RNA only, no normal | FUNDAMENTAL — needs new sequencing |
| Alignment | STAR 2-pass / HISAT2 + GTF, full genome | STAR chr1 pilot + HISAT2 full (44.7%) | pilot scope; full BAMs retired |
| Variant calling | Mutect2+Strelka+VarScan+Pindel ensemble, gnomAD AF filter, bam-readcount VAF (pVACtools) | FreeBayes + GATK-HC, rsID proxy filter | no ensemble/Mutect2; proxy AF |
| Annotation | VEP + VAtools (coverage/VAF/expression/phasing), MANE/canonical/TSL | VEP REST gene map (canonical), TSL unused | TSL/MANE not enforced |
| Expression | StringTie/Kallisto TPM + RNA VAF; allele expr = TPM×VAF ≥ 2.5 (pVACview) | StringTie TPM ≥ 1 (gene-level) | no RNA-VAF allele expr |
| HLA | Clinical PCR typing (Keskin); else OptiType (best class-I RNA, 99.2%) + arcasHLA/HLA-HD for class II (benchmarks) | OptiType class I only, RNA, 20–29 reads | no class II at all; low depth |
| Peptides | 8–11mers class I + 12–18mers class II, anchor analysis (pVACtools) | 9mers class I only | length + class II gap |
| Binding | Ensemble (NetMHCpan + MHCflurry + MixMHCpred + PRIME), allele-specific cutoffs, %rank ≤ 2 | MHCflurry + NetMHCpan BA, 50/500 + TESLA 34 nM audit | ensemble of 2; no PRIME/MixMHCpred |
| Recognition | Agretopicity (DAI, WT/MT), foreignness, stability (TESLA); reference-proteome filter (pVACtools) | k-mer proteome attribution only; DAI/foreignness/stability PENDING | open, scripted as next step |
| Prioritisation | Clonality (PyClone), multi-allele binding, ML tiering (pVACview Pass tiers) | Affinity + TPM + novelty rank | no clonality/ML tier |
| Validation | ELISpot/multimer, TILs, scTCR-seq (Keskin); 90% response real-world (Nat Commun 2024) | None (computational only) | wet lab required |
| Delivery context | Long peptides ≤20 (Keskin), mRNA ≤34 neoantigens (mRNA-4157), +PD-1 combos | 9mer list only | no construct design |

## Recommended upgrade path (in order)
1. Tumor–normal WES → Mutect2 ensemble + VEP/VAtools (closes germline + clonality).
2. arcasHLA (RNA, class I+II) cross-check + clinical HLA before any ranking is acted on.
3. 8–11mer + class-II screen, ensemble predictors, allele-specific cutoffs, anchor filter.
4. Agretopicity/foreignness/stability scoring (TESLA recognition half).
5. Long-peptide/mRNA construct design (≤20–34 epitopes) + PD-1 combo rationale.
6. ELISpot/multimer + TIL/scTCR validation (the Nature-level requirement).

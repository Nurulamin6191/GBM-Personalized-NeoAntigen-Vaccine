# Benchmark vs Nature/Cell-level literature — where this project stands

## Reference set
- **Keskin et al., Nature 2019** (GBM NeoVax, NCT02287428): tumor+normal WES, tumor RNA,
  clinical HLA, NetMHCpan 2.4 IC50<500, ≤20 long peptides, CD4+CD8 responses, TIL
  infiltration, pembro combination cohorts.
- **Hilf et al., Nature 2019** (GAPVAC-101): warehouse TAAs + neoepitopes,
  HLA-A*02:01/A*24:02, transcriptome + immunopeptidome, 30–50 GBM mutations, CD8/CD4 responses.
- **2024 GBM field**: 173-patient real-world peptide vaccines (Nat Commun, 90% response,
  multi-response 53 mo OS); CVGBM mRNA; DNA+PD-1 NCT05743595; O.R.A.N/TOFU mRNA platform.
- **Wells et al., Cell 2020 (TESLA)**: presentation (affinity, abundance, stability) +
  recognition (agretopicity, foreignness); 34 nM filter; 98% non-immunogenic removed at
  0.70 precision; 6% validation rate baseline.
- **HLA benchmarks 2022–2023**: OptiType best class-I RNA (99.2%); arcasHLA best overall
  RNA + class II; HLA-HD class II.

## Element-by-element verdict

| Element | Literature bar | This project | Verdict |
|---|---|---|---|
| Data pairing | Tumor+normal WES + RNA | Tumor RNA only | GAP — needs sequencing |
| Replicate design | Multi-patient cohorts | 2 runs, 1 BioSample, cross-checked | Adequate for methods paper |
| Alignment QC | Full-genome, GTF-aware | chr1 pilot + full HISAT2 w/ rates | MEETS (as pilot + audit) |
| Variant evidence | Ensemble callers, VAF filters | Two callers, 19% overlap disclosed | APPROACHES (honest discordance) |
| Germline control | Matched normal | rsID proxy, 49.9% removed | PARTIAL (proxy, stated) |
| HLA | Clinical typing | OptiType RNA, benchmark-backed, instability disclosed | APPROACHES for class I; class II MISSING |
| Binding | Ensemble + allele cutoffs | MHCflurry + NetMHCpan API, 34 nM TESLA audit 10/10 | MEETS computationally |
| Recognition | DAI/foreignness/stability | PENDING (scripted next) | GAP |
| Expression | TPM + RNA VAF allele expr | TPM gate 1,047/1,083 | APPROACHES (no allele expr) |
| Negative results | Rarely published | SDF4 triple-strike exclusion documented | EXCEEDS (audit-grade honesty) |
| Reproducibility | Code + env (pVACtools standard) | 11 figures regenerable, env pinned, audit log | MEETS |
| Immune validation | ELISpot/multimer/TILs/scTCR | None | GAP — wet lab required |

## What would make this Nature/Cell-level
1. Tumor–normal WES on the same PDX (closes germline + clonality + VAF).
2. Class-II HLA (arcasHLA now; clinical typing before action) + 8–18mer screen.
3. Recognition half: agretopicity/foreignness/stability per TESLA.
4. ELISpot + TIL/scTCR on Table 1 v2 (the non-negotiable clinical-journal requirement).
5. Framing as a **methods-and-reproducibility study with an audited negative result**
   (SDF4) is itself publishable now at the bioinformatics-journal tier
   (e.g. BMC Bioinformatics, Briefings in Bioinformatics, Frontiers) — the Nature/Cell
   tier additionally requires 1–4.

## Positioning statement (for cover letter)
Unlike trial reports, this work contributes an end-to-end reproducible audit: it
reproduces a prior claim three ways, refutes it four ways (QUAL, caller, germline,
binding ×2 predictors), replaces it with a germline-filtered, expression-gated,
TESLA-audited Table 1, and states every remaining boundary. That audit discipline —
not a new vaccine — is the manuscript's claim.

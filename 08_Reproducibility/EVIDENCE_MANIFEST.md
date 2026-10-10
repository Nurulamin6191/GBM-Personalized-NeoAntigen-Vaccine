# Evidence manifest — curated publication repo (numbered stages)

Provenance note: the original `Last/` work folder was removed during consolidation;
its curated evidence below is what survives, plus public sources (SRA PRJNA1161693,
ENSEMBL r111, UniProt UP000005640, Ensembl REST) that regenerate the rest.

| Repo path | Origin | Note |
|---|---|---|
| 01_Data/SAMPLES.csv + WHY_DATASET.md + fetch.sh | ENA + SRA toolkit | 33M spots total |
| 02_QC/*_trimming_report.txt | local Trim Galore paired runs (run_qc.sh) | KB | adapter rates feed Fig6 |
| 03_Alignment/SRR_starLog.* + markdup + Log.final + flagstats + hisat2 rates | Galaxy STAR + local STAR/HISAT2 re-runs | 13.9% chr1 ×3; 44.7% full-genome |
| 04_Variants/variant_counts.tsv + SDF4 excerpt + pileup + PASS.stats + missense_dbsnp.tsv | FreeBayes VCF (165,818) + local GATK + VEP REST | 4,224 missense; SDF4 triple-strike |
| 04_Variants/missense_dbsnp.tsv | VEP REST rsID filter script (filter_gnomad.py) | 4,265 rows | reconciled 2126/1815/324 |
| 05_HLA/result.tsv + coverage PDFs + HLA PNGs + optitype_rerun/ | prior OptiType + local OptiType 1.5.0 rerun | 899 exact; 902 unstable (documented) |
| 06_Neoantigens/candidates + batch + high_affinity + focused | Galaxy Python pipeline | 36,203 9mers; 1,953 binders |
| 06_Neoantigens/IEDB/ | IEDB NetMHCpan EL/BA runs | orthogonal validation |
| 06_Neoantigens/per-sample + per_sample_A0206/ | per-sample runs re-scored vs A*02:06 | A*02:01-default bug fixed |
| 06_Neoantigens/sdf4_validation_*.csv + proteome_attribution.tsv | fresh MHCflurry 2.1.2 + UniProt | SDF4 tiling negative; ALAAAVAAV synthetic |
| 06_Neoantigens/Table1*.csv + build/gate/fetch scripts | germline filter + StringTie TPM + VEP gene map | 1,083 rare; v2 expressed |
| 06_Neoantigens/rare_genes.tsv | VEP REST canonical mapping (fetch_rare_genes.py) | 1,364 genes | feeds both Table builds |
| 06_Neoantigens/fetch_iedb.py + IEDB_api/ | script-driven NetMHCpan via IEDB API | replaces web-only tables |
| 07_Figures/fig01-10 + generators | all regenerated from checked-in tables | verified reproducible |
| 08_Reproducibility/config.yaml + gbm-pub.yml + METHODS/AUDIT/VALIDATION/ORGANIZATION | pipeline records | |
| 09_Manuscript/Report_Final.docx + manuscript_outline.md | author drafts | wet-lab confirmations listed |

## EXCLUDED (unpushable / irrelevant)
- *.bam/*.bai, *.fastq.gz, *.sra, STAR/HISAT2 indexes, VEP cache tar, .idx/.tbi,
  82M raw gene counts (superseded by StringTie TPM), 155M file1 duplicates,
  PI_request.txt (process advice), submission zips, personal SOP odt, 0-byte .Rhistory.
- Regenerate via stage scripts (01_Data/fetch.sh … 06_Neoantigens/run_expression.sh).

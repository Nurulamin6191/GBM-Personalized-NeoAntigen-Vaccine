# Evidence manifest — previous work (Last/) → curated (evidence (now 01-06 folders)/)

Source is READ-ONLY provenance; this repo carries only pushable curation (<10 MB).
Large originals stay in Last/ and are referenced here, excluded from git.

| Curated path | Source in Last/ | Size | Note |
|---|---|---|---|
| 01_QC/*/trimming_report.txt + *_fastqc.html | Pipeline_outcome_Data/*/Raw_qc/ | 7 MB | zips excluded (regenerable) |
| 02_Align/SRR_starLog.* + marked_dup_metrics | .../SRR30677899/Alignment_records + .txt | KB | SJ table included |
| 03_Variants/variant_counts.tsv + SDF4 excerpt | Galaxy_Pipelines/SRR_annotated.vcf (165,818) | bytes | counts script in VALIDATION |
| 04_HLA/result.tsv + coverage_plot.pdf | .../hla_result/2026_01_19_14_55_21/ | 108 KB | proves C*07:01 homozygous |
| 05_Predict/candidates + batch + high_affinity + focused | Galaxy_Pipelines/*.csv + galaxy/final_results.csv | 1.5 MB | |
| scripts/*.py | Galaxy_Pipelines/*.py | KB | |

## EXCLUDED as irrelevant / unpushable (kept in Last/, NOT in repo)
- PI_request.txt (AI process advice, not data)
- file1.csv / file1.txt (155 MB duplicates of candidate.vcf)
- *.bam / *.bai (~GBs), *.fastq / *.fq.gz (GBs), star_index/ binaries (2 GB), tmp/vep tar
- sorted.gtf.gz(.tbi), *.idx, _STARtmp/, SRR_annotated.vcf 8B empty stub + 101M .vcf.vcf raw
- gene_counts.txt 82M raw (needs summarising to TPM table — TODO)
- Data/ fastq duplicates of rebuild/data/fastq

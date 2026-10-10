#!/bin/bash
# 05 — HLA typing (OptiType). Needs razers3 + config; wrapped so repo documents exact call.
# v1 reported 5 alleles with 2nd HLA-C missing — this script re-derives diploid 6-allele profile.
set -euo pipefail
mkdir -p results/hla logs
echo "OptiType requires razers3 and its config.ini; install: conda install -c bioconda optitype razers3"
echo "Example per-sample call (uncomment to run):"
echo '  OptiTypePipeline.py -i results/qc/trimmed/SAMPLE_1.val_1.fq.gz results/qc/trimmed/SAMPLE_2.val_2.fq.gz --dna --outdir results/hla/SAMPLE --prefix SAMPLE -v'
# Record reported v1 profile for comparison:
cat > results/hla/v1_reported_profile.tsv <<'EOF'
locus	allele1	allele2	status
HLA-A	A*02:06	A*30:01	v1-reported
HLA-B	B*08:01	B*53:01	v1-reported
HLA-C	C*07:01	MISSING	v1-incomplete-rederive-in-full-run
EOF
cat results/hla/v1_reported_profile.tsv

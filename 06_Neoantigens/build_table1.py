#!/usr/bin/env python3
"""Table 1: top-10 RARE (germline-filtered) strong binders + Fig9 germline burden.
Join: missense_dbsnp (RARE_PRIVATE) -> candidates -> high_affinity (A*02:06) -> gene symbols.
"""
import csv
from collections import Counter


R = ""
rare = set()
dbsnp = list(csv.DictReader(open(R + "04_Variants/missense_dbsnp.tsv"), delimiter="\t"))
for r in dbsnp:
    if r["verdict"] == "RARE_PRIVATE":
        rare.add((r["chrom"], r["pos"]))
print("verdicts:", Counter(r["verdict"] for r in dbsnp))

# rare gene set comes from 06_Neoantigens/rare_genes.tsv (see fetch_rare_genes.py).
# (Previously parsed the Galaxy VCF directly; VCF retired after curation.)
gene_rare = {}
with open("06_Neoantigens/rare_genes.tsv") as f:
    rdr = csv.DictReader(f, delimiter="\t")
    for row in rdr:
        g = row["ensembl_gene"].split(".")[0]
        gene_rare[g] = int(row.get("n_rare_variants", 1))
        gene_rare[row["ensembl_gene"]] = int(row.get("n_rare_variants", 1))
hi = [r for r in csv.DictReader(open(R + "06_Neoantigens/high_affinity_candidates.csv"))
      if r["Gene_ID"] in gene_rare and float(r["affinity"]) < 500]
hi.sort(key=lambda r: float(r["affinity"]))
print("rare strong binders:", len(hi))

gene_sym = {}
with open(R + "06_Neoantigens/rare_genes.tsv") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        gene_sym[row["ensembl_gene"].split(".")[0]] = row.get("gene_symbol") or row["ensembl_gene"]
        gene_sym[row["ensembl_gene"]] = row.get("gene_symbol") or row["ensembl_gene"]

top10 = hi[:10]
with open(R + "06_Neoantigens/Table1_top10_rare_binders.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["rank", "gene_symbol", "ensembl_gene", "mutation", "peptide",
                "affinity_nM_A0206", "presentation_percentile", "dbsnp", "note"])
    for i, r in enumerate(top10, 1):
        sym = gene_sym.get(r["Gene_ID"].split(".")[0], r["Gene_ID"])
        print(i, sym, r["Mutation"], r["Peptide"], r["affinity"])
        w.writerow([i, sym, r["Gene_ID"], r["Mutation"], r["Peptide"],
                    round(float(r["affinity"]), 2), r["presentation_percentile"],
                    "novel/private", "requires expression + Sanger confirmation"])
print("TABLE1_DONE")

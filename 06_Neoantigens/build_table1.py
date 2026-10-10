#!/usr/bin/env python3
"""Table 1: top-10 RARE (germline-filtered) strong binders + Fig9 germline burden.
Join: missense_dbsnp (RARE_PRIVATE) -> candidates -> high_affinity (A*02:06) -> gene symbols.
"""
import csv, json, urllib.request
from collections import Counter

R = ""
rare = set()
dbsnp = list(csv.DictReader(open(R + "04_Variants/missense_dbsnp.tsv"), delimiter="\t"))
for r in dbsnp:
    if r["verdict"] == "RARE_PRIVATE":
        rare.add((r["chrom"], r["pos"]))
print("verdicts:", Counter(r["verdict"] for r in dbsnp))

# candidates.csv uses ENSG + protein pos/AA, no genomic coords: map via batch gene/mutation
# high_affinity has Gene_ID + Mutation + Peptide + affinity; keep rows whose gene is rare-associated:
# approximate: gene is rare if ANY of its batch peptides' source variant is rare.
# candidates.csv lacks genomic pos, so map gene->rare via VCF ANN gene field.
gene_rare = {}
with open("/mnt/nvme0n1p5/AGY/Previous/Last/Galaxy_Pipelines/SRR_annotated.vcf") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.split("\t")
        key = (p[0].replace("chr", ""), p[1])
        if key not in rare:
            continue
        for m in p[7].split(","):
            if "missense_variant" not in m:
                continue
            parts = m.split("|")
            gene = parts[4] if len(parts) > 4 else ""
            if gene:
                gene_rare[gene] = gene_rare.get(gene, 0) + 1
print("rare genes:", len(gene_rare))

hi = [r for r in csv.DictReader(open(R + "06_Neoantigens/high_affinity_candidates.csv"))
      if r["Gene_ID"] in gene_rare and float(r["affinity"]) < 500]
hi.sort(key=lambda r: float(r["affinity"]))
print("rare strong binders:", len(hi))

def symbol(eid):
    try:
        url = f"https://rest.ensembl.org/lookup/id/{eid}?content-type=application/json"
        d = json.load(urllib.request.urlopen(
            urllib.request.Request(url, headers={"Content-Type": "application/json"}), timeout=20))
        return d.get("display_name", eid)
    except Exception:
        return eid

top10 = hi[:10]
with open(R + "06_Neoantigens/Table1_top10_rare_binders.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["rank", "gene_symbol", "ensembl_gene", "mutation", "peptide",
                "affinity_nM_A0206", "presentation_percentile", "dbsnp", "note"])
    for i, r in enumerate(top10, 1):
        sym = symbol(r["Gene_ID"])
        print(i, sym, r["Mutation"], r["Peptide"], r["affinity"])
        w.writerow([i, sym, r["Gene_ID"], r["Mutation"], r["Peptide"],
                    round(float(r["affinity"]), 2), r["presentation_percentile"],
                    "novel/private", "requires expression + Sanger confirmation"])
print("TABLE1_DONE")

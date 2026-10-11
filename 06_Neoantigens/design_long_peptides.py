#!/usr/bin/env python3
"""Keskin-style long-peptide design: 25mers centered on each Table 1 v2 mutation.
Locates the 9mer inside the Ensembl protein, expands to 25aa (or protein bounds).
Out: 06_Neoantigens/long_peptides_25mer.csv
(Delivery-format step toward trial-grade design; binding of longmers not predicted here.)
"""
import csv
import json
import time
import urllib.request

IN = "06_Neoantigens/Table1v2_expressed.csv"
OUT = "06_Neoantigens/long_peptides_25mer.csv"
HALF = 12


def protein(eid):
    url = f"https://rest.ensembl.org/sequence/id/{eid}?type=protein;multiple_sequences=1;content-type=application/json"
    for _ in range(4):
        try:
            d = json.load(urllib.request.urlopen(
                urllib.request.Request(url, headers={"Content-Type": "application/json"}),
                timeout=60))
            d.sort(key=lambda x: len(x["seq"]), reverse=True)
            return d[0]["seq"]
        except Exception as e:
            print("retry", eid, e, flush=True)
            time.sleep(10)
    return None


rows = list(csv.DictReader(open(IN)))
out = []
for r in rows:
    seq = protein(r["ensembl_gene"])
    if not seq:
        out.append((r["rank"], r["gene_symbol"], r["mutation"], r["peptide"], "PROTEIN_FETCH_FAIL", ""))
        continue
    i = seq.find(r["peptide"])
    if i < 0:
        out.append((r["rank"], r["gene_symbol"], r["mutation"], r["peptide"], "9MER_NOT_FOUND", ""))
        continue
    s = max(0, i + 4 - HALF)
    longmer = seq[s:s + 25]
    out.append((r["rank"], r["gene_symbol"], r["mutation"], r["peptide"], longmer, len(longmer)))
    print(r["gene_symbol"], i, longmer, flush=True)
    time.sleep(1)

with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["rank", "gene_symbol", "mutation", "core_9mer", "long_25mer", "length"])
    w.writerows(out)
print(f"LONGPEP_DONE -> {OUT}")

#!/usr/bin/env python3
"""Keskin-style long-peptide design: 25mers centered on each Table 1 v2 mutation.
Mutant 9mers never occur in the reference proteome, so each is mapped by reverting
the mutation (L34F etc.) at every register and searching the Ensembl protein for the
wild-type 9mer; the 25mer is then built around that site with the mutant AA restored.
Out: 06_Neoantigens/long_peptides_25mer.csv
"""
import csv
import json
import re
import time
import urllib.request

IN = "06_Neoantigens/Table1v2_expressed.csv"
OUT = "06_Neoantigens/long_peptides_25mer.csv"
HALF = 12


def proteins(eid):
    url = (f"https://rest.ensembl.org/sequence/id/{eid}?type=protein;"
           "multiple_sequences=1;content-type=application/json")
    for _ in range(4):
        try:
            d = json.load(urllib.request.urlopen(
                urllib.request.Request(url, headers={"Content-Type": "application/json"}),
                timeout=60))
            d.sort(key=lambda x: len(x["seq"]), reverse=True)
            return [(x.get("id", "?"), x["seq"]) for x in d]
        except Exception as e:
            print("retry", eid, e, flush=True)
            time.sleep(10)
    return []


rows = list(csv.DictReader(open(IN)))
tilings = {}
with open("06_Neoantigens/batch_peptides.csv") as f:
    for r in csv.DictReader(f):
        tilings.setdefault(r["Gene_ID"].split(".")[0], []).append(r["Peptide"])


def mutant_contig(peps):
    """Assemble maximal tiling contig from overlapping 9mers (each shifted by 1)."""
    remaining = list(dict.fromkeys(peps))
    contig = remaining.pop(0)
    changed = True
    while changed and remaining:
        changed = False
        for p in list(remaining):
            if p[:8] == contig[-8:]:
                contig += p[-1]
                remaining.remove(p)
                changed = True
            elif p[-8:] == contig[:8]:
                contig = p[0] + contig
                remaining.remove(p)
                changed = True
    return contig


out = []
for r in rows:
    m = re.match(r"([A-Z])(\d+)([A-Z])", r["mutation"])
    isos = proteins(r["ensembl_gene"]) if m else []
    hit = None
    if m:
        ref, alt = m.group(1), m.group(3)
        mut9 = r["peptide"]
        for tid, seq in isos:
            for i, a in enumerate(mut9):
                if a != alt:
                    continue
                wt9 = mut9[:i] + ref + mut9[i + 1:]
                j = seq.find(wt9)
                if j >= 0:
                    hit = (tid, seq, j, i, wt9)
                    break
            if hit:
                break
    if not hit and m:
        # fallback: assemble the mutant tiling contig, revert each alt
        # occurrence, and anchor the wild-type contig in any isoform
        ref, alt = m.group(1), m.group(3)
        contig = mutant_contig(tilings.get(r["ensembl_gene"].split(".")[0], []))
        for tid, seq in isos:
            for k, a in enumerate(contig):
                if a != alt:
                    continue
                wtc = contig[:k] + ref + contig[k + 1:]
                j = seq.find(wtc)
                if j >= 0:
                    c9 = contig.find(r["peptide"])
                    if c9 < 0 or k < c9 or k >= c9 + 9:
                        continue
                    hit = (tid, seq, j + c9, k - c9, wtc)
                    break
            if hit:
                break
    if not hit:
        out.append((r["rank"], r["gene_symbol"], r["mutation"], r["peptide"],
                    "UNMAPPED", "", ""))
        print(r["gene_symbol"], "UNMAPPED", flush=True)
        continue
    j, i, wt9 = hit[2], hit[3], hit[4]
    seq = hit[1]
    center = j + i  # mutant residue index in protein
    s = max(0, center - HALF)
    ctx = seq[s:s + 25]
    mut25 = ctx[:center - s] + mut9[i] + ctx[center - s + 1:]
    out.append((r["rank"], r["gene_symbol"], r["mutation"], r["peptide"], mut25, len(mut25), hit[0]))
    print(r["gene_symbol"], f"WT9@{j}", mut25, flush=True)
    time.sleep(1)

with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["rank", "gene_symbol", "mutation", "core_9mer", "long_25mer", "length", "isoform"])
    w.writerows(out)
print(f"LONGPEP_DONE -> {OUT}")

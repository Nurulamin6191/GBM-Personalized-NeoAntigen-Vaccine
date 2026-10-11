#!/usr/bin/env python3
"""TESLA-style audit of Table 1 v2 using only checked-in data (no predictor env needed).
TESLA (Wells et al., Cell 2020): immunogenic pMHC pass PRESENTATION filters
(affinity <= 34 nM intervention; tumor abundance; stability) and RECOGNITION filters
(low agretopicity = WT binds much worse, or high foreignness).
What we can score offline: affinity<=34, TPM (abundance), dbSNP-novelty, NetMHCpan IC50<=34.
What we cannot (no MHCflurry/NetMHC MFC env, no WT contexts): agretopicity, foreignness,
stability -> marked PENDING, listed as the exact upgrade path.
Out: 06_Neoantigens/Table1_TESLA_audit.csv
"""
import csv

t1 = list(csv.DictReader(open("06_Neoantigens/Table1v2_expressed.csv")))
iedb = {}
with open("06_Neoantigens/IEDB_api/netmhcpan_ba_A0206.tsv") as f:
    rdr = csv.reader(f, delimiter="\t")
    next(rdr)
    for row in rdr:
        # columns: set, gene, peptide, <10 IEDB fields ending ic50, percentile_rank>
        iedb[row[2]] = float(row[-2])

rows = []
for r in t1:
    aff = float(r["affinity_nM_A0206"])
    tpm = float(r["max_TPM"])
    ic50 = iedb.get(r["peptide"])
    rows.append({
        "rank": r["rank"], "gene_symbol": r["gene_symbol"], "mutation": r["mutation"],
        "peptide": r["peptide"],
        "P_affinity_le34": aff <= 34,
        "P_abundance_TPM": tpm,
        "P_netmhcpan_le34": (ic50 <= 34) if ic50 is not None else "NA",
        "R_agretopicity": "PENDING (needs WT affinities)",
        "R_foreignness": "PENDING (needs proteome-wide scoring)",
        "stability": "PENDING (needs NetMHCstabpan)",
        "TESLA_presentation_pass": (aff <= 34) and (tpm >= 1) and
                                   (ic50 is not None and ic50 <= 34),
    })

with open("06_Neoantigens/Table1_TESLA_audit.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
n = sum(1 for r in rows if r["TESLA_presentation_pass"])
print(f"TESLA presentation pass: {n}/{len(rows)}")
for r in rows:
    print(r["rank"], r["gene_symbol"], r["peptide"], "affOK" if r["P_affinity_le34"] else "AFFAIL",
          "netMHC" if r["P_netmhcpan_le34"] is True else r["P_netmhcpan_le34"])

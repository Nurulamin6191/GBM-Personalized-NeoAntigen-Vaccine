#!/usr/bin/env python3
"""Script-reproducible NetMHCpan validation via the IEDB API (POST, TSV out).
Replaces reliance on point-and-click IEDB web tables.
Input : 06_Neoantigens/Table1v2_expressed.csv + SDF4/controls
Output: 06_Neoantigens/IEDB_api/netmhcpan_ba_A0206.tsv
Method codes: netmhcpan_ba (binding affinity IC50) and netmhcpan_el (elution ligand).
"""
import csv
import time
import urllib.parse
import urllib.request

API = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
ALLELE = "HLA-A*02:06"
OUT = "06_Neoantigens/IEDB_api/netmhcpan_ba_A0206.tsv"


def predict(peptide, method="netmhcpan_ba"):
    data = urllib.parse.urlencode({
        "method": method, "sequence_text": peptide,
        "allele": ALLELE, "length": len(peptide)}).encode()
    req = urllib.request.Request(API, data=data)
    for _ in range(3):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                lines = r.read().decode().strip().splitlines()
            if len(lines) >= 2:
                return lines[1]
        except Exception:
            time.sleep(5)
    return None


def main():
    peps = []
    for r in csv.DictReader(open("06_Neoantigens/Table1v2_expressed.csv")):
        peps.append(("Table1v2", r["peptide"], r["gene_symbol"]))
    for p in ["SAKEMQRWV", "EMQRWVMEK", "ALAAAVAAV", "LLFGYPVYV", "SYFPEITHI"]:
        peps.append(("control", p, "control"))
    seen, rows = set(), []
    for src, pep, gene in peps:
        if pep in seen:
            continue
        seen.add(pep)
        line = predict(pep)
        print(pep, "->", (line.split("\t")[8] if line else "FAIL"), flush=True)
        rows.append((src, gene, pep, line))
        time.sleep(1)
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["set", "gene", "peptide", "iedb_allele", "seq_num", "start",
                    "end", "length", "core", "icore", "ic50", "percentile_rank"])
        for src, gene, pep, line in rows:
            parts = line.split("\t") if line else [""] * 10
            w.writerow([src, gene, pep] + parts[:10])
    print(f"IEDB_API_DONE n={len(rows)} -> {OUT}")


if __name__ == "__main__":
    main()

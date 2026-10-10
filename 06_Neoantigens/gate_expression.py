#!/usr/bin/env python3
"""Gate Table 1 by expression: keep genes with TPM>=1 in either replicate.
Inputs: results/expr/SRR*.tab (StringTie gene abundance), 06_Neoantigens/Table1_top10_rare_binders.csv
        + full rare-binder list rebuild: join all 1083 via high_affinity table.
Out: 06_Neoantigens/Table1v2_expressed.csv
"""
import csv, glob, json, urllib.request


def symbol(eid):
    try:
        url = f"https://rest.ensembl.org/lookup/id/{eid}?content-type=application/json"
        d = json.load(urllib.request.urlopen(
            urllib.request.Request(url, headers={"Content-Type": "application/json"}), timeout=20))
        return d.get("display_name", eid)
    except Exception:
        return eid

tpm = {}
for f in glob.glob("results/expr/*.tab"):
    with open(f) as fh:
        for line in fh:
            if line.startswith("Gene ID"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            gene, cov, fpkm, tpmv = p[0], p[6], p[7], p[7]
            try:
                t = float(tpmv)
            except ValueError:
                continue
            tpm[gene] = max(tpm.get(gene, 0.0), t)
print(f"genes with TPM: {len(tpm)}")

hi = list(csv.DictReader(open("06_Neoantigens/high_affinity_candidates.csv")))
rare_genes = set()
with open("06_Neoantigens/Table1_top10_rare_binders.csv") as f:
    pass
# full rare set: genes from build_table1 logic
import sys
sys.path.insert(0, "06_Neoantigens")
rare = set()
dbsnp = list(csv.DictReader(open("04_Variants/missense_dbsnp.tsv"), delimiter="\t"))
# map via VCF ANN like build_table1
with open("/mnt/nvme0n1p5/AGY/Previous/Last/Galaxy_Pipelines/SRR_annotated.vcf") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.split("\t")
        key = (p[0].replace("chr", ""), p[1])
        hit = [r for r in dbsnp if (r["chrom"], r["pos"]) == key and r["verdict"] == "RARE_PRIVATE"]
        if not hit:
            continue
        for m in p[7].split(","):
            if "missense_variant" in m:
                parts = m.split("|")
                if len(parts) > 4 and parts[4]:
                    rare.add(parts[4])
rows = [r for r in hi if r["Gene_ID"] in rare and float(r["affinity"]) < 500]
rows.sort(key=lambda r: float(r["affinity"]))
kept, dropped = [], 0
for r in rows:
    # StringTie gene IDs are versionless ENSG; strip version
    gid = r["Gene_ID"].split(".")[0]
    t = tpm.get(gid, tpm.get(r["Gene_ID"], 0.0))
    if t >= 1.0:
        r["max_TPM"] = round(t, 2)
        kept.append(r)
    else:
        dropped += 1
print(f"rare binders: {len(rows)}, expressed TPM>=1: {len(kept)}, dropped: {dropped}")
with open("06_Neoantigens/Table1v2_expressed.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["rank", "gene_symbol", "ensembl_gene", "mutation", "peptide",
                "affinity_nM_A0206", "max_TPM"])
    for i, r in enumerate(kept[:10], 1):
        w.writerow([i, symbol(r["Gene_ID"]), r["Gene_ID"], r["Mutation"], r["Peptide"],
                    round(float(r["affinity"]), 2), r["max_TPM"]])
print("TABLE1V2_DONE")

#!/usr/bin/env python3
"""dbSNP germline filter via Ensembl VEP REST (threaded batch POST, 4 workers).
Input : 04_Variants/missense_dbsnp.tsv (chrom/pos/ref/alt; re-runnable, no VCF needed)
Output: same file refreshed with rsid + verdict columns.
Verdict: COMMON_GERMLINE_PROXY if a colocated dbSNP rsID exists, else RARE_PRIVATE.
(Caveat: rare germline variants are also novel; tumor-only design requires this proxy.
VEP REST returns no allele frequencies in this Ensembl version — verified with rs1799983.)
"""
import csv, json, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

INOUT = "04_Variants/missense_dbsnp.tsv"
CHUNK = 200
WORKERS = 4


def post(chunk):
    body = json.dumps({"variants": chunk}).encode()
    req = urllib.request.Request(
        "https://rest.ensembl.org/vep/homo_sapiens/region", data=body,
        headers={"Content-Type": "application/json", "Accept": "application/json"})
    for _ in range(3):
        try:
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.load(r)
        except Exception:
            time.sleep(5)
    return []


def verdict_for(v):
    for col in v.get("colocated_variants", []):
        if str(col.get("id", "")).startswith("rs"):
            return col["id"], "COMMON_GERMLINE_PROXY"
    return "NA", "RARE_PRIVATE"


def main():
    rows = list(csv.DictReader(open(INOUT), delimiter="\t"))
    print(f"records: {len(rows)}", flush=True)
    chunks = [[f"{r['chrom']} {r['pos']} {r['pos']} {r['ref']}/{r['alt']} +"
               for r in rows[i:i + CHUNK]]
              for i in range(0, len(rows), CHUNK)]
    out_rows = []
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for n, (chunk_rows, res) in enumerate(
                zip([rows[i:i + CHUNK] for i in range(0, len(rows), CHUNK)],
                    ex.map(post, chunks)), 1):
            got = res if isinstance(res, list) else []
            for j, r in enumerate(chunk_rows):
                if j < len(got):
                    rs, verdict = verdict_for(got[j])
                else:
                    rs, verdict = "NA", "UNRESOLVED"
                out_rows.append((r["chrom"], r["pos"], r["ref"], r["alt"], rs, verdict))
            print(f"chunk {n}/{len(chunks)}", flush=True)
    with open(INOUT, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom", "pos", "ref", "alt", "rsid", "verdict"])
        w.writerows(out_rows)
    from collections import Counter
    print(Counter(v for *_, v in out_rows), "DBSNP_DONE")


if __name__ == "__main__":
    main()

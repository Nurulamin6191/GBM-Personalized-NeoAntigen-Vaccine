#!/usr/bin/env python3
"""Map RARE_PRIVATE missense positions -> Ensembl gene IDs via VEP REST (threaded).
Out: 06_Neoantigens/rare_genes.tsv (ensembl_gene, n_rare_variants, gene_symbol)
Replaces the deleted-VCF parsing previously done inside gate_expression.py.
"""
import csv, json, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

MISS = "04_Variants/missense_dbsnp.tsv"
OUT = "06_Neoantigens/rare_genes.tsv"
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


def main():
    rare = [(r["chrom"], r["pos"], r["ref"], r["alt"])
            for r in csv.DictReader(open(MISS), delimiter="\t")
            if r["verdict"] == "RARE_PRIVATE"]
    print(f"rare positions: {len(rare)}", flush=True)
    chunks = [[f"{c} {p} {p} {r}/{a} +" for c, p, r, a in rare[i:i + CHUNK]]
              for i in range(0, len(rare), CHUNK)]
    genes = {}
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for n, res in enumerate(ex.map(post, chunks), 1):
            for v in res if isinstance(res, list) else []:
                best = None
                for tc in v.get("transcript_consequences", []):
                    if "missense_variant" in tc.get("consequence_terms", []):
                        if tc.get("canonical") == 1:
                            best = tc
                            break
                        best = best or tc
                if best and best.get("gene_id"):
                    g = best["gene_id"]
                    genes[g] = (genes.get(g, (0, best.get("gene_symbol", g)))[0] + 1,
                                best.get("gene_symbol", g))
            print(f"chunk {n}/{len(chunks)} genes={len(genes)}", flush=True)
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["ensembl_gene", "n_rare_variants", "gene_symbol"])
        for g, (cnt, sym) in sorted(genes.items()):
            w.writerow([g, cnt, sym])
    print(f"RARE_GENES_DONE n={len(genes)} -> {OUT}")


if __name__ == "__main__":
    main()

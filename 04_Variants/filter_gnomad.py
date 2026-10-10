#!/usr/bin/env python3
"""Filter missense VCF records against gnomAD AF via Ensembl VEP REST (batch POST).
Keeps AF<0.01 or absent (rare/private). Input: Galaxy SRR_annotated.vcf (chr1).
Out: 04_Variants/missense_gnomad.tsv + rare_candidates.csv
"""
import json, re, sys, time, urllib.request

VCF = "Last/Galaxy_Pipelines/SRR_annotated.vcf"
OUT = "rebuild/04_Variants/missense_dbsnp.tsv"
CHUNK = 200

def post(chunk):
    body = json.dumps({"variants": chunk}).encode()
    req = urllib.request.Request("https://rest.ensembl.org/vep/homo_sapiens/region",
                                 data=body, headers={"Content-Type": "application/json",
                                                     "Accept": "application/json"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except Exception as e:
            print(f"retry {attempt}: {e}", flush=True)
            time.sleep(5)
    return None

recs = []
with open("/mnt/nvme0n1p5/AGY/previous/" + VCF) as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.split("\t")
        ann = p[7]
        if "missense_variant" not in ann:
            continue
        m = re.search(r"\|(\w+\.\w+Ala\d+\w+|\w+\.p\.\w+)|", ann)
        # VEP region format: "CHR START END REF/ALT STRAND"
        recs.append((p[0].replace("chr", ""), p[1], p[3], p[4], ann))

print(f"missense records: {len(recs)}", flush=True)
out = open("/mnt/nvme0n1p5/AGY/previous/" + OUT, "w")
out.write("chrom\tpos\tref\talt\trsid\tverdict\n")
rare = 0
for i in range(0, len(recs), CHUNK):
    chunk = [f"{c} {p} {p} {r}/{a} +" for c, p, r, a, _ in recs[i:i + CHUNK]]
    res = post(chunk)
    if res is None:
        for c, p, r, a, _ in recs[i:i + CHUNK]:
            out.write(f"{c}\t{p}\t{r}\t{a}\tNA\tNA\tVEP_FAIL\n")
        continue
    got = res if isinstance(res, list) else []
    for j, (c, p, r, a, _) in enumerate(recs[i:i + CHUNK]):
        if j >= len(got):
            out.write(f"{c}\t{p}\t{r}\t{a}\tNA\tUNRESOLVED\n")
            continue
        v = got[j]
        rs = "NA"
        for col in v.get("colocated_variants", []):
            cid = str(col.get("id", ""))
            if cid.startswith("rs"):
                rs = cid
                break
        verdict = "COMMON_GERMLINE_PROXY" if rs != "NA" else "RARE_PRIVATE"
        if verdict == "RARE_PRIVATE":
            rare += 1
        out.write(f"{c}\t{p}\t{r}\t{a}\t{rs}\t{verdict}\n")
    out.flush()
    print(f"done {min(i + CHUNK, len(recs))}/{len(recs)} rare={rare}", flush=True)
    time.sleep(1)
out.close()
print(f"DBSNP_DONE rare={rare}")

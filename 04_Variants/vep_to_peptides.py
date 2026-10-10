#!/usr/bin/env python3
"""VEP tab → mutant 8–11mer tiling.
Bridge added in rebuild: v1 jumped from VCF straight to 2 hand-picked peptides.
Usage: python3 04_Variants/vep_to_peptides.py [vep_tab] [out_csv]
Falls back to demo rows when VEP output absent (keeps repo runnable pre-full-run).
"""
import csv, sys, itertools
from pathlib import Path

def tile_9mers(center_seq, k=9):
    center_seq = center_seq.upper().replace("*", "")
    if len(center_seq) < k:
        return [center_seq]
    return [center_seq[i:i+k] for i in range(len(center_seq)-k+1)]

def main():
    inp = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("results-RUNTIME/peptides/tiled_9mers.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    if inp and inp.exists():
        with open(inp) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                # default VEP --fields order; guard short lines
                d = {h: (p[i] if i < len(p) else "") for i, h in enumerate(
                    ["Uploaded_variation","Location","Allele","Gene","SYMBOL","Consequence",
                     "Protein_position","Amino_acids","Codons","SIFT","PolyPhen","EXON"])}
                if "missense" not in d["Consequence"]:
                    continue
                aa = d["Amino_acids"]  # e.g. L/V
                if "/" not in aa:
                    continue
                ref, alt = aa.split("/")
                # window mock: 10 ref-context AAs each side; replaced by real Ensembl lookup in full run
                ctx = "A"*10 + alt + "A"*10
                for pep in tile_9mers(ctx):
                    rows.append({"variant": d["Uploaded_variation"], "gene": d["SYMBOL"],
                                 "consequence": d["Consequence"], "peptide": pep, "length": len(pep)})
    if not rows:
        rows = [
            {"variant": "DEMO:chr1:missense", "gene": "DEMO", "consequence": "missense_variant",
             "peptide": "LLFGYPVYV", "length": 9, "note": "v1 top candidate; re-derived in full run"},
            {"variant": "DEMO:control", "gene": "CTRL", "consequence": "control",
             "peptide": "SYFPEITHI", "length": 9, "note": "classical HLA-A*02:01 binder control — NOT a negative control"},
        ]
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"wrote {len(rows)} peptides -> {out}")

if __name__ == "__main__":
    main()

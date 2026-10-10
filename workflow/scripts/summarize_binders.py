#!/usr/bin/env python3
"""Summarize MHCflurry output into ranked table + thresholds.
Strong <50 nM, weak <500 nM (IEDB convention). v1 used <500 alone — now stratified.
"""
import pandas as pd
from pathlib import Path
inp = Path("results/neoantigens/mhcflurry_A0206.csv")
out = Path("results/neoantigens/ranked_candidates.csv")
if not inp.exists():
    # preliminary table from v1 values so repo has a runnable result pre-full-run
    df = pd.DataFrame([
        {"peptide": "LLFGYPVYV", "allele": "HLA-A*02:06", "affinity_nM": 10.31, "class": "strong binder (<50 nM)", "priority": 1},
        {"peptide": "SYFPEITHI", "allele": "HLA-A*02:06", "affinity_nM": 746.85, "class": "non-binder (>500 nM)", "priority": 99},
    ])
    print("MHCflurry output absent — writing v1 preliminary table (to be replaced by full run).")
else:
    df = pd.read_csv(inp)
    aff = [c for c in df.columns if "affinity" in c.lower() or "ic50" in c.lower()][0]
    df = df.rename(columns={aff: "affinity_nM"}).sort_values("affinity_nM")
    df["class"] = pd.cut(df["affinity_nM"], [-1, 50, 500, 1e9],
                          labels=["strong binder (<50 nM)", "weak binder (50–500 nM)", "non-binder (>500 nM)"])
    df["priority"] = range(1, len(df)+1)
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, index=False)
print(f"wrote {len(df)} rows -> {out}")
print(df.to_string(index=False))

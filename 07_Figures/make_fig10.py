#!/usr/bin/env python3
"""Fig10: Table 1 v2 — affinity vs expression.
Run from repo root: python3 07_Figures/make_fig10.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

FIG = Path("07_Figures")
FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 11, "figure.dpi": 300,
                     "axes.titlesize": 12, "axes.titleweight": "bold"})

rows = pd.read_csv("06_Neoantigens/Table1v2_expressed.csv")
fig, ax = plt.subplots(figsize=(8, 5))
ax.scatter(rows["affinity_nM_A0206"], rows["max_TPM"], c="#16a34a", s=70)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("MHCflurry affinity, nM vs HLA-A*02:06 (log)")
ax.set_ylabel("max TPM across replicates (log)")
ax.set_title("Table 1 candidates are all expressed (TPM 2.6-95)")
for _, r in rows.iterrows():
    ax.text(r["affinity_nM_A0206"] * 1.04, r["max_TPM"], r["gene_symbol"], fontsize=10)
fig.tight_layout()
fig.savefig(FIG / "fig10_table1v2_expression.png", bbox_inches="tight")
print("wrote fig10")

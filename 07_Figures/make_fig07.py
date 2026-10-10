#!/usr/bin/env python3
"""Fig7: SDF4 tiling 9mers across the full diploid HLA profile (MHCflurry 2.1.2).
Run from repo root: python3 07_Figures/make_fig07.py
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

a0206 = pd.read_csv("06_Neoantigens/sdf4_validation_out.csv")
other = pd.read_csv("06_Neoantigens/sdf4_validation_other_alleles.csv")
df = pd.concat([a0206[["allele", "peptide", "mhcflurry_affinity"]],
                other[["allele", "peptide", "mhcflurry_affinity"]]], ignore_index=True)
sdf4 = ["SAKEMQRWV", "AKEMQRWVM", "KEMQRWVME", "EMQRWVMEK", "MQRWVMEKT",
        "QRWVMEKTA", "RWVMEKTAE", "WVMEKTAEH", "VMEKTAEHF"]
piv = df[df.peptide.isin(sdf4)].pivot(index="peptide", columns="allele",
                                      values="mhcflurry_affinity")
piv = piv.reindex(sdf4)

fig, ax = plt.subplots(figsize=(11, 5.2))
x = range(len(sdf4))
w = 0.15
alleles = ["HLA-A*02:06", "HLA-A*30:01", "HLA-B*08:01", "HLA-B*53:01", "HLA-C*07:01"]
for i, al in enumerate(alleles):
    vals = [min(piv.loc[p, al], 35000) for p in sdf4]
    ax.bar([xi + (i - 2) * w for xi in x], vals, width=w, label=al)
ax.set_xticks(list(x))
ax.set_xticklabels(sdf4, rotation=25, ha="right", fontsize=10)
ax.axhline(50, color="#16a34a", ls="--", lw=1.2)
ax.text(8.6, 65, "strong < 50 nM", fontsize=10, color="#16a34a")
ax.axhline(500, color="#d97706", ls="--", lw=1.2)
ax.text(8.6, 650, "weak < 500 nM", fontsize=10, color="#d97706")
ax.set_yscale("log")
ax.set_ylabel("affinity, nM (log scale, capped at 35k)")
ax.set_title("SDF4 tiling fails on every patient allele (best: 192 nM, weak)")
ax.legend(fontsize=10, ncol=3, loc="upper left")
fig.tight_layout()
fig.savefig(FIG / "fig07_sdf4_full_allele.png", bbox_inches="tight")
print("wrote fig07")

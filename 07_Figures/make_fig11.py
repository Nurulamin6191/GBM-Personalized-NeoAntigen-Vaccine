#!/usr/bin/env python3
"""Fig11: MHCflurry vs NetMHCpan (IEDB API) concordance on Table 1 + controls.
Run from repo root: python3 07_Figures/make_fig11.py
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

ied = pd.read_csv("06_Neoantigens/IEDB_api/netmhcpan_ba_A0206.tsv", sep="\t")
ied["ic50"] = ied["ic50"].astype(float)
mbf = pd.read_csv("06_Neoantigens/sdf4_validation_out.csv")
mbf = dict(zip(mbf["peptide"], mbf["mhcflurry_affinity"]))
t1 = pd.read_csv("06_Neoantigens/Table1v2_expressed.csv")
t1 = dict(zip(t1["peptide"], t1["affinity_nM_A0206"]))

xs, ys, labels, colors = [], [], [], []
for _, r in ied.iterrows():
    pep = r["peptide"]
    m = mbf.get(pep, t1.get(pep))
    if m is None or pd.isna(m):
        continue
    xs.append(float(m))
    ys.append(float(r["ic50"]))
    labels.append(pep if r["set"] == "control" else r["gene"])
    colors.append("#d97706" if r["set"] == "control" else "#2563eb")

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(xs, ys, c=colors, s=70)
ax.set_xscale("log")
ax.set_yscale("log")
ax.plot([1, 40000], [1, 40000], color="#6b7280", ls="--", lw=1.2, label="perfect agreement")
ax.axhline(50, color="#16a34a", ls=":", lw=1.2)
ax.axvline(50, color="#16a34a", ls=":", lw=1.2)
ax.text(60, 1.6, "strong < 50 nM", fontsize=10, color="#16a34a")
ax.set_xlabel("MHCflurry affinity, nM (log)")
ax.set_ylabel("NetMHCpan BA IC50, nM (log)")
ax.set_title("Two predictors agree: Table 1 strong, SDF4 tiling not")
for x, y, lab in zip(xs, ys, labels):
    ax.text(x * 1.08, y, lab, fontsize=9)
ax.legend(fontsize=10)
fig.tight_layout()
fig.savefig(FIG / "fig11_predictor_concordance.png", bbox_inches="tight")
print("wrote fig11")

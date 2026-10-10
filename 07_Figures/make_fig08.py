#!/usr/bin/env python3
"""Fig8: replicate consistency + independent NetMHCpan (IEDB) validation.
Run from repo root: python3 07_Figures/make_fig08.py
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
BLUE, GREEN = "#2563eb", "#16a34a"

fig, ax = plt.subplots(1, 2, figsize=(12, 5))

# (a) replicate variant + peptide counts
cats = ["FreeBayes\nvariants", "Strong-binder\npeptides"]
v899 = [45456, 453]
v902 = [51435, 430]
x = range(len(cats))
w = 0.35
ax[0].bar([i - w / 2 for i in x], v899, width=w, label="SRR30677899", color=BLUE)
ax[0].bar([i + w / 2 for i in x], v902, width=w, label="SRR30677902", color="#d97706")
ax[0].set_xticks(list(x))
ax[0].set_xticklabels(cats, fontsize=11)
ax[0].set_title("(a) Replicates agree (same BioSample)")
ax[0].legend(fontsize=10)
for i, (a, b) in enumerate(zip(v899, v902)):
    ax[0].text(i - w / 2, a * 1.04, f"{a:,}", ha="center", fontsize=10)
    ax[0].text(i + w / 2, b * 1.04, f"{b:,}", ha="center", fontsize=10)
ax[0].text(1, 30000, "155 peptides shared", ha="center", fontsize=11, weight="bold",
           bbox=dict(fc="#fef3c7", alpha=0.9))

# (b) IEDB NetMHCpan BA IC50 (independent validation set)
ied = []
for s in ["SRR30677899", "SRR30677902"]:
    d = pd.read_csv(f"06_Neoantigens/IEDB/IEDB_peptide_table_{s}.csv")
    ied.append(d)
ied = pd.concat(ied, ignore_index=True)
ax[1].bar(range(len(ied)), ied["netmhcpan_ba ic50"],
          color=[BLUE] * 9 + ["#d97706"] * 10)
ax[1].axhline(50, color=GREEN, ls="--", lw=1.2)
ax[1].text(0, 58, "strong binder cutoff 50 nM", fontsize=10, color=GREEN)
ax[1].set_xticks(range(len(ied)))
ax[1].set_xticklabels(ied["peptide"], rotation=45, ha="right", fontsize=9)
ax[1].set_ylabel("NetMHCpan BA IC50, nM")
ax[1].set_title("(b) Independent NetMHCpan check: all strong")
fig.suptitle("Replicates agree; an independent predictor confirms strong binding",
             fontsize=13, weight="bold")
fig.tight_layout()
fig.savefig(FIG / "fig08_replicates_IEDB.png", bbox_inches="tight")
print("wrote fig08")

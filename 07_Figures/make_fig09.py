#!/usr/bin/env python3
"""Fig9: germline burden + Table 1 top-10 rare binders.
Run from repo root: python3 07_Figures/make_fig09.py
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

fig, ax = plt.subplots(1, 2, figsize=(12, 5))

ax[0].pie([2126, 1815, 324],
          labels=["Common germline\n2,126 (49.9%)", "Rare / private\n1,815 (42.6%)",
                  "Unresolved\n324 (7.6%)"],
          colors=["#dc2626", "#16a34a", "#6b7280"],
          textprops={"fontsize": 11})
ax[0].set_title("(a) Half the missense pool is germline")

top = pd.read_csv("06_Neoantigens/Table1_top10_rare_binders.csv")
ax[1].barh(top["peptide"][::-1], top["affinity_nM_A0206"][::-1], color="#16a34a")
ax[1].set_xlabel("MHCflurry affinity, nM vs HLA-A*02:06")
ax[1].set_title("(b) Top-10 rare binders (all under 13 nM)")
for i, (pep, v) in enumerate(zip(top["peptide"][::-1], top["affinity_nM_A0206"][::-1])):
    ax[1].text(v + 0.1, i, f"{v:.1f}", va="center", fontsize=9)
fig.suptitle("Germline filtering halves the pool; rare top-10 all bind strongly",
             fontsize=13, weight="bold")
fig.tight_layout()
fig.savefig(FIG / "fig09_germline_table1.png", bbox_inches="tight")
print("wrote fig09")

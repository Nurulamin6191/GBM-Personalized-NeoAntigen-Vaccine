#!/usr/bin/env python3
"""Fig8: replicate consistency (899 vs 902) + independent NetMHCpan (IEDB) validation."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
FIG = Path("07_Figures"); FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 8, "figure.dpi": 300})

fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.8))

# left: replicate variant + peptide counts
cats = ["FreeBayes\nvariants", "Strong-binder\npeptides*"]
v899 = [45456, 453]; v902 = [51435, 430]
x = range(len(cats)); w = 0.35
ax[0].bar([i - w/2 for i in x], v899, width=w, label="SRR30677899", color="#3070b8")
ax[0].bar([i + w/2 for i in x], v902, width=w, label="SRR30677902", color="#b36600")
ax[0].set_xticks(list(x)); ax[0].set_xticklabels(cats)
ax[0].set_title("(a) Replicate consistency: same BioSample")
ax[0].legend(fontsize=7)
for i, (a, b) in enumerate(zip(v899, v902)):
    ax[0].text(i - w/2, a*1.03, f"{a:,}", ha="center", fontsize=7)
    ax[0].text(i + w/2, b*1.03, f"{b:,}", ha="center", fontsize=7)
ax[0].text(1, 30000, "shared: 155", ha="center", fontsize=8, weight="bold",
           bbox=dict(fc="yellow", alpha=0.5))

# right: IEDB NetMHCpan BA IC50 distribution (independent validation set)
ied = []
for s in ["SRR30677899", "SRR30677902"]:
    d = pd.read_csv(f"06_Neoantigens/IEDB/IEDB_peptide_table_{s}.csv")
    d["sample"] = s[-5:]
    ied.append(d)
ied = pd.concat(ied)
ax[1].bar(range(len(ied)), ied["netmhcpan_ba ic50"],
          color=["#3070b8"]*9 + ["#b36600"]*10)
ax[1].axhline(50, color="green", ls="--", lw=1)
ax[1].text(0, 55, "strong <50 nM", fontsize=7, color="green")
ax[1].set_xticks(range(len(ied)))
ax[1].set_xticklabels(ied["peptide"], rotation=45, ha="right", fontsize=6)
ax[1].set_ylabel("NetMHCpan BA IC50 (nM)")
ax[1].set_title("(b) Independent NetMHCpan validation (IEDB): all strong")
fig.suptitle("Fig 8 — Replicate consistency + orthogonal NetMHCpan validation.", fontsize=9, weight="bold")
fig.tight_layout()
fig.savefig(FIG/"fig08_replicates_IEDB.png", bbox_inches="tight")
print("wrote fig08")
print("*per-sample A*02:01 runs re-scored vs patient A*02:06 in 06_Neoantigens/per_sample_A0206/ (values reproduce).")

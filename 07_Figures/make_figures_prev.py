#!/usr/bin/env python3
"""Fig4-6: top binders, consequence spectrum, adapter burden.
Run from repo root: python3 07_Figures/make_figures_prev.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import re
from pathlib import Path

FIG = Path("07_Figures")
FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 11, "figure.dpi": 300,
                     "axes.titlesize": 12, "axes.titleweight": "bold"})
BLUE, GREEN, AMBER, RED = "#2563eb", "#16a34a", "#d97706", "#dc2626"

# ---- Fig4: top-20 batch binders ----
hi = pd.read_csv("06_Neoantigens/high_affinity_candidates.csv")
hi = hi.sort_values("affinity").head(20).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(9, 7.5))
ax.barh(hi["Peptide"], hi["affinity"], color=GREEN)
ax.invert_yaxis()
ax.set_xlabel("MHCflurry affinity, nM vs HLA-A*02:06 (lower = stronger)")
ax.set_title("Top-20 batch-screen binders (all under 13 nM)")
ax.axvline(50, color=GREEN, ls="--", lw=1.2)
ax.text(51, 0.4, "strong binder cutoff 50 nM", fontsize=10, color=GREEN)
for i, v in enumerate(hi["affinity"]):
    ax.text(v + 0.12, i, f"{v:.1f}", va="center", fontsize=9)
fig.tight_layout()
fig.savefig(FIG / "fig04_top20_affinity.png", bbox_inches="tight")
print("wrote fig04", len(hi))

# ---- Fig5: consequence spectrum ----
vc = pd.read_csv("04_Variants/variant_counts.tsv", sep="\t",
                 header=None, names=["kind", "count"])
vc = vc[vc.kind != "total_variants"].head(9)
fig, ax = plt.subplots(figsize=(9, 4.6))
ax.bar(range(len(vc)), vc["count"], color=BLUE)
ax.set_xticks(range(len(vc)))
ax.set_xticklabels(vc["kind"], rotation=25, ha="right", fontsize=10)
ax.set_ylabel("variant count (n = 165,818)")
ax.set_title("FreeBayes consequence spectrum: 4,224 missense feed the screen")
for x, v in enumerate(vc["count"]):
    ax.text(x, v * 1.02, f"{v:,}", ha="center", fontsize=10)
fig.tight_layout()
fig.savefig(FIG / "fig05_consequences.png", bbox_inches="tight")
print("wrote fig05")

# ---- Fig6: adapter burden (local Trim Galore paired reports) ----
rows = []
for s in ["SRR30677899", "SRR30677902"]:
    for r in ["1", "2"]:
        p = Path(f"02_QC/{s}_{r}.fastq.gz_trimming_report.txt")
        if not p.exists():
            continue
        t = p.read_text()
        m = re.search(r"Reads with adapters:\s+([\d,]+)", t)
        rows.append({"label": f"{s[-5:]} R{r}",
                     "adapters": int(m.group(1).replace(",", "")) if m else -1})
df = pd.DataFrame(rows)
fig, ax = plt.subplots(figsize=(7, 4.4))
ax.bar(df["label"], df["adapters"], color=[AMBER, BLUE, AMBER, BLUE])
ax.set_ylabel("reads with adapters trimmed")
ax.set_title("Adapter burden: R2 higher in both replicates")
for i, v in enumerate(df["adapters"]):
    ax.text(i, v * 1.02, f"{v:,}", ha="center", fontsize=10)
fig.tight_layout()
fig.savefig(FIG / "fig06_adapter_burden.png", bbox_inches="tight")
print("wrote fig06")

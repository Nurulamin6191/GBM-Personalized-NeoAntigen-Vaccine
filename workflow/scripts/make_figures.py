#!/usr/bin/env python3
"""Publication figures (preliminary; regenerated with real counts after full run)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
FIG = Path("rebuild/figures")
FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 9, "figure.dpi": 300})

# Fig1: pipeline with missing bridge highlighted
fig, ax = plt.subplots(figsize=(10, 3.2))
ax.axis("off")
steps = ["SRA\n15–18M reads", "QC\nQ20 + MultiQC", "Align\nHISAT2 full\n+ STAR-chr1 pilot",
         "Variants\nGATK RNA-seq", "VEP\nmissense", "TILE\n8–11mers\n★ NEW BRIDGE",
         "HLA\nOptiType\n6 alleles", "MHCflurry\n+ NetMHCpan\n50/500 nM"]
for i, s in enumerate(steps):
    fc = "#ffe08a" if "TILE" in s else ("#cfe8cf" if i < 5 else "#cfe0f1")
    ec = "red" if "TILE" in s else "black"
    lw = 2.2 if "TILE" in s else 1.0
    ax.add_patch(mpatches.FancyBboxPatch((i*1.22+0.05, 0.35), 1.12, 1.1, boxstyle="round,pad=0.02",
                                         fc=fc, ec=ec, lw=lw))
    ax.text(i*1.22+0.61, 0.9, s, ha="center", va="center", fontsize=7.5,
            weight="bold" if "TILE" in s else "normal")
    if i < len(steps)-1:
        ax.annotate("", xy=(i*1.22+1.17, 0.9), xytext=(i*1.22+1.27, 0.9),
                    arrowprops=dict(arrowstyle="->", lw=1.2))
ax.text(5*1.22+0.61, 0.28, "v1 skipped this → tested 2 hand-picked peptides",
        ha="center", fontsize=7, color="red", style="italic")
ax.set_xlim(0, len(steps)*1.22)
ax.set_ylim(0, 1.7)
ax.set_title("Fig 1 — End-to-end pipeline. Yellow = bridge added in rebuild (VEP→peptide tiling).", fontsize=9, weight="bold")
fig.tight_layout()
fig.savefig(FIG/"fig01_pipeline.png", bbox_inches="tight")
print("wrote fig01")

# Fig2: affinity barplot (v1 values, IEDB thresholds)
fig, ax = plt.subplots(figsize=(5.2, 3.4))
peps = ["LLFGYPVYV\n(top candidate)", "SYFPEITHI\n(binder control*)"]
vals = [10.31, 746.85]
cols = ["#1a7f37", "#808080"]
bars = ax.bar(peps, vals, color=cols)
ax.axhspan(0, 50, color="green", alpha=0.12, label="strong <50 nM")
ax.axhspan(50, 500, color="orange", alpha=0.12, label="weak 50–500 nM")
ax.axhspan(500, 1200, color="red", alpha=0.08, label="non-binder >500 nM")
ax.set_yscale("log")
ax.set_ylabel("MHCflurry affinity (nM, log scale)")
ax.set_title("Fig 2 — Preliminary binding (HLA-A*02:06, v1 values).\n*SYFPEITHI is a known A*02:01 binder, not a negative control.", fontsize=8)
for b, v in zip(bars, vals):
    ax.text(b.get_x()+b.get_width()/2, v*1.25, f"{v:.2f} nM", ha="center", fontsize=8, weight="bold")
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(FIG/"fig02_affinity.png", bbox_inches="tight")
print("wrote fig02")

# Fig3: chr1 pilot mapping explained
fig, ax = plt.subplots(figsize=(5.2, 3.2))
ax.barh(["chr1 pilot\nuniquely mapped", "other chromosomes\n(unmapped in pilot)"],
        [13, 87], color=["#3070b8", "#d3d3d3"])
ax.set_xlabel("% of RNA-seq reads")
ax.set_title("Fig 3 — Why ~13% is expected for chr1-only alignment\n(chr1 ≈ 8% of genome + expression bias), not 'Excellent'.", fontsize=8)
for i, v in enumerate([13, 87]):
    ax.text(v+1, i, f"{v}%", va="center", fontsize=9, weight="bold")
fig.tight_layout()
fig.savefig(FIG/"fig03_mapping.png", bbox_inches="tight")
print("wrote fig03")
print("FIGURES DONE")

#!/usr/bin/env python3
"""Fig1-3: pipeline overview, controls vs SDF4 tiling, mapping rates.
Run from repo root: python3 07_Figures/make_figures.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

FIG = Path("07_Figures")
FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 11, "figure.dpi": 300,
                     "axes.titlesize": 12, "axes.titleweight": "bold"})
BLUE, GREEN, AMBER, RED, GRAY = "#2563eb", "#16a34a", "#d97706", "#dc2626", "#6b7280"

# ---- Fig1: pipeline ----
fig, ax = plt.subplots(figsize=(13, 4.2))
ax.axis("off")
steps = ["SRA\n15-18M reads", "QC\nQ20 + MultiQC", "Align\nHISAT2 full\n+ STAR chr1",
         "Variants\nGATK + FreeBayes", "Annotate\nmissense", "TILE\n8-11mers",
         "HLA\nOptiType", "MHCflurry\n+ NetMHCpan"]
for i, s in enumerate(steps):
    fc = "#fef3c7" if "TILE" in s else ("#dcfce7" if i < 5 else "#dbeafe")
    ec = RED if "TILE" in s else "black"
    lw = 2.4 if "TILE" in s else 1.2
    ax.add_patch(mpatches.FancyBboxPatch((i * 1.6 + 0.05, 0.5), 1.45, 1.5,
                                         boxstyle="round,pad=0.02", fc=fc, ec=ec, lw=lw))
    ax.text(i * 1.6 + 0.775, 1.25, s, ha="center", va="center", fontsize=10,
            weight="bold" if "TILE" in s else "normal")
    if i < len(steps) - 1:
        ax.annotate("", xy=(i * 1.6 + 1.5, 1.25), xytext=(i * 1.6 + 1.65, 1.25),
                    arrowprops=dict(arrowstyle="->", lw=1.4))
ax.text(5 * 1.6 + 0.775, 0.32, "bridge v1 skipped (it tested 2 hand-picked peptides)",
        ha="center", fontsize=10, color=RED, style="italic")
ax.set_xlim(0, len(steps) * 1.6)
ax.set_ylim(0, 2.2)
ax.set_title("End-to-end neoantigen pipeline (yellow = added tiling bridge)")
fig.tight_layout()
fig.savefig(FIG / "fig01_pipeline.png", bbox_inches="tight")
print("wrote fig01")

# ---- Fig2: controls vs SDF4 tiling (HLA-A*02:06, MHCflurry 2.1.2) ----
fig, ax = plt.subplots(figsize=(8, 5))
labels = ["LLFGYPVYV\ncontrol binder", "ALAAAVAAV\nfocus (unattributed)",
          "SAKEMQRWV\nbest SDF4 9mer"]
vals = [10.31, 16.94, 8426.37]
cols = [GREEN, AMBER, RED]
bars = ax.bar(labels, vals, color=cols)
ax.axhspan(0.1, 50, color="green", alpha=0.10)
ax.axhspan(50, 500, color="orange", alpha=0.10)
ax.text(2.4, 12, "strong < 50 nM", fontsize=10, color=GREEN)
ax.text(2.4, 120, "weak 50-500 nM", fontsize=10, color="#b45309")
ax.text(2.4, 1500, "non-binder > 500 nM", fontsize=10, color=RED)
ax.set_yscale("log")
ax.set_ylabel("MHCflurry affinity, nM (log scale)")
ax.set_title("Positive control binds; SDF4 tiling does not (HLA-A*02:06)")
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width() / 2, v * 1.3, f"{v:,.2f} nM",
            ha="center", fontsize=11, weight="bold")
fig.tight_layout()
fig.savefig(FIG / "fig02_affinity.png", bbox_inches="tight")
print("wrote fig02")

# ---- Fig3: mapping rates ----
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
runs = ["Galaxy", "Local 899", "Local 902"]
uniq = [13.90, 13.91, 13.83]
ax[0].bar(runs, uniq, color=BLUE)
ax[0].set_ylabel("uniquely mapped reads (%)")
ax[0].set_title("(a) chr1 pilot: ~14% three times")
for i, v in enumerate(uniq):
    ax[0].text(i, v + 0.3, f"{v:.2f}%", ha="center", fontsize=11, weight="bold")
full = [44.72, 44.82]
ax[1].bar(["899 full", "902 full"], full, color=GREEN)
ax[1].set_ylabel("overall alignment rate (%)")
ax[1].set_title("(b) Full genome: ~45% (only 1.4% concordantly paired)")
for i, v in enumerate(full):
    ax[1].text(i, v + 0.8, f"{v:.2f}%", ha="center", fontsize=11, weight="bold")
fig.suptitle("chr1-only mapping is expected near 14%, not a quality failure",
             fontsize=12, weight="bold")
fig.tight_layout()
fig.savefig(FIG / "fig03_mapping.png", bbox_inches="tight")
print("wrote fig03")
print("FIGURES DONE")

#!/usr/bin/env python3
"""Fig4-6 from curated previous-work tables (no large downloads needed)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
FIG = Path("07_Figures"); FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 9, "figure.dpi": 300})

# Fig4: top-20 affinities + SDF4 focused marker
hi = pd.read_csv("06_Neoantigens/high_affinity_candidates.csv")
hi = hi.sort_values("affinity").head(20).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(8.5, 4.2))
ax.barh(hi["Peptide"], hi["affinity"], color="#1a7f37")
ax.invert_yaxis()
ax.set_xlabel("MHCflurry affinity (nM) vs HLA-A*02:06 — lower = stronger")
ax.set_title("Fig 4 — Top-20 batch-screen binders (all <13 nM, strong). SDF4 focus peptide ALAAAVAAV = 16.94 nM\n(not top-affinity; biologically prioritized, provenance unconfirmed — see VALIDATION).", fontsize=8)
ax.axvline(50, color="green", ls="--", lw=1); ax.text(51, 0.5, "strong <50", fontsize=7, color="green")
ax.axvline(16.94, color="red", ls=":", lw=1.5); ax.text(17.5, 10, "SDF4 focus 16.94", fontsize=7, color="red")
for i, v in enumerate(hi["affinity"]):
    ax.text(v+0.15, i, f"{v:.1f}", va="center", fontsize=6.5)
fig.tight_layout(); fig.savefig(FIG/"fig04_top20_affinity.png", bbox_inches="tight")
print("wrote fig04", len(hi))

# Fig5: variant consequence spectrum (165,818 FreeBayes raw)
vc = pd.read_csv("04_Variants/variant_counts.tsv", sep="\t", header=None, names=["k","v"])
vc = vc[vc.k != "total_variants"].head(9)
fig, ax = plt.subplots(figsize=(7.5, 3.6))
ax.bar(vc["k"], vc["v"], color="#3070b8")
ax.set_xticklabels(vc["k"], rotation=28, ha="right", fontsize=7)
ax.set_ylabel("count (n=165,818)")
ax.set_title("Fig 5 — FreeBayes chr1 raw consequence spectrum. Missense = 4,224 (~2.5%) → candidate pool.", fontsize=9)
for x, v in zip(range(len(vc)), vc["v"]):
    ax.text(x, v*1.02, f"{v:,}", ha="center", fontsize=7)
fig.tight_layout(); fig.savefig(FIG/"fig05_consequences.png", bbox_inches="tight")
print("wrote fig05")

# Fig6: QC adapter burden (Galaxy Trim Galore reports, both samples)
import re
rows = []
for s in ["SRR30677899", "SRR30677902"]:
    for r in ["1", "2"]:
        p = Path(f"02_QC/{s}/{s}_{r}.fastq.gz_trimming_report.txt")
        if p.exists():
            t = p.read_text()
            m = re.search(r"Reads with adapters:\s+([\d,]+)", t)
            rows.append({"label": f"{s[-5:]} R{r}",
                         "adapters": int(m.group(1).replace(",", "")) if m else -1})
df = pd.DataFrame(rows)
fig, ax = plt.subplots(figsize=(5.5, 3.4))
ax.bar(df["label"], df["adapters"], color=["#b36600","#3070b8","#b36600","#3070b8"])
ax.set_ylabel("reads with adapters (Galaxy Trim Galore, single-end mode)")
ax.set_title("Fig 6 — Adapter burden per read (Galaxy reports). R2 consistently higher.", fontsize=9)
for i, v in enumerate(df["adapters"]):
    ax.text(i, v*1.02, f"{v:,}", ha="center", fontsize=7)
fig.tight_layout(); fig.savefig(FIG/"fig06_adapter_burden.png", bbox_inches="tight")
print("wrote fig06")

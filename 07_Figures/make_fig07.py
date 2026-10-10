#!/usr/bin/env python3
"""Fig7: SDF4 tiling 9mers + focus/control across full diploid HLA profile (fresh MHCflurry 2.1.2)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
FIG = Path("07_Figures"); FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 8, "figure.dpi": 300})

a0206 = pd.read_csv("06_Neoantigens/sdf4_validation_out.csv")
other = pd.read_csv("06_Neoantigens/sdf4_validation_other_alleles.csv")
df = pd.concat([a0206[["allele", "peptide", "mhcflurry_affinity"]],
                other[["allele", "peptide", "mhcflurry_affinity"]]], ignore_index=True)
sdf4 = ["SAKEMQRWV","AKEMQRWVM","KEMQRWVME","EMQRWVMEK","MQRWVMEKT",
        "QRWVMEKTA","RWVMEKTAE","WVMEKTAEH","VMEKTAEHF"]
piv = df[df.peptide.isin(sdf4)].pivot(index="peptide", columns="allele", values="mhcflurry_affinity")
piv = piv.reindex(sdf4)

fig, ax = plt.subplots(figsize=(9, 4.4))
x = range(len(sdf4)); w = 0.15
for i, al in enumerate(["HLA-A*02:06","HLA-A*30:01","HLA-B*08:01","HLA-B*53:01","HLA-C*07:01"]):
    vals = [min(piv.loc[p, al], 35000) for p in sdf4]
    ax.bar([xi + (i-2)*w for xi in x], vals, width=w, label=al)
ax.set_xticks(list(x)); ax.set_xticklabels(sdf4, rotation=25, ha="right", fontsize=7)
ax.axhline(50, color="green", ls="--", lw=1); ax.text(8.6, 60, "strong<50", fontsize=7, color="green")
ax.axhline(500, color="orange", ls="--", lw=1); ax.text(8.6, 600, "weak<500", fontsize=7, color="orange")
ax.set_yscale("log"); ax.set_ylabel("affinity (nM, log; capped 35k)")
ax.set_title("Fig 7 — SDF4 A120V 9mer tiling × full diploid HLA (MHCflurry 2.1.2, fresh run).\n"
             "All A*02:06 >6,800 nM (non-binders). Only EMQRWVMEK/A*30:01 = 192 nM (weak). "
             "Focus ALAAAVAAV (16.94, A*02:06-only) is NOT SDF4-derived.", fontsize=8)
ax.legend(fontsize=7, ncol=5, loc="upper left")
fig.tight_layout(); fig.savefig(FIG/"fig07_sdf4_full_allele.png", bbox_inches="tight")
print("wrote fig07")

# print exact table for VALIDATION
print(df[df.peptide.isin(sdf4+["ALAAAVAAV","LLFGYPVYV"])].sort_values(["peptide","allele"]).to_string(index=False))

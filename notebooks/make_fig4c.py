# Figure 4C — candidate evidence matrix across tiers
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = os.path.join("..", "..", "results")
FIG = os.path.join("..", "..", "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"figure.dpi": 140, "savefig.dpi": 300, "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.spines.left": False, "axes.spines.bottom": False})

ev = pd.read_csv(os.path.join(RES, "j45_candidate_evidence_matrix.csv"))

genes = ["GLIS3", "TP73", "EOMES", "EGR3", "FOXE1", "KLF5", "SETBP1"]
NA_CELLS = {("SETBP1", 0), ("SETBP1", 1), ("SETBP1", 4)}
REPCEN = {"GLIS3": 0.47, "KLF5": 0.55, "TP73": 0.34, "EOMES": 0.36,
          "EGR3": 0.29, "FOXE1": 0.55, "SETBP1": np.nan}
cols = [("Hs27 perm\nFDR < 0.05", "hs27_fdr_perm", lambda v: v < 0.05),
        ("RPE1 direction\nconsistent", "dir", None),
        ("Coloc\nPP_H4 > 0.9", "coloc_PP_H4_max", lambda v: v > 0.9),
        ("Open Targets\nscore", "ot_score", lambda v: not pd.isna(v)),
        ("Hs27 rep_cen\n(guide split)", "repcen", lambda v: v > 0.1)]

rows = []
for g in genes:
    r = ev[ev.gene == g].iloc[0]
    hits = []
    for label, key, test in cols:
        if key == "repcen":
            v = REPCEN.get(g, np.nan)
            h = False if pd.isna(v) else bool(test(v))
            hits.append((h, v if not pd.isna(v) else None))
            continue
        if key == "dir":
            v = r.get("rpe1_score", np.nan)
            h = (not pd.isna(v)) and np.sign(v) == np.sign(r.get("hs27_score", 0))
            hits.append((h, "+" if h else None))
        else:
            v = r.get(key, np.nan)
            h = False if pd.isna(v) else bool(test(v))
            hits.append((h, v))
    rows.append(hits)

fig, ax = plt.subplots(figsize=(6.4, 3.4))
cmap = {True: "#1d4ed8", False: "#e5e7eb"}
for i, hits in enumerate(rows):
    for j, (h, v) in enumerate(hits):
        if (genes[i], j) in NA_CELLS:
            ax.scatter(j, len(genes)-1-i, s=560, marker="s",
                       c="#f8fafc", edgecolors="0.6", linewidths=0.5,
                       hatch="///", zorder=2)
            ax.text(j, len(genes)-1-i, "N/A", ha="center", va="center",
                    fontsize=5.5, color="0.45", zorder=3)
            continue
        ax.scatter(j, len(genes)-1-i, s=560, marker="s",
                   c=cmap[h], edgecolors="0.6", linewidths=0.5, zorder=2)
        if h and not pd.isna(v):
            txt = v if isinstance(v, str) else ("<0.01" if v < 0.01 else (f"{v:.3f}" if abs(v) < 10 else f"{v:.1f}"))
            ax.text(j, len(genes)-1-i, txt, ha="center", va="center",
                    fontsize=6.5, color="white", zorder=3)
ax.set_xlim(-0.6, len(cols)-0.4)
ax.set_ylim(-0.6, len(genes)-0.4)
ax.set_xticks(range(len(cols)), [c[0] for c in cols], fontsize=7)
ax.set_yticks(range(len(genes)), genes[::-1], fontsize=8.5)
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
ax.set_title("Candidate evidence tiers (screen permutation vs orthogonal genetic annotation)", fontsize=9)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor="#1d4ed8", edgecolor="0.6", linewidth=0.5, label="tier met"),
                   Patch(facecolor="#e5e7eb", edgecolor="0.6", linewidth=0.5, label="not met"),
                   Patch(facecolor="#f8fafc", edgecolor="0.6", linewidth=0.5, hatch="///",
                         label="not assayed (N/A)")],
          loc="upper right", bbox_to_anchor=(1.42, 1.02), frameon=False, fontsize=7,
          handlelength=1.2, handleheight=1.2)
fig.savefig(os.path.join(FIG, "fig4c_evidence_matrix.png"), bbox_inches="tight")


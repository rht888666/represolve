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
cols = [("Hs27 perm\nFDR<0.05", "hs27_fdr_perm", lambda v: v < 0.05),
        ("RPE1 direction\nconsistent", "dir", None),
        ("Coloc\nPP_H4>0.9", "coloc_PP_H4_max", lambda v: v > 0.9),
        ("AG2 DNase\nemp. sig.", "ag2_emp_q_dn", lambda v: v < 0.05),
        ("Open Targets\nscore", "ot_score", lambda v: not pd.isna(v))]

rows = []
for g in genes:
    r = ev[ev.gene == g].iloc[0]
    hits = []
    for label, key, test in cols:
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
        ax.scatter(j, len(genes)-1-i, s=560, marker="s",
                   c=cmap[h], edgecolors="0.6", linewidths=0.5, zorder=2)
        if h and not pd.isna(v):
            txt = v if isinstance(v, str) else (f"{v:.2f}" if abs(v) < 10 else f"{v:.1f}")
            ax.text(j, len(genes)-1-i, txt, ha="center", va="center",
                    fontsize=6.5, color="white", zorder=3)
ax.set_xlim(-0.6, len(cols)-0.4)
ax.set_ylim(-0.6, len(genes)-0.4)
ax.set_xticks(range(len(cols)), [c[0] for c in cols], fontsize=7)
ax.set_yticks(range(len(genes)), genes[::-1], fontsize=8.5)
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
ax.set_title("Candidate evidence tiers (screen permutation vs genetic)", fontsize=9)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0],[0], marker="s", color="w", markerfacecolor="#1d4ed8",
                          markersize=9, label="tier met"),
                   Line2D([0],[0], marker="s", color="w", markerfacecolor="#e5e7eb",
                          markeredgecolor="0.6", markersize=9, label="not met")],
          loc="upper right", bbox_to_anchor=(1.42, 1.02), frameon=False, fontsize=7)
fig.savefig(os.path.join(FIG, "fig4c_evidence_matrix.png"), bbox_inches="tight")
plt.show()

# -*- coding: utf-8 -*-
# Figure audit round-1 fixes:
#  fig1b add DAE; fig2a add DAE column; fig2c display names;
#  fig4a de-overlap labels; fig4b clean contrast labels; fig5b move legend.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = os.path.join("..", "..", "results")
FIG = os.path.join("..", "..", "figures")
plt.rcParams.update({"figure.dpi": 140, "savefig.dpi": 300, "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})

MODELS = ["lw", "svd20", "moppca5", "nn10", "nn50", "dae"]
MODEL_LABELS = {"lw": "LW shrinkage", "svd20": "SVD-20", "moppca5": "moPPCA-5",
                "nn10": "centroid-10", "nn50": "centroid-50", "dae": "DAE"}
CONTEXTS = [
    ("K562 essential", "j22_k562_foldfrozen_centroid.csv", "j47_k562_dae_foldfrozen.csv"),
    ("K562 GWPS", "j22_k562_gwps_foldfrozen_centroid.csv", "j47_gwps_dae_foldfrozen.csv"),
    ("Norman", "j49_norman_foldfrozen.csv", "j47_norman_dae.csv"),
    ("RPE1", "j49_rpe1_foldfrozen.csv", "j47_rpe1_dae.csv"),
    ("sci-Plex3", "j49_sciplex3.csv", "j47_sciplex3_dae.csv"),
]
CTX_DISPLAY = {"K562_essential": "K562 essential", "K562_GWPS": "K562 GWPS",
               "Norman": "Norman", "RPE1_Replogle": "RPE1",
               "sciPlex3_chemical": "sci-Plex3"}

def load_ctx(cname, f5, fdae):
    df = pd.read_csv(os.path.join(RES, f5))
    dae = pd.read_csv(os.path.join(RES, fdae))
    key = "pert" if "pert" in df.columns and "pert" in dae.columns else None
    if key:
        df = df.merge(dae[["pert", "dae_A_cen", "dae_B_cen"]], on="pert", how="left")
    else:
        df["dae_A_cen"] = dae["dae_A_cen"].values[: len(df)]
        df["dae_B_cen"] = dae["dae_B_cen"].values[: len(df)]
    return df

# ---------- fig1b: same vs independent half, now incl. DAE ----------
df, dae = (pd.read_csv(os.path.join(RES, CONTEXTS[0][1])),
           pd.read_csv(os.path.join(RES, CONTEXTS[0][2])))
rows = []
for m in MODELS:
    col_a, col_b = f"{m}_A_cen", f"{m}_B_cen"
    if m == "dae":
        rows.append({"model": "DAE", "same half": dae.dae_A_cen.mean(),
                     "independent half": dae.dae_B_cen.mean()})
    else:
        rows.append({"model": MODEL_LABELS[m],
                     "same half": df[col_a].mean(),
                     "independent half": df[col_b].mean()})
t = pd.DataFrame(rows).set_index("model")
ax = t.plot.bar(figsize=(7.6, 3.2), color=["#93c5fd", "#1d4ed8"])
ax.set_ylabel("centered correlation")
ax.set_title("Same-half scores exceed independent-half scores for every model\n"
             "(K562 essential, n=%d perturbations)" % len(df))
ax.legend(frameon=False)
plt.xticks(rotation=0)
ax.get_figure().savefig(os.path.join(FIG, "fig1b_same_vs_indep.png"),
                        bbox_inches="tight")
plt.close()

# ---------- fig2a: heatmap incl. DAE ----------
mat = np.zeros((len(CONTEXTS), len(MODELS)))
for i, (cname, f5, fd) in enumerate(CONTEXTS):
    d = load_ctx(cname, f5, fd)
    for j, m in enumerate(MODELS):
        mat[i, j] = d[f"{m}_B_cen"].mean()
fig, ax = plt.subplots(figsize=(7.0, 3.2))
im = ax.imshow(mat, cmap="viridis", aspect="auto")
ax.set_xticks(range(len(MODELS)), [MODEL_LABELS[m] for m in MODELS])
ax.set_yticks(range(len(CONTEXTS)), [c[0] for c in CONTEXTS])
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        ax.text(j, i, f"{mat[i,j]:.3f}", ha="center", va="center",
                color="white" if mat[i, j] < mat.max() * 0.6 else "black",
                fontsize=8)
ax.set_title("Independent-half centered score: five contexts, five orderings")
fig.colorbar(im, ax=ax, shrink=0.8)
fig.savefig(os.path.join(FIG, "fig2a_ranking_heatmap.png"), bbox_inches="tight")
plt.close()

# ---------- fig2c: display context names ----------
ci = pd.read_csv(os.path.join(RES, "j53_pairwise_ranking_ci.csv"))
lead = ci.groupby("context").first().reset_index()
lead["disp"] = lead.context.map(CTX_DISPLAY).fillna(lead.context)
fig, ax = plt.subplots(figsize=(6.8, 2.8))
ypos = np.arange(len(lead))
for i, r in lead.iterrows():
    ax.plot([r.ci95_low, r.ci95_high], [i, i], color="0.3", lw=2)
    ax.plot(r.mean_diff, i, "o", color="#b91c1c")
    ax.text(r.ci95_high + 0.004, i, f"{r.better}>{r.worse}",
            va="center", fontsize=8)
ax.axvline(0, color="0.6", ls=":")
ax.set_yticks(ypos, lead["disp"])
ax.set_xlabel("mean paired difference (95% bootstrap CI)")
ax.set_title("Winner vs rank-2: significant in all five contexts")
ax.set_ylim(-0.6, len(lead) - 0.2)
fig.savefig(os.path.join(FIG, "fig2c_pairwise_ci.png"), bbox_inches="tight")
plt.close()

# ---------- fig4a: de-overlapped hit labels ----------
pn = pd.read_csv(os.path.join(RES, "j32_keloid_perm_null.csv"))
hs = pn[(pn.cell == "hs27") & (pn["variant"] == "full")].copy()
hs["mlog10"] = -np.log10(hs.fdr_perm.clip(lower=1e-6))
hits = hs[hs.fdr_perm < 0.05].sort_values("score")
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.scatter(hs.score, hs.mlog10, s=4, c="0.7", label="all TFs")
ax.scatter(hits.score, hits.mlog10, s=18, c="#b91c1c", label="FDR<0.05")
offsets = {"GLIS3": (-6, 8), "TP73": (4, 6), "EOMES": (-14, 8),
           "EGR3": (4, -4), "FOXE1": (4, 2), "KLF5": (4, 2)}
for _, r in hits.iterrows():
    dx, dy = offsets.get(r.tf, (4, 2))
    ax.annotate(r.tf, (r.score, r.mlog10), fontsize=7,
                xytext=(dx, dy), textcoords="offset points")
ax.axhline(-np.log10(0.05), color="0.4", ls=":")
ax.set_xlabel("permutation-null TF score (Hs27)")
ax.set_ylabel("-log10(FDR)")
ax.set_title(f"6 of {len(hs)} Hs27 TFs pass the fixed-size permutation null")
ax.legend(frameon=False)
fig.savefig(os.path.join(FIG, "fig4a_perm_volcano.png"), bbox_inches="tight")
plt.close()

# ---------- fig4b: clean contrast labels ----------
au = pd.read_csv(os.path.join(RES, "j28_signature_external_auc.csv"))
pretty = {"Keloid_vs_Control": "vs control",
          "Keloid_vs_SurroundingSkin": "vs surrounding skin",
          "Keloid_vs_HealthySkin": "vs healthy skin",
          "Keloid_vs_Regression": "vs regressing",
          "Keloid_vs_ALL_NONKELOID": "vs all non-keloid",
          "Keloid_vs_HypertrophicScar": "vs hypertrophic scar",
          "Keloid_vs_NormotrophicScar": "vs normotrophic scar"}
au["label"] = au.contrast.map(pretty).fillna(au.contrast) + \
    "\n(" + au.cohort + ", n=" + au.n_pos.astype(str) + "/" + \
    au.n_neg.astype(str) + ")"
fig, ax = plt.subplots(figsize=(8.0, 3.4))
colors = ["#1d4ed8" if a > 0.9 else ("#b91c1c" if a < 0.6 else "#f59e0b")
          for a in au.auroc]
ax.bar(range(len(au)), au.auroc, color=colors)
ax.set_xticks(range(len(au)), au["label"], rotation=35, ha="right", fontsize=6.5)
ax.axhline(0.5, color="0.5", ls=":")
ax.set_ylabel("AUROC")
ax.set_title("Signature separates keloid from non-keloid — except hypertrophic scar")
fig.savefig(os.path.join(FIG, "fig4b_external_auc.png"), bbox_inches="tight")
plt.close()
print("contrast names:", au.contrast.tolist())

# ---------- fig5b: legend to upper-left, named points annotated ----------
dn = pd.read_csv(os.path.join(RES, "ag2_14_empirical_null_DNASE.csv"))
order = dn.sort_values("maxabs").reset_index(drop=True)
fig, ax = plt.subplots(figsize=(6.4, 3.6))
sig = order.q_emp < 0.05
ax.scatter(order.index[~sig], order.maxabs[~sig], s=8, c="0.7",
           label="credible-set variants")
ax.scatter(order.index[sig], order.maxabs[sig], s=22, c="#b91c1c",
           label="empirical adj. p<0.05", zorder=3)
named = {"rs72905594": "SETBP1", "rs2254759": "ASIP/NEDD4 locus",
         "rs140716753": "JAG1"}
for _, r in dn[dn.rsid.isin(named)].iterrows():
    i = order.index[order.rsid == r.rsid][0]
    ax.annotate(f"{r.rsid}\n({named[r.rsid]})", (i, r.maxabs), fontsize=7,
                xytext=(-70, -4), textcoords="offset points",
                arrowprops=dict(arrowstyle="-", color="0.4", lw=0.6))
ax.set_xlabel("rank by |DNase effect| (246 credible-set variants)")
ax.set_ylabel("max |AG2 DNase log2 effect|")
ax.set_title("AG2 empirical null: rs72905594 (SETBP1) exceeds all 560 matched nulls")
ax.legend(frameon=False, fontsize=7.5, loc="upper left")
fig.savefig(os.path.join(FIG, "fig5b_ag2_null.png"), bbox_inches="tight")
plt.close()
print("done")

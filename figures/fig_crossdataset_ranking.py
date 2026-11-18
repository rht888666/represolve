"""Core benchmark figure: model ranking reverses across replicate-
reliability regimes and perturbation densities.

Reads results/j51_crossdataset_model_ranking.csv (frozen merge of
j22/j47/j49 outputs) and renders a heatmap of independent-half
centered scores, annotated with per-context median replicate agreement
and perturbation counts.

Usage: python figures/fig_crossdataset_ranking.py <j51.csv> <out.png>
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main(csv_path: str, out_path: str) -> None:
    df = pd.read_csv(csv_path, index_col=0)
    models = ["lw", "svd20", "moppca5", "nn10", "nn50", "dae"]
    labels = {
        "lw": "Ledoit-Wolf", "svd20": "SVD-20", "moppca5": "moPPCA-5",
        "nn10": "nn-10", "nn50": "nn-50", "dae": "DAE",
    }
    mat = df[models].T  # models x contexts
    mat.columns = [c.replace("_", "\n") for c in mat.columns]

    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [3.2, 1]},
        constrained_layout=True)

    im = ax.imshow(mat.values, cmap="viridis", aspect="auto",
                   vmin=0, vmax=mat.values.max())
    ax.set_xticks(range(len(mat.columns)), mat.columns, fontsize=8)
    ax.set_yticks(range(len(models)), [labels[m] for m in models], fontsize=9)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.3f}", ha="center", va="center",
                        fontsize=7, color="w" if v > 0.45 else "black")
    ax.set_title("Independent-half centered score (B_cen)")
    fig.colorbar(im, ax=ax, fraction=0.046)

    # regime annotation: median replicate agreement per context
    order = df.index
    rep = df.loc[order, "rep_cen_median"].values
    nperts = df.loc[order, "n_perts"].values
    xpos = np.arange(len(order))
    ax2.barh(xpos, rep, color="steelblue")
    ax2.set_yticks(xpos, [c.replace("_", "\n") for c in order], fontsize=8)
    ax2.set_xlim(0, 1.05)
    ax2.set_xlabel("median rep_cen")
    for i, (r, n) in enumerate(zip(rep, nperts)):
        ax2.text(r + 0.02, i, f"{r:.2f} (n={int(n)})", va="center", fontsize=7)
    ax2.set_title("Replicate agreement")
    ax2.invert_yaxis()

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print("wrote", out_path)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

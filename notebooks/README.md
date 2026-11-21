# Figure notebooks

One notebook per main figure, named after the figure number. Each notebook:

1. loads the frozen result CSVs listed in `paper_map/FIGURE_TO_CODE.md`,
2. regenerates the figure exactly as it appears in the manuscript,
3. is committed **with executed output** so reviewers can verify without
   re-running.

Shipped executed notebooks (matching `SUPPLEMENT_IMETA.md` S/F1–S/F6):

- `fig1_protocol.ipynb` — protocol schematic + same-half vs independent-half inflation (`figures/fig1a/1b`)
- `fig2_ranking_regimes.ipynb` — cross-context ranking heatmap, `rep_cen` ECDFs, and paired bootstrap CIs (`figures/fig2a–c`)
- `fig3_reversal_dae.ipynb` — within-context decile reversal, masked-MLP, and 5-fold scGPT audits (`figures/fig3_decile_reversal`, `fig3b`)
- `fig4_keloid.ipynb` (`+ make_fig4c.py`) — keloid CRISPRa permutation volcano, external cohort AUROCs, and 5-tier evidence matrix (`figures/fig4a/4b/4c`)
- `fig5_genetics_ag2.ipynb` — multi-cohort colocalization and AlphaGenome empirical-null ranking (`figures/fig5a/5b`)
- `fig6_decision_schematic.ipynb` — 2D routing phase map (`obs_agreement` vs `d_cos10`) and pre-registered external validation (`figures/fig6`)

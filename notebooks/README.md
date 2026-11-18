# Figure notebooks

One notebook per main figure, named after the figure number. Each notebook:

1. loads the frozen result CSVs listed in `paper_map/FIGURE_TO_CODE.md`,
2. regenerates the figure exactly as it appears in the manuscript,
3. is committed **with executed output** so reviewers can verify without
   re-running.

Planned notebooks (created as figures are finalized):

- `fig1_protocol_and_leakage.ipynb` — evaluation task diagram + same-vs-
  independent-half inflation by response regime (j11/j22 outputs)
- `fig2_reliability_and_ranking.ipynb` — per-TF centered scores, Hs27/RPE1
  responsive strata, uncertainty across perturbations
- `fig3_k562_gwps.ipynb` — essential vs GWPS ordering reversal, decile
  ranking, mask/k-choice/partition sensitivity panels
- `fig4_keloid_application.ipynb` — multi-cohort signature, screen
  permutation-null hits, orthogonal validation summary
- `fig5_genetics_ag2.ipynb` — evidence matrix, SETBP1 convergent chain,
  AG2 empirical null, state-conditioned deltas
- `tableS1_evidence_matrix.ipynb` — candidate × evidence-layer master table

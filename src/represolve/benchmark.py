"""Fold-frozen masked-completion benchmark.

The core evaluation loop used in the paper. Given two replicate measurement
halves of the same perturbations (A and B), perturbations are split into
K folds; gene selection and all model fitting happen inside each training
fold, and each held-out perturbation is predicted from a random observed-gene
mask, then scored against BOTH halves. Comparing the same- half (A) and
independent-half (B) scores quantifies how much of a benchmark number is
shared-response recovery versus reproducible target-specific prediction.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold

from .metrics import pearson, replicate_agreement
from .models import default_models


def masked_completion_benchmark(
    half_a: np.ndarray,
    half_b: np.ndarray,
    pert_names: np.ndarray,
    genes: np.ndarray,
    n_folds: int = 5,
    n_genes_selected: int = 5000,
    mask_fraction: float = 0.1,
    mask_seed_base: int = 1000,
    models: list | None = None,
    random_state: int = 0,
) -> pd.DataFrame:
    """Run the fold-frozen masked-completion benchmark.

    Parameters
    ----------
    half_a, half_b : (n_perturbations, n_genes) arrays
        Mean expression profiles of the same perturbations measured in two
        independent partitions (e.g. guide halves or cell-pool halves).
    pert_names, genes : (n_perturbations,), (n_genes,) arrays
    n_genes_selected : int
        Number of top-variance genes selected inside each training fold.
    mask_fraction : float
        Fraction of selected genes treated as observed for each test
        perturbation (the paper uses 0.1).
    mask_seed_base : int
        Per-perturbation mask seed = mask_seed_base + perturbation index,
        making masks deterministic and comparable across runs.
    models : list of BaseCompletion, optional
        Defaults to `default_models()` (svd20, lw, moppca5, nn10, nn50).

    Returns
    -------
    DataFrame with one row per (perturbation, fold): rep_cen plus
    {model}_A_cen and {model}_B_cen centered scores.
    """
    if models is None:
        models = default_models()
    half_a = np.asarray(half_a, dtype=np.float64)
    half_b = np.asarray(half_b, dtype=np.float64)
    mean_all = (half_a + half_b) / 2.0

    rows = []
    for fold, (train, test) in enumerate(
        KFold(n_folds, shuffle=True, random_state=random_state).split(pert_names)
    ):
        selected = np.argsort(-mean_all[train].var(axis=0), kind="stable")[
            :n_genes_selected
        ]
        T = mean_all[np.ix_(train, selected)]
        mu = T.mean(axis=0)
        for m in models:
            m.fit(T)

        for idx in test:
            observed = (
                np.random.RandomState(mask_seed_base + idx).rand(len(selected))
                < mask_fraction
            )
            hidden = ~observed
            hid_pos = np.flatnonzero(hidden)
            aa = half_a[idx, selected]
            bb = half_b[idx, selected]
            row = {
                "pert": pert_names[idx],
                "fold": fold,
                "n_observed": int(observed.sum()),
                "rep_cen": replicate_agreement(aa[hid_pos], bb[hid_pos], mu[hid_pos]),
            }
            for m in models:
                pred = m.predict(aa, observed)
                row[f"{m.name}_A_cen"] = pearson(
                    pred[hid_pos] - mu[hid_pos], aa[hid_pos] - mu[hid_pos]
                )
                row[f"{m.name}_B_cen"] = pearson(
                    pred[hid_pos] - mu[hid_pos], bb[hid_pos] - mu[hid_pos]
                )
            rows.append(row)
    return pd.DataFrame(rows)


def summarize(df: pd.DataFrame, models: list | None = None) -> pd.DataFrame:
    """Per-model mean A-half and B-half centered scores and the A-B gap."""
    model_names = models or [
        c[: -len("_A_cen")] for c in df.columns if c.endswith("_A_cen")
    ]
    out = []
    for m in model_names:
        a, b = df[f"{m}_A_cen"], df[f"{m}_B_cen"]
        out.append(
            {
                "model": m,
                "same_half": a.mean(),
                "independent_half": b.mean(),
                "inflation": a.mean() - b.mean(),
                "n": a.notna().sum(),
            }
        )
    return pd.DataFrame(out).sort_values("independent_half", ascending=False)

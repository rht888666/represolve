"""Prediction-time regime diagnostics for routing.

`observed_agreement` measures between-half centered agreement on the
observed-gene mask only — no hidden-gene information — so the screen's
reliability regime is declarable before hidden-gene scoring.

`cosine_knn_density` is a scale-invariant perturbation-density diagnostic:
the median cosine distance from each perturbation to its k-th nearest
neighbour on L2-normalised profiles. It is invariant to signal magnitude
and preprocessing scale, unlike raw Euclidean kNN distance.

`route` maps the two diagnostics to the first-generation routing map
(Table 3 of the manuscript). It is a falsifiable hypothesis, not a fitted
model: new screens re-score on frozen partitions and update the map.
"""
from __future__ import annotations

import numpy as np

from .metrics import pearson


def observed_agreement(half_a_obs: np.ndarray, half_b_obs: np.ndarray,
                       train_mean: np.ndarray) -> float:
    """Centered correlation between the two halves on observed genes."""
    return float(pearson(half_a_obs - train_mean, half_b_obs - train_mean))


def cosine_knn_density(train_profiles: np.ndarray, k: int = 10,
                       center: bool = True) -> float:
    """Median cosine distance to the k-th nearest neighbour.

    train_profiles: (n_perturbations, n_genes) matrix (e.g. the pooled
    two-half mean on the selected genes). Returns a value in [0, 2];
    smaller means a denser perturbation grid.
    """
    X = np.asarray(train_profiles, dtype=np.float64)
    if center:
        X = X - X.mean(axis=0)
    nrm = np.linalg.norm(X, axis=1, keepdims=True)
    nrm[nrm == 0] = 1.0
    U = X / nrm
    G = U @ U.T
    np.fill_diagonal(G, -np.inf)
    # distance to the k-th nearest neighbour: k-th largest cosine similarity
    # (1-indexed), i.e. element k-1 of the distance partition. argpartition
    # does not sort, so the k-th element must be taken by position, not
    # by slicing an unordered top-k set.
    return float(np.median(np.partition(1.0 - G, k - 1, axis=1)[:, k - 1]))


def route(obs_agreement: float, cos_density: float | None = None,
          n_perts: int | None = None) -> str:
    """First-generation routing hypothesis (Table 3), fit on the same
    nine contexts it summarizes — in-sample, falsifiable.

    Returns one of: 'flag_near_floor', 'centroid_nn50', 'centroid_nn10',
    'nn10_dense_grid', 'lw_shrinkage', 'inherit_design_priors'.
    """
    if obs_agreement is None:
        return "inherit_design_priors"
    if obs_agreement < 0.05:
        return "flag_near_floor"
    if obs_agreement < 0.15:
        return "centroid_nn50"
    if obs_agreement < 0.8:
        return "centroid_nn10"
    if cos_density is not None and cos_density < 0.05:
        return "nn10_dense_grid"
    return "lw_shrinkage"

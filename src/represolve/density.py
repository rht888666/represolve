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
    dense_ok = cos_density is not None and cos_density < 0.05
    # Table 3's dense-grid branch additionally requires a large perturbation
    # panel (n ≳ 10^3): a small dense screen is routed to shrinkage instead.
    if dense_ok and (n_perts is None or n_perts >= 1000):
        return "nn10_dense_grid"
    return "lw_shrinkage"


# Soft-regime router: continuous experts in (obs_agreement, d_cos10) space.
# Prototypes are the Table-3 regime centroids (not fitted to held-out scores);
# LOCO analyses in the paper re-estimate prototypes from training contexts.
ROUTE_SOFT_PROTOTYPES = {
    "flag_near_floor": np.array([-0.02, 0.90]),
    "centroid_nn50": np.array([0.07, 0.85]),
    "centroid_nn10": np.array([0.35, 0.55]),
    "nn10_dense_grid": np.array([0.95, 0.01]),
    "lw_shrinkage": np.array([0.90, 0.25]),
}

_ROUTE_SOFT_TO_MODEL = {
    "flag_near_floor": "flag",
    "centroid_nn50": "nn50",
    "centroid_nn10": "nn10",
    "nn10_dense_grid": "nn10",
    "lw_shrinkage": "lw",
}


def _soft_feature(obs_agreement: float, cos_density: float) -> np.ndarray:
    """z = (obs_agreement, d_cos10) with mild clipping for numerical stability."""
    oa = float(np.clip(obs_agreement, -0.2, 1.05))
    dc = float(np.clip(cos_density, 0.0, 2.0))
    return np.array([oa, dc], dtype=np.float64)


def soft_route_weights(
    obs_agreement: float,
    cos_density: float,
    temperature: float = 0.08,
    prototypes: dict[str, np.ndarray] | None = None,
    scale: tuple[float, float] = (1.0, 1.2),
) -> dict[str, float]:
    """RBF soft weights over regime experts in (obs_agreement, d_cos10).

    Higher temperature → softer mixture. ``scale`` stretches the two axes
    so density and agreement contribute on comparable length scales.
    Returns a dict keyed by hard-route action labels (sums to 1).
    """
    if temperature <= 0:
        raise ValueError("temperature must be positive")
    protos = prototypes or ROUTE_SOFT_PROTOTYPES
    z = _soft_feature(obs_agreement, cos_density)
    s = np.asarray(scale, dtype=np.float64)
    logits = []
    keys = list(protos.keys())
    for k in keys:
        d = (z - np.asarray(protos[k], dtype=np.float64)) / s
        logits.append(-0.5 * float(d @ d) / temperature)
    logits = np.asarray(logits, dtype=np.float64)
    logits -= logits.max()
    w = np.exp(logits)
    w /= w.sum()
    return {k: float(wi) for k, wi in zip(keys, w)}


def route_soft(
    obs_agreement: float,
    cos_density: float,
    temperature: float = 0.08,
    prototypes: dict[str, np.ndarray] | None = None,
) -> str:
    """Argmax of soft_route_weights — continuous analogue of Table 3."""
    w = soft_route_weights(
        obs_agreement, cos_density, temperature=temperature, prototypes=prototypes
    )
    return max(w, key=w.get)


def route_soft_model(
    obs_agreement: float,
    cos_density: float,
    temperature: float = 0.08,
    prototypes: dict[str, np.ndarray] | None = None,
) -> str:
    """Map soft route action to an estimator name (or ``flag``)."""
    action = route_soft(
        obs_agreement, cos_density, temperature=temperature, prototypes=prototypes
    )
    return _ROUTE_SOFT_TO_MODEL[action]


def adaptive_nn_k(
    obs_agreement: float,
    k_low: int = 50,
    k_mid: int = 10,
    lo: float = 0.05,
    hi: float = 0.80,
) -> float:
    """Continuous neighborhood size along the mid-reliability band.

    Maps obs_agreement in [lo, hi] from ``k_low`` (near floor) to ``k_mid``
    (mid band). Outside the band, returns the endpoint k. Useful for
    blending nn10/nn50 scores without a hard cut.
    """
    oa = float(obs_agreement)
    if oa <= lo:
        return float(k_low)
    if oa >= hi:
        return float(k_mid)
    t = (oa - lo) / (hi - lo)
    return float(k_low + t * (k_mid - k_low))


def blend_nn_scores(
    nn10_score: float,
    nn50_score: float,
    obs_agreement: float,
    lo: float = 0.05,
    hi: float = 0.80,
) -> float:
    """Convex blend of nn10/nn50 guided by adaptive_nn_k."""
    k = adaptive_nn_k(obs_agreement, lo=lo, hi=hi)
    # k=50 → all nn50; k=10 → all nn10
    w10 = (50.0 - k) / 40.0
    w10 = float(np.clip(w10, 0.0, 1.0))
    return w10 * float(nn10_score) + (1.0 - w10) * float(nn50_score)

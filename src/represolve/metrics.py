"""Metrics for replicate-resolved perturbation-prediction evaluation.

All scores are *per-perturbation* centered Pearson correlations: both the
prediction and the measurement are expressed as deviations from the training
mean profile before correlating, so a predictor is credited only for
perturbation-specific signal, not for recovering shared variation.
"""
from __future__ import annotations

import numpy as np


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    """Safe Pearson correlation; returns nan on degenerate input."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    sx = x.std()
    sy = y.std()
    if sx == 0 or sy == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def centered_score(pred: np.ndarray, measured: np.ndarray, mu: np.ndarray) -> float:
    """Centered correlation between a prediction and a measured profile.

    Parameters
    ----------
    pred, measured : (n_genes,) arrays
    mu : (n_genes,) array
        Reference mean profile (e.g. the training-fold mean). Both vectors
        are centered on ``mu`` before correlating.
    """
    return pearson(pred - mu, measured - mu)


def replicate_agreement(half_a: np.ndarray, half_b: np.ndarray, mu: np.ndarray) -> float:
    """rep_cen: centered correlation between two replicate measurement halves.

    This is the empirical reliability axis of the benchmark — the agreement
    between two independent measurement partitions of the *same* perturbation.
    It is a descriptive statistic, not a theoretical prediction ceiling.
    """
    return pearson(half_a - mu, half_b - mu)


def neighbour_similarity(
    train_profiles: np.ndarray,
    train_mu: np.ndarray,
    obs_pos: np.ndarray,
    test_obs: np.ndarray,
) -> np.ndarray:
    """Centered correlation of every training perturbation to the test
    perturbation's observed-gene effect. Used by the nearest-centroid model
    and by coverage-aware analyses.

    Returns one similarity per training perturbation (-inf where undefined).
    """
    train_eff = train_profiles[:, obs_pos] - train_mu[obs_pos]
    train_eff = train_eff - train_eff.mean(axis=1, keepdims=True)
    test_eff = test_obs - train_mu[obs_pos]
    test_eff = test_eff - test_eff.mean()
    tn = np.linalg.norm(train_eff, axis=1)
    den = tn * np.linalg.norm(test_eff)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = (train_eff @ test_eff) / den
    out[~np.isfinite(out)] = -np.inf
    return out

"""Permutation nulls and type-I calibration.

Two nulls used in the paper:
1. Label-permutation null for signature-projection screen scores — permute
   the signature labels to calibrate which screen hits exceed chance.
2. Synthetic no-signal screens — generate screens where no true signal exists
   and count FDR<0.05 false positives, verifying the permutation procedure
   itself is calibrated (paper: ~6.5% at alpha=0.05 across simulations).
"""
from __future__ import annotations

import numpy as np
from scipy.stats import false_discovery_control  # type: ignore[attr-defined]


def _bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values."""
    try:
        return false_discovery_control(p)
    except Exception:
        order = np.argsort(p)
        q = np.empty_like(p, dtype=float)
        n = len(p)
        ranked = p[order] * n / np.arange(1, n + 1)
        q[order] = np.minimum.accumulate(ranked[::-1])[::-1]
        return np.clip(q, 0, 1)


def permutation_null(observed: np.ndarray, permuted: np.ndarray,
                     alternative: str = "two-sided") -> np.ndarray:
    """Empirical p-values from a permutation distribution.

    observed : (n_tests,) observed statistics
    permuted : (n_perms, n_tests) or (n_perms,) permutation draws
    """
    observed = np.asarray(observed, dtype=float)
    permuted = np.asarray(permuted, dtype=float)
    if permuted.ndim == 1:
        permuted = permuted[:, None]
    if alternative == "two-sided":
        extreme = np.abs(permuted) >= np.abs(observed)[None, :]
    elif alternative == "greater":
        extreme = permuted >= observed[None, :]
    else:
        extreme = permuted <= observed[None, :]
    return (extreme.sum(axis=0) + 1) / (permuted.shape[0] + 1)


def bh_qvalues(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR q-values (alias of _bh for public API)."""
    return _bh(np.asarray(p, dtype=float))


def type1_check(simulate_screen, n_sims: int = 200, alpha: float = 0.05,
                seed: int = 0) -> dict:
    """Run `simulate_screen(rng)` under the null and count FDR<alpha hits.

    `simulate_screen` must return a 1-D array of empirical p-values for a
    no-signal screen. Returns the fraction of simulations with any
    FDR-significant hit and the mean number of false positives.
    """
    rng = np.random.RandomState(seed)
    any_fp = 0
    total_fp = 0
    for _ in range(n_sims):
        p = np.asarray(simulate_screen(rng), dtype=float)
        fp = int((_bh(p) < alpha).sum())
        any_fp += fp > 0
        total_fp += fp
    return {
        "n_sims": n_sims,
        "alpha": alpha,
        "any_hit_rate": any_fp / n_sims,
        "mean_false_positives": total_fp / n_sims,
    }

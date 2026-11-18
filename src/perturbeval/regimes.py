"""Regime-stratified analysis: model ranking by replicate-agreement decile.

The paper's headline analysis: rank perturbations by replicate agreement
(rep_cen) into deciles, then re-compute per-model independent-half scores
within each decile. The central finding is that model ranking reverses
across the detectable-response regime — there is no universally optimal
model.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def assign_deciles(df: pd.DataFrame, score_col: str = "rep_cen", n: int = 10,
                   by_fold: bool = True) -> pd.Series:
    """Assign each perturbation to a rep_cen decile (within fold if requested)."""
    if by_fold:
        return df.groupby("fold")[score_col].transform(
            lambda s: pd.qcut(s.rank(method="first"), n, labels=False)
        )
    return pd.qcut(df[score_col].rank(method="first"), n, labels=False)


def decile_ranking(df: pd.DataFrame, score_suffix: str = "_B_cen",
                   n_deciles: int = 10) -> pd.DataFrame:
    """Mean independent-half score per model per agreement decile.

    Returns a long DataFrame: decile, model, mean, sem, n.
    """
    df = df.copy()
    df["decile"] = assign_deciles(df, n=n_deciles)
    model_cols = [c for c in df.columns if c.endswith(score_suffix)]
    rows = []
    for d, sub in df.groupby("decile"):
        for c in model_cols:
            v = sub[c].dropna()
            rows.append(
                {
                    "decile": int(d),
                    "model": c[: -len(score_suffix)],
                    "mean": v.mean(),
                    "sem": v.std() / np.sqrt(max(len(v), 1)),
                    "n": len(v),
                }
            )
    return pd.DataFrame(rows)


def bootstrap_ci(df: pd.DataFrame, col_a: str, col_b: str, n_boot: int = 10000,
                 seed: int = 0) -> dict:
    """Paired perturbation bootstrap for the mean difference col_a - col_b."""
    rng = np.random.RandomState(seed)
    d = (df[col_a] - df[col_b]).dropna().to_numpy()
    idx = rng.randint(0, len(d), size=(n_boot, len(d)))
    boot = d[idx].mean(axis=1)
    return {
        "diff": d.mean(),
        "ci_low": np.percentile(boot, 2.5),
        "ci_high": np.percentile(boot, 97.5),
        "n": len(d),
    }

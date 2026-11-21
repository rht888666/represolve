"""Interactive terminal walkthrough for the represolve iMeta video tutorial.

Runs in ~20 seconds on synthetic replicate halves with known ground truth,
demonstrating the exact five-step API used across the benchmark:
  1. Inspect frozen datasets and build/load replicate halves (A | B)
  2. Compute prediction-time diagnostics (obs_agreement, d_cos10) — no hidden genes
  3. Route the screen with Table 3 (`represolve.route`)
  4. Run fold-frozen masked-gene completion (`masked_completion_benchmark`)
  5. Separate same-half inflation from independent-half reproducible signal

Usage:
    python demo/run_video_demo.py            # fast mode (no pauses)
    python demo/run_video_demo.py --present  # paced pauses for live screen recording
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import represolve as rv  # noqa: E402


def step(title: str, pause: float) -> None:
    print(f"\n{'=' * 68}\n  {title}\n{'=' * 68}", flush=True)
    if pause > 0:
        time.sleep(pause)


def make_regime_halves(n_perts: int, n_genes: int, noise_sd: float, seed: int):
    """Synthetic replicate halves with low-rank biological signal + half-specific noise."""
    rng = np.random.default_rng(seed)
    factors = rng.normal(0, 1.0, (n_perts, 8))
    loadings = rng.normal(0, 0.6, (8, n_genes))
    truth = factors @ loadings
    half_a = truth + rng.normal(0, noise_sd, truth.shape)
    half_b = truth + rng.normal(0, noise_sd, truth.shape)
    perts = np.array([f"PERT_{i:03d}" for i in range(n_perts)])
    genes = np.array([f"GENE_{j:04d}" for j in range(n_genes)])
    return half_a, half_b, perts, genes


def inspect_regime(label: str, ha: np.ndarray, hb: np.ndarray,
                   perts: np.ndarray, genes: np.ndarray, pause: float) -> None:
    pooled = 0.5 * (ha + hb)
    mu = pooled.mean(axis=0)
    obs_scores = []
    for i in range(len(perts)):
        m = np.random.RandomState(1000 + i).rand(len(genes)) < 0.1
        obs_scores.append(rv.observed_agreement(ha[i, m], hb[i, m], mu[m]))
    obs_med = float(np.median(obs_scores))
    d_cos10 = rv.cosine_knn_density(pooled, k=10)
    decision = rv.route(obs_med, d_cos10, n_perts=len(perts))
    print(f"  [{label}]")
    print(f"    Prediction-time obs_agreement median : {obs_med:+.3f}")
    print(f"    Perturbation density d_cos10         : {d_cos10:.3f}")
    print(f"    Table 3 pre-scoring route            : {decision}")
    if pause > 0:
        time.sleep(pause)
    df = rv.masked_completion_benchmark(
        ha, hb, perts, genes, n_folds=3, n_genes_selected=300,
    )
    summary = rv.summarize(df)
    print(summary.to_string(index=False, float_format=lambda x: f"{x:+.3f}"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--present", action="store_true",
                    help="Add short pauses between sections for OBS/PowerPoint recording")
    args = ap.parse_args()
    p = 1.5 if args.present else 0.0

    step("Step 1 · represolve — Replicate-resolved perturbation benchmarking", p)
    print(f"  Package version     : {rv.__version__}")
    print(f"  Frozen benchmark DB : {len(rv.list_datasets())} archived contexts with SHA-256")
    for k in rv.list_datasets()[:5]:
        print(f"    • {k}")
    print("    • ... (9 ranking contexts + held-out external validation contexts)")

    step("Step 2 · High-reliability regime (low measurement noise)", p)
    ha_hi, hb_hi, perts, genes = make_regime_halves(120, 300, noise_sd=0.4, seed=1)
    inspect_regime("High-agreement screen", ha_hi, hb_hi, perts, genes, p)

    step("Step 3 · Low-reliability regime (high measurement noise)", p)
    ha_lo, hb_lo, perts, genes = make_regime_halves(120, 300, noise_sd=2.5, seed=2)
    inspect_regime("Low-agreement screen", ha_lo, hb_lo, perts, genes, p)

    step("Step 4 · Takeaway — No universally optimal model", p)
    print("  • Same-half scores (A_cen) inflate every estimator by shared-half noise.")
    print("  • Independent-half scoring (B_cen) reveals regime-dependent model rankings.")
    print("  • obs_agreement + d_cos10 declare the regime BEFORE hidden-gene scoring.")


if __name__ == "__main__":
    main()

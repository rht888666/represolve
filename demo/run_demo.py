"""Minimal runnable demo (<1 min): synthetic two-replicate perturbation data
through the full fold-frozen benchmark, then decile ranking and bootstrap.

Synthetic design: perturbations combine a small set of shared response
programs plus an idiosyncratic component, so hidden genes are predictable
from observed genes — mimicking real Perturb-seq structure.
"""
import numpy as np

from perturbeval import (
    masked_completion_benchmark,
    summarize,
    decile_ranking,
    bootstrap_ci,
)

rng = np.random.RandomState(0)
n_pert, n_genes, n_programs = 80, 400, 8

# shared response programs + per-perturbation loadings
programs = rng.randn(n_programs, n_genes)
weights = rng.randn(n_pert, n_programs) * 0.6
latent = weights @ programs + rng.randn(n_pert, n_genes) * 0.3

# two replicate measurement halves: shared latent effect + independent noise
half_a = latent + rng.randn(n_pert, n_genes) * 0.5
half_b = latent + rng.randn(n_pert, n_genes) * 0.5

perts = np.array([f"pert{i}" for i in range(n_pert)])
genes = np.array([f"g{j}" for j in range(n_genes)])

df = masked_completion_benchmark(
    half_a, half_b, perts, genes, n_folds=4, n_genes_selected=300
)
print("\n=== Model summary (sorted by independent-half score) ===")
print(summarize(df).to_string(index=False))
print(f"\nMean replicate agreement (rep_cen): {df['rep_cen'].mean():.3f}")

print("\n=== Decile ranking (independent half) ===")
print(
    decile_ranking(df, n_deciles=5)
    .pivot(index="decile", columns="model", values="mean")
    .to_string()
)

print("\n=== nn10 vs svd20 bootstrap CI ===")
print(bootstrap_ci(df, "nn10_B_cen", "svd20_B_cen", n_boot=2000))

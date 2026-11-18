"""Synthetic smoke test: construct two correlated replicate halves of
perturbation profiles and check the benchmark runs end-to-end."""
import numpy as np

from perturbeval import (
    masked_completion_benchmark,
    summarize,
    decile_ranking,
    bootstrap_ci,
    bh_qvalues,
    permutation_null,
    type1_check,
)


def make_toy(n_pert=60, n_genes=300, seed=0):
    rng = np.random.RandomState(seed)
    latent = rng.randn(n_pert, n_genes) * 0.3
    a = latent + rng.randn(n_pert, n_genes) * 0.7
    b = latent + rng.randn(n_pert, n_genes) * 0.7
    perts = np.array([f"pert{i}" for i in range(n_pert)])
    genes = np.array([f"g{j}" for j in range(n_genes)])
    return a, b, perts, genes


def test_end_to_end():
    a, b, perts, genes = make_toy()
    df = masked_completion_benchmark(
        a, b, perts, genes, n_folds=3, n_genes_selected=200
    )
    assert len(df) == len(perts)
    s = summarize(df)
    assert set(s.columns) == {"model", "same_half", "independent_half", "inflation", "n"}
    assert "rep_cen" in df.columns


def test_decile_and_bootstrap():
    a, b, perts, genes = make_toy()
    df = masked_completion_benchmark(a, b, perts, genes, n_folds=3, n_genes_selected=200)
    dr = decile_ranking(df, n_deciles=5)
    assert {"decile", "model", "mean", "sem", "n"} <= set(dr.columns)
    ci = bootstrap_ci(df, "nn10_B_cen", "svd20_B_cen", n_boot=200)
    assert ci["ci_low"] <= ci["diff"] <= ci["ci_high"]


def test_nulls():
    rng = np.random.RandomState(0)
    obs = rng.randn(50)
    perm = rng.randn(1000, 50)
    p = permutation_null(obs, perm)
    assert p.shape == (50,)
    q = bh_qvalues(p)
    assert q.shape == (50,)

    def sim(r):
        return r.uniform(0, 1, size=30)

    res = type1_check(sim, n_sims=20)
    assert 0 <= res["any_hit_rate"] <= 1


def test_load_halves(tmp_path):
    import numpy as np
    from perturbeval.data import load_halves, list_datasets

    a, b, perts, genes = make_toy(n_pert=10, n_genes=20)
    path = tmp_path / "toy_halves.npz"
    np.savez(path, S=np.stack([a, b]).astype(np.float32),
             N=np.array([[10] * 10, [10] * 10]),
             perts=perts, genes=genes)
    ha, hb, hp, hg = load_halves(path)
    assert ha.shape == a.shape and np.allclose(ha, a)
    assert list(hp) == list(perts)
    names = [d.name for d in list_datasets()]
    assert "sciplex3" in names and "k562_gwps" in names

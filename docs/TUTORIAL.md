# perturbeval tutorial — reproduce the paper's core result

This walkthrough reproduces the K562-essential masked-completion
benchmark on the frozen replicate halves (~5 min on a laptop) and
shows how the five context results assemble into the cross-dataset
ranking figure.

## 0. Install and get the data

```bash
pip install perturbeval
# download frozen halves (Zenodo DOI: TBD) into ./data/
```

`data/` should contain the five `.npz` files listed in
`src/perturbeval/DATA_MANIFEST.csv` (verify checksums first):

```bash
sha256sum data/*.npz    # compare against DATA_MANIFEST.csv
```

## 1. Run one context

```python
import perturbeval as pv

ha, hb, perts, genes = pv.load_halves("data/replogle_k562_essential_halves.npz")
df = pv.masked_completion_benchmark(
    ha, hb, perts, genes,
    n_folds=5, n_genes_selected=5000, mask_fraction=0.1,
)
print(pv.summarize(df))
```

Expected output (K562 essential): independent-half centered scores
approximately `nn10 0.431 > nn50 0.420 > svd20 0.404 > moppca5 0.402 >
lw 0.355` with median `rep_cen ≈ 0.279`.

## 2. Regime-stratified ranking

```python
dr = pv.decile_ranking(df)          # model ranking per rep_cen decile
ci = pv.bootstrap_ci(df, "nn10_B_cen", "lw_B_cen")
print(dr.head(15)); print(ci)
```

In this low-reliability context the nearest-centroid retrieval leads and
Ledoit-Wolf shrinkage loses — the reverse of what the Norman and RPE1
contexts produce. Compare `data/norman_singles_gemhalves.npz`
(median rep_cen 0.953) where `lw` wins at every decile.

## 3. Chemical-modality context

```python
ha, hb, perts, genes = pv.load_halves("data/sciplex3_rep_halves.npz")
df = pv.masked_completion_benchmark(ha, hb, perts, genes)
print(pv.summarize(df))
```

sci-Plex3 halves are the dataset's own rep1/rep2 biological replicates
— no synthetic partition. Median replicate agreement is ~0.99 and
`nn10` leads anyway because drug×dose×cell-line units form a dense
grid where neighbour interpolation is cheap.

## 4. Cross-dataset figure

```bash
python figures/fig_crossdataset_ranking.py \
    ../results/j51_crossdataset_model_ranking.csv fig_ranking.png
```

## 5. Reproducing every manuscript number

`paper_map/FIGURE_TO_CODE.md` maps each figure/table to the exact
script (`j22`/`j47`/`j49` series in the analysis repo) and output CSV.
The masked-completion runs in this tutorial use the same fold seeds,
gene-selection and mask seeds as those scripts.

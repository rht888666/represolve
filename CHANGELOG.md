# Changelog

## 1.0.4 — 2026-11-24

### Fixed
- `route()` (Python and R) now applies the Table 3 large-panel requirement on the dense-grid branch: when `n_perts` is supplied and is below 1,000, a dense high-agreement screen is routed to `lw_shrinkage` instead of `nn10_dense_grid`. Routes for the nine benchmark and four external contexts are unchanged.
- `DATA_MANIFEST.csv` is now shipped inside the wheel (`package-data`).

### Changed
- `environment.yml` adds `matplotlib` and `jupyter` (needed for the notebooks and figures).
- Package description, README and maintainer e-mail aligned with the submitted manuscript title.

## 1.0.3 — 2026-11-21

### Added
- Native R package (`represolve` for R >= 4.1.0, base `stats` and `utils` only, installable via `remotes::install_github("rht888666/represolve")`)
- Soft Regime Router (`soft_route_weights`, `route_soft`, `route_soft_model`) and continuous adaptive neighbourhood sizing (`adaptive_nn_k`, `blend_nn_scores`) in both Python and R
- GitHub Actions CI testing both Python and R packages (`R CMD build` and `R CMD check`)

### Changed
- `load_halves` auto-detects raw-sum `.npz` archives, applies cell-count filtering (`min_cells = 25`) and optional `non-targeting` exclusion, and divides by `N` in-place
- `LinearCompletion` and `NearestCentroid` accelerated via exact analytical Gram/Woodbury push-through identities with zero numerical drift

## 1.0.0 — 2026-11-21

### Added
- Prediction-time regime helpers: `observed_agreement`, `cosine_knn_density`, `route`
- `MaskedMLPCompletion` as the canonical nonlinear baseline name
- Frozen dataset registry for 9 benchmark + 4 external contexts (`list_datasets` / `load_halves`)
- CI workflow running smoke tests + `demo/run_demo.py`
- Publication-oriented README, CITATION.cff, MIT LICENSE

### Changed
- Package version bumped from 0.1.0 → 1.0.0 for first public release
- README reframed around masked-gene completion (not generic perturbation prediction)

### Notes
- Large frozen `.npz` arrays ship via Zenodo (DOI filled after deposit)
- `DAECompletion` retained as compatibility alias of `MaskedMLPCompletion`

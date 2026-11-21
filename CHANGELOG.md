# Changelog

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

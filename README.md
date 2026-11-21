# represolve

[![PyPI](https://img.shields.io/pypi/v/represolve.svg)](https://pypi.org/project/represolve/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24409221.svg)](https://doi.org/10.5281/zenodo.24409221)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Replicate-resolved masked-gene completion benchmarking.**

Installable companion to the manuscript:

> *represolve: Replicate-resolved masked-completion benchmarking reveals
> regime-dependent model rankings and calibrates keloid fibroblast
> candidates*

`represolve` scores masked-gene completion against an **independent
replicate measurement half** (held-out guides, cell pools, batches, or
biological replicates), not against the same partition used as model input.
It ships fold-frozen splits, SHA-256 manifests, six baselines, prediction-time
regime routing (`obs_agreement` + `d_cos10`), and a keloid case-study toolkit.

## Install

**Python (PyPI):**

```bash
pip install represolve          # https://pypi.org/project/represolve/
# or from a local checkout:
pip install -e .
# or conda:
conda env create -f environment.yml && conda activate represolve
```

Requires Python ≥ 3.10. Core dependencies: NumPy, pandas, SciPy, scikit-learn.

**R (GitHub):**

```r
# install.packages("remotes")
remotes::install_github("rht888666/represolve")
```

Requires R ≥ 4.1.0 (base `stats` and `utils` only; no compilation required). A pre-built source tarball (`represolve_1.0.0.tar.gz`) is also archived in `dist/` and on Zenodo.

Docker (optional):

```bash
docker build -t represolve .
docker run --rm represolve
```

## 30-second demo

```bash
python demo/run_demo.py
```

```python
import numpy as np
from represolve import (
    masked_completion_benchmark, summarize, decile_ranking,
    observed_agreement, cosine_knn_density, route,
)

# half_a, half_b: (n_perturbations, n_genes) mean profiles of the SAME
# perturbations measured in two independent partitions.
df = masked_completion_benchmark(half_a, half_b, pert_names, genes)
print(summarize(df))          # same-half vs independent-half scores
print(decile_ranking(df))     # ranking by replicate-agreement decile

# Prediction-time regime (before hidden-gene scoring)
obs = observed_agreement(half_a, half_b, mask_frac=0.1, seed=0)
dcos = cosine_knn_density(half_a)
print(route(obs_agreement=float(np.nanmedian(obs)), d_cos10=float(dcos)))
```

In R:

```r
library(represolve)

df <- masked_completion_benchmark(half_a, half_b, pert_names, genes)
summarize_benchmark(df)
decile_ranking(df)

route(obs_agreement = 0.29, cos_density = 0.65)
soft_route_weights(obs_agreement = 0.29, cos_density = 0.65)
```

## Core API

| Symbol | Purpose |
|---|---|
| `masked_completion_benchmark` | Fold-frozen masked-gene completion vs two replicate halves |
| `summarize` / `decile_ranking` / `bootstrap_ci` | Tables and paired bootstrap CIs |
| `replicate_agreement` / `observed_agreement` | Hidden-half vs prediction-time agreement |
| `cosine_knn_density` / `route` | Scale-invariant density + first-generation regime router |
| `SVDCompletion`, `LinearCompletion`, `MoPPCACompletion`, `NearestCentroid`, `MaskedMLPCompletion` | Baselines (`DAECompletion` is kept as an alias of `MaskedMLPCompletion`) |
| `list_datasets` / `load_halves` | Frozen 9-benchmark + 4-external contexts (arrays archived at Zenodo) |
| `permutation_null`, `bh_qvalues`, `type1_check` | Screen-null calibration utilities |

## Frozen benchmark inputs

| Dataset | Split axis | Modality | Perturbations |
|---|---|---|---|
| `k562_essential` | gem-group | CRISPRi | 1,777 |
| `k562_gwps` | gem-group | CRISPRi GWPS | 9,211 |
| `norman` | gem-group | CRISPRa | 105 |
| `rpe1` / `rpe1_e150` | batch / gem-group | CRISPRi | 1,487 / 134 |
| `sciplex3` | biological replicate | chemical | 2,224 |
| `frangieh` | guide split | Cas9 KO | 233 |
| `tian2021i` | batch | CRISPRi | 183 |
| `tahoe_plate3` | sample split | chemical | 145 |
| External ×4 | gem / wells / batch | HepG2, Jurkat, sci-Plex4, Tian CRISPRa | 96–1,704 |

Checksums live in `src/represolve/DATA_MANIFEST.csv`. Large `.npz` arrays are
archived at Zenodo ([DOI: 10.5281/zenodo.24409221](https://doi.org/10.5281/zenodo.24409221)).

## Add your own predictor

```python
from represolve.models import BaseCompletion  # or duck-type fit/predict

class MyModel:
    def fit(self, X_train, genes): ...
    def predict(self, x_obs, observed_mask): ...
```

Plug into `masked_completion_benchmark(..., models={"mine": MyModel()})`.

## Tests & CI

```bash
python -m pytest tests/ -q   # or: python -c "import runpy; runpy.run_path('tests/test_smoke.py')"
python demo/run_demo.py
R CMD build . && R CMD check --no-manual represolve_1.0.0.tar.gz
```

GitHub Actions workflow: `.github/workflows/ci.yml`.

## Scope

- Replicate halves are **measurement partitions**, not necessarily donor biology.
- `rep_cen` is a descriptive reliability statistic, not a theoretical ceiling.
- This package is an evaluation protocol + baselines, not a new architecture.
- The keloid analysis is a **calibrated case study**, not a wet-lab-validated target claim.

## Citation

See `CITATION.cff`. After Zenodo deposition, update the DOI there and in the manuscript Data Availability statement.

## License

MIT — see `LICENSE`.

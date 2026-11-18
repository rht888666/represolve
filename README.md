# perturbeval

**Replicate-resolved evaluation for perturbation-response models.**

Companion analysis framework to the manuscript:

> *Replicate-aware evaluation exposes context-dependent perturbation
> predictability and re-prioritizes keloid fibroblast programs*
> (manuscript in preparation; preprint and archived release to follow).

`perturbeval` implements the evaluation protocol used in the paper: it scores
perturbation predictions against an **independent replicate measurement
half** (held-out guides or cell pools), not against the same measurements used
as model inputs, and stratifies model ranking by the detectable-response
regime. This separates reproducible target-specific prediction from recovery
of shared variation.

## Why

Same-sample evaluation conflates a model's ability to recover systematic or
shared variation with its ability to predict reproducible perturbation
effects. Scoring against a held-out replicate half provides an empirical
reliability axis (`rep_cen`) on which:

- model inflation (same-half minus independent-half score) is concentrated in
  weakly responsive perturbations, and
- **model rankings reverse across the detectable-response regime** — no
  universally optimal model exists (essential-scale screens favour a
  k-nearest-perturbation-centroid baseline; genome-wide screens reorder the
  field entirely).

## Install

```bash
conda env create -f environment.yml && conda activate perturbeval
# or:
pip install -e .
```

## Quick start (<1 min demo)

```bash
python demo/run_demo.py
```

```python
import numpy as np
from perturbeval import masked_completion_benchmark, summarize, decile_ranking

# half_a, half_b: (n_perturbations, n_genes) mean profiles of the SAME
# perturbations measured in two independent partitions.
df = masked_completion_benchmark(half_a, half_b, pert_names, genes)
print(summarize(df))          # per-model same-half vs independent-half scores
print(decile_ranking(df))     # model ranking by replicate-agreement decile
```

## Core API

| Function | Purpose |
|---|---|
| `masked_completion_benchmark` | Fold-frozen masked-gene completion benchmark scored against two replicate halves |
| `summarize` | Per-model same-half / independent-half / inflation table |
| `decile_ranking` | Model ranking stratified by replicate-agreement (`rep_cen`) decile |
| `bootstrap_ci` | Paired perturbation-bootstrap CI for model differences |
| `replicate_agreement` | `rep_cen` — centered correlation between replicate halves |
| `permutation_null`, `bh_qvalues`, `type1_check` | Empirical nulls, FDR, synthetic type-I calibration |
| `SVDCompletion`, `LinearCompletion`, `MoPPCACompletion`, `NearestCentroid`, `DAECompletion` | Baseline predictors (fit-once/predict-per-perturbation interface); `DAECompletion` is the nonlinear denoising-autoencoder baseline |
| `list_datasets`, `load_halves`, `DATASETS` | Frozen replicate-half benchmark inputs (five contexts) with provenance |

### Frozen benchmark inputs

`perturbeval.data` ships the metadata and loader for the five frozen
replicate-half datasets used in the paper (`DATA_MANIFEST.csv` holds
SHA-256 checksums):

| Dataset | Split axis | Modality | Perturbations |
|---|---|---|---|
| `k562_essential` | gem-group cell pools | CRISPRi | 1,777 |
| `k562_gwps` | gem-group cell pools | CRISPRi genome-scale | 9,211 |
| `norman` | gem-group | CRISPRa (single genes) | 105 |
| `rpe1` | batch odd/even | CRISPRi | 1,487 |
| `sciplex3` | **biological replicate (rep1/rep2)** | chemical | 2,224 |

The `.npz` files themselves are archived at Zenodo (DOI TBD); raw sources
are the Replogle et al. 2022 screens and the scPerturb collection
(Zenodo record 7041849).

Note on metric equivalence: our centered score subtracts the **training-perturbation mean profile** before correlating — i.e. it is a delta-vs-perturbed-centroid metric in the sense of Systema (Nat Biotechnol 2025), applied per replicate partition. What differs is the scoring axis (held-out measurement replicate) and the regime-stratified ranking, not the reference choice.

## Scope and terminology

- Replicate halves are **measurement partitions** (guide halves, cell-pool
  halves), not donor-level biological replication — the same terminology
  caution the manuscript makes.
- `rep_cen` is a descriptive reliability statistic, not a theoretical
  prediction ceiling.
- This package contains the evaluation protocol and baselines; it is not a
  new predictive architecture.

## Relationship to prior work

Systema (Nat Biotechnol 2025) changes the reference point; AmbiModBench
(arXiv 2025/2026) audits absolute metrics with matched references and
permutation screens. `perturbeval` complements both by changing the
**scoring partition** — held-out biological-measurement replicates — and by
stratifying model ranking on the resulting reliability axis.

## Paper reproduction

Analysis scripts that regenerate the manuscript's results live in the parent
repository (`../j*.py`); the mapping from each figure/result to scripts and
outputs is in `paper_map/`. Processed benchmark inputs will be deposited at
Zenodo with a DOI (link TBD).

## License

MIT (see `LICENSE`). Citation information in `CITATION.cff`.

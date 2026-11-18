"""perturbeval: replicate-resolved evaluation for perturbation-response models.

Companion framework to the manuscript "Replicate-aware evaluation exposes
context-dependent perturbation predictability and re-prioritizes keloid
fibroblast programs". Provides the fold-frozen masked-completion benchmark,
nearest-perturbation-centroid and shrinkage baselines, regime-stratified
ranking, and the permutation/type-I calibration utilities used in the paper.
"""
from .metrics import centered_score, neighbour_similarity, pearson, replicate_agreement
from .models import (
    LinearCompletion,
    MoPPCACompletion,
    NearestCentroid,
    SVDCompletion,
    default_models,
)
from .benchmark import masked_completion_benchmark, summarize
from .regimes import assign_deciles, bootstrap_ci, decile_ranking
from .null import bh_qvalues, permutation_null, type1_check

__version__ = "0.1.0"
__all__ = [
    "pearson",
    "centered_score",
    "replicate_agreement",
    "neighbour_similarity",
    "SVDCompletion",
    "LinearCompletion",
    "MoPPCACompletion",
    "NearestCentroid",
    "default_models",
    "masked_completion_benchmark",
    "summarize",
    "assign_deciles",
    "decile_ranking",
    "bootstrap_ci",
    "permutation_null",
    "bh_qvalues",
    "type1_check",
]

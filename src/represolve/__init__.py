"""represolve: replicate-resolved evaluation for perturbation-response models.

Companion framework to the manuscript "represolve: Replicate-resolved
masked-completion benchmarking reveals regime-dependent model rankings and
calibrates keloid fibroblast candidates". Provides the fold-frozen masked-completion benchmark,
nearest-perturbation-centroid and shrinkage baselines, regime-stratified
ranking, and the permutation/type-I calibration utilities used in the paper.
"""
from .metrics import centered_score, neighbour_similarity, pearson, replicate_agreement
from .models import (
    DAECompletion,
    LinearCompletion,
    MaskedMLPCompletion,
    MoPPCACompletion,
    NearestCentroid,
    SVDCompletion,
    default_models,
)
from .density import (
    adaptive_nn_k,
    blend_nn_scores,
    cosine_knn_density,
    observed_agreement,
    route,
    route_soft,
    route_soft_model,
    soft_route_weights,
)
from .benchmark import masked_completion_benchmark, summarize
from .regimes import assign_deciles, bootstrap_ci, decile_ranking
from .null import bh_qvalues, permutation_null, type1_check
from .data import DATASETS, HalfDataset, list_datasets, load_halves

__version__ = "1.0.3"
__all__ = [
    "pearson",
    "centered_score",
    "replicate_agreement",
    "neighbour_similarity",
    "SVDCompletion",
    "DAECompletion",
    "MaskedMLPCompletion",
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
    "DATASETS",
    "HalfDataset",
    "list_datasets",
    "load_halves",
    "cosine_knn_density",
    "observed_agreement",
    "route",
    "soft_route_weights",
    "route_soft",
    "route_soft_model",
    "adaptive_nn_k",
    "blend_nn_scores",
]

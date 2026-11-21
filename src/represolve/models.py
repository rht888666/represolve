"""Masked-completion predictors.

A predictor is fit once per training fold on the training perturbation
profiles (perturbations x genes). Given the observed subset of genes for a
held-out perturbation, it predicts the hidden genes. All predictors return a
full-length prediction vector; the caller slices hidden positions for scoring.
"""
from __future__ import annotations

import numpy as np
from sklearn.covariance import LedoitWolf
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture

from .metrics import neighbour_similarity


class BaseCompletion:
    """Fit-once / predict-per-perturbation interface."""

    name = "base"

    def fit(self, train_profiles: np.ndarray) -> "BaseCompletion":
        raise NotImplementedError

    def predict(self, observed: np.ndarray, obs_mask: np.ndarray) -> np.ndarray:
        """observed: full-length profile (only obs_mask positions used).
        obs_mask: boolean array, True = observed gene."""
        raise NotImplementedError


class SVDCompletion(BaseCompletion):
    """Low-rank completion via PCA loadings."""

    def __init__(self, n_components: int = 20, random_state: int = 0):
        self.n_components = n_components
        self.random_state = random_state
        self.name = f"svd{n_components}"

    def fit(self, T: np.ndarray) -> "SVDCompletion":
        self.mu_ = T.mean(axis=0)
        k = min(self.n_components, min(T.shape) - 1)
        self.W_ = PCA(k, random_state=self.random_state).fit(T).components_.T
        return self

    def predict(self, observed: np.ndarray, obs_mask: np.ndarray) -> np.ndarray:
        coef = np.linalg.lstsq(
            self.W_[obs_mask], observed[obs_mask] - self.mu_[obs_mask], rcond=None
        )[0]
        return self.mu_ + self.W_ @ coef


class LinearCompletion(BaseCompletion):
    """Conditional Gaussian completion with a Ledoit-Wolf covariance."""

    def __init__(self):
        self.name = "lw"

    def fit(self, T: np.ndarray) -> "LinearCompletion":
        self.mu_ = T.mean(axis=0)
        self.cov_ = LedoitWolf().fit(T).covariance_
        return self

    def predict(self, observed: np.ndarray, obs_mask: np.ndarray) -> np.ndarray:
        hid = ~obs_mask
        pred = self.mu_.copy()
        pred[hid] = self.mu_[hid] + self.cov_[np.ix_(hid, obs_mask)] @ np.linalg.solve(
            self.cov_[np.ix_(obs_mask, obs_mask)],
            observed[obs_mask] - self.mu_[obs_mask],
        )
        return pred


class MoPPCACompletion(BaseCompletion):
    """Mixture of probabilistic PCA experts (moPPCA)."""

    def __init__(self, n_components: int = 5, latent_dim: int = 50, expert_dim: int = 10,
                 random_state: int = 0):
        self.n_components = n_components
        self.latent_dim = latent_dim
        self.expert_dim = expert_dim
        self.random_state = random_state
        self.name = f"moppca{n_components}"

    def fit(self, T: np.ndarray) -> "MoPPCACompletion":
        self.mu_ = T.mean(axis=0)
        k = min(self.latent_dim, min(T.shape) - 1)
        self.W50_ = PCA(k, random_state=self.random_state).fit(T).components_.T
        latent = T @ self.W50_
        self.mix_ = GaussianMixture(
            self.n_components, covariance_type="diag", reg_covar=1e-4,
            random_state=self.random_state,
        ).fit(latent)
        labels = self.mix_.predict(latent)
        self.experts_ = {
            c: PCA(self.expert_dim, random_state=self.random_state).fit(T[labels == c])
            for c in range(self.n_components)
            if (labels == c).sum() >= 15
        }
        return self

    def predict(self, observed: np.ndarray, obs_mask: np.ndarray) -> np.ndarray:
        query = np.linalg.lstsq(
            self.W50_[obs_mask], observed[obs_mask] - self.mu_[obs_mask], rcond=None
        )[0]
        weights = self.mix_.predict_proba(query[None])[0]
        weighted = np.zeros_like(self.mu_)
        total = 0.0
        for c, expert in self.experts_.items():
            if weights[c] < 1e-6:
                continue
            W = expert.components_.T
            estimate = expert.mean_ + W @ np.linalg.lstsq(
                W[obs_mask], observed[obs_mask] - expert.mean_[obs_mask], rcond=None
            )[0]
            weighted += weights[c] * estimate
            total += weights[c]
        return weighted / total if total else self.mu_.copy()


class MaskedMLPCompletion(BaseCompletion):
    """Nonlinear baseline: masked-input multilayer perceptron (masked-MLP).

    Trained on augmented training perturbations: each sample is masked with
    n_masks random masks, hidden positions are filled with the training mean,
    and a mask indicator channel is appended. A 256-unit hidden layer, no
    bottleneck — the deep-learning analogue of the completion task.
    """

    def __init__(self, hidden_layer_sizes=(256,), n_masks: int = 6,
                 mask_fraction: float = 0.1, max_iter: int = 40,
                 random_state: int = 0):
        self.hidden_layer_sizes = hidden_layer_sizes
        self.n_masks = n_masks
        self.mask_fraction = mask_fraction
        self.max_iter = max_iter
        self.random_state = random_state
        self.name = "masked_mlp"

    def fit(self, T: np.ndarray) -> "MaskedMLPCompletion":
        from sklearn.neural_network import MLPRegressor

        rng = np.random.RandomState(self.random_state)
        self.mu_ = T.mean(axis=0)
        G = T.shape[1]
        reps = np.repeat(T, self.n_masks, axis=0)
        masks = rng.rand(reps.shape[0], G) < self.mask_fraction
        # hidden positions -> training mean
        masked = reps.copy()
        masked[~masks] = np.tile(self.mu_, (reps.shape[0], 1))[~masks]
        X_in = np.concatenate([masked, masks.astype(np.float32)], axis=1).astype(np.float32)
        y = reps.astype(np.float32)
        self.net_ = MLPRegressor(
            hidden_layer_sizes=self.hidden_layer_sizes,
            activation="relu",
            max_iter=self.max_iter,
            early_stopping=True,
            random_state=self.random_state,
        )
        self.net_.fit(X_in, y)
        return self

    def predict(self, observed: np.ndarray, obs_mask: np.ndarray) -> np.ndarray:
        masked = np.where(obs_mask, observed, self.mu_)
        x = np.concatenate([masked, obs_mask.astype(np.float32)])[None].astype(np.float32)
        return self.net_.predict(x)[0]


class NearestCentroid(BaseCompletion):
    """Nearest-perturbation-centroid baseline.

    For each held-out perturbation, find the k training perturbations whose
    observed-gene effects are most (centered-)correlated with the test
    observed effect, then predict hidden genes as their centroid.
    """

    def __init__(self, k: int = 10):
        self.k = k
        self.name = f"nn{k}"

    def fit(self, T: np.ndarray) -> "NearestCentroid":
        self.T_ = T
        self.mu_ = T.mean(axis=0)
        return self

    def predict(self, observed: np.ndarray, obs_mask: np.ndarray) -> np.ndarray:
        obs_pos = np.flatnonzero(obs_mask)
        sims = neighbour_similarity(self.T_, self.mu_, obs_pos, observed[obs_mask])
        nn = np.argsort(-sims, kind="stable")[: self.k]
        pred = self.mu_.copy()
        pred[~obs_mask] = self.T_[nn][:, ~obs_mask].mean(axis=0)
        pred[obs_mask] = observed[obs_mask]
        return pred


def default_models(nn_ks=(10, 50)) -> list:
    """The model set used in the paper's fold-frozen benchmarks."""
    models = [
        SVDCompletion(20),
        LinearCompletion(),
        MoPPCACompletion(5),
    ]
    models += [NearestCentroid(k) for k in nn_ks]
    return models


# Backward-compatible alias for pre-rename notebooks/scripts.
DAECompletion = MaskedMLPCompletion

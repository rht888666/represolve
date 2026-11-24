.ROUTE_SOFT_PROTOTYPES <- list(
    flag_near_floor = c(-0.02, 0.90),
    centroid_nn50 = c(0.07, 0.85),
    centroid_nn10 = c(0.35, 0.55),
    nn10_dense_grid = c(0.95, 0.01),
    lw_shrinkage = c(0.90, 0.25)
)

.ROUTE_SOFT_TO_MODEL <- c(
    flag_near_floor = "flag",
    centroid_nn50 = "nn50",
    centroid_nn10 = "nn10",
    nn10_dense_grid = "nn10",
    lw_shrinkage = "lw"
)

#' Observed-gene agreement between replicate halves
#'
#' @param half_a_obs Numeric vector (or matrix) on observed genes in half A.
#' @param half_b_obs Numeric vector (or matrix) on observed genes in half B.
#' @param train_mean Numeric vector of training means on observed genes.
#' @return Scalar numeric centered Pearson correlation.
observed_agreement <- function(half_a_obs, half_b_obs, train_mean) {
    pearson(as.numeric(half_a_obs) - as.numeric(train_mean),
            as.numeric(half_b_obs) - as.numeric(train_mean))
}

#' Scale-invariant cosine k-NN density diagnostic (d_cos10)
#'
#' Computes the median cosine distance from each training perturbation to its
#' k-th nearest training neighbour on L2-normalized profiles.
#'
#' @param train_profiles Numeric matrix of shape (n_perturbations, n_genes).
#' @param k Integer neighbour rank (default 10).
#' @param center Logical; whether to column-center before L2 normalization.
#' @return Scalar median cosine distance in \code{[0, 2]}.
cosine_knn_density <- function(train_profiles, k = 10L, center = TRUE) {
    X <- as.matrix(train_profiles)
    storage.mode(X) <- "double"
    if (center) {
        X <- sweep(X, 2L, colMeans(X), "-")
    }
    nrm <- sqrt(rowSums(X^2))
    nrm[nrm == 0] <- 1.0
    U <- X / nrm
    G <- tcrossprod(U)
    diag(G) <- -Inf
    k <- as.integer(min(k, nrow(G) - 1L))
    d_mat <- 1.0 - G
    kth <- vapply(seq_len(nrow(d_mat)), function(i) {
        sort(d_mat[i, ], partial = k)[k]
    }, numeric(1L))
    as.numeric(stats::median(kth))
}

#' First-generation hard regime router (Table 3)
#'
#' @param obs_agreement Scalar prediction-time observed-gene agreement (or NULL).
#' @param cos_density Scalar cosine k-NN density \code{d_cos10} (optional).
#' @param n_perts Optional perturbation count; the dense-grid branch requires
#'   a large panel (\code{n_perts >= 1000}) when supplied.
#' @return Character routing decision.
route <- function(obs_agreement = NULL, cos_density = NULL, n_perts = NULL) {
    if (is.null(obs_agreement) || is.na(obs_agreement)) {
        return("inherit_design_priors")
    }
    if (obs_agreement < 0.05) {
        return("flag_near_floor")
    }
    if (obs_agreement < 0.15) {
        return("centroid_nn50")
    }
    if (obs_agreement < 0.8) {
        return("centroid_nn10")
    }
    dense_ok <- !is.null(cos_density) && !is.na(cos_density) && cos_density < 0.05
    if (dense_ok && (is.null(n_perts) || is.na(n_perts) || n_perts >= 1000)) {
        return("nn10_dense_grid")
    }
    "lw_shrinkage"
}

.soft_feature <- function(obs_agreement, cos_density) {
    oa <- max(-0.2, min(1.05, as.numeric(obs_agreement)))
    dc <- max(0.0, min(2.0, as.numeric(cos_density)))
    c(oa, dc)
}

#' RBF Soft Regime Router weights
#'
#' @param obs_agreement Scalar observed-gene agreement.
#' @param cos_density Scalar cosine k-NN density \code{d_cos10}.
#' @param temperature Positive RBF kernel temperature (default 0.08).
#' @param prototypes Optional named list of 2-element regime prototypes.
#' @param scale Numeric length-2 axis scaling vector (default \code{c(1.0, 1.2)}).
#' @return Named numeric vector of regime weights summing to 1.
soft_route_weights <- function(obs_agreement, cos_density,
                               temperature = 0.08,
                               prototypes = NULL,
                               scale = c(1.0, 1.2)) {
    if (temperature <= 0) {
        stop("temperature must be positive")
    }
    protos <- if (is.null(prototypes)) .ROUTE_SOFT_PROTOTYPES else prototypes
    z <- .soft_feature(obs_agreement, cos_density)
    s <- as.numeric(scale)
    keys <- names(protos)
    logits <- vapply(keys, function(k) {
        d <- (z - as.numeric(protos[[k]])) / s
        -0.5 * sum(d^2) / temperature
    }, numeric(1L))
    logits <- logits - max(logits)
    w <- exp(logits)
    w <- w / sum(w)
    stats::setNames(as.numeric(w), keys)
}

#' Soft Regime Router argmax action
#'
#' @inheritParams soft_route_weights
#' @return Character action label.
route_soft <- function(obs_agreement, cos_density,
                       temperature = 0.08,
                       prototypes = NULL) {
    w <- soft_route_weights(obs_agreement, cos_density,
                            temperature = temperature,
                            prototypes = prototypes)
    names(w)[which.max(w)]
}

#' Soft Regime Router estimator label
#'
#' @inheritParams soft_route_weights
#' @return Character estimator name (or \code{"flag"}).
route_soft_model <- function(obs_agreement, cos_density,
                             temperature = 0.08,
                             prototypes = NULL) {
    action <- route_soft(obs_agreement, cos_density,
                         temperature = temperature,
                         prototypes = prototypes)
    unname(.ROUTE_SOFT_TO_MODEL[[action]])
}

#' Continuous neighbourhood size along the mid-reliability band
#'
#' @param obs_agreement Scalar observed-gene agreement.
#' @param k_low Neighbourhood size near floor (default 50).
#' @param k_mid Neighbourhood size at mid-band upper end (default 10).
#' @param lo Lower agreement cut (default 0.05).
#' @param hi Upper agreement cut (default 0.80).
#' @return Scalar numeric neighbourhood size.
adaptive_nn_k <- function(obs_agreement, k_low = 50L, k_mid = 10L,
                          lo = 0.05, hi = 0.80) {
    oa <- as.numeric(obs_agreement)
    if (oa <= lo) {
        return(as.numeric(k_low))
    }
    if (oa >= hi) {
        return(as.numeric(k_mid))
    }
    t <- (oa - lo) / (hi - lo)
    as.numeric(k_low + t * (k_mid - k_low))
}

#' Convex blend of nn10 and nn50 scores guided by adaptive_nn_k
#'
#' @param nn10_score Scalar score for nn10.
#' @param nn50_score Scalar score for nn50.
#' @param obs_agreement Scalar observed-gene agreement.
#' @param lo Lower agreement cut (default 0.05).
#' @param hi Upper agreement cut (default 0.80).
#' @return Scalar blended score.
blend_nn_scores <- function(nn10_score, nn50_score, obs_agreement,
                            lo = 0.05, hi = 0.80) {
    k <- adaptive_nn_k(obs_agreement, lo = lo, hi = hi)
    w10 <- max(0.0, min(1.0, (50.0 - k) / 40.0))
    w10 * as.numeric(nn10_score) + (1.0 - w10) * as.numeric(nn50_score)
}

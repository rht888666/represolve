#' Low-rank SVD completion predictor
#'
#' @param n_components Integer rank of truncated SVD (default 20).
#' @return A model descriptor with \code{$name} and \code{$fit(T_mat)}.
svd_completion <- function(n_components = 20L) {
    n_components <- as.integer(n_components)
    list(
        name = paste0("svd", n_components),
        fit = function(T_mat) {
            T_mat <- as.matrix(T_mat)
            mu <- colMeans(T_mat)
            X <- sweep(T_mat, 2L, mu, "-")
            k <- min(n_components, min(dim(T_mat)) - 1L)
            W <- svd(X, nu = 0L, nv = k)$v
            list(
                name = paste0("svd", n_components),
                mu = mu,
                W = W,
                predict = function(observed, obs_mask) {
                    obs_mask <- as.logical(obs_mask)
                    y <- as.numeric(observed[obs_mask]) - mu[obs_mask]
                    coef <- qr.solve(W[obs_mask, , drop = FALSE], y)
                    as.numeric(mu + W %*% coef)
                }
            )
        }
    )
}

#' Conditional Gaussian completion with analytical Ledoit-Wolf shrinkage
#'
#' Uses exact Gram shrinkage and the Woodbury push-through identity when
#' \code{n <= p_obs} so high-dimensional screens fit and predict quickly.
#'
#' @return A model descriptor with \code{$name = "lw"} and \code{$fit(T_mat)}.
linear_completion <- function() {
    list(
        name = "lw",
        fit = function(T_mat) {
            T_mat <- as.matrix(T_mat)
            storage.mode(T_mat) <- "double"
            n <- nrow(T_mat)
            p <- ncol(T_mat)
            mu <- colMeans(T_mat)
            X <- sweep(T_mat, 2L, mu, "-")
            X2 <- X * X
            tr <- sum(X2) / n
            mu_trace <- tr / p
            if (n <= p) {
                G <- tcrossprod(X)
                delta_ <- sum(G * G) / (n * n)
            } else {
                S <- crossprod(X) / n
                delta_ <- sum(S * S)
            }
            row_s <- rowSums(X2)
            beta_ <- sum(row_s * row_s)
            beta <- (beta_ / n - delta_) / (p * n)
            delta <- (delta_ - p * mu_trace * mu_trace) / p
            beta <- min(beta, delta)
            shrinkage <- if (beta <= 0 || delta <= 0) 0.0 else beta / delta
            alpha <- (1.0 - shrinkage) / n
            lam <- shrinkage * mu_trace
            cov_mat <- NULL
            if (n > 500L && p <= 10000L) {
                cov_mat <- alpha * crossprod(X)
                diag(cov_mat) <- diag(cov_mat) + lam
            }
            list(
                name = "lw",
                mu = mu,
                shrinkage = shrinkage,
                predict = function(observed, obs_mask) {
                    obs_mask <- as.logical(obs_mask)
                    hid <- !obs_mask
                    pred <- mu
                    y_obs <- as.numeric(observed[obs_mask]) - mu[obs_mask]
                    if (!is.null(cov_mat)) {
                        rows <- cov_mat[obs_mask, , drop = FALSE]
                        v <- solve(rows[, obs_mask, drop = FALSE], y_obs)
                        pred[hid] <- mu[hid] + as.numeric(crossprod(v, rows))[hid]
                        return(pred)
                    }
                    Xo <- X[, obs_mask, drop = FALSE]
                    Xh <- X[, hid, drop = FALSE]
                    p_obs <- ncol(Xo)
                    if (lam <= 0) {
                        cov_oo <- alpha * crossprod(Xo)
                        cov_ho <- alpha * crossprod(Xh, Xo)
                        pred[hid] <- mu[hid] + as.numeric(cov_ho %*% qr.solve(cov_oo, y_obs))
                    } else if (n <= p_obs) {
                        K <- alpha * tcrossprod(Xo)
                        diag(K) <- diag(K) + lam
                        w <- solve(K, alpha * as.numeric(Xo %*% y_obs))
                        pred[hid] <- mu[hid] + as.numeric(crossprod(Xh, w))
                    } else {
                        C <- alpha * crossprod(Xo)
                        diag(C) <- diag(C) + lam
                        v <- solve(C, y_obs)
                        pred[hid] <- mu[hid] + alpha * as.numeric(crossprod(Xh, Xo %*% v))
                    }
                    pred
                }
            )
        }
    )
}

#' Mixture of probabilistic PCA experts (moPPCA) completion predictor
#'
#' @param n_components Number of mixture experts (default 5).
#' @param latent_dim Rank of global latent subspace (default 50).
#' @param expert_dim Rank of each local expert PCA (default 10).
#' @param seed Integer RNG seed for k-means clustering (default 0).
#' @return A model descriptor with \code{$name} and \code{$fit(T_mat)}.
moppca_completion <- function(n_components = 5L, latent_dim = 50L,
                              expert_dim = 10L, seed = 0L) {
    n_components <- as.integer(n_components)
    latent_dim <- as.integer(latent_dim)
    expert_dim <- as.integer(expert_dim)
    seed <- as.integer(seed)
    mname <- paste0("moppca", n_components)
    list(
        name = mname,
        fit = function(T_mat) {
            T_mat <- as.matrix(T_mat)
            mu <- colMeans(T_mat)
            Xc <- sweep(T_mat, 2L, mu, "-")
            k_lat <- min(latent_dim, min(dim(T_mat)) - 1L)
            W_lat <- svd(Xc, nu = 0L, nv = k_lat)$v
            latent <- Xc %*% W_lat
            k_cl <- min(n_components, max(1L, nrow(latent) - 1L))
            set.seed(seed)
            km <- stats::kmeans(latent, centers = k_cl, nstart = 5L, iter.max = 50L)
            labels <- km$cluster
            centers <- km$centers
            vars <- matrix(1.0, nrow = k_cl, ncol = ncol(latent))
            experts <- vector("list", k_cl)
            for (c in seq_len(k_cl)) {
                idx <- which(labels == c)
                if (length(idx) > 1L) {
                    sub_lat <- latent[idx, , drop = FALSE]
                    vars[c, ] <- pmax(apply(sub_lat, 2L, stats::var), 1e-4)
                }
                if (length(idx) >= 15L) {
                    sub_T <- T_mat[idx, , drop = FALSE]
                    e_mu <- colMeans(sub_T)
                    e_Xc <- sweep(sub_T, 2L, e_mu, "-")
                    ke <- min(expert_dim, min(dim(sub_T)) - 1L)
                    experts[[c]] <- list(
                        mu = e_mu,
                        W = svd(e_Xc, nu = 0L, nv = ke)$v
                    )
                }
            }
            list(
                name = mname,
                mu = mu,
                predict = function(observed, obs_mask) {
                    obs_mask <- as.logical(obs_mask)
                    query <- as.numeric(
                        qr.solve(W_lat[obs_mask, , drop = FALSE],
                                 as.numeric(observed[obs_mask]) - mu[obs_mask])
                    )
                    logp <- vapply(seq_len(k_cl), function(c) {
                        -0.5 * sum(((query - centers[c, ])^2) / vars[c, ] +
                                   log(vars[c, ]))
                    }, numeric(1L))
                    logp <- logp - max(logp)
                    w <- exp(logp)
                    w <- w / sum(w)
                    weighted <- numeric(length(mu))
                    total <- 0.0
                    for (c in seq_len(k_cl)) {
                        ex <- experts[[c]]
                        if (is.null(ex) || w[c] < 1e-6) {
                            next
                        }
                        coef <- qr.solve(ex$W[obs_mask, , drop = FALSE],
                                         as.numeric(observed[obs_mask]) - ex$mu[obs_mask])
                        est <- ex$mu + as.numeric(ex$W %*% coef)
                        weighted <- weighted + w[c] * est
                        total <- total + w[c]
                    }
                    if (total > 0) weighted / total else mu
                }
            )
        }
    )
}

#' Nearest-perturbation-centroid predictor
#'
#' @param k Integer number of nearest training neighbours (default 10).
#' @return A model descriptor with \code{$name} and \code{$fit(T_mat)}.
nearest_centroid <- function(k = 10L) {
    k <- as.integer(k)
    mname <- paste0("nn", k)
    list(
        name = mname,
        fit = function(T_mat) {
            T_mat <- as.matrix(T_mat)
            storage.mode(T_mat) <- "double"
            mu <- colMeans(T_mat)
            list(
                name = mname,
                mu = mu,
                predict = function(observed, obs_mask) {
                    obs_mask <- as.logical(obs_mask)
                    obs_pos <- which(obs_mask)
                    if (length(obs_pos) == 0L) {
                        return(mu)
                    }
                    sims <- neighbour_similarity(T_mat, mu, obs_pos, observed[obs_pos])
                    nn <- order(-sims)[seq_len(min(k, length(sims)))]
                    pred <- colMeans(T_mat[nn, , drop = FALSE])
                    pred[obs_mask] <- observed[obs_mask]
                    as.numeric(pred)
                }
            )
        }
    )
}

#' Default estimator panel used in the paper
#'
#' @param nn_ks Integer vector of neighbourhood sizes (default \code{c(10L, 50L)}).
#' @return List of model descriptors (\code{svd20}, \code{lw}, \code{moppca5}, \code{nn10}, \code{nn50}).
default_models <- function(nn_ks = c(10L, 50L)) {
    c(
        list(svd_completion(20L), linear_completion(), moppca_completion(5L)),
        lapply(nn_ks, nearest_centroid)
    )
}

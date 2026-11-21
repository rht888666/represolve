#' Safe Pearson correlation
#'
#' Computes the Pearson correlation between two numeric vectors, returning
#' \code{NA_real_} when either vector has zero standard deviation or fewer
#' than two finite observations.
#'
#' @param x Numeric vector.
#' @param y Numeric vector.
#' @return Scalar numeric correlation or \code{NA_real_}.
pearson <- function(x, y) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    if (length(x) < 2L) {
        return(NA_real_)
    }
    sx <- stats::sd(x)
    sy <- stats::sd(y)
    if (!is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) {
        return(NA_real_)
    }
    as.numeric(stats::cor(x, y))
}

#' Centered correlation between prediction and measurement
#'
#' Centers both \code{pred} and \code{measured} on reference mean \code{mu}
#' before computing Pearson correlation.
#'
#' @param pred Numeric vector of predicted values.
#' @param measured Numeric vector of measured values.
#' @param mu Numeric vector of reference training-fold gene means.
#' @return Scalar numeric correlation.
centered_score <- function(pred, measured, mu) {
    pearson(as.numeric(pred) - as.numeric(mu), as.numeric(measured) - as.numeric(mu))
}

#' Replicate agreement (rep_cen) between independent halves
#'
#' @param half_a Numeric vector for half A.
#' @param half_b Numeric vector for half B.
#' @param mu Numeric vector of reference training-fold gene means.
#' @return Scalar numeric correlation.
replicate_agreement <- function(half_a, half_b, mu) {
    pearson(as.numeric(half_a) - as.numeric(mu), as.numeric(half_b) - as.numeric(mu))
}

#' Centered correlation of each training perturbation to a test observed profile
#'
#' @param train_profiles Numeric matrix of shape (n_train, n_genes).
#' @param train_mu Numeric vector of length n_genes.
#' @param obs_pos Integer indices of observed genes (1-based).
#' @param test_obs Numeric vector of observed values at \code{obs_pos}.
#' @return Numeric vector of similarities of length n_train (-Inf where undefined).
neighbour_similarity <- function(train_profiles, train_mu, obs_pos, test_obs) {
    train_eff <- sweep(
        train_profiles[, obs_pos, drop = FALSE],
        2L,
        train_mu[obs_pos],
        "-"
    )
    train_eff <- train_eff - rowMeans(train_eff)
    test_eff <- as.numeric(test_obs) - as.numeric(train_mu[obs_pos])
    test_eff <- test_eff - mean(test_eff)
    tn <- sqrt(rowSums(train_eff^2))
    den <- tn * sqrt(sum(test_eff^2))
    out <- as.numeric(train_eff %*% test_eff) / den
    out[!is.finite(out)] <- -Inf
    out
}

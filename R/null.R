#' Empirical p-values from a permutation distribution
#'
#' Uses \code{(r + 1) / (B + 1)} calibration across \code{B} permutation draws.
#'
#' @param observed Numeric vector of observed test statistics.
#' @param permuted Numeric matrix (n_perms x n_tests) or vector of permutation draws.
#' @param alternative One of \code{"two-sided"}, \code{"greater"}, or \code{"less"}.
#' @return Numeric vector of empirical p-values.
permutation_null <- function(observed, permuted, alternative = "two-sided") {
    observed <- as.numeric(observed)
    if (is.null(dim(permuted))) {
        permuted <- matrix(as.numeric(permuted), ncol = 1L)
    } else {
        permuted <- as.matrix(permuted)
    }
    B <- nrow(permuted)
    n_tests <- length(observed)
    p <- numeric(n_tests)
    for (j in seq_len(n_tests)) {
        col_j <- permuted[, min(j, ncol(permuted))]
        if (alternative == "two-sided") {
            r <- sum(abs(col_j) >= abs(observed[j]))
        } else if (alternative == "greater") {
            r <- sum(col_j >= observed[j])
        } else {
            r <- sum(col_j <= observed[j])
        }
        p[j] <- (r + 1.0) / (B + 1.0)
    }
    p
}

#' Benjamini-Hochberg FDR q-values
#'
#' @param p Numeric vector of p-values.
#' @return Numeric vector of BH-adjusted q-values.
bh_qvalues <- function(p) {
    as.numeric(stats::p.adjust(as.numeric(p), method = "BH"))
}

#' Empirical type-I error rate check on simulated null screens
#'
#' @param simulate_screen Function taking no arguments (or integer sim index) and
#'   returning a numeric vector of null p-values for one simulated screen.
#' @param n_sims Integer number of simulated null screens (default 200).
#' @param alpha Significance threshold on BH q-values (default 0.05).
#' @param seed Integer RNG seed (default 0).
#' @return Named list with \code{n_sims}, \code{alpha}, \code{any_hit_rate},
#'   and \code{mean_false_positives}.
type1_check <- function(simulate_screen, n_sims = 200L, alpha = 0.05, seed = 0L) {
    set.seed(as.integer(seed))
    n_sims <- as.integer(n_sims)
    any_fp <- 0L
    total_fp <- 0L
    has_arg <- length(formals(simulate_screen)) > 0L
    for (i in seq_len(n_sims)) {
        p <- if (has_arg) simulate_screen(i) else simulate_screen()
        fp <- sum(bh_qvalues(p) < alpha, na.rm = TRUE)
        if (fp > 0L) {
            any_fp <- any_fp + 1L
        }
        total_fp <- total_fp + fp
    }
    list(
        n_sims = n_sims,
        alpha = alpha,
        any_hit_rate = any_fp / n_sims,
        mean_false_positives = total_fp / n_sims
    )
}

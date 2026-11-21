.qcut_rank <- function(x, n) {
    r <- rank(x, ties.method = "first", na.last = "keep")
    br <- stats::quantile(r, probs = seq(0, 1, length.out = n + 1L), na.rm = TRUE)
    as.integer(cut(r, breaks = br, include.lowest = TRUE, labels = FALSE)) - 1L
}

#' Assign perturbations to replicate-agreement deciles
#'
#' @param df Benchmark output data frame.
#' @param score_col Column name containing replicate agreement (default \code{"rep_cen"}).
#' @param n Integer number of bins (default 10).
#' @param by_fold Logical; whether to bin within each fold.
#' @return Integer vector of 0-indexed bin assignments.
assign_deciles <- function(df, score_col = "rep_cen", n = 10L, by_fold = TRUE) {
    n <- as.integer(n)
    if (!by_fold) {
        return(.qcut_rank(df[[score_col]], n))
    }
    out <- integer(nrow(df))
    for (f in unique(df$fold)) {
        idx <- which(df$fold == f)
        out[idx] <- .qcut_rank(df[[score_col]][idx], n)
    }
    out
}

#' Mean independent-half score per model per replicate-agreement decile
#'
#' @param df Benchmark output data frame.
#' @param score_suffix Suffix identifying score columns (default \code{"_B_cen"}).
#' @param n_deciles Integer number of quantile bins (default 10).
#' @return Long data frame with columns \code{decile}, \code{model}, \code{mean}, \code{sem}, \code{n}.
decile_ranking <- function(df, score_suffix = "_B_cen", n_deciles = 10L) {
    df$decile <- assign_deciles(df, n = n_deciles)
    pat <- paste0(score_suffix, "$")
    model_cols <- grep(pat, names(df), value = TRUE)
    rows <- list()
    k <- 1L
    for (d in sort(unique(df$decile[!is.na(df$decile)]))) {
        sub <- df[df$decile == d, , drop = FALSE]
        for (c in model_cols) {
            v <- sub[[c]]
            v <- v[is.finite(v)]
            rows[[k]] <- data.frame(
                decile = as.integer(d),
                model = sub(pat, "", c),
                mean = if (length(v)) mean(v) else NA_real_,
                sem = if (length(v)) stats::sd(v) / sqrt(max(length(v), 1L)) else NA_real_,
                n = length(v),
                stringsAsFactors = FALSE
            )
            k <- k + 1L
        }
    }
    do.call(rbind, rows)
}

#' Paired perturbation bootstrap 95 percent interval for mean score difference
#'
#' @param df Benchmark output data frame.
#' @param col_a First score column name.
#' @param col_b Second score column name.
#' @param n_boot Integer number of bootstrap resamples (default 10000).
#' @param seed Integer RNG seed (default 0).
#' @return Named list with \code{diff}, \code{ci_low}, \code{ci_high}, and \code{n}.
bootstrap_ci <- function(df, col_a, col_b, n_boot = 10000L, seed = 0L) {
    d <- df[[col_a]] - df[[col_b]]
    d <- d[is.finite(d)]
    n <- length(d)
    set.seed(as.integer(seed))
    idx <- matrix(sample.int(n, size = as.integer(n_boot) * n, replace = TRUE),
                  nrow = as.integer(n_boot), ncol = n)
    boot <- rowMeans(matrix(d[idx], nrow = as.integer(n_boot), ncol = n))
    qs <- stats::quantile(boot, probs = c(0.025, 0.975), names = FALSE)
    list(diff = mean(d), ci_low = qs[1L], ci_high = qs[2L], n = n)
}

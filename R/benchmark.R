#' Fold-frozen masked-completion benchmark
#'
#' Evaluates masked-gene completion across \code{n_folds} perturbation folds,
#' selecting top-variance genes and fitting all models strictly inside each
#' training fold and scoring predictions against both same-half (\code{A}) and
#' independent-half (\code{B}) measurements.
#'
#' @param half_a Numeric matrix (n_perturbations x n_genes) for partition A.
#' @param half_b Numeric matrix (n_perturbations x n_genes) for partition B.
#' @param pert_names Character vector of perturbation identifiers.
#' @param genes Character vector of gene names.
#' @param n_folds Integer number of cross-validation folds (default 5).
#' @param n_genes_selected Number of top-variance genes selected per training fold.
#' @param mask_fraction Fraction of selected genes observed per test perturbation.
#' @param seed Integer RNG seed for fold assignment and deterministic masks.
#' @param models Optional list of model descriptors (defaults to \code{default_models()}).
#' @return Data frame with one row per perturbation containing \code{pert},
#'   \code{fold}, \code{n_observed}, \code{rep_cen}, and
#'   \code{<model>_A_cen} / \code{<model>_B_cen}.
masked_completion_benchmark <- function(half_a, half_b, pert_names, genes,
                                        n_folds = 5L, n_genes_selected = 5000L,
                                        mask_fraction = 0.1, seed = 0L,
                                        models = NULL) {
    if (is.null(models)) {
        models <- default_models()
    }
    half_a <- as.matrix(half_a)
    half_b <- as.matrix(half_b)
    storage.mode(half_a) <- "double"
    storage.mode(half_b) <- "double"
    mean_all <- (half_a + half_b) / 2.0
    n_pert <- nrow(half_a)
    n_folds <- as.integer(n_folds)
    set.seed(as.integer(seed))
    perm <- sample.int(n_pert)
    fold_ids <- integer(n_pert)
    fold_ids[perm] <- rep(seq_len(n_folds), length.out = n_pert)

    n_sel <- min(as.integer(n_genes_selected), ncol(mean_all))
    rows <- list()
    row_i <- 1L
    for (fold in seq_len(n_folds)) {
        train <- which(fold_ids != fold)
        test <- which(fold_ids == fold)
        tr_sub <- mean_all[train, , drop = FALSE]
        vars <- apply(tr_sub, 2L, stats::var)
        selected <- order(-vars)[seq_len(n_sel)]
        T_mat <- mean_all[train, selected, drop = FALSE]
        mu <- colMeans(T_mat)
        fitted_models <- lapply(models, function(m) m$fit(T_mat))

        for (idx in test) {
            set.seed(as.integer(seed) + 1000L + idx)
            observed <- stats::runif(n_sel) < mask_fraction
            if (!any(observed)) {
                observed[1L] <- TRUE
            }
            if (all(observed)) {
                observed[n_sel] <- FALSE
            }
            hid_pos <- which(!observed)
            aa <- half_a[idx, selected]
            bb <- half_b[idx, selected]
            row_entry <- list(
                pert = as.character(pert_names[idx]),
                fold = fold - 1L,
                n_observed = sum(observed),
                rep_cen = replicate_agreement(aa[hid_pos], bb[hid_pos], mu[hid_pos])
            )
            for (fm in fitted_models) {
                pred <- fm$predict(aa, observed)
                row_entry[[paste0(fm$name, "_A_cen")]] <- pearson(
                    pred[hid_pos] - mu[hid_pos], aa[hid_pos] - mu[hid_pos]
                )
                row_entry[[paste0(fm$name, "_B_cen")]] <- pearson(
                    pred[hid_pos] - mu[hid_pos], bb[hid_pos] - mu[hid_pos]
                )
            }
            rows[[row_i]] <- as.data.frame(row_entry, stringsAsFactors = FALSE)
            row_i <- row_i + 1L
        }
    }
    do.call(rbind, rows)
}

#' Summarize benchmark same-half vs independent-half scores
#'
#' @param df Data frame returned by \code{masked_completion_benchmark}.
#' @param models Optional character vector of model names.
#' @return Data frame sorted by \code{independent_half} descending.
summarize_benchmark <- function(df, models = NULL) {
    if (is.null(models)) {
        a_cols <- grep("_A_cen$", names(df), value = TRUE)
        models <- sub("_A_cen$", "", a_cols)
    }
    out <- lapply(models, function(m) {
        a <- df[[paste0(m, "_A_cen")]]
        b <- df[[paste0(m, "_B_cen")]]
        ma <- mean(a, na.rm = TRUE)
        mb <- mean(b, na.rm = TRUE)
        data.frame(
            model = m,
            same_half = ma,
            independent_half = mb,
            inflation = ma - mb,
            n = sum(is.finite(a)),
            stringsAsFactors = FALSE
        )
    })
    res <- do.call(rbind, out)
    res[order(-res$independent_half), , drop = FALSE]
}

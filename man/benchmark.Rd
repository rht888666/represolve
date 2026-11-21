\name{benchmark}
\alias{masked_completion_benchmark}
\alias{summarize_benchmark}
\title{Fold-frozen masked-completion benchmark}
\description{
Runs fold-frozen K-fold cross-validation where top-variance genes and all
estimators are fit strictly inside each training fold and scored against both
same-half (\code{A_cen}) and independent-half (\code{B_cen}) profiles.
}
\usage{
masked_completion_benchmark(half_a, half_b, pert_names, genes, n_folds = 5L,
                            n_genes_selected = 5000L, mask_fraction = 0.1,
                            seed = 0L, models = NULL)
summarize_benchmark(df, models = NULL)
}
\arguments{
\item{half_a, half_b}{Numeric matrices of shape \code{(n_perturbations, n_genes)}.}
\item{pert_names, genes}{Character vectors of perturbation and gene names.}
\item{n_folds}{Integer number of cross-validation folds (default \code{5L}).}
\item{n_genes_selected}{Integer number of top-variance genes selected per fold.}
\item{mask_fraction}{Fraction of selected genes treated as observed.}
\item{seed}{Integer RNG seed for fold splits and deterministic masks.}
\item{models}{Optional list of model descriptors (or character model names
for \code{summarize_benchmark}).}
\item{df}{Data frame returned by \code{masked_completion_benchmark}.}
}
\value{
Per-perturbation score data frame from \code{masked_completion_benchmark}, or a
per-model summary data frame from \code{summarize_benchmark}.
}
\examples{
set.seed(1)
latent <- matrix(rnorm(40 * 80, sd = 0.4), nrow = 40)
a <- latent + matrix(rnorm(40 * 80, sd = 0.6), nrow = 40)
b <- latent + matrix(rnorm(40 * 80, sd = 0.6), nrow = 40)
df <- masked_completion_benchmark(
    a, b, paste0("p", 1:40), paste0("g", 1:80),
    n_folds = 3L, n_genes_selected = 60L,
    models = list(svd_completion(5L), linear_completion(), nearest_centroid(5L))
)
summarize_benchmark(df)
}

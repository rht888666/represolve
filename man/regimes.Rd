\name{regimes}
\alias{assign_deciles}
\alias{decile_ranking}
\alias{bootstrap_ci}
\title{Regime stratification and paired perturbation bootstrap}
\description{
Stratifies perturbations by replicate agreement (\code{rep_cen}) into quantile
bins and computes paired perturbation-bootstrap 95 percent intervals for
mean model score differences.
}
\usage{
assign_deciles(df, score_col = "rep_cen", n = 10L, by_fold = TRUE)
decile_ranking(df, score_suffix = "_B_cen", n_deciles = 10L)
bootstrap_ci(df, col_a, col_b, n_boot = 10000L, seed = 0L)
}
\arguments{
\item{df}{Benchmark output data frame.}
\item{score_col}{Column name containing replicate agreement.}
\item{n, n_deciles}{Integer number of quantile bins.}
\item{by_fold}{Logical; whether to bin within each fold.}
\item{score_suffix}{Suffix identifying model score columns.}
\item{col_a, col_b}{Column names of the two model score vectors to compare.}
\item{n_boot}{Integer number of bootstrap resamples.}
\item{seed}{Integer RNG seed.}
}
\value{
Integer vector of 0-based bin assignments for \code{assign_deciles}; a long
data frame for \code{decile_ranking}; a named list with \code{diff},
\code{ci_low}, \code{ci_high}, and \code{n} for \code{bootstrap_ci}.
}
\examples{
df <- data.frame(fold = rep(0:1, each = 10), rep_cen = seq(0.1, 0.9, length.out = 20),
                 lw_B_cen = rnorm(20, 0.5, 0.1), svd20_B_cen = rnorm(20, 0.4, 0.1))
decile_ranking(df, n_deciles = 5L)
bootstrap_ci(df, "lw_B_cen", "svd20_B_cen", n_boot = 200L)
}

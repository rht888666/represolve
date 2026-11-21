\name{null}
\alias{permutation_null}
\alias{bh_qvalues}
\alias{type1_check}
\title{Permutation-null calibration and type-I error check}
\description{
Computes exact \code{(r + 1)/(B + 1)} empirical p-values against permutation
draws, Benjamini-Hochberg FDR q-values, and empirical type-I error rates across
simulated null screens.
}
\usage{
permutation_null(observed, permuted, alternative = "two-sided")
bh_qvalues(p)
type1_check(simulate_screen, n_sims = 200L, alpha = 0.05, seed = 0L)
}
\arguments{
\item{observed}{Numeric vector of observed test statistics.}
\item{permuted}{Numeric matrix or vector of permutation draws.}
\item{alternative}{One of \code{"two-sided"}, \code{"greater"}, or \code{"less"}.}
\item{p}{Numeric vector of p-values.}
\item{simulate_screen}{Function returning a numeric vector of null p-values.}
\item{n_sims}{Integer number of simulated null screens.}
\item{alpha}{FDR threshold.}
\item{seed}{Integer RNG seed.}
}
\value{
Numeric vectors for \code{permutation_null} and \code{bh_qvalues}; a named
list with \code{n_sims}, \code{alpha}, \code{any_hit_rate}, and
\code{mean_false_positives} for \code{type1_check}.
}
\examples{
set.seed(0)
p <- permutation_null(c(2.1, 0.2), matrix(rnorm(200), nrow = 100, ncol = 2))
bh_qvalues(p)
}

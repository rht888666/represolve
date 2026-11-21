\name{metrics}
\alias{pearson}
\alias{centered_score}
\alias{replicate_agreement}
\alias{neighbour_similarity}
\title{Centered correlation and neighbour similarity metrics}
\description{
Computes safe Pearson correlation, training-mean-centered completion scores,
between-half replicate agreement (\code{rep_cen}), and observed-gene centroid
similarities.
}
\usage{
pearson(x, y)
centered_score(pred, measured, mu)
replicate_agreement(half_a, half_b, mu)
neighbour_similarity(train_profiles, train_mu, obs_pos, test_obs)
}
\arguments{
\item{x, y}{Numeric vectors.}
\item{pred, measured}{Numeric vectors of predicted and measured expression.}
\item{mu}{Numeric vector of reference training-fold gene means.}
\item{half_a, half_b}{Numeric vectors for independent replicate halves.}
\item{train_profiles}{Numeric matrix of shape \code{(n_train, n_genes)}.}
\item{train_mu}{Numeric vector of length \code{n_genes}.}
\item{obs_pos}{Integer indices of observed genes (1-based).}
\item{test_obs}{Numeric vector of observed values at \code{obs_pos}.}
}
\value{
A scalar numeric correlation (\code{NA_real_} when degenerate) for
\code{pearson}, \code{centered_score}, and \code{replicate_agreement}; a numeric
vector of similarities (\code{-Inf} where undefined) for
\code{neighbour_similarity}.
}
\examples{
mu <- c(0, 0, 0)
centered_score(c(1, 2, 3), c(1, 2, 4), mu)
}

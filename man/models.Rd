\name{models}
\alias{svd_completion}
\alias{linear_completion}
\alias{moppca_completion}
\alias{nearest_centroid}
\alias{default_models}
\title{Masked-completion baseline estimators}
\description{
Constructors for the baseline masked-completion estimators evaluated in the
\code{represolve} benchmark: truncated SVD, analytical Ledoit-Wolf shrinkage,
mixture of probabilistic PCA experts, and nearest-perturbation-centroid
retrieval.
}
\usage{
svd_completion(n_components = 20L)
linear_completion()
moppca_completion(n_components = 5L, latent_dim = 50L, expert_dim = 10L,
                  seed = 0L)
nearest_centroid(k = 10L)
default_models(nn_ks = c(10L, 50L))
}
\arguments{
\item{n_components}{Integer rank (for \code{svd_completion}) or number of
mixture experts (for \code{moppca_completion}).}
\item{latent_dim}{Rank of the global latent subspace for \code{moppca_completion}.}
\item{expert_dim}{Rank of each local expert PCA in \code{moppca_completion}.}
\item{seed}{Integer RNG seed for k-means clustering.}
\item{k}{Integer number of nearest training neighbours.}
\item{nn_ks}{Integer vector of neighbourhood sizes for \code{default_models}.}
}
\value{
A model descriptor list with \code{name} and \code{fit(T_mat)} (or a list of
model descriptors for \code{default_models}).
}
\examples{
mods <- default_models()
length(mods)
}

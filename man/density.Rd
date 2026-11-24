\name{density}
\alias{observed_agreement}
\alias{cosine_knn_density}
\alias{route}
\alias{soft_route_weights}
\alias{route_soft}
\alias{route_soft_model}
\alias{adaptive_nn_k}
\alias{blend_nn_scores}
\title{Prediction-time regime diagnostics and hard/soft routing}
\description{
Computes prediction-time observed-gene replicate agreement
(\code{obs_agreement}), scale-invariant cosine k-NN perturbation density
(\code{d_cos10}), hard regime routing (Table 3), RBF Soft Regime Router weights,
and adaptive neighbourhood interpolation along the mid-reliability band.
}
\usage{
observed_agreement(half_a_obs, half_b_obs, train_mean)
cosine_knn_density(train_profiles, k = 10L, center = TRUE)
route(obs_agreement = NULL, cos_density = NULL, n_perts = NULL)
soft_route_weights(obs_agreement, cos_density, temperature = 0.08,
                   prototypes = NULL, scale = c(1.0, 1.2))
route_soft(obs_agreement, cos_density, temperature = 0.08, prototypes = NULL)
route_soft_model(obs_agreement, cos_density, temperature = 0.08,
                 prototypes = NULL)
adaptive_nn_k(obs_agreement, k_low = 50L, k_mid = 10L, lo = 0.05, hi = 0.80)
blend_nn_scores(nn10_score, nn50_score, obs_agreement, lo = 0.05, hi = 0.80)
}
\arguments{
\item{half_a_obs, half_b_obs}{Numeric vectors on observed genes in half A and half B.}
\item{train_mean}{Numeric vector of training means on observed genes.}
\item{train_profiles}{Numeric matrix of shape \code{(n_perturbations, n_genes)}.}
\item{k}{Integer neighbour rank (default \code{10L}).}
\item{center}{Logical; whether to column-center before L2 normalization.}
\item{obs_agreement}{Scalar observed-gene agreement (or \code{NULL} for \code{route}).}
\item{cos_density}{Scalar cosine k-NN density \code{d_cos10}.}
\item{n_perts}{Optional perturbation count; the dense-grid branch requires a large panel (\code{n_perts >= 1000}) when supplied.}
\item{temperature}{Positive RBF kernel temperature (default \code{0.08}).}
\item{prototypes}{Optional named list of 2-element regime prototypes.}
\item{scale}{Numeric length-2 axis scaling vector (default \code{c(1.0, 1.2)}).}
\item{k_low, k_mid}{Neighbourhood sizes at lower and upper mid-band endpoints.}
\item{lo, hi}{Agreement interval bounds for adaptive neighbourhood sizing.}
\item{nn10_score, nn50_score}{Scalar scores for \code{nn10} and \code{nn50}.}
}
\value{
Numeric scalars for \code{observed_agreement}, \code{cosine_knn_density},
\code{adaptive_nn_k}, and \code{blend_nn_scores}; character routing labels for
\code{route}, \code{route_soft}, and \code{route_soft_model}; named numeric
vector summing to 1 for \code{soft_route_weights}.
}
\examples{
route(0.29, 0.65)
soft_route_weights(0.29, 0.65)
route_soft_model(0.96, 0.20)
}

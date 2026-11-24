library(represolve)

set.seed(0)
n_pert <- 60L
n_genes <- 200L
latent <- matrix(rnorm(n_pert * n_genes, sd = 0.35), nrow = n_pert)
half_a <- latent + matrix(rnorm(n_pert * n_genes, sd = 0.65), nrow = n_pert)
half_b <- latent + matrix(rnorm(n_pert * n_genes, sd = 0.65), nrow = n_pert)
perts <- paste0("pert", seq_len(n_pert))
genes <- paste0("g", seq_len(n_genes))

mods <- list(
    svd_completion(10L),
    linear_completion(),
    moppca_completion(3L, latent_dim = 20L, expert_dim = 5L),
    nearest_centroid(10L)
)
df <- masked_completion_benchmark(
    half_a, half_b, perts, genes,
    n_folds = 3L, n_genes_selected = 150L, models = mods
)
stopifnot(nrow(df) == n_pert)
stopifnot("rep_cen" %in% names(df))

s <- summarize_benchmark(df)
stopifnot(all(c("model", "same_half", "independent_half", "inflation", "n") %in% names(s)))

dr <- decile_ranking(df, n_deciles = 5L)
stopifnot(all(c("decile", "model", "mean", "sem", "n") %in% names(dr)))

ci <- bootstrap_ci(df, "lw_B_cen", "svd10_B_cen", n_boot = 200L)
stopifnot(ci$ci_low <= ci$diff && ci$diff <= ci$ci_high)

stopifnot(route(0.29, 0.65) == "centroid_nn10")
stopifnot(route(0.99, 0.001, n_perts = 2224) == "nn10_dense_grid")
stopifnot(route(0.99, 0.001, n_perts = 198) == "lw_shrinkage")
w <- soft_route_weights(0.29, 0.65, temperature = 0.08)
stopifnot(abs(sum(w) - 1.0) < 1e-6)
stopifnot(route_soft_model(0.96, 0.20) == "lw")
stopifnot(route_soft_model(0.07, 0.85) == "nn50")
stopifnot(route_soft(-0.02, 0.90) == "flag_near_floor")
stopifnot(adaptive_nn_k(0.05) == 50.0)
stopifnot(adaptive_nn_k(0.80) == 10.0)
mid <- adaptive_nn_k(0.425)
stopifnot(mid > 10.0 && mid < 50.0)
b <- blend_nn_scores(0.4, 0.2, 0.425)
stopifnot(b >= 0.2 && b <= 0.4)

dcos <- cosine_knn_density(half_a, k = 10L)
stopifnot(is.finite(dcos) && dcos >= 0 && dcos <= 2)
oa <- observed_agreement(half_a[1, 1:20], half_b[1, 1:20], colMeans(half_a)[1:20])
stopifnot(is.finite(oa))

p <- permutation_null(rnorm(30), matrix(rnorm(500 * 30), nrow = 500, ncol = 30))
q <- bh_qvalues(p)
stopifnot(length(p) == 30L && length(q) == 30L)
t1 <- type1_check(function() runif(25), n_sims = 20L)
stopifnot(t1$any_hit_rate >= 0 && t1$any_hit_rate <= 1)

ds <- list_datasets()
stopifnot(nrow(ds) == 9L && "sciplex3" %in% ds$dataset)
cat("ALL R TESTS PASSED\n")

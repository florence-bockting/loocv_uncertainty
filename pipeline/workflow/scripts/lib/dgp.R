# Data generating process of Sivula et al. (2025), Section 4.
# Mirrors ProblemRun.make_data in simulated/general_setup.py:
#   y = X beta + eps,  X[, -1] ~ N(0, 1),  eps ~ N(mu_d, sigma2_d),
# where the first n_obs_out observations get the outlier mean
# out_dev * sqrt(sigma2_d + sum(beta^2)).

# beta = (intercept, other effects..., beta_t); beta_t is only in model B.
make_beta <- function(beta_t, n_dim = 3, beta_other = 1, beta_intercept = 0) {
  beta <- c(rep(beta_other, n_dim - 1), beta_t)
  beta[1] <- beta_intercept
  beta
}

# Draws n_sets data sets of size n_obs_max and keeps the first n_obs rows,
# so that cells with different n_obs share the same random numbers.
# Returns X as an [n_obs, n_dim, n_sets] array and y as an [n_obs, n_sets]
# matrix.
make_data <- function(n_sets, n_obs, n_obs_max, beta, out_dev,
                      n_obs_out = 1, sigma2_d = 1) {
  n_dim <- length(beta)
  x <- array(rnorm(n_obs_max * (n_dim - 1) * n_sets),
             c(n_obs_max, n_dim - 1, n_sets))
  eps <- matrix(rnorm(n_obs_max * n_sets, sd = sqrt(sigma2_d)),
                n_obs_max, n_sets)
  mu_d <- numeric(n_obs_max)
  mu_d[seq_len(n_obs_out)] <- out_dev * sqrt(sigma2_d + sum(beta^2))
  eps <- eps + mu_d

  keep <- seq_len(n_obs)
  X <- array(1, c(n_obs, n_dim, n_sets))
  X[, -1, ] <- x[keep, , , drop = FALSE]
  y <- eps[keep, , drop = FALSE]
  for (k in seq_len(n_dim)) {
    y <- y + beta[k] * matrix(X[, k, ], n_obs, n_sets)
  }
  list(X = X, y = y)
}

# Stacks all data sets into one design matrix and one response vector.
pool_sets <- function(data) {
  dims <- dim(data$X)
  X <- aperm(data$X, c(1, 3, 2))
  dim(X) <- c(dims[1] * dims[3], dims[2])
  list(X = X, y = as.vector(data$y))
}

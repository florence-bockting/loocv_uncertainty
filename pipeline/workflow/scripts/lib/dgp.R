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

# Binomial variant: y_i ~ Bernoulli(plogis(X_i beta)). The outlier analogue
# of mu_d shifts the linear predictor of the first n_obs_out observations
# away from their own sign, so that a large out_dev makes those
# observations almost surely the class the covariates argue against.
# As in make_data, n_obs_max rows are drawn and the first n_obs are kept,
# so that cells with different n_obs share the same random numbers.
make_data_binomial <- function(n_sets, n_obs, n_obs_max, beta, out_dev,
                               n_obs_out = 1) {
  n_dim <- length(beta)
  x <- array(rnorm(n_obs_max * (n_dim - 1) * n_sets),
             c(n_obs_max, n_dim - 1, n_sets))
  u <- matrix(runif(n_obs_max * n_sets), n_obs_max, n_sets)

  keep <- seq_len(n_obs)
  X <- array(1, c(n_obs, n_dim, n_sets))
  X[, -1, ] <- x[keep, , , drop = FALSE]
  eta <- matrix(0, n_obs, n_sets)
  for (k in seq_len(n_dim)) {
    eta <- eta + beta[k] * matrix(X[, k, ], n_obs, n_sets)
  }
  out <- seq_len(min(n_obs_out, n_obs))
  eta[out, ] <- eta[out, , drop = FALSE] -
    out_dev * sign_pos(eta[out, , drop = FALSE])
  y <- (u[keep, , drop = FALSE] < plogis(eta)) * 1
  list(X = X, y = y)
}

# sign() with ties broken towards +1, so that the outlier shift never
# vanishes.
sign_pos <- function(x) ifelse(x < 0, -1, 1)

# Draws the data of one cell, for any family.
make_data_family <- function(family, ...) {
  switch(family,
    gaussian = make_data(...),
    binomial = make_data_binomial(...),
    poisson = make_data_poisson(...),
    stop("unknown family: ", family)
  )
}

# Poisson variant: y_i ~ Poisson(exp(eta_i)). The rate grows as exp(eta),
# so the coefficients of make_beta() are scaled down and an intercept is
# added: with POIS_SCALE and POIS_INTERCEPT the rate stays near e and the
# counts stay small enough for a count grid. beta_t = 1 is still a clear
# difference between the models, as in the other families.
#
# The outlier shifts the linear predictor of the first n_obs_out
# observations away from their own side of the intercept, as in
# make_data_binomial, so out_dev is on the log scale.
POIS_SCALE <- 0.5
POIS_INTERCEPT <- 1

make_data_poisson <- function(n_sets, n_obs, n_obs_max, beta, out_dev,
                              n_obs_out = 1) {
  n_dim <- length(beta)
  b <- beta * POIS_SCALE
  b[1] <- POIS_INTERCEPT
  x <- array(rnorm(n_obs_max * (n_dim - 1) * n_sets),
             c(n_obs_max, n_dim - 1, n_sets))

  keep <- seq_len(n_obs)
  X <- array(1, c(n_obs, n_dim, n_sets))
  X[, -1, ] <- x[keep, , , drop = FALSE]
  eta <- matrix(0, n_obs, n_sets)
  for (k in seq_len(n_dim)) {
    eta <- eta + b[k] * matrix(X[, k, ], n_obs, n_sets)
  }
  out <- seq_len(min(n_obs_out, n_obs))
  eta[out, ] <- eta[out, , drop = FALSE] -
    out_dev * sign_pos(eta[out, , drop = FALSE] - POIS_INTERCEPT)
  y <- matrix(rpois(length(eta), exp(eta)), n_obs, n_sets)
  list(X = X, y = y)
}

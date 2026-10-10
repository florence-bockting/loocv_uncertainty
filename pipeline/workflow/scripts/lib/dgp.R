# Data generating process of Sivula et al. (2025), Section 4.
# Mirrors ProblemRun.make_data in simulated/general_setup.py:
#   y = X beta + eps,  X[, -1] ~ N(0, 1),  eps ~ N(mu_d, sigma2_d),
# where the first n_obs_out observations get the outlier mean
# out_dev * sqrt(sigma2_d + sum(beta^2)).
#
# Extension (gaussian family only): x_dist sets the distribution of the
# covariates X[, -1]. The outlier mean keeps the formula above.
#   "normal"   N(0, 1), the paper.
#   "t<df>"    standard t with df degrees of freedom, from Cauchy (t1)
#              towards the normal. Not rescaled, because its variance is
#              infinite for df <= 2.
#   "skewt<df>" Azzalini skew-t with df degrees of freedom and slant
#              SKEWT_ALPHA, standardised to mean 0 and variance 1. Its
#              fourth moment is finite for df > 4.
#   "gamma<k>" gamma with shape k, standardised to mean 0 and variance 1.

# Slant of the skew-t covariate.
SKEWT_ALPHA <- 5

# Splits a covariate label such as "t3" into its distribution and parameter.
parse_x_dist <- function(x_dist) {
  pattern <- "^(normal|t|skewt|gamma)([0-9.]*)$"
  m <- regmatches(x_dist, regexec(pattern, x_dist))[[1]]
  if (length(m) == 0 || (m[2] == "normal") != (m[3] == "")) {
    stop("unknown covariate distribution: ", x_dist)
  }
  list(dist = m[2], par = as.numeric(m[3]))
}

# Covariates of distribution x_dist from the standard normals z and the
# uniforms u. The skew-t draws one more standard normal per covariate.
covariate_draw <- function(z, u, x_dist) {
  d <- parse_x_dist(x_dist)
  switch(d$dist,
    t = z / sqrt(qchisq(u, d$par) / d$par),
    skewt = {
      nu <- d$par
      if (nu <= 2) stop("skewt needs df > 2 for a finite variance")
      delta <- SKEWT_ALPHA / sqrt(1 + SKEWT_ALPHA^2)
      z0 <- rnorm(length(z))
      x <- (delta * abs(z0) + sqrt(1 - delta^2) * z) /
        sqrt(qchisq(u, nu) / nu)
      mu <- delta * sqrt(nu / pi) * exp(lgamma((nu - 1) / 2) - lgamma(nu / 2))
      (x - mu) / sqrt(nu / (nu - 2) - mu^2)
    },
    gamma = (qgamma(u, d$par) - d$par) / sqrt(d$par)
  )
}

# Draws n_sets data sets of size n_obs_max and keeps the first n_obs rows,
# so that cells with different n_obs share the same random numbers.
# Returns X as an [n_obs, n_dim, n_sets] array and y as an [n_obs, n_sets]
# matrix.
#
# z and eps are drawn as for the normal covariate. A covariate other than
# the normal then draws u, and the skew-t draws z0 last. So the normal cells
# keep their random numbers, and all other cells share the same z, eps and
# u: they differ only in x_dist.
make_data_gaussian <- function(n_sets, n_obs, n_obs_max, beta, out_dev,
                      n_obs_out = 1, sigma2_d = 1, x_dist = "normal") {
  n_dim <- length(beta)
  x <- array(rnorm(n_obs_max * (n_dim - 1) * n_sets),
             c(n_obs_max, n_dim - 1, n_sets))
  eps <- matrix(rnorm(n_obs_max * n_sets, sd = sqrt(sigma2_d)),
                n_obs_max, n_sets)
  if (x_dist != "normal") {
    x[] <- covariate_draw(x, runif(length(x)), x_dist)
  }
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

# Binomial variant: y_i ~ Bernoulli(plogis(X_i beta)). The outlier analogue
# of mu_d shifts the linear predictor of the first n_obs_out observations
# away from their own sign, so that a large out_dev makes those
# observations almost surely the class the covariates argue against.
# As in make_data_gaussian, n_obs_max rows are drawn and the first n_obs
# are kept, so that cells with different n_obs share the same random numbers.
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

# Poisson variant: y_i ~ Poisson(exp(eta_i)). The rate grows as exp(eta),
# so the coefficients of make_beta_coeff() are scaled down and an intercept is
# added: with POIS_SCALE and POIS_INTERCEPT the median rate is e and the
# counts stay small enough for a count grid. beta_t = 1 is still a clear
# difference between the models, as in the other families.
#
# The outlier adds out_dev to the linear predictor of the first n_obs_out
# observations, so out_dev is on the log scale. The shift is always
# upwards: a downward shift gives counts near 0, which are common at a rate
# near e, so those outliers would be much weaker than the upward ones.
#
# The counts are drawn by inversion from n_obs_max uniforms, as in
# make_data_binomial, so that cells with different n_obs share them.
POIS_SCALE <- 0.5
POIS_INTERCEPT <- 1

make_data_poisson <- function(n_sets, n_obs, n_obs_max, beta, out_dev,
                              n_obs_out = 1) {
  n_dim <- length(beta)
  b <- beta * POIS_SCALE
  b[1] <- POIS_INTERCEPT
  x <- array(rnorm(n_obs_max * (n_dim - 1) * n_sets),
             c(n_obs_max, n_dim - 1, n_sets))
  u <- matrix(runif(n_obs_max * n_sets), n_obs_max, n_sets)

  keep <- seq_len(n_obs)
  X <- array(1, c(n_obs, n_dim, n_sets))
  X[, -1, ] <- x[keep, , , drop = FALSE]
  eta <- matrix(0, n_obs, n_sets)
  for (k in seq_len(n_dim)) {
    eta <- eta + b[k] * matrix(X[, k, ], n_obs, n_sets)
  }
  out <- seq_len(min(n_obs_out, n_obs))
  eta[out, ] <- eta[out, , drop = FALSE] + out_dev
  y <- matrix(qpois(u[keep, , drop = FALSE], exp(eta)), n_obs, n_sets)
  list(X = X, y = y)
}

# Draws the data of one cell, for any family. Only the gaussian family has
# other covariates than the normal, so x_dist must be "normal" for the others.
make_data_family <- function(family, ..., x_dist = "normal") {
  if (family != "gaussian" && x_dist != "normal") {
    stop("x_dist is for the gaussian family only, not ", family)
  }
  switch(family,
    gaussian = make_data_gaussian(..., x_dist = x_dist),
    binomial = make_data_binomial(...),
    poisson = make_data_poisson(...),
    stop("unknown family: ", family)
  )
}

# sign() with ties broken towards +1, so that the outlier shift never
# vanishes.
sign_pos <- function(x) ifelse(x < 0, -1, 1)

# beta = (intercept, other effects..., beta_t); beta_t is only in model B.
make_beta_coeff <- function(beta_t, n_dim = 3, beta_other = 1,
                            beta_intercept = 0) {
  beta <- c(rep(beta_other, n_dim - 1), beta_t)
  beta[1] <- beta_intercept
  beta
}

# Stacks all data sets into one design matrix and one response vector.
pool_sets <- function(data) {
  dims <- dim(data$X)
  X <- aperm(data$X, c(1, 3, 2))
  dim(X) <- c(dims[1] * dims[3], dims[2])
  list(X = X, y = as.vector(data$y))
}

# Logistic regression with a normal prior N(0, prior_sd^2 I) on the
# coefficients, its Laplace posterior, and the exact LOO predictive.
#
# The normal linear models of Sivula et al. (2025) use a flat prior, but a
# flat prior leaves the logistic posterior improper under separation, which
# is common at the small n of the grid. prior_sd is therefore weakly
# informative rather than flat.
#
# A binary predictive is list(p), the posterior predictive probability of
# y = 1, so that lib/scores.R can score it without knowing how p was made.

PRIOR_SD <- 5

# Penalised IRLS. Returns the posterior mode and the Laplace covariance
# V = (X'WX + P)^-1, with P = I / prior_sd^2.
fit_logistic <- function(X, y, prior_sd = PRIOR_SD, start = NULL,
                         tol = 1e-10, max_iter = 100) {
  d <- ncol(X)
  P <- diag(1 / prior_sd^2, d)
  beta <- if (is.null(start)) numeric(d) else start
  for (iter in seq_len(max_iter)) {
    eta <- drop(X %*% beta)
    p <- plogis(eta)
    w <- pmax(p * (1 - p), 1e-10)
    H <- crossprod(X * w, X) + P
    step <- solve(H, crossprod(X, y - p) - P %*% beta)
    beta <- beta + drop(step)
    if (max(abs(step)) < tol) break
  }
  eta <- drop(X %*% beta)
  p <- plogis(eta)
  w <- pmax(p * (1 - p), 1e-10)
  H <- crossprod(X * w, X) + P
  list(beta_hat = beta, V = chol2inv(chol(H)), iter = iter)
}

# Step and half range of the quadrature, in standard deviations of eta.
QUAD_STEP <- 0.5
QUAD_RANGE <- 9
# Largest quadrature matrix held at one time, in elements.
QUAD_MAX_CELLS <- 2e6

# E[plogis(eta)] for eta ~ N(m, s^2), for vectors m and s. The predictive
# depends on the coefficients only through eta, so this one dimensional
# integral is the whole posterior predictive.
#
# The rule is the trapezoid rule over the standardised eta, which is
# spectrally accurate for an analytic integrand. Two terms set the error.
# The normal density needs a step below its own width, which the standard
# scale fixes at 1. plogis(m + s t) has its nearest poles at a distance
# pi / s from the real axis, and the error falls as exp(-2 pi^2 / (s h)),
# so the step shrinks as s grows.
predictive_prob <- function(m, s) {
  h <- QUAD_STEP / max(1, max(s))
  t <- seq(-QUAD_RANGE, QUAD_RANGE, by = h)
  w <- h * dnorm(t)
  # The rule needs a length(m) by length(t) matrix, which is too large for
  # the pooled test set of a big cell, so it runs in blocks.
  block <- max(1, floor(QUAD_MAX_CELLS / length(t)))
  p <- numeric(length(m))
  for (from in seq(1, length(m), by = block)) {
    idx <- seq(from, min(from + block - 1, length(m)))
    p[idx] <- drop(plogis(outer(s[idx], t) + m[idx]) %*% w)
  }
  p
}

# Exact LOO predictive: refit without each observation in turn, warm
# started at the full-data mode.
loo_predictive_binary <- function(X, y, prior_sd = PRIOR_SD, start = NULL) {
  n <- nrow(X)
  if (is.null(start)) start <- fit_logistic(X, y, prior_sd)$beta_hat
  m <- numeric(n)
  s <- numeric(n)
  for (i in seq_len(n)) {
    fit_i <- fit_logistic(X[-i, , drop = FALSE], y[-i], prior_sd, start = start)
    x_i <- X[i, ]
    m[i] <- sum(x_i * fit_i$beta_hat)
    s[i] <- sqrt(sum(x_i * drop(fit_i$V %*% x_i)))
  }
  list(p = predictive_prob(m, s))
}

# Fit to all of (X, y): the Laplace posterior, and the LOO predictive.
fit_model_binary <- function(X, y, prior_sd = PRIOR_SD) {
  fit <- fit_logistic(X, y, prior_sd)
  fit$prior_sd <- prior_sd
  fit$loo <- loo_predictive_binary(X, y, prior_sd, start = fit$beta_hat)
  fit
}

# Predictive of a full-data fit for new points X_new.
predict_new_binary <- function(fit, X_new) {
  m <- drop(X_new %*% fit$beta_hat)
  s <- sqrt(rowSums((X_new %*% fit$V) * X_new))
  list(p = predictive_prob(m, s))
}


# ---------------------------------------------------------------------------
# Poisson regression with the same normal prior, its Laplace posterior, and
# the exact LOO predictive.
#
# A count predictive is list(m, s): the rate is lambda = exp(eta) with
# eta ~ N(m, s^2). lib/scores.R turns it into a pmf over the counts, so it
# never needs to know how m and s were made.

# Largest linear predictor kept, so that exp(eta) cannot overflow.
ETA_MAX <- 30

# Penalised IRLS with the log link. Returns the posterior mode and the
# Laplace covariance V = (X'WX + P)^-1, with W = diag(exp(eta)) and
# P = I / prior_sd^2.
fit_poisson <- function(X, y, prior_sd = PRIOR_SD, start = NULL,
                        tol = 1e-10, max_iter = 100) {
  d <- ncol(X)
  P <- diag(1 / prior_sd^2, d)
  beta <- if (is.null(start)) numeric(d) else start
  mu <- NULL
  for (iter in seq_len(max_iter)) {
    mu <- exp(pmin(drop(X %*% beta), ETA_MAX))
    w <- pmax(mu, 1e-10)
    H <- crossprod(X * w, X) + P
    step <- solve(H, crossprod(X, y - mu) - P %*% beta)
    beta <- beta + drop(step)
    if (max(abs(step)) < tol) break
  }
  mu <- exp(pmin(drop(X %*% beta), ETA_MAX))
  H <- crossprod(X * pmax(mu, 1e-10), X) + P
  list(beta_hat = beta, V = chol2inv(chol(H)), iter = iter)
}

# Half range of the quadrature over eta, in standard deviations. The rate
# grows as exp(eta), so a range as wide as the logistic one would push the
# count grid out of reach. A range of 5 leaves 6e-7 of the normal mass out.
QUAD_RANGE_POIS <- 5

# Largest count in the predictive grid, and the mixture tail the grid must
# leave outside. The tail cannot be smaller than the mass that the range
# above already drops.
POIS_K_MAX <- 20000
POIS_TAIL <- 1e-6

# Trapezoid nodes and weights over the standardised eta. As a function of
# t, dpois(k, exp(m + s t)) peaks at t = (log k - m) / s with a width of
# about 1 / (s sqrt(k)). The step must stay below that width, so it shrinks
# with s and with `kref`, the largest count the rule has to resolve.
pois_nodes <- function(s, kref) {
  h <- QUAD_STEP / max(1, max(s) * sqrt(1 + kref))
  t <- seq(-QUAD_RANGE_POIS, QUAD_RANGE_POIS, by = h)
  w <- stats::dnorm(t)
  list(t = t, w = w / sum(w))
}

# Largest count the grid needs, over all of m and s. The envelope
# predictive, with the largest m and the largest s, bounds every row, and
# its tail is halved by doubling k until it falls below POIS_TAIL.
pois_k_max <- function(m, s) {
  nd <- pois_nodes(max(s), 0)
  lambda <- exp(pmin(max(m) + max(s) * nd$t, ETA_MAX))
  k <- 8
  while (k < POIS_K_MAX &&
         sum(nd$w * stats::ppois(k, lambda, lower.tail = FALSE)) > POIS_TAIL) {
    k <- 2 * k
  }
  min(k, POIS_K_MAX)
}

# Posterior predictive pmf of the counts 0:k_max, one row per element of m,
# renormalised so that the cut at k_max cannot bias a score.
predictive_pmf <- function(m, s, k_max = pois_k_max(m, s),
                           nodes = pois_nodes(s, k_max)) {
  k <- 0:k_max
  lg <- lgamma(k + 1)
  pmf <- matrix(0, length(m), length(k))
  for (j in seq_along(nodes$t)) {
    lambda <- exp(pmin(m + s * nodes$t[j], ETA_MAX))
    pmf <- pmf + nodes$w[j] * exp(outer(log(lambda), k) - lambda -
                                    rep(lg, each = length(m)))
  }
  list(k = k, pmf = pmf / rowSums(pmf))
}

# Predictive probability of the observed count of each row. This needs no
# count grid, so it stays exact however heavy the tail of the predictive is.
predictive_count_dens <- function(y, m, s, nodes = pois_nodes(s, max(y))) {
  p <- numeric(length(y))
  for (j in seq_along(nodes$t)) {
    p <- p + nodes$w[j] * stats::dpois(y, exp(pmin(m + s * nodes$t[j],
                                                   ETA_MAX)))
  }
  p
}

# Exact LOO predictive: refit without each observation in turn, warm
# started at the full-data mode.
loo_predictive_poisson <- function(X, y, prior_sd = PRIOR_SD, start = NULL) {
  n <- nrow(X)
  if (is.null(start)) start <- fit_poisson(X, y, prior_sd)$beta_hat
  m <- numeric(n)
  s <- numeric(n)
  for (i in seq_len(n)) {
    fit_i <- fit_poisson(X[-i, , drop = FALSE], y[-i], prior_sd, start = start)
    x_i <- X[i, ]
    m[i] <- sum(x_i * fit_i$beta_hat)
    s[i] <- sqrt(sum(x_i * drop(fit_i$V %*% x_i)))
  }
  list(m = m, s = s)
}

# Fit to all of (X, y): the Laplace posterior, and the LOO predictive.
fit_model_poisson <- function(X, y, prior_sd = PRIOR_SD) {
  fit <- fit_poisson(X, y, prior_sd)
  fit$prior_sd <- prior_sd
  fit$loo <- loo_predictive_poisson(X, y, prior_sd, start = fit$beta_hat)
  fit
}

# Predictive of a full-data fit for new points X_new.
predict_new_poisson <- function(fit, X_new) {
  list(m = drop(X_new %*% fit$beta_hat),
       s = sqrt(rowSums((X_new %*% fit$V) * X_new)))
}

# Predictive of a full-data fit, for either generalised linear family.
predict_new_family <- function(family, fit, X_new) {
  switch(family,
    gaussian = predict_new(fit, X_new),
    binomial = predict_new_binary(fit, X_new),
    poisson = predict_new_poisson(fit, X_new),
    stop("unknown family: ", family)
  )
}

# Fit of one training set, for any family. tau2 is used by the normal model
# only.
fit_model_family <- function(family, X, y, tau2 = NULL) {
  switch(family,
    gaussian = fit_model(X, y, tau2),
    binomial = fit_model_binary(X, y),
    poisson = fit_model_poisson(X, y),
    stop("unknown family: ", family)
  )
}

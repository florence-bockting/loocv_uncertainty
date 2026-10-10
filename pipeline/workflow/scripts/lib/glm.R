# Posterior predictive of the normal linear model with unknown noise
# variance tau2 and a uniform prior on (beta, log(tau2)).
# Mirrors calc_loo_ti and calc_elpd_tl in simulated/general_setup.py.
# A predictive is list(mu, scale, df), a Student-t location-scale predictive.

# Model A drops the last covariate, model B uses all covariates.
model_columns <- function(model, n_dim) {
  switch(model,
    A = seq_len(n_dim - 1),
    B = seq_len(n_dim),
    stop("unknown model: ", model)
  )
}

# Exact LOO predictive from one fit, with the hat-matrix identities
#   y_i - mu_{-i}              = e_i / (1 - h_i)
#   x_i' (X_{-i}'X_{-i})^-1 x_i = h_i / (1 - h_i)
#   RSS_{-i}                   = RSS - e_i^2 / (1 - h_i).
loo_predictive_gaussian <- function(X, y) {
  n <- nrow(X)
  d <- ncol(X)
  qr_X <- qr(X)
  h <- rowSums(qr.Q(qr_X)^2)
  e <- qr.resid(qr_X, y)
  mu <- y - e / (1 - h)
  s2 <- (sum(e^2) - e^2 / (1 - h)) / (n - 1 - d)
  list(mu = mu, scale = sqrt(s2 / (1 - h)), df = n - 1 - d)
}

# Fit to all of (X, y): the posterior of beta.
fit_model_gaussian <- function(X, y) {
  n <- nrow(X)
  d <- ncol(X)
  V <- chol2inv(chol(crossprod(X)))
  beta_hat <- drop(V %*% crossprod(X, y))
  s2 <- sum((y - X %*% beta_hat)^2) / (n - d)
  list(beta_hat = beta_hat, V = V, s2 = s2, df = n - d)
}

# Predictive of a full-data fit for new points X_new.
predict_new_gaussian <- function(fit, X_new) {
  mu <- drop(X_new %*% fit$beta_hat)
  xVx <- rowSums((X_new %*% fit$V) * X_new)
  list(mu = mu, scale = sqrt((xVx + 1) * fit$s2), df = fit$df)
}


# ---------------------------------------------------------------------------
# Logistic regression with a normal prior N(0, prior_sd^2 I) on the
# coefficients, its Laplace posterior, and the exact LOO predictive.
#
# The normal linear models of Sivula et al. (2025) use a flat prior, but a
# flat prior leaves the logistic posterior improper under separation, which
# is common at the small n of the grid. prior_sd is therefore weakly
# informative rather than flat.
#
# A binary predictive is list(p), the posterior predictive probability of
# y = 1, so that lib/measures.R can score it without knowing how p was made.

PRIOR_SD <- 5

# Largest linear predictor kept, so that exp(eta) cannot overflow.
ETA_MAX <- 30

# Mean mu and IRLS weight w at the linear predictor eta, for the logit link
# (binomial) and the log link (poisson).
mean_weight <- function(eta, family) {
  switch(family,
    binomial = {
      mu <- plogis(eta)
      list(mu = mu, w = pmax(mu * (1 - mu), 1e-10))
    },
    poisson = {
      mu <- exp(pmin(eta, ETA_MAX))
      list(mu = mu, w = pmax(mu, 1e-10))
    },
    stop("unknown family: ", family)
  )
}

# Log likelihood at the linear predictor eta. The poisson rate is not
# capped here, so that a step into a huge eta gets a low value.
log_lik <- function(eta, y, family) {
  switch(family,
    binomial = sum(y * eta + plogis(-eta, log.p = TRUE)),
    poisson = sum(y * eta - exp(eta)),
    stop("unknown family: ", family)
  )
}

# Largest number of step halvings in one Newton iteration.
MAX_HALVING <- 30

# Laplace approximation of the posterior by penalised IRLS. Returns the
# posterior mode and the covariance V = (X'WX + P)^-1, with
# W = diag(w) and P = I / prior_sd^2. Without convergence it warns and
# returns NA, so that the trial is dropped and the job does not fail.
laplace_approx <- function(X, y, family, prior_sd = PRIOR_SD, start = NULL,
                           tol = 1e-10, max_iter = 100) {
  d <- ncol(X)
  P <- diag(1 / prior_sd^2, d)
  log_post <- function(beta) {
    log_lik(drop(X %*% beta), y, family) - sum(beta^2) / (2 * prior_sd^2)
  }
  beta <- if (is.null(start) || anyNA(start)) numeric(d) else start
  lp <- log_post(beta)
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    mw <- mean_weight(drop(X %*% beta), family)
    H <- crossprod(X * mw$w, X) + P
    step <- drop(solve(H, crossprod(X, y - mw$mu) - P %*% beta))
    if (max(abs(step)) < tol) {
      beta <- beta + step
      converged <- TRUE
      break
    }
    # A full step can overshoot far when one count is large, so the step
    # is halved until the log posterior does not decrease. A decrease of
    # the size of the rounding error does not count.
    lp_new <- log_post(beta + step)
    halvings <- 0
    while (lp_new < lp - 1e-10 * (1 + abs(lp)) && halvings < MAX_HALVING) {
      step <- step / 2
      lp_new <- log_post(beta + step)
      halvings <- halvings + 1
    }
    beta <- beta + step
    lp <- lp_new
  }
  if (!converged) {
    warning("laplace_approx: no convergence in ", max_iter, " iterations")
    return(list(beta_hat = rep(NA_real_, d), V = matrix(NA_real_, d, d),
                iter = iter))
  }
  mw <- mean_weight(drop(X %*% beta), family)
  H <- crossprod(X * mw$w, X) + P
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
  # A fit that did not converge gives NA; its trial is dropped.
  if (anyNA(m) || anyNA(s)) return(rep(NA_real_, length(m)))
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
loo_predictive_binomial <- function(X, y, prior_sd = PRIOR_SD, start = NULL) {
  n <- nrow(X)
  if (is.null(start)) {
    start <- laplace_approx(X, y, "binomial", prior_sd)$beta_hat
  }
  m <- numeric(n)
  s <- numeric(n)
  for (i in seq_len(n)) {
    fit_i <- laplace_approx(X[-i, , drop = FALSE], y[-i], "binomial", prior_sd,
                            start = start)
    x_i <- X[i, ]
    m[i] <- sum(x_i * fit_i$beta_hat)
    s[i] <- sqrt(sum(x_i * drop(fit_i$V %*% x_i)))
  }
  list(p = predictive_prob(m, s))
}

# Fit to all of (X, y): the Laplace posterior.
fit_model_binomial <- function(X, y, prior_sd = PRIOR_SD) {
  fit <- laplace_approx(X, y, "binomial", prior_sd)
  fit$prior_sd <- prior_sd
  fit
}

# Predictive of a full-data fit for new points X_new.
predict_new_binomial <- function(fit, X_new) {
  m <- drop(X_new %*% fit$beta_hat)
  s <- sqrt(rowSums((X_new %*% fit$V) * X_new))
  list(p = predictive_prob(m, s))
}


# ---------------------------------------------------------------------------
# Poisson regression with the same normal prior, its Laplace posterior, and
# the exact LOO predictive.
#
# A count predictive is list(m, s): the rate is lambda = exp(eta) with
# eta ~ N(m, s^2). lib/measures.R turns it into a pmf over the counts, so it
# never needs to know how m and s were made.

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
  if (is.null(start)) {
    start <- laplace_approx(X, y, "poisson", prior_sd)$beta_hat
  }
  m <- numeric(n)
  s <- numeric(n)
  for (i in seq_len(n)) {
    fit_i <- laplace_approx(X[-i, , drop = FALSE], y[-i], "poisson", prior_sd,
                            start = start)
    x_i <- X[i, ]
    m[i] <- sum(x_i * fit_i$beta_hat)
    s[i] <- sqrt(sum(x_i * drop(fit_i$V %*% x_i)))
  }
  list(m = m, s = s)
}

# Fit to all of (X, y): the Laplace posterior.
fit_model_poisson <- function(X, y, prior_sd = PRIOR_SD) {
  fit <- laplace_approx(X, y, "poisson", prior_sd)
  fit$prior_sd <- prior_sd
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
    gaussian = predict_new_gaussian(fit, X_new),
    binomial = predict_new_binomial(fit, X_new),
    poisson = predict_new_poisson(fit, X_new),
    stop("unknown family: ", family)
  )
}

# Fit of one training set, for any family.
fit_model_family <- function(family, X, y) {
  switch(family,
    gaussian = fit_model_gaussian(X, y),
    binomial = fit_model_binomial(X, y),
    poisson = fit_model_poisson(X, y),
    stop("unknown family: ", family)
  )
}

# LOO predictive of one training set, for any family. start warm starts the
# GLM refits, usually at the full-data mode; the Gaussian ignores it.
loo_predictive_family <- function(family, X, y, start = NULL) {
  switch(family,
    gaussian = loo_predictive_gaussian(X, y),
    binomial = loo_predictive_binomial(X, y, start = start),
    poisson = loo_predictive_poisson(X, y, start = start),
    stop("unknown family: ", family)
  )
}

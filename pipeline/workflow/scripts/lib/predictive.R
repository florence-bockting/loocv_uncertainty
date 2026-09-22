# Posterior predictive of the normal linear model with a uniform prior on
# beta, and on log(tau2) when tau2 is unknown (tau2 = NULL).
# Mirrors calc_loo_ti and calc_elpd_tl in simulated/general_setup.py.
# A predictive is list(mu, scale, df); df = Inf is a normal predictive,
# otherwise a Student-t location-scale predictive.

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
loo_predictive <- function(X, y, tau2 = NULL) {
  n <- nrow(X)
  d <- ncol(X)
  qr_X <- qr(X)
  h <- rowSums(qr.Q(qr_X)^2)
  e <- qr.resid(qr_X, y)
  mu <- y - e / (1 - h)
  if (is.null(tau2)) {
    s2 <- (sum(e^2) - e^2 / (1 - h)) / (n - 1 - d)
    df <- n - 1 - d
  } else {
    s2 <- tau2
    df <- Inf
  }
  list(mu = mu, scale = sqrt(s2 / (1 - h)), df = df)
}

# Fit to all of (X, y): the posterior of beta, and the LOO predictive.
fit_model <- function(X, y, tau2 = NULL) {
  n <- nrow(X)
  d <- ncol(X)
  V <- chol2inv(chol(crossprod(X)))
  beta_hat <- drop(V %*% crossprod(X, y))
  if (is.null(tau2)) {
    s2 <- sum((y - X %*% beta_hat)^2) / (n - d)
    df <- n - d
  } else {
    s2 <- tau2
    df <- Inf
  }
  list(beta_hat = beta_hat, V = V, s2 = s2, df = df,
       loo = loo_predictive(X, y, tau2))
}

# Predictive of a full-data fit for new points X_new.
predict_new <- function(fit, X_new) {
  mu <- drop(X_new %*% fit$beta_hat)
  xVx <- rowSums((X_new %*% fit$V) * X_new)
  list(mu = mu, scale = sqrt((xVx + 1) * fit$s2), df = fit$df)
}

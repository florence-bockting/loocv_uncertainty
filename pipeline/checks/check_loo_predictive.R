# Checks that loo_predictive_gaussian (one fit, hat-matrix identities) gives the
# same LOO predictive as n explicit refits, as in calc_loo_ti in
# simulated/general_setup.py.
source("pipeline/workflow/scripts/lib/dgp.R")
source("pipeline/workflow/scripts/lib/glm.R")

# Explicit refit without observation i, predicting observation i.
loo_refit <- function(X, y) {
  n <- nrow(X)
  d <- ncol(X)
  mu <- scale <- numeric(n)
  for (i in seq_len(n)) {
    X_i <- X[-i, , drop = FALSE]
    y_i <- y[-i]
    V <- solve(crossprod(X_i))
    beta_hat <- V %*% crossprod(X_i, y_i)
    x <- X[i, ]
    mu[i] <- sum(x * beta_hat)
    s2 <- sum((y_i - X_i %*% beta_hat)^2) / (n - 1 - d)
    scale[i] <- sqrt((drop(x %*% V %*% x) + 1) * s2)
  }
  list(mu = mu, scale = scale)
}

set.seed(2)
for (n in c(16, 128)) {
  d <- make_data_gaussian(1, n, n, make_beta_coeff(0.5), 20)
  X <- d$X[, , 1]
  y <- d$y[, 1]
  fast <- loo_predictive_gaussian(X, y)
  slow <- loo_refit(X, y)
  cat(sprintf("n = %4d  max |diff| mu = %.1e, scale = %.1e\n", n,
              max(abs(fast$mu - slow$mu)),
              max(abs(fast$scale - slow$scale))))
}

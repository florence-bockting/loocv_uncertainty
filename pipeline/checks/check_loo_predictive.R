# Checks that loo_predictive (one fit, hat-matrix identities) gives the
# same LOO predictive as n explicit refits, as in calc_loo_ti in
# simulated/general_setup.py.
source("pipeline/workflow/scripts/lib/dgp.R")
source("pipeline/workflow/scripts/lib/predictive.R")

# Explicit refit without observation i, predicting observation i.
loo_refit <- function(X, y, tau2 = NULL) {
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
    s2 <- if (is.null(tau2)) {
      sum((y_i - X_i %*% beta_hat)^2) / (n - 1 - d)
    } else {
      tau2
    }
    scale[i] <- sqrt((drop(x %*% V %*% x) + 1) * s2)
  }
  list(mu = mu, scale = scale)
}

set.seed(2)
for (n in c(16, 128)) {
  d <- make_data(1, n, n, make_beta(0.5), 20)
  X <- d$X[, , 1]
  y <- d$y[, 1]
  for (tau2 in list(NULL, 1)) {
    fast <- loo_predictive(X, y, tau2)
    slow <- loo_refit(X, y, tau2)
    cat(sprintf("n = %4d, tau2 = %-4s  max |diff| mu = %.1e, scale = %.1e\n",
                n, if (is.null(tau2)) "NULL" else tau2,
                max(abs(fast$mu - slow$mu)),
                max(abs(fast$scale - slow$scale))))
  }
}

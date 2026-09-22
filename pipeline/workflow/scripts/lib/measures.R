# Model differences Ma - Mb on the utility scale, computed from the
# pointwise scores. Equation numbers refer to
# Predictive_measures_and_model_comparison_in_the_loo_package.pdf.

# The score that each measure is computed from.
MEASURE_SCORE <- c(elpd = "elpd", mse = "sqe", rmse = "sqe", r2 = "sqe",
                   mae = "ae", rps = "rps", srps = "srps")

# LOO estimate of Ma - Mb and its SE, from pointwise LOO scores.
# elpd is on the sum scale (2.12, 2.14), the other additive measures on the
# mean scale (2.11, 2.13). rmse uses (3.3), r2 uses (3.4, 3.5).
loo_difference <- function(measure, loo_a, loo_b, y) {
  n <- length(y)
  mean_se <- function(d) c(mean(d), sqrt(var(d) / n))
  a <- loo_a[[MEASURE_SCORE[[measure]]]]
  b <- loo_b[[MEASURE_SCORE[[measure]]]]
  switch(measure,
    elpd = n * mean_se(a - b),
    mse = , mae = , rps = mean_se(b - a),
    srps = mean_se(a - b),
    rmse = {
      rmse_a <- sqrt(mean(a))
      rmse_b <- sqrt(mean(b))
      z <- a / rmse_a - b / rmse_b
      c(rmse_b - rmse_a, sqrt(var(z) / n) / 2)
    },
    r2 = {
      sq_dev <- (y - mean(y))^2
      mse_y <- mean(sq_dev)
      delta <- a - b
      c_delta <- mean(delta) / mse_y
      c(-c_delta, sqrt(var(delta - c_delta * sq_dev) / n) / mse_y)
    },
    stop("unknown measure: ", measure)
  )
}

# Target Ma - Mb from the test-set means of the scores (1.5). For rmse and
# r2 the function g is applied to the means. elpd is scaled to the sum over
# n_obs points, as elpd_t in simulated/run_all.py.
target_difference <- function(measure, test_a, test_b, n_obs, y_test_var) {
  a <- test_a[[MEASURE_SCORE[[measure]]]]
  b <- test_b[[MEASURE_SCORE[[measure]]]]
  switch(measure,
    elpd = n_obs * (a - b),
    mse = , mae = , rps = b - a,
    srps = a - b,
    rmse = sqrt(b) - sqrt(a),
    r2 = (b - a) / y_test_var,
    stop("unknown measure: ", measure)
  )
}

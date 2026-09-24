# Model differences Ma - Mb on the utility scale, computed from the
# pointwise scores. Equation numbers refer to
# Predictive_measures_and_model_comparison_in_the_loo_package.pdf.

# The score that each measure is computed from.
MEASURE_SCORE <- c(elpd = "elpd", mlpd = "elpd", ic = "elpd",
                   mse = "sqe", rmse = "sqe", r2 = "sqe", brier = "sqe",
                   mae = "ae", rps = "rps", srps = "srps",
                   acc = "acc", bacc = "acc")

# Measures that need a binary outcome. The others apply to both families.
MEASURES_BINARY <- c("brier", "acc", "bacc")

# The measures of `measures` that the family supports.
measures_for_family <- function(measures, family) {
  if (family == "binomial") measures else setdiff(measures, MEASURES_BINARY)
}

# An SE of zero carries no information: the z-score and the PIT of the
# trial are then undefined. acc and bacc reach it whenever both models
# classify every observation alike, which happens for about a tenth of the
# trials when the models are similar, at every n. Return NA for the SE and
# keep the estimate, so that the trial is dropped, not silently counted.
finite_se <- function(est, se) {
  if (!is.finite(se) || se == 0) c(est, NA_real_) else c(est, se)
}

# LOO estimate of Ma - Mb and its SE, from pointwise LOO scores.
# elpd and ic are on the sum scale (2.12, 2.14), the other additive
# measures on the mean scale (2.11, 2.13). rmse uses (3.3), r2 uses
# (3.4, 3.5). bacc averages two class-wise means, so its SE adds the
# variances of the two strata, as in loo:::.se_diff_bacc().
loo_difference <- function(measure, loo_a, loo_b, y) {
  n <- length(y)
  mean_se <- function(d) finite_se(mean(d), sqrt(var(d) / n))
  a <- loo_a[[MEASURE_SCORE[[measure]]]]
  b <- loo_b[[MEASURE_SCORE[[measure]]]]
  switch(measure,
    elpd = n * mean_se(a - b),
    mlpd = mean_se(a - b),
    ic = 2 * n * mean_se(a - b),
    brier = mean_se(b - a),
    acc = mean_se(a - b),
    bacc = stratified_mean_se(a - b, y),
    mse = , mae = , rps = mean_se(b - a),
    srps = mean_se(a - b),
    rmse = {
      rmse_a <- sqrt(mean(a))
      rmse_b <- sqrt(mean(b))
      z <- a / rmse_a - b / rmse_b
      finite_se(rmse_b - rmse_a, sqrt(var(z) / n) / 2)
    },
    r2 = {
      sq_dev <- (y - mean(y))^2
      mse_y <- mean(sq_dev)
      delta <- a - b
      c_delta <- mean(delta) / mse_y
      finite_se(-c_delta, sqrt(var(delta - c_delta * sq_dev) / n) / mse_y)
    },
    stop("unknown measure: ", measure)
  )
}

# Mean of the class-wise means of d, and its SE. The classes are disjoint
# sets of observations, so their variances add; within a class both models
# score the same observations, so d is already paired. Both are NA unless
# every class holds two observations or more. A class that the training set
# does not contain changes the estimand: bacc is then the accuracy of the
# one class present, but the target still averages the classes of the test
# set. A class of one observation contributes no variance, so the SE is too
# small. Drop these trials instead.
stratified_mean_se <- function(d, y, n_classes = 2) {
  strata <- split(d, y)
  if (length(strata) != n_classes || any(lengths(strata) < 2)) {
    return(c(NA_real_, NA_real_))
  }
  est <- mean(vapply(strata, mean, numeric(1)))
  var_c <- vapply(strata, function(x) var(x) / length(x), numeric(1))
  finite_se(est, sqrt(sum(var_c)) / n_classes)
}

# Target Ma - Mb from the test-set means of the scores (1.5). For rmse and
# r2 the function g is applied to the means. elpd and ic are scaled to the
# sum over n_obs points, as elpd_t in simulated/run_all.py. bacc uses the
# class-wise test means in `by_class`.
target_difference <- function(measure, test_a, test_b, n_obs, y_test_var,
                              by_class_a = NULL, by_class_b = NULL) {
  if (measure == "bacc") {
    a <- by_class_a[[MEASURE_SCORE[[measure]]]]
    b <- by_class_b[[MEASURE_SCORE[[measure]]]]
    return(mean(a) - mean(b))
  }
  a <- test_a[[MEASURE_SCORE[[measure]]]]
  b <- test_b[[MEASURE_SCORE[[measure]]]]
  switch(measure,
    elpd = n_obs * (a - b),
    mlpd = a - b,
    ic = 2 * n_obs * (a - b),
    acc = a - b,
    mse = , mae = , rps = , brier = b - a,
    srps = a - b,
    rmse = sqrt(b) - sqrt(a),
    r2 = (b - a) / y_test_var,
    stop("unknown measure: ", measure)
  )
}

# Sign of the pointwise term of an additive measure, relative to the score
# difference a - b, and the factor 2 of ic. The mean of the terms is the LOO
# estimate divided by measure_scale(). rmse and r2 have no pointwise terms:
# g is applied to the mean, not to each point. bacc has one mean per class.
POINTWISE_SIGN <- c(elpd = 1, mlpd = 1, ic = 2, srps = 1, acc = 1,
                    mse = -1, mae = -1, rps = -1, brier = -1)

# The measures of `measures` that have pointwise terms.
measures_pointwise <- function(measures) {
  intersect(measures, names(POINTWISE_SIGN))
}

# Pointwise terms of a measure, one per observation.
pointwise_terms <- function(measure, loo_a, loo_b) {
  s <- MEASURE_SCORE[[measure]]
  POINTWISE_SIGN[[measure]] * (loo_a[[s]] - loo_b[[s]])
}

# Factor between the scale of the measure and the per-observation scale.
# elpd and ic are sums over the n observations, the others are means.
measure_scale <- function(measure, n_obs) {
  if (measure %in% c("elpd", "ic")) n_obs else 1
}

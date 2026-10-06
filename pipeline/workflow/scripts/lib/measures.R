# Pointwise scores and the model differences Ma - Mb built from them.
# Equation numbers refer to
# Predictive_measures_and_model_comparison_in_the_loo_package.pdf.

# ---------------------------------------------------------------------------
# Pointwise scores of a Student-t predictive, on their
# natural scale: elpd and srps are utilities; sqe, ae and rps are losses.
# ---------------------------------------------------------------------------

SCORES <- c("elpd", "sqe", "ae", "rps", "srps")

# The point prediction is the predictive mean (2.29). rps and srps use the
# closed forms of E|X - y| and E|X - X'| (2.20, 2.21).
pointwise_scores <- function(y, pred, scores = SCORES) {
  z <- (y - pred$mu) / pred$scale
  nu <- pred$df
  out <- list()
  if ("elpd" %in% scores) {
    out$elpd <- dt(z, nu, log = TRUE) - log(pred$scale)
  }
  if ("sqe" %in% scores) out$sqe <- (y - pred$mu)^2
  if ("ae" %in% scores) out$ae <- abs(y - pred$mu)
  if (any(c("rps", "srps") %in% scores)) {
    m <- abs_moments(z, nu)
    abs_dev <- pred$scale * m$abs_dev
    spread <- pred$scale * m$spread
    if ("rps" %in% scores) out$rps <- abs_dev - spread / 2
    if ("srps" %in% scores) out$srps <- -abs_dev / spread - log(spread) / 2
  }
  out[scores]
}

# E|Z - z| and E|Z - Z'| for Z, Z' iid t(nu).
abs_moments <- function(z, nu) {
  list(abs_dev = z * (2 * pt(z, nu) - 1) +
         2 * dt(z, nu) * (nu + z^2) / (nu - 1),
       spread = 4 * sqrt(nu) / (nu - 1) *
         exp(lbeta(0.5, nu - 0.5) - 2 * lbeta(0.5, nu / 2)))
}

# Scores that a binomial cell adds to SCORES.
SCORES_BINOMIAL <- c(SCORES, "acc")

# Pointwise scores of a Bernoulli predictive, list(p). The point prediction
# is again the predictive mean, here p itself, so sqe is the Brier score and
# ae is |y - p|. For X, X' ~ Bernoulli(p) the two absolute moments are
# E|X - y| = y (1 - p) + (1 - y) p and E|X - X'| = 2 p (1 - p).
pointwise_scores_binary <- function(y, pred, scores = SCORES_BINOMIAL) {
  p <- pred$p
  out <- list()
  if ("elpd" %in% scores) {
    out$elpd <- ifelse(y == 1, log(p), log1p(-p))
  }
  if ("sqe" %in% scores) out$sqe <- (y - p)^2
  if ("ae" %in% scores) out$ae <- abs(y - p)
  if (any(c("rps", "srps") %in% scores)) {
    abs_dev <- ifelse(y == 1, 1 - p, p)
    spread <- 2 * p * (1 - p)
    if ("rps" %in% scores) out$rps <- abs_dev - spread / 2
    if ("srps" %in% scores) out$srps <- -abs_dev / spread - log(spread) / 2
  }
  # The classification is the majority class of the predictive; a tie at
  # p = 0.5 counts as class 1, as in loo::measure_acc().
  if ("acc" %in% scores) out$acc <- as.numeric((p >= 0.5) == (y == 1))
  out[scores]
}

# The scores of `scores` that the family supports. Only the binomial family
# classifies, so only it scores the accuracy.
scores_for_family <- function(scores, family) {
  if (family == "binomial") scores else setdiff(scores, "acc")
}

# Pointwise scores of any family. A predictive with NA, from a fit that did
# not converge, gives NA for every score.
pointwise_scores_family <- function(family, y, pred, scores) {
  if (anyNA(unlist(pred))) {
    return(sapply(scores, function(s) rep(NA_real_, length(y)),
                  simplify = FALSE))
  }
  switch(family,
    gaussian = pointwise_scores(y, pred, scores),
    binomial = pointwise_scores_binary(y, pred, scores),
    poisson = pointwise_scores_count(y, pred, scores),
    stop("unknown family: ", family)
  )
}


# Scores of a count predictive, list(m, s), as lib/glm.R builds it. The
# definitions are the ones of pointwise_scores(): the point prediction is
# the predictive mean, rps is E|X - y| - E|X - X'| / 2, and srps its scaled
# form. On an ordered count support that rps is the ranked probability
# score, which is why a count family makes rps a measure of its own rather
# than a copy of the Brier score.
#
# Only rps and srps need the pmf of the counts. The elpd uses the
# predictive density of the observed count alone, and the point prediction
# is the lognormal mean exp(m + s^2 / 2), so neither of them depends on
# where the count grid ends.
SCORES_POISSON <- SCORES

# Largest pmf matrix held at one time, in elements.
POIS_MAX_CELLS <- 2e6

pointwise_scores_count <- function(y, pred, scores = SCORES_POISSON) {
  out <- list()
  if ("elpd" %in% scores) {
    out$elpd <- log(predictive_count_dens(y, pred$m, pred$s))
  }
  # The rate is lognormal, so the predictive mean is exact.
  mu <- exp(pred$m + pred$s^2 / 2)
  if ("sqe" %in% scores) out$sqe <- (y - mu)^2
  if ("ae" %in% scores) out$ae <- abs(y - mu)
  if (any(c("rps", "srps") %in% scores)) {
    m <- count_abs_moments(y, pred)
    if ("rps" %in% scores) out$rps <- m$abs_dev - m$spread / 2
    if ("srps" %in% scores) {
      out$srps <- -m$abs_dev / m$spread - log(m$spread) / 2
    }
  }
  out[scores]
}

# E|X - y| and E|X - X'| for X, X' from the predictive of each observation.
# The spread uses the identity E|X - X'| = 2 sum_k F_k (1 - F_k), which
# holds for a distribution on the integers. The pmf is built in blocks, so
# that a heavy-tailed predictive cannot exhaust the memory.
count_abs_moments <- function(y, pred) {
  n <- length(y)
  k_max <- pois_k_max(pred$m, pred$s)
  nodes <- pois_nodes(pred$s, k_max)
  block <- max(1, floor(POIS_MAX_CELLS / (k_max + 1)))
  abs_dev <- numeric(n)
  spread <- numeric(n)
  for (from in seq(1, n, by = block)) {
    idx <- seq(from, min(from + block - 1, n))
    p <- predictive_pmf(pred$m[idx], pred$s[idx], k_max, nodes)
    cdf <- p$pmf
    for (j in seq_along(p$k)[-1]) cdf[, j] <- cdf[, j - 1] + cdf[, j]
    abs_dev[idx] <- rowSums(p$pmf * abs(outer(y[idx], p$k, "-")))
    spread[idx] <- 2 * rowSums(cdf * (1 - cdf))
  }
  list(abs_dev = abs_dev, spread = spread)
}

# ---------------------------------------------------------------------------
# Model differences Ma - Mb on the utility scale, computed from the
# pointwise scores.
# ---------------------------------------------------------------------------

# Each measure as in loo's .measure_spec (R/pred_measure-builtin.R, branch
# pred_measure of the loo package): the score it is computed from, whether
# the measure is a loss, and how a difference Ma - Mb is aggregated.
# "sum" and "mean" aggregate the paired pointwise differences: elpd and ic
# are on the sum scale (2.12, 2.14), the others on the mean scale (2.11,
# 2.13). A "measure_specific" measure is not a mean of pointwise terms: it
# has an estimate for one model, and SE_DIFF_FUNS gives the SE of the
# difference. weight turns the score into the pointwise value of the
# measure: ic is -2 elpd, as elpd_transform in loo.
MEASURE_SPEC <- list(
  elpd = list(score = "elpd", loss = FALSE, diff_method = "sum"),
  mlpd = list(score = "elpd", loss = FALSE, diff_method = "mean"),
  ic = list(score = "elpd", loss = TRUE, diff_method = "sum", weight = -2),
  mse = list(score = "sqe", loss = TRUE, diff_method = "mean"),
  rmse = list(score = "sqe", loss = TRUE, diff_method = "measure_specific",
              estimate = function(s, y) sqrt(mean(s)),
              se_diff_fun = "rmse"),
  r2 = list(score = "sqe", loss = FALSE, diff_method = "measure_specific",
            estimate = function(s, y) 1 - mean(s) / mean((y - mean(y))^2),
            se_diff_fun = "r2"),
  brier = list(score = "sqe", loss = TRUE, diff_method = "mean"),
  mae = list(score = "ae", loss = TRUE, diff_method = "mean"),
  rps = list(score = "rps", loss = TRUE, diff_method = "mean"),
  srps = list(score = "srps", loss = FALSE, diff_method = "mean"),
  acc = list(score = "acc", loss = FALSE, diff_method = "mean"),
  bacc = list(score = "acc", loss = FALSE, diff_method = "measure_specific",
              estimate = function(s, y) {
                if (!bacc_strata_ok(y)) return(NA_real_)
                mean(vapply(split(s, y), mean, numeric(1)))
              },
              se_diff_fun = "bacc")
)

# The score that each measure is computed from.
MEASURE_SCORE <- vapply(MEASURE_SPEC, `[[`, character(1), "score")

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
# loo reports an SE of zero instead.
finite_se <- function(est, se) {
  if (!is.finite(se) || se == 0) c(est, NA_real_) else c(est, se)
}

# LOO estimate of Ma - Mb on the utility scale and its SE, from pointwise
# LOO scores, as loo:::.pair_measure_stats() with Ma as cmp and Mb as ref.
loo_difference <- function(measure, loo_a, loo_b, y) {
  spec <- MEASURE_SPEC[[measure]]
  if (is.null(spec)) stop("unknown measure: ", measure)
  if (spec$diff_method == "measure_specific") {
    a <- loo_a[[spec$score]]
    b <- loo_b[[spec$score]]
    sign <- if (spec$loss) -1 else 1
    est <- sign * (spec$estimate(a, y) - spec$estimate(b, y))
    return(finite_se(est, SE_DIFF_FUNS[[spec$se_diff_fun]](a, b, y)))
  }
  d <- pointwise_terms(measure, loo_a, loo_b)
  n <- length(d)
  est_se <- finite_se(mean(d), sqrt(var(d) / n))
  if (spec$diff_method == "sum") n * est_se else est_se
}

# The SE of a difference of measure_specific measures, from the scores a
# and b of the two models, copied from loo (R/pred_measure-builtin.R).

# rmse: delta method from the MSE scale (3.3), as loo:::.se_diff_rmse().
# The paired contrast z gives the two variances and the covariance in one
# variance, so two models with the same scores cancel exactly.
se_diff_rmse <- function(a, b, y) {
  n <- length(a)
  mse_a <- mean(a)
  mse_b <- mean(b)
  if (n <= 1L || mse_a <= 0 || mse_b <= 0) return(0)
  z <- a / sqrt(mse_a) - b / sqrt(mse_b)
  0.5 * sqrt(var(z) / n)
}

# r2: delta method (3.4, 3.5), as loo:::.se_diff_r2(). The difference is
# -mean(a - b) / MSE(y), so it is the single-model expansion with the
# squared errors replaced by their pointwise differences.
se_diff_r2 <- function(a, b, y) {
  se_r2_delta(a - b, (y - mean(y))^2)
}

# loo:::.se_r2_delta().
se_r2_delta <- function(sqe, mse_y_i) {
  mse_y_hat <- mean(mse_y_i)
  scaled <- sqe - (mean(sqe) / mse_y_hat) * mse_y_i
  sqrt(var(scaled) / length(sqe)) / mse_y_hat
}

# bacc: a difference of stratified means, as loo:::.se_diff_bacc(). The
# classes are disjoint sets of observations, so their variances add; within
# a class both models score the same observations, so a - b is paired.
se_diff_bacc <- function(a, b, y) {
  if (!bacc_strata_ok(y)) return(NA_real_)
  var_c <- vapply(split(a - b, y), function(d) var(d) / length(d),
                  numeric(1))
  sqrt(sum(var_c)) / length(var_c)
}

SE_DIFF_FUNS <- list(rmse = se_diff_rmse, r2 = se_diff_r2,
                     bacc = se_diff_bacc)

# Whether y supports bacc: every class holds two observations or more. A
# class that the training set does not contain changes the estimand: bacc
# is then the accuracy of the one class present, but the target still
# averages the classes of the test set. A class of one observation
# contributes no variance, so the SE is too small. bacc and its SE are NA
# for these trials, so that they are dropped. loo keeps them instead.
bacc_strata_ok <- function(y, n_classes = 2) {
  n_c <- table(y)
  length(n_c) == n_classes && all(n_c >= 2)
}

# The measures of `measures` that have pointwise terms. rmse, r2 and bacc
# have none: g is applied to the means, not to each point, and bacc has
# one mean per class.
measures_pointwise <- function(measures) {
  is_pointwise <- vapply(MEASURE_SPEC[measures], function(s) {
    s$diff_method != "measure_specific"
  }, logical(1))
  measures[is_pointwise]
}

# Pointwise terms of a measure, one per observation, on the utility scale.
# Their mean is the LOO estimate divided by measure_scale().
pointwise_terms <- function(measure, loo_a, loo_b) {
  spec <- MEASURE_SPEC[[measure]]
  weight <- if (is.null(spec$weight)) 1 else spec$weight
  sign <- if (spec$loss) -weight else weight
  sign * (loo_a[[spec$score]] - loo_b[[spec$score]])
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

# Factor between the scale of the measure and the per-observation scale.
# elpd and ic are sums over the n observations, the others are means.
measure_scale <- function(measure, n_obs) {
  if (measure %in% c("elpd", "ic")) n_obs else 1
}

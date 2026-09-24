# Pointwise scores of a normal (df = Inf) or Student-t predictive, on their
# natural scale: elpd and srps are utilities; sqe, ae and rps are losses.
# Equation numbers refer to
# Predictive_measures_and_model_comparison_in_the_loo_package.pdf.

SCORES <- c("elpd", "sqe", "ae", "rps", "srps")

# The point prediction is the predictive mean (2.29). rps and srps use the
# closed forms of E|X - y| and E|X - X'| (2.20, 2.21).
pointwise_scores <- function(y, pred, scores = SCORES) {
  z <- (y - pred$mu) / pred$scale
  nu <- pred$df
  out <- list()
  if ("elpd" %in% scores) {
    log_dens <- if (is.infinite(nu)) dnorm(z, log = TRUE) else dt(z, nu, log = TRUE)
    out$elpd <- log_dens - log(pred$scale)
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

# E|Z - z| and E|Z - Z'| for Z, Z' iid standard normal or t(nu).
abs_moments <- function(z, nu) {
  if (is.infinite(nu)) {
    list(abs_dev = z * (2 * pnorm(z) - 1) + 2 * dnorm(z),
         spread = 2 / sqrt(pi))
  } else {
    list(abs_dev = z * (2 * pt(z, nu) - 1) +
           2 * dt(z, nu) * (nu + z^2) / (nu - 1),
         spread = 4 * sqrt(nu) / (nu - 1) *
           exp(lbeta(0.5, nu - 0.5) - 2 * lbeta(0.5, nu / 2)))
  }
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

# Pointwise scores of any family.
pointwise_scores_family <- function(family, ...) {
  switch(family,
    gaussian = pointwise_scores(...),
    binomial = pointwise_scores_binary(...),
    poisson = pointwise_scores_count(...),
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

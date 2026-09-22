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

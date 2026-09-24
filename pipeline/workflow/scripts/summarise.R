# Rule summarise: moments of the estimate, the target and the error over
# the trials of each cell and measure, with Bayesian-bootstrap intervals.
# Mirrors bb_mean_sd_skew in simulated/general_setup.py (Fig. 9).
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])

trials <- readRDS(snakemake@input[[1]])
n_bb <- snakemake@params[["n_bb"]]
set.seed(snakemake@params[["seed"]])

# Dirichlet(1, ..., 1) weights, one row per bootstrap draw.
dirichlet_weights <- function(n_bb, n) {
  w <- matrix(rexp(n_bb * n), n_bb, n)
  w / rowSums(w)
}

# Weighted mean/SD and skewness for each bootstrap draw.
bb_moments <- function(x, w) {
  sum_w2 <- rowSums(w^2)
  sum_w3 <- rowSums(w^3)
  m <- drop(w %*% x)
  m2 <- drop(w %*% x^2) - m^2
  m3 <- drop(w %*% x^3) - 3 * m * drop(w %*% x^2) + 2 * m^3
  sd <- sqrt(m2 / (1 - sum_w2))
  skew <- m3 / sd^3 / (1 - 3 * sum_w2 + 2 * sum_w3)
  list(`mean/SD` = m / sd, skewness = skew)
}

# A trial without a finite SE has no z-score and no PIT, so it cannot enter
# the moments. acc and bacc produce such trials; see finite_se() in
# lib/measures.R. Drop them here and count them for each cell and measure.
usable <- function(g) {
  g[is.finite(g$estimate) & is.finite(g$se) & is.finite(g$target), ]
}

# Fewer usable trials than this leave the bootstrap moments meaningless.
MIN_TRIALS <- 10

groups <- split(trials, list(trials$cell, trials$measure), drop = TRUE)
weights <- list()
rows <- lapply(groups, function(g) {
  label <- paste(g$cell[1], g$measure[1], sep = " / ")
  n_all <- nrow(g)
  g <- usable(g)
  n_dropped <- n_all - nrow(g)
  if (n_dropped > 0) {
    message(sprintf("%s: dropped %d of %d trials without a finite SE",
                    label, n_dropped, n_all))
  }
  if (nrow(g) < MIN_TRIALS) {
    message(sprintf("%s: skipped, only %d usable trials", label, nrow(g)))
    return(NULL)
  }
  key <- as.character(nrow(g))
  if (is.null(weights[[key]])) weights[[key]] <<- dirichlet_weights(n_bb, nrow(g))
  w <- weights[[key]]
  quantities <- list(estimate = g$estimate, target = g$target,
                     error = g$estimate - g$target)
  do.call(rbind, lapply(names(quantities), function(q) {
    mom <- bb_moments(quantities[[q]], w)
    do.call(rbind, lapply(names(mom), function(s) {
      qs <- stats::quantile(mom[[s]], c(0.025, 0.5, 0.975), na.rm = TRUE)
      data.frame(cell = g$cell[1], family = g$family[1],
                 n_obs = g$n_obs[1], beta_t = g$beta_t[1],
                 out_dev = g$out_dev[1], tau2 = g$tau2[1],
                 measure = g$measure[1], quantity = q, stat = s,
                 q025 = qs[[1]], q500 = qs[[2]], q975 = qs[[3]],
                 n_trial = nrow(g), n_dropped = n_dropped)
    }))
  }))
})
res <- do.call(rbind, rows)
rownames(res) <- NULL
saveRDS(res, snakemake@output[[1]])

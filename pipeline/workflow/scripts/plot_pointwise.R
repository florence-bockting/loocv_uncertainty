# Rule plot_pointwise: distribution of the pointwise utility terms of one
# cell, with the expected utility (target) and the LOO estimate as vertical
# lines. One panel per measure. Page 1 shows one trial, page 2 all trials.
# Everything is per observation: elpd and ic are divided by n_obs.
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])
snakemake@source("lib/plot.R")
snakemake@source("lib/measures.R")

a <- lapply(snakemake@input[["a"]], readRDS)
b <- lapply(snakemake@input[["b"]], readRDS)
cmp <- do.call(rbind, lapply(snakemake@input[["compare"]], readRDS))
cell <- snakemake@params[["cell"]]
cell_name <- snakemake@wildcards[["cell"]]
measures <- measures_pointwise(
  measures_for_family(unlist(snakemake@params[["measures"]]), cell$family))

trials <- unlist(lapply(a, function(x) x$trials))
n_obs <- cell$n_obs
# Pointwise terms of every trial: one n_obs x n_trial matrix per measure.
terms <- sapply(measures, function(m) {
  do.call(cbind, lapply(seq_along(a), function(k) {
    pointwise_terms(m, a[[k]]$loo, b[[k]]$loo)
  }))
}, simplify = FALSE)

# The terms and the compare output must agree: the mean of the terms is the
# LOO estimate on the per-observation scale.
for (m in measures) {
  rows <- cmp[cmp$measure == m, ]
  rows <- rows[match(trials, rows$trial), ]
  est <- rows$estimate / measure_scale(m, n_obs)
  stopifnot(max(abs(colMeans(terms[[m]]) - est)) < 1e-8)
}

# Estimate, target and SE of one trial, or the mean over the trials, per
# obs. The SE is the one of the LOO estimate, sd(terms) / sqrt(n_obs).
# A dropped acc trial has no SE, so the means skip the missing values.
per_obs <- function(m, trial = NULL) {
  rows <- cmp[cmp$measure == m, ]
  if (!is.null(trial)) rows <- rows[rows$trial == trial, ]
  s <- measure_scale(m, n_obs)
  c(estimate = mean(rows$estimate, na.rm = TRUE) / s,
    target = mean(rows$target, na.rm = TRUE) / s,
    se = mean(rows$se, na.rm = TRUE) / s)
}

# Number of quantile dots per panel, as `quantiles` in arviz.plot_dot.
N_DOT <- 50

# The x limits of a measure: a robust window around the median of all its
# terms, plus the two lines. The tails of the outlier cells reach far, so a
# 1 %-99 % window would leave every dot in one stack. Both pages use the
# same limits, so the two pages are comparable.
limits <- sapply(measures, function(m) {
  x <- terms[[m]]
  spread <- stats::IQR(x)
  if (spread == 0) spread <- diff(stats::quantile(x, c(0.05, 0.95)))
  at <- c("estimate", "target")
  v <- range(c(stats::median(x) + c(-4, 4) * spread,
               per_obs(m)[at], per_obs(m, trials[1])[at]))
  v + c(-1, 1) * 0.04 * diff(v)
}, simplify = FALSE)

# N_DOT dots of equal probability mass, as arviz.plot_dot. Fewer terms than
# dots: every term keeps its own dot.
quantile_dots <- function(x) {
  if (length(x) <= N_DOT) sort(x)
  else stats::quantile(x, (seq_len(N_DOT) - 0.5) / N_DOT, names = FALSE)
}

# Wilkinson stacking of sorted dots: a dot opens a new stack when it lies
# more than one bin width from the first dot of the current stack.
wilkinson <- function(q, bw) {
  stack <- integer(length(q))
  k <- 1L
  first <- q[1]
  for (i in seq_along(q)) {
    if (q[i] - first > bw) {
      k <- k + 1L
      first <- q[i]
    }
    stack[i] <- k
  }
  stack
}

# The largest bin width whose tallest stack still fits the panel. A dot is
# as wide as the bin, so a wider bin also makes the dots taller.
fit_binwidth <- function(q, lims, pin) {
  for (bw in diff(lims) / seq(6, 200, by = 2)) {
    height <- max(tabulate(wilkinson(q, bw)))
    if (height * bw * pin[1] / diff(lims) <= pin[2]) return(bw)
  }
  diff(lims) / 200
}

# The dots of a measure: the position and the stack height of each dot.
dots <- function(m, x) {
  lims <- limits[[m]]
  # A dot outside the limits stays on the border, as a ring.
  q <- quantile_dots(x)
  out <- q < lims[1] | q > lims[2]
  q <- pmin(pmax(q, lims[1]), lims[2])
  bw <- fit_binwidth(q, lims, par("pin"))
  stack <- wilkinson(q, bw)
  list(at = stats::ave(q, stack, FUN = function(v) mean(range(v))),
       height = stats::ave(stack, stack, FUN = seq_along),
       out = out, bw = bw)
}

# Dot diameter in inches, for the cex of pch 16.
DOT_INCH <- function(cex) 0.75 * par("cin")[1] * cex

DOT_COL <- grDevices::adjustcolor("grey45", alpha.f = 0.35)
RING_COL <- grDevices::adjustcolor(MPL[4], alpha.f = 0.7)

panel <- function(m, d, cex, lines, label, se_name) {
  lims <- limits[[m]]
  # One dot row per y unit, so that the dots touch but do not overlap.
  rows <- par("pin")[2] / DOT_INCH(cex)
  plot(NA, xlim = lims, ylim = c(0, rows), xlab = label, ylab = "",
       yaxt = "n", yaxs = "i")
  points(d$at, d$height - 0.5, pch = ifelse(d$out, 1, 16), cex = cex,
         col = ifelse(d$out, RING_COL, DOT_COL))
  abline(v = lines[["target"]], col = MPL[3], lwd = 2)
  abline(v = lines[["estimate"]], col = MPL[2], lwd = 2, lty = 2)
  mtext(m, side = 3, line = 0.3, font = 2, cex = 0.9, adj = 0)
  mtext(sprintf("%s = %.3g", se_name, lines[["se"]]), side = 3, line = 0.3,
        cex = 0.75, adj = 1)
}

page <- function(trial, title) {
  ncol <- ceiling(sqrt(length(measures)))
  nrow <- ceiling(length(measures) / ncol)
  par(mfrow = c(nrow, ncol), mar = c(3.2, 1.2, 1.8, 1.2), oma = c(0, 0, 3.5, 0),
      mgp = c(1.9, 0.5, 0), tcl = -0.25, bty = "n")
  # Every panel uses the dot size of the panel with the narrowest bins.
  d <- sapply(measures, function(m) {
    x <- if (is.null(trial)) as.vector(terms[[m]]) else terms[[m]][, trials == trial]
    dots(m, x)
  }, simplify = FALSE)
  cex <- min(vapply(measures, function(m) {
    d[[m]]$bw * par("pin")[1] / diff(limits[[m]]) / (0.75 * par("cin")[1])
  }, numeric(1)))
  for (m in measures) {
    panel(m, d[[m]], cex, per_obs(m, trial), "pointwise term (per obs.)",
          if (is.null(trial)) "mean SE" else "SE")
  }
  # The legend goes into the free space of the grid, or on the last panel.
  if (nrow * ncol > length(measures)) plot.new()
  n_term <- if (is.null(trial)) length(terms[[1]]) else n_obs
  legend("center", bty = "n", lwd = c(2, 2, NA, NA), lty = c(1, 2, NA, NA),
         pch = c(NA, NA, 16, 1), col = c(MPL[c(3, 2)], DOT_COL, RING_COL),
         legend = c("expected utility (target)", "LOO estimate",
                    if (n_term <= N_DOT) "dot = 1 pointwise term"
                    else sprintf("dot = %g %% of the terms", 100 / N_DOT),
                    "dot outside the limits"))
  mtext(title, side = 3, outer = TRUE, line = 1.2, font = 2)
}

open_pdf(snakemake@output[[1]], 720, 576)
first <- trials[1]
page(first, sprintf("%s: pointwise terms of trial %d (n = %d)",
                    cell_name, first, n_obs))
page(NULL, sprintf("%s: pointwise terms of %d trials pooled (n = %d)",
                   cell_name, length(trials), n_obs))
dev.off()

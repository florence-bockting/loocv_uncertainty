# Rule plot_overview: one page per family and measure, with the six
# diagnostics of the paper side by side. Cells without an outlier only.
#
#   1 calibration   PIT of the normal approximation, minus the uniform
#   2 variance      the stated SE over the actual SD of the error
#   3 bias          mean of the error over its SD
#   4 skewness      skewness of the error
#   5 joint         the estimate against the target, hardest cell
#   6 error shape   the error over its own SD, against the standard normal
#
# Panels 1, 2, 5 and 6 read the trials, panels 3 and 4 the moments. None of
# them assumes that the measure decomposes into pointwise terms.
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])
snakemake@source("lib/plot.R")

trials <- readRDS(snakemake@input[["trials"]])
mom <- readRDS(snakemake@input[["moments"]])
measure <- snakemake@wildcards[["measure"]]
family <- snakemake@wildcards[["family"]]
sel <- snakemake@params[["sel"]]
n_sel <- unlist_num(sel$n_obs)
beta_sel <- unlist_num(sel$beta_t)
labels <- measure_labels(measure)

tr <- trials[trials$measure == measure & trials$family == family &
               trials$out_dev == 0 & trials$tau2 == sel$tau2, ]
mo <- mom[mom$measure == measure & mom$family == family &
            mom$out_dev == 0 & mom$tau2 == sel$tau2, ]
stopifnot(nrow(tr) > 0)
n_all <- sort(unique(tr$n_obs))
n_dropped <- sum(!is.finite(tr$se) | !is.finite(tr$estimate) |
                   !is.finite(tr$target))

# The trials of one cell of the no-outlier grid.
cell_of <- function(n, b) {
  tr[tr$n_obs == n & abs(tr$beta_t - b) < 1e-12, ]
}

# Line type per n, colour per beta_t.
lty_n <- function(n) match(n, n_sel)

panel <- function(main) title(main = main, font.main = 1, cex.main = 1)

open_pdf(snakemake@output[[1]], 864, 540)
par(mfrow = c(2, 3), mar = c(3.4, 3.6, 2.2, 1), oma = c(0, 0, 2.6, 0),
    mgp = c(2.1, 0.5, 0), tcl = -0.25, bty = "l")

# 1. calibration ---------------------------------------------------------
# The PIT of a calibrated estimate is uniform. The band is the 95 % KS band
# for the largest cell, so a line outside it is not a sampling artefact.
plot(NA, xlim = c(0, 1), ylim = c(-0.35, 0.35), xlab = "PIT",
     ylab = "ECDF - uniform")
abline(h = 0, col = "grey60")
grid_p <- seq(0, 1, length.out = 201)
for (b_i in seq_along(beta_sel)) {
  for (n in n_sel) {
    p <- pit(cell_of(n, beta_sel[b_i]))
    p <- p[!is.na(p)]
    if (length(p) < 10) next
    lines(grid_p, stats::ecdf(p)(grid_p) - grid_p, col = MPL[b_i],
          lty = lty_n(n), lwd = 1.3)
  }
}
band <- 1.358 / sqrt(max(table(tr$n_obs, tr$beta_t)))
polygon(c(0, 1, 1, 0), c(band, band, -band, -band),
        col = grDevices::adjustcolor("grey40", 0.15), border = NA)
panel("calibration of the stated SE")

# 2. variance ratio ------------------------------------------------------
plot(NA, xlim = range(n_all), ylim = c(0, 1.6), log = "x", xaxt = "n",
     xlab = "n", ylab = expression(widehat(SE)[LOO] / SE(err[LOO])))
axis(1, at = n_all, labels = n_all)
abline(h = 1, col = adjust_lightness("red", 1.3))
for (b_i in seq_along(beta_sel)) {
  q <- vapply(n_all, function(n) {
    rows <- cell_of(n, beta_sel[b_i])
    r <- rows$se / stats::sd(loo_error(rows), na.rm = TRUE)
    stats::quantile(r, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  }, numeric(3))
  polygon(c(n_all, rev(n_all)), c(q[1, ], rev(q[3, ])),
          col = grDevices::adjustcolor(MPL[b_i], 0.18), border = NA)
  lines(n_all, q[2, ], col = MPL[b_i], lwd = 1.6)
}
panel("is the uncertainty understated?")

# 3 and 4. moments of the error ------------------------------------------
for (s in c("mean/SD", "skewness")) {
  d <- mo[mo$quantity == "error" & mo$stat == s & mo$beta_t %in% beta_sel, ]
  ylim <- if (s == "mean/SD") c(-0.35, 0.35) else
    range(c(d$q025, d$q975, 0), finite = TRUE)
  plot(NA, xlim = range(n_all), ylim = ylim, log = "x", xaxt = "n",
       xlab = "n", ylab = if (s == "mean/SD") "mean(err) / SD(err)" else
         "skewness of err")
  axis(1, at = n_all, labels = n_all)
  abline(h = 0, col = "grey60")
  for (b_i in seq_along(beta_sel)) {
    db <- d[abs(d$beta_t - beta_sel[b_i]) < 1e-12, ]
    db <- db[order(db$n_obs), ]
    if (nrow(db) == 0) next
    polygon(c(db$n_obs, rev(db$n_obs)), c(db$q025, rev(db$q975)),
            col = grDevices::adjustcolor(MPL[b_i], 0.18), border = NA)
    lines(db$n_obs, db$q500, col = MPL[b_i], lwd = 1.6)
  }
  panel(if (s == "mean/SD") "bias of the estimate" else
          "skewness of the error")
}

# 5. joint, the hardest cell ---------------------------------------------
# Smallest n and the smallest beta_t: the two models predict alike, so the
# comparison is least informative.
rows <- cell_of(min(n_sel), min(beta_sel))
# Each axis gets its own zoom, as in plot_joint.R: the estimate spreads far
# wider than the target, so a common range hides the cloud.
zoom <- function(v) {
  lims <- stats::quantile(v, c(0.01, 0.99), na.rm = TRUE)
  lims + c(-1, 1) * 0.25 * diff(lims)
}
xlim <- zoom(rows$estimate)
ylim <- zoom(rows$target)
plot(rows$estimate, rows$target, xlim = xlim, ylim = ylim, pch = 16,
     cex = 0.7, col = grDevices::adjustcolor(MPL[1], 0.45),
     xlab = labels$estimate, ylab = labels$target)
abline(h = 0, v = 0, col = "grey80")
lo <- max(xlim[1], ylim[1])
hi <- min(xlim[2], ylim[2])
if (lo < hi) segments(lo, lo, hi, hi, col = MPL[3], lwd = 1.5)
panel(sprintf("n = %d, beta = %g", min(n_sel), min(beta_sel)))

# 6. shape of the error --------------------------------------------------
# Against the standard normal, which the normal approximation assumes.
plot(NA, xlim = c(-4, 4), ylim = c(0, 0.55), xlab = "err / SD(err)",
     ylab = "density")
# The normal sits behind the densities, as a wide solid grey reference.
x <- seq(-4, 4, length.out = 200)
lines(x, stats::dnorm(x), col = grDevices::adjustcolor("grey45", 0.55),
      lwd = 4)
for (b_i in seq_along(beta_sel)) {
  for (n in c(min(n_sel), max(n_sel))) {
    rows <- cell_of(n, beta_sel[b_i])
    e <- loo_error(rows)
    e <- e[is.finite(e)] / stats::sd(e, na.rm = TRUE)
    if (length(e) < 10) next
    lines(stats::density(e), col = MPL[b_i], lty = lty_n(n), lwd = 1.3)
  }
}
legend("topleft", bty = "n", cex = 0.75, ncol = 2,
       legend = c(sprintf("beta = %g", beta_sel), sprintf("n = %d", n_sel)),
       col = c(MPL[seq_along(beta_sel)], rep("grey40", length(n_sel))),
       lty = c(rep(1, length(beta_sel)), lty_n(n_sel)), lwd = 1.3)
panel("shape of the error")

mtext(sprintf("%s, %s, no outlier%s", measure, family,
              if (n_dropped > 0) sprintf(", %d trials dropped", n_dropped)
              else ""),
      side = 3, outer = TRUE, line = 0.6, font = 2, cex = 1.1)
dev.off()

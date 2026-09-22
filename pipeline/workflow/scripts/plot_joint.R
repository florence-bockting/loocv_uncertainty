# Rule plot_joint: joint distribution of the LOO estimate and the target,
# as Fig. 5 (simulated/plot_all__joint.py): 2D histogram in greys, KDE
# contours in copper, the diagonal in green and the zero lines in grey.
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])
snakemake@source("lib/plot.R")

trials <- readRDS(snakemake@input[[1]])
measure <- snakemake@wildcards[["measure"]]
sel <- snakemake@params[["sel"]]
n_sel <- unlist_num(sel$n_obs)
beta_sel <- unlist_num(sel$beta_t)
labels <- measure_labels(measure)

# Zoom as in the script: 1 %-99 % range plus 25 % on each side.
zoom <- function(v) {
  lims <- stats::quantile(v, c(0.01, 0.99))
  lims + c(-1, 1) * 0.25 * diff(lims)
}

open_pdf(snakemake@output[[1]], 648, 576)
par(mfrow = c(length(beta_sel), length(n_sel)), mar = c(3, 3, 1, 0.8),
    oma = c(0, 7, 3, 0), mgp = c(1.8, 0.5, 0), tcl = -0.25, bty = "l")
for (b_i in seq_along(beta_sel)) {
  for (n_i in seq_along(n_sel)) {
    rows <- cell_rows(trials, measure, n_sel[n_i], beta_sel[b_i], 0, sel$tau2)
    xlim <- zoom(rows$estimate)
    ylim <- zoom(rows$target)
    keep <- rows$estimate > xlim[1] & rows$estimate < xlim[2] &
      rows$target > ylim[1] & rows$target < ylim[2]
    x <- rows$estimate[keep]
    y <- rows$target[keep]

    plot(NA, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i",
         xlab = if (b_i == length(beta_sel)) labels$estimate else "",
         ylab = if (n_i == 1) labels$target else "")
    if (xlim[1] < 0 && xlim[2] > 0) abline(v = 0, col = "grey80", lwd = 0.8)
    if (ylim[1] < 0 && ylim[2] > 0) abline(h = 0, col = "grey80", lwd = 0.8)

    # 2D histogram in place of hexbin (gridsize 45), empty cells white.
    bx <- seq(xlim[1], xlim[2], length.out = 46)
    by <- seq(ylim[1], ylim[2], length.out = 46)
    counts <- table(cut(x, bx), cut(y, by))
    counts[counts == 0] <- NA
    image(bx, by, matrix(counts, 45, 45), col = greys(64), add = TRUE)

    # scipy gaussian_kde(bw_method = 0.4): kernel SD = 0.4 * data SD.
    # MASS::kde2d uses h = 4 * kernel SD. Contours at fractions of the range.
    kde <- MASS::kde2d(x, y, h = 4 * 0.4 * c(stats::sd(x), stats::sd(y)),
                       n = 100, lims = c(xlim, ylim))
    levels <- min(kde$z) + diff(range(kde$z)) * c(0.1, 0.3, 0.6, 0.9)
    contour(kde, levels = levels, col = copper(4), drawlabels = FALSE,
            add = TRUE)

    lo <- max(xlim[1], ylim[1])
    hi <- min(xlim[2], ylim[2])
    if (lo < hi) segments(lo, lo, hi, hi, col = MPL[3], lwd = 1.5)

    if (b_i == 1) mtext(bquote(n == .(n_sel[n_i])), side = 3, line = 1)
    if (n_i == 1) mtext(bquote(beta[Delta] == .(beta_sel[b_i])), side = 2,
                        line = 4.5, las = 1, adj = 1)
  }
}
dev.off()

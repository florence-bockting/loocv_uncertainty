# Rule plot_joint_grouped: the joint distribution of plot_joint.R as points,
# coloured by a statistic of the training design (rule design_stats).
#   var: quantile of mean(x2^2) within the cell
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])
snakemake@source("lib/plot.R")

trials <- readRDS(snakemake@input[["trials"]])
trials <- select_x_dist(trials, snakemake@params[["x_dist"]])
trials <- select_tau2(trials, snakemake@params[["tau2"]])
trials <- merge(trials, readRDS(snakemake@input[["design"]]),
                by = c("cell", "trial"))
measure <- snakemake@wildcards[["measure"]]
family <- snakemake@wildcards[["family"]]
sel <- snakemake@params[["sel"]]
n_sel <- unlist_num(sel$n_obs)
beta_sel <- unlist_num(sel$beta_t)
labels <- measure_labels(measure)

# Yellow (low bin) to blue (high bin): reversed viridis without its purple end.
yl_bu <- function(k) grDevices::hcl.colors(k + 1, "viridis", rev = TRUE)[1:k]
g <- switch(snakemake@wildcards[["group"]],
  var = list(col = "mx2", probs = c(0, 0.5, 0.8, 0.95, 1),
             labels = c("< 50 %", "50-80 %", "80-95 %", "> 95 %"),
             title = expression("quantile of training mean(" * x[2]^2 *
                                  ") within the cell"),
             pal = yl_bu(4)))

# Zoom as in plot_joint.R: 1 %-99 % range plus 25 % on each side.
zoom <- function(v) {
  lims <- stats::quantile(v, c(0.01, 0.99), na.rm = TRUE)
  lims + c(-1, 1) * 0.25 * diff(lims)
}

open_pdf(snakemake@output[[1]], 648, 600)
par(mfrow = c(length(beta_sel), length(n_sel)), mar = c(3, 3, 1, 0.8),
    oma = c(0, 7, 5, 0), mgp = c(1.8, 0.5, 0), tcl = -0.25, bty = "l")
for (b_i in seq_along(beta_sel)) {
  for (n_i in seq_along(n_sel)) {
    rows <- cell_rows(trials, measure, n_sel[n_i], beta_sel[b_i], 0, family)
    xlim <- zoom(rows$estimate)
    ylim <- zoom(rows$target)
    inside <- rows$estimate > xlim[1] & rows$estimate < xlim[2] &
      rows$target > ylim[1] & rows$target < ylim[2]
    v <- rows[[g$col]]
    # Quantile bins are widened at both ends so every trial falls in a bin.
    breaks <- stats::quantile(v, g$probs) +
      c(-1, rep(0, length(g$probs) - 2), 1)
    bin <- cut(v, breaks, right = FALSE)
    ord <- order(bin)  # high bins drawn last, on top

    plot(NA, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i",
         xlab = if (b_i == length(beta_sel)) labels$estimate else "",
         ylab = if (n_i == 1) labels$target else "")
    if (xlim[1] < 0 && xlim[2] > 0) abline(v = 0, col = "grey90", lwd = 0.8)
    if (ylim[1] < 0 && ylim[2] > 0) abline(h = 0, col = "grey90", lwd = 0.8)
    lo <- max(xlim[1], ylim[1])
    hi <- min(xlim[2], ylim[2])
    if (lo < hi) segments(lo, lo, hi, hi, col = MPL[4], lwd = 1.5)
    points(rows$estimate[ord], rows$target[ord], pch = 16, cex = 0.45,
           col = g$pal[as.integer(bin)][ord])
    if (any(!inside))
      mtext(sprintf("%d outside", sum(!inside)), side = 3, line = -1,
            adj = 0.98, cex = 0.6, col = "grey40")

    if (b_i == 1) mtext(bquote(n == .(n_sel[n_i])), side = 3, line = 1)
    if (n_i == 1) mtext(bquote(beta[Delta] == .(beta_sel[b_i])), side = 2,
                        line = 4.5, las = 1, adj = 1)
  }
}
# Legend across the top of the page.
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("top", legend = g$labels, col = g$pal, pch = 16, horiz = TRUE,
       bty = "n", title = g$title, cex = 0.85, pt.cex = 1.2)
dev.off()

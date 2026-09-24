# Rules plot_err and plot_errdirection: the LOO error in units of its own SD,
# as simulated/plot_all__err_n_b_both.py and
# simulated/plot_all__errdirection_n_b_both.py. Rows: beta_t. Columns: the
# outlier level. One box per n, plus the mean as a short line.
#
# plot_errdirection multiplies the error by the sign of the target. A
# positive value then means that the estimate favours the better model too
# strongly, whichever model that is.
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])
snakemake@source("lib/plot.R")

trials <- readRDS(snakemake@input[[1]])
measure <- snakemake@wildcards[["measure"]]
family <- snakemake@wildcards[["family"]]
sel <- snakemake@params[["sel"]]
signed <- isTRUE(snakemake@params[["signed"]])
beta_sel <- unlist_num(sel$beta_t)
out_sel <- unlist_num(sel$out_dev)

tr <- trials[trials$measure == measure & trials$family == family &
               trials$tau2 == sel$tau2, ]
n_sel <- sort(unique(tr$n_obs))
stopifnot(nrow(tr) > 0)

# The error of one cell, divided by its own SD over the trials.
scaled_error <- function(b, o, n) {
  rows <- tr[abs(tr$beta_t - b) < 1e-12 & abs(tr$out_dev - o) < 1e-12 &
               tr$n_obs == n, ]
  err <- loo_error(rows)
  if (signed) err <- sign(rows$target) * err
  err / stats::sd(loo_error(rows), na.rm = TRUE)
}

# One entry per beta_t, out_dev and n, so that a row can share its y axis.
vals <- lapply(beta_sel, function(b) {
  lapply(out_sel, function(o) lapply(n_sel, function(n) scaled_error(b, o, n)))
})

open_pdf(snakemake@output[[1]], 648, 576)
par(mfcol = c(length(beta_sel), length(out_sel)),
    mar = c(0.4, 0.8, 0.3, 0.8), oma = c(3.6, 6, 2.4, 0.4),
    mgp = c(2.4, 0.5, 0), tcl = -0.25)
for (o_i in seq_along(out_sel)) {
  for (b_i in seq_along(beta_sel)) {
    # The y axis of a row ignores the points outside the outermost box.
    ylim <- range(vapply(vals[[b_i]], function(per_n) {
      range(vapply(per_n, boxen_range, numeric(2)), na.rm = TRUE)
    }, numeric(2)), na.rm = TRUE)
    plot(NA, xlim = c(0.4, length(n_sel) + 0.6), ylim = ylim, axes = FALSE,
         xlab = "", ylab = "")
    for (n_i in seq_along(n_sel)) {
      x <- vals[[b_i]][[o_i]][[n_i]]
      boxen(x, n_i)
      segments(n_i - 0.21, mean(x, na.rm = TRUE),
               n_i + 0.21, mean(x, na.rm = TRUE), col = MPL[2], lwd = 1.5)
    }
    # The reference line goes on top, so that a wide box cannot hide it.
    abline(h = 0, col = adjust_lightness("gray", 1.3))
    if (o_i == 1) axis(2, las = 1) else axis(2, labels = FALSE)
    axis(1, at = seq_along(n_sel),
         labels = if (b_i == length(beta_sel)) n_sel else FALSE)
    box(bty = "l")
    if (o_i == 1) mtext(bquote(beta[Delta] == .(beta_sel[b_i])), side = 2,
                        line = 2.8, cex = 0.9)
    if (b_i == 1) mtext(if (out_sel[o_i] > 0) "outlier" else "no outlier",
                        side = 3, line = 0.6, cex = 0.9)
  }
}
mtext("n", side = 1, outer = TRUE, line = 2.2)
ylab <- if (signed) {
  bquote(sign(.(as.name(measure))) ~ err[LOO] ~ "/" ~ SD(err[LOO]))
} else {
  bquote(err[LOO] ~ "/" ~ SD(err[LOO]))
}
mtext(ylab, side = 2, outer = TRUE, line = 4.4)
dev.off()

# Rule plot_var_ratio: the stated SE divided by the actual SD of the LOO
# error over the trials of a cell, as simulated/plot_all__var_ratioi_n_b.py.
# One row per beta_t, one box per n. A ratio below 1 means that the SE
# understates the error. The quantity is scale free, so it applies to every
# measure without an additivity assumption.
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])
snakemake@source("lib/plot.R")

trials <- readRDS(snakemake@input[[1]])
measure <- snakemake@wildcards[["measure"]]
family <- snakemake@wildcards[["family"]]
sel <- snakemake@params[["sel"]]
beta_sel <- unlist_num(sel$beta_t)

tr <- trials[trials$measure == measure & trials$family == family &
               trials$out_dev == 0 & trials$tau2 == sel$tau2, ]
# The paper figure uses the whole n grid, not the plot selection.
n_sel <- sort(unique(tr$n_obs))
stopifnot(nrow(tr) > 0)

open_pdf(snakemake@output[[1]], 504, 374)
par(mfrow = c(length(beta_sel), 1), mar = c(0.4, 1.2, 0.3, 0.8),
    oma = c(3.6, 5.6, 0.6, 0.4), mgp = c(2.4, 0.5, 0), tcl = -0.25)
for (b_i in seq_along(beta_sel)) {
  b <- beta_sel[b_i]
  plot(NA, xlim = c(0.4, length(n_sel) + 0.6), ylim = c(0, 1.75),
       axes = FALSE, xlab = "", ylab = "")
  for (n_i in seq_along(n_sel)) {
    rows <- tr[abs(tr$beta_t - b) < 1e-12 & tr$n_obs == n_sel[n_i], ]
    boxen(rows$se / stats::sd(loo_error(rows), na.rm = TRUE), n_i)
  }
  # The reference line goes on top, so that a wide box cannot hide it.
  abline(h = 1, col = adjust_lightness("red", 1.3))
  axis(2, at = c(0, 0.5, 1, 1.5), las = 1)
  axis(1, at = seq_along(n_sel),
       labels = if (b_i == length(beta_sel)) n_sel else FALSE)
  box(bty = "l")
  mtext(bquote(beta[Delta] == .(b)), side = 2, line = 2.4, cex = 0.9)
}
mtext("n", side = 1, outer = TRUE, line = 2.2)
mtext(bquote(widehat(SE)[LOO] ~ "/" ~ SE(err[LOO])), side = 2, outer = TRUE,
      line = 4.2)
dev.off()

# Rule plot_calibration: PIT histograms of the normal approximation, as
# Fig. 6 (simulated/plot_all__calibration_both.py, normal column only).
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])
snakemake@source("lib/plot.R")

trials <- readRDS(snakemake@input[[1]])
measure <- snakemake@wildcards[["measure"]]
family <- snakemake@wildcards[["family"]]
sel <- snakemake@params[["sel"]]
n_sel <- unlist_num(sel$n_obs)
beta_sel <- unlist_num(sel$beta_t)
out_sel <- unlist_num(sel$out_dev)
n_bins <- 7

open_pdf(snakemake@output[[1]], 748.8 / 2 + 180, 648)
# One column per n; per beta_t one block of rows, one row per outlier level.
n_rows <- length(beta_sel) * length(out_sel)
par(mfrow = c(n_rows, length(n_sel)), mar = c(0.6, 0.8, 0.4, 0.8),
    oma = c(4, 9, 4, 1), mgp = c(2, 0.4, 0), tcl = -0.25)
for (b_i in seq_along(beta_sel)) {
  for (o_i in seq_along(out_sel)) {
    for (n_i in seq_along(n_sel)) {
      rows <- cell_rows(trials, measure, n_sel[n_i], beta_sel[b_i],
                        out_sel[o_i], sel$tau2, family)
      p <- pit(rows)
      n_na <- sum(is.na(p))
      p <- p[!is.na(p)]
      counts <- graphics::hist(p, breaks = seq(0, 1, length.out = n_bins + 1),
                               plot = FALSE)$counts
      # The envelope counts the trials in the histogram, not the dropped ones.
      q005 <- stats::qbinom(0.005, length(p), 1 / n_bins)
      q995 <- stats::qbinom(0.995, length(p), 1 / n_bins)
      ymax <- max(counts, q995) * 1.05
      plot(NA, xlim = c(0, 1), ylim = c(0, ymax), axes = FALSE, xlab = "",
           ylab = "", xaxs = "i", yaxs = "i")
      rect(seq(0, 1 - 1 / n_bins, length.out = n_bins), 0,
           seq(1 / n_bins, 1, length.out = n_bins), counts,
           col = adjust_lightness(MPL[1], 1.6), border = NA)
      rect(0, q005, 1, q995, col = grDevices::adjustcolor(MPL[2], 0.3),
           border = NA)
      last_row <- b_i == length(beta_sel) && o_i == length(out_sel)
      axis(1, at = c(0, 0.5, 1), labels = if (last_row) c("0", "0.5", "1") else FALSE)
      if (n_na > 0) mtext(sprintf("dropped: %d", n_na), side = 3, line = -1,
                          cex = 0.6, adj = 1)
      if (n_i == 1) {
        mtext(if (out_sel[o_i] > 0) "outlier" else "no outlier", side = 2,
              line = 0.5, las = 1, adj = 1, cex = 0.8)
        if (o_i == 1) mtext(bquote(beta[Delta] == .(beta_sel[b_i])), side = 2,
                            line = 5.5, las = 1, adj = 1, cex = 1)
      }
      if (b_i == 1 && o_i == 1) mtext(bquote(n == .(n_sel[n_i])), side = 3,
                                      line = 1.2, cex = 1)
    }
  }
}
mtext(bquote(p(widehat(.(as.name(measure))) < .(as.name(measure)))),
      side = 1, outer = TRUE, line = 2.5)
mtext(sprintf("%s: normal approximation", measure), side = 3, outer = TRUE,
      line = 2.5, font = 2)
dev.off()

# Rule plot_moments: mean/SD and skewness of the estimate, the target and
# the error against n, as Fig. 9 (simulated/plot_all__moments_n_b.py).
# Lines: Bayesian-bootstrap median; bands: 95 % interval.
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])
snakemake@source("lib/plot.R")

mom <- readRDS(snakemake@input[[1]])
measure <- snakemake@wildcards[["measure"]]
sel <- snakemake@params[["sel"]]
beta_sel <- unlist_num(sel$moments_beta_t)
labels <- measure_labels(measure)
mom <- mom[mom$measure == measure & mom$out_dev == 0 & mom$tau2 == sel$tau2, ]

quantities <- c("estimate", "target", "error")
stats_names <- c("mean/SD", "skewness")

open_pdf(snakemake@output[[1]], 720, 432)
layout(matrix(c(1:6, 7, 7), 2, 4), widths = c(1, 1, 1, 0.35))
par(mar = c(3.2, 4, 2.2, 0.8), mgp = c(2, 0.5, 0), tcl = -0.25, bty = "l")
for (q in quantities) {
  for (s in stats_names) {
    d <- mom[mom$quantity == q & mom$stat == s &
               mom$beta_t %in% beta_sel, ]
    ylim <- if (q == "error" && s == "mean/SD") c(-0.31, 0.31) else
      range(c(d$q025, d$q975, 0), finite = TRUE)
    plot(NA, xlim = range(d$n_obs), ylim = ylim,
         xlab = if (s == "skewness") "n" else "",
         ylab = if (q == "estimate") s else "")
    abline(h = 0, col = "grey50", lwd = 1)
    for (b_i in seq_along(beta_sel)) {
      db <- d[abs(d$beta_t - beta_sel[b_i]) < 1e-12, ]
      db <- db[order(db$n_obs), ]
      if (nrow(db) == 0) next
      polygon(c(db$n_obs, rev(db$n_obs)), c(db$q025, rev(db$q975)),
              col = grDevices::adjustcolor(MPL[b_i], 0.2), border = NA)
      lines(db$n_obs, db$q500, col = MPL[b_i], lwd = 1.5)
    }
    if (s == "mean/SD") title(main = labels[[q]], font.main = 1)
  }
}
par(mar = c(0, 0, 0, 0))
plot.new()
legend("left", legend = beta_sel, col = MPL[seq_along(beta_sel)], lwd = 1.5,
       bty = "n", title = expression(beta[Delta]))
dev.off()

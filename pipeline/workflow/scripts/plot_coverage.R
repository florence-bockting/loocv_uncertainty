# Rule plot_coverage: coverage of the 95 % normal interval against the
# observed |diff|, in quantile bins, for the cells without an outlier.
# The cutoff of a measure is the lower edge of the first bin after the
# coverage minimum whose coverage reaches the target.
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])
snakemake@source("lib/plot.R")

n_bins <- snakemake@params[["n_bins"]]
target_cov <- snakemake@params[["target_cov"]]
families <- unlist(snakemake@params[["families"]])
fam_col <- stats::setNames(MPL[seq_along(families)], families)

trials <- readRDS(snakemake@input[[1]])
trials <- trials[trials$out_dev == 0 & is.finite(trials$se) &
                   is.finite(trials$estimate), ]
trials$covered <- abs(trials$estimate - trials$target) <=
  stats::qnorm(0.975) * trials$se
measures <- unique(trials$measure)

# Coverage per quantile bin of `ad`, for each family and measure.
bin_coverage <- function(ad) {
  groups <- split(seq_len(nrow(trials)),
                  list(trials$family, trials$measure), drop = TRUE)
  do.call(rbind, lapply(groups, function(i) {
    br <- unique(stats::quantile(ad[i], seq(0, 1, length.out = n_bins + 1)))
    k <- cut(ad[i], br, include.lowest = TRUE)
    data.frame(family = trials$family[i[1]], measure = trials$measure[i[1]],
               bin = seq_along(levels(k)), lo = utils::head(br, -1),
               mid = as.vector(tapply(ad[i], k, stats::median)),
               cov = as.vector(tapply(trials$covered[i], k, mean)),
               n = as.vector(table(k)))
  }))
}

cutoffs <- function(b) {
  groups <- split(b, list(b$family, b$measure), drop = TRUE)
  out <- do.call(rbind, lapply(groups, function(d) {
    m <- which.min(d$cov)
    i <- which(d$cov >= target_cov & seq_len(nrow(d)) > m)[1]
    above <- if (is.na(i)) NA else d$cov[i:nrow(d)]
    data.frame(family = d$family[1], measure = d$measure[1],
               min_cov = round(d$cov[m], 3), at_diff = signif(d$mid[m], 3),
               thr_diff = if (is.na(i)) NA else signif(d$lo[i], 3),
               thr_quantile = if (is.na(i)) NA else round((i - 1) / n_bins, 2),
               max_cov_above = round(max(above), 3),
               min_cov_above = round(min(above), 3))
  }))
  rownames(out) <- NULL
  out[order(out$measure, out$family), ]
}

# One panel per measure; x on a log scale, free per panel.
coverage_panels <- function(path, b, fams, xlab, title, cut = NULL) {
  b <- b[b$mid > 0 & b$family %in% fams, ]
  ms <- intersect(measures, b$measure)
  n_col <- if (is.null(cut)) 4 else 3
  n_row <- ceiling(length(ms) / n_col)
  open_pdf(path, 720, 110 + 150 * n_row)
  par(mfrow = c(n_row, n_col), mar = c(2.5, 3, 1.6, 0.8),
      oma = c(2, 1, 4, 0), mgp = c(1.8, 0.5, 0), tcl = -0.25, bty = "l")
  for (m in ms) {
    d <- b[b$measure == m, ]
    c_m <- if (is.null(cut)) NULL else cut$thr_diff[cut$measure == m]
    main <- if (length(c_m) && !is.na(c_m)) {
      sprintf("%s  (cutoff %g)", m, round(c_m))
    } else m
    plot(d$mid, d$cov, type = "n", log = "x", ylim = c(0.4, 1),
         xlab = "", ylab = "coverage", main = main, font.main = 1,
         cex.main = 1)
    abline(h = target_cov, col = "grey50", lwd = 0.8)
    abline(h = 0.95, col = "grey60", lty = 2, lwd = 0.8)
    if (length(c_m) && !is.na(c_m)) abline(v = round(c_m), col = MPL[2])
    for (f in fams) {
      e <- d[d$family == f, ]
      col <- if (length(fams) > 1) fam_col[f] else MPL[1]
      lines(e$mid, e$cov, col = col, type = "o", pch = 16, cex = 0.5)
    }
    if (m == ms[1] && length(fams) > 1) {
      legend("bottomright", fams, col = fam_col[fams], lty = 1, pch = 16,
             pt.cex = 0.5, bty = "n", cex = 0.8)
    }
  }
  mtext(xlab, side = 1, outer = TRUE, line = 0.5, cex = 0.8)
  mtext(title, side = 3, outer = TRUE, line = 2.2)
  legend_txt <- sprintf("solid grey: %.2f target   dashed grey: 0.95 nominal",
                        target_cov)
  if (!is.null(cut)) legend_txt <- paste(legend_txt, "  orange: cutoff")
  mtext(legend_txt, side = 3, outer = TRUE, line = 0.6, cex = 0.7)
  invisible(dev.off())
}

# |diff| on its own scale.
raw <- bin_coverage(abs(trials$estimate))
coverage_panels(snakemake@output[["raw"]], raw, families,
                "|diff|, bin median, log scale",
                sprintf("Coverage by |diff|, out_dev = 0, %d quantile bins",
                        n_bins))

# |diff| x n_obs for the mean-scale measures; elpd and ic are sums already.
per_obs <- ifelse(trials$measure %in% c("elpd", "ic"), 1, trials$n_obs)
scaled <- bin_coverage(abs(trials$estimate) * per_obs)
cut <- cutoffs(scaled)
utils::write.csv(cut, snakemake@output[["cutoffs"]], row.names = FALSE)
print(cut, row.names = FALSE)
xlab_n <- "|diff| x n_obs (elpd, ic: |diff|), bin median, log scale"
coverage_panels(snakemake@output[["scaled"]], scaled, families, xlab_n,
                sprintf("Coverage by |diff| x n_obs, out_dev = 0, %d quantile bins",
                        n_bins))
for (f in families) {
  coverage_panels(snakemake@output[[paste0("cutoff_", f)]], scaled, f, xlab_n,
                  sprintf("%s: coverage by |diff| x n_obs, with cutoff", f),
                  cut = cut[cut$family == f, ])
}

# Plot helpers that copy the matplotlib style of simulated/plot_all__*.py.

# matplotlib default colour cycle C0, C1, ...
MPL <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")

# adjust_lightness in the plot scripts: scales the HLS lightness.
adjust_lightness <- function(col, amount) {
  rgb <- grDevices::col2rgb(col)[, 1] / 255
  mx <- max(rgb)
  mn <- min(rgb)
  l <- (mx + mn) / 2
  s <- if (mx == mn) 0 else if (l <= 0.5) (mx - mn) / (mx + mn) else (mx - mn) / (2 - mx - mn)
  h <- if (mx == mn) 0 else {
    d <- mx - mn
    hh <- if (mx == rgb[1]) (rgb[2] - rgb[3]) / d else if (mx == rgb[2]) 2 + (rgb[3] - rgb[1]) / d else 4 + (rgb[1] - rgb[2]) / d
    (hh / 6) %% 1
  }
  l <- min(1, max(0, amount * l))
  q <- if (l < 0.5) l * (1 + s) else l + s - l * s
  p <- 2 * l - q
  hue <- function(t) {
    t <- t %% 1
    if (t < 1 / 6) p + (q - p) * 6 * t else if (t < 1 / 2) q else if (t < 2 / 3) p + (q - p) * (2 / 3 - t) * 6 else p
  }
  grDevices::rgb(hue(h + 1 / 3), hue(h), hue(h - 1 / 3))
}

# matplotlib "copper" colormap, truncated to [0.3, 1] as kde_cmap.
copper <- function(k) {
  x <- seq(0.3, 1, length.out = k)
  grDevices::rgb(pmin(1, 1.25 * x), 0.7812 * x, 0.4975 * x)
}

# matplotlib "Greys" colormap, truncated to [0.3, 1] as the hexbin cmap.
greys <- function(k) {
  grDevices::grey(seq(1 - 0.3, 0, length.out = k) * 0.95)
}

# Trials of one cell, selected by its parameters.
cell_rows <- function(trials, measure, n_obs, beta_t, out_dev, tau2) {
  trials[trials$measure == measure & trials$n_obs == n_obs &
           abs(trials$beta_t - beta_t) < 1e-12 &
           abs(trials$out_dev - out_dev) < 1e-12 &
           trials$tau2 == tau2, ]
}

# Probability integral transform of the target under N(estimate, SE).
pit <- function(rows) {
  stats::pnorm(rows$target, mean = rows$estimate, sd = rows$se)
}

# Axis names of a measure: the estimate and the target.
measure_labels <- function(measure) {
  list(estimate = bquote(widehat(.(as.name(measure)))[LOO]),
       target = as.name(measure),
       error = bquote(err[LOO]))
}

# Opens a pdf of the given size in points, as the matplotlib figures.
open_pdf <- function(path, width_pt, height_pt) {
  grDevices::pdf(path, width = width_pt / 72, height = height_pt / 72)
}

unlist_num <- function(x) as.numeric(unlist(x))

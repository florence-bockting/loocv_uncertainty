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

# Rows of one covariate distribution x_dist (see lib/dgp.R). A table
# without an x_dist column comes from an earlier run: with an x_df column,
# Inf is the normal and a finite df the t covariate; without one, it holds
# normal cells only.
select_x_dist <- function(d, x_dist) {
  if (is.null(d$x_dist)) {
    x_df <- if (is.null(d$x_df)) rep(Inf, nrow(d)) else d$x_df
    d$x_dist <- ifelse(is.finite(x_df), paste0("t", x_df), "normal")
  }
  d[d$x_dist == x_dist, ]
}

# Title of a covariate distribution, such as "t(3)" for "t3". The slant of
# the skew-t mirrors SKEWT_ALPHA in lib/dgp.R.
x_dist_title <- function(x_dist) {
  dist <- sub("[0-9.]+$", "", x_dist)
  par <- sub("^[a-z]+", "", x_dist)
  switch(dist,
    normal = "N(0, 1)",
    t = sprintf("t(%s)", par),
    skewt = sprintf("skew-t(%s, alpha = 5), standardised", par),
    gamma = sprintf("gamma(%s), standardised", par),
    stop("unknown covariate distribution: ", x_dist)
  )
}

# Trials of one cell, selected by its parameters.
cell_rows <- function(trials, measure, n_obs, beta_t, out_dev,
                      family = "gaussian") {
  trials[trials$measure == measure & trials$family == family &
           trials$n_obs == n_obs &
           abs(trials$beta_t - beta_t) < 1e-12 &
           abs(trials$out_dev - out_dev) < 1e-12, ]
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

# Letter-value ("boxen") plot of x at the position `at`, as the boxenplot of
# seaborn with its default depth. The innermost box spans the quartiles;
# each box further out halves the tail probability, and is narrower and
# lighter. The outer boxes are drawn first, so the inner ones stay visible.
boxen <- function(x, at, width = 0.8, col = MPL[1]) {
  x <- x[is.finite(x)]
  if (length(x) < 8) return(invisible(NULL))
  k <- max(1, floor(log2(length(x))) - 3)
  for (i in rev(seq_len(k))) {
    p <- 2^-(i + 1)
    q <- stats::quantile(x, c(p, 1 - p), names = FALSE)
    w <- width / 2 * 0.5^((i - 1) / 2)
    graphics::rect(at - w, q[1], at + w, q[2],
                   col = adjust_lightness(col, 1.3 + 0.1 * (i - 1)),
                   border = "black", lwd = 0.5)
  }
  graphics::segments(at - width / 2, stats::median(x),
                     at + width / 2, stats::median(x), lwd = 1.2)
}

# The deepest quantiles that boxen() draws, for an axis range that ignores
# the points outside the outermost box.
boxen_range <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 8) return(c(NA_real_, NA_real_))
  p <- 2^-(max(1, floor(log2(length(x))) - 3) + 1)
  stats::quantile(x, c(p, 1 - p), names = FALSE)
}

# LOO error of each trial: the estimate minus the target.
loo_error <- function(rows) rows$estimate - rows$target

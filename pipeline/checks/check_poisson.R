# Checks of the poisson cells, without Snakemake. Run from the repository
# root:
#     Rscript pipeline/checks/check_poisson.R
#
# 1. the penalised IRLS mode and its Laplace covariance against optim();
# 2. the predictive pmf against a brute-force integration, and its mean
#    against the lognormal mean exp(m + s^2 / 2);
# 3. the LOO predictive against a brute-force refit;
# 4. the pointwise scores against their definitions, and rps against the
#    ranked probability score written out as sum_k (F_k - 1{y <= k})^2;
# 5. one whole cell: data -> fit -> score -> compare.

source("pipeline/workflow/scripts/lib/dgp.R")
source("pipeline/workflow/scripts/lib/predictive.R")
source("pipeline/workflow/scripts/lib/glm.R")
source("pipeline/workflow/scripts/lib/scores.R")
source("pipeline/workflow/scripts/lib/measures.R")

ok <- function(label, value, tol) {
  pass <- is.finite(value) && value <= tol
  cat(sprintf("%-46s %9.2e  %s\n", label, value, if (pass) "ok" else "FAIL"))
  if (!pass) stop(label, " failed")
}

set.seed(20250924)
n <- 50
X <- cbind(1, matrix(rnorm(n * 2), n, 2))
beta <- c(1, 0.5, 0.25)
y <- rpois(n, exp(X %*% beta))

# 1. mode and covariance ------------------------------------------------
neg_log_post <- function(b, X, y) {
  eta <- drop(X %*% b)
  -sum(y * eta - exp(eta)) + sum(b^2) / (2 * PRIOR_SD^2)
}
fit <- fit_poisson(X, y)
opt <- optim(numeric(3), neg_log_post, X = X, y = y, method = "BFGS",
             hessian = TRUE, control = list(reltol = 1e-15))
ok("mode vs optim()", max(abs(fit$beta_hat - opt$par)), 1e-5)
ok("covariance vs optim() hessian",
   max(abs(fit$V - solve(opt$hessian))), 1e-4)

# 2. predictive pmf ------------------------------------------------------
# integrate() over the whole line misses the peak of dpois(k, exp(m + s t))
# once k is large, because that peak sits far from zero and is narrow. The
# interval is therefore centred on the peak, at t = (log k - m) / s, and
# spans many peak widths of 1 / (s sqrt(k)) on each side.
exact_pmf <- function(m, s, k) {
  vapply(k, function(kk) {
    f <- function(t) dpois(kk, exp(m + s * t)) * dnorm(t)
    t0 <- (log(max(kk, 0.5)) - m) / s
    w <- 1 + 12 / (s * sqrt(kk + 1))
    integrate(f, max(-12, t0 - w), min(12, t0 + w), rel.tol = 1e-12,
              subdivisions = 1000)$value
  }, numeric(1))
}
# s reaches about 1.5 at the smallest n, when a test point lies far from
# the training covariates.
grid <- expand.grid(m = c(-1, 0, 1, 2.5), s = c(0.05, 0.5, 1.5))
err <- mapply(function(m, s) {
  p <- predictive_pmf(m, s)
  max(abs(p$pmf[1, ] - exact_pmf(m, s, p$k)))
}, grid$m, grid$s)
ok("predictive pmf vs integrate()", max(err), 1e-6)

# E[X] = E[lambda] = exp(m + s^2 / 2), because the rate is lognormal. The
# grid holds all but POIS_TAIL of the mass, so the summed mean agrees only
# to that order. pointwise_scores_count() therefore uses the exact mean.
err <- mapply(function(m, s) {
  p <- predictive_pmf(m, s)
  abs(drop(p$pmf %*% p$k) / exp(m + s^2 / 2) - 1)
}, grid$m, grid$s)
ok("summed mean vs lognormal mean, relative", max(err), 1e-3)
ok("each row of the pmf sums to one",
   max(abs(rowSums(predictive_pmf(grid$m, grid$s)$pmf) - 1)), 1e-12)

# The density of the observed count needs no count grid, so only the range
# of the quadrature limits it. That range drops 6e-7 of the normal mass.
err <- mapply(function(m, s) {
  y_k <- 0:20
  max(abs(predictive_count_dens(y_k, rep(m, 21), rep(s, 21)) -
            exact_pmf(m, s, y_k)))
}, grid$m, grid$s)
ok("density of the observed count vs integrate()", max(err), 1e-6)

# 3. LOO predictive -----------------------------------------------------
got <- loo_predictive_poisson(X, y)
brute <- vapply(seq_len(n), function(i) {
  o <- optim(numeric(3), neg_log_post, X = X[-i, ], y = y[-i], method = "BFGS",
             hessian = TRUE, control = list(reltol = 1e-15))
  c(sum(X[i, ] * o$par),
    sqrt(drop(X[i, ] %*% solve(o$hessian) %*% X[i, ])))
}, numeric(2))
ok("LOO predictive mean vs brute force", max(abs(got$m - brute[1, ])), 1e-5)
ok("LOO predictive sd vs brute force", max(abs(got$s - brute[2, ])), 1e-5)

# 4. pointwise scores ---------------------------------------------------
s <- pointwise_scores_count(y, got)
p <- predictive_pmf(got$m, got$s)
mu <- exp(got$m + got$s^2 / 2)
# The score takes the density without a grid, the pmf is renormalised
# over one, so the two agree only to the mass that the grid drops.
ok("elpd score", max(abs(s$elpd - log(p$pmf[cbind(seq_len(n), y + 1)]))), 1e-6)
ok("sqe score", max(abs(s$sqe - (y - mu)^2)), 1e-10)
ok("ae score", max(abs(s$ae - abs(y - mu))), 1e-10)

# The ranked probability score of an ordered count support, written out.
cdf <- t(apply(p$pmf, 1, cumsum))
rps_direct <- rowSums((cdf - outer(y, p$k, "<="))^2)
ok("rps equals the ranked probability score",
   max(abs(s$rps - rps_direct)), 1e-8)

abs_dev <- rowSums(p$pmf * abs(outer(y, p$k, "-")))
spread <- 2 * rowSums(cdf * (1 - cdf))
ok("srps score",
   max(abs(s$srps - (-abs_dev / spread - log(spread) / 2))), 1e-10)

# 5. one cell -----------------------------------------------------------
cell <- list(family = "poisson", n_obs = 40, beta_t = 0.5, out_dev = 2)
n_trial <- 20
set.seed(1)
train <- make_data_family(cell$family, n_trial, cell$n_obs, cell$n_obs,
                          make_beta(cell$beta_t), cell$out_dev)
test <- pool_sets(make_data_family(cell$family, 50, cell$n_obs, cell$n_obs,
                                   make_beta(cell$beta_t), cell$out_dev))
cat(sprintf("%-46s %9.1f\n", "mean count of the test set", mean(test$y)))
cat(sprintf("%-46s %9d\n", "largest count of the test set", max(test$y)))
scores <- scores_for_family(SCORES_POISSON, cell$family)
measures <- measures_for_family(names(MEASURE_SCORE), cell$family)

start <- Sys.time()
rows <- lapply(seq_len(n_trial), function(t) {
  y_t <- train$y[, t]
  per_model <- lapply(c(A = "A", B = "B"), function(model) {
    cols <- model_columns(model, dim(train$X)[2])
    f <- fit_model_family(cell$family, train$X[, cols, t], y_t)
    test_i <- pointwise_scores_family(
      cell$family, test$y,
      predict_new_family(cell$family, f, test$X[, cols]), scores)
    list(loo = pointwise_scores_family(cell$family, y_t, f$loo, scores),
         test = lapply(test_i, mean))
  })
  est <- vapply(measures, loo_difference, numeric(2),
                loo_a = per_model$A$loo, loo_b = per_model$B$loo, y = y_t)
  tgt <- vapply(measures, target_difference, numeric(1),
                test_a = per_model$A$test, test_b = per_model$B$test,
                n_obs = cell$n_obs,
                y_test_var = mean((test$y - mean(test$y))^2))
  data.frame(trial = t, measure = measures, estimate = est[1, ],
             se = est[2, ], target = tgt, row.names = NULL)
})
res <- do.call(rbind, rows)
cat(sprintf("%-46s %9.1f\n", "seconds for 20 trials of both models",
            as.numeric(Sys.time() - start, units = "secs")))

ok("every measure is finite", sum(!is.finite(as.matrix(res[, -(1:2)]))), 0)
ok("every SE is positive", sum(res$se <= 0), 0)
ok("brier, acc and bacc are not poisson measures",
   sum(measures %in% MEASURES_BINARY), 0)
# ic is 2 elpd on the same scale, so both must agree up to the factor.
wide <- reshape(res[res$measure %in% c("elpd", "ic", "mlpd"), ],
                direction = "wide", idvar = "trial", timevar = "measure")
ok("ic equals 2 elpd",
   max(abs(wide$estimate.ic - 2 * wide$estimate.elpd)), 1e-10)
ok("mlpd equals elpd / n",
   max(abs(wide$estimate.mlpd - wide$estimate.elpd / cell$n_obs)), 1e-10)
# On a count support rps is its own measure, not a copy of mse.
wide2 <- reshape(res[res$measure %in% c("mse", "rps"), ],
                 direction = "wide", idvar = "trial", timevar = "measure")
ok("rps differs from mse",
   -max(abs(wide2$estimate.mse - wide2$estimate.rps)), -1e-6)

cat("\nall checks passed\n")
print(aggregate(cbind(estimate, se, target) ~ measure, res, mean))

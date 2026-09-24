# Checks of the binomial cells, without Snakemake. Run from the repository
# root:
#     Rscript pipeline/checks/check_binomial.R
#
# 1. the penalised IRLS mode and its Laplace covariance against optim();
# 2. the quadrature for E[plogis(eta)] against integrate();
# 3. the LOO predictive against a brute-force refit and integration;
# 4. the pointwise scores against their definitions;
# 5. the bacc difference and its SE against a direct computation, and the
#    guards that drop a degenerate trial;
# 6. one whole cell: data -> fit -> score -> compare.

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

set.seed(20250923)
n <- 50
X <- cbind(1, matrix(rnorm(n * 2), n, 2))
beta <- c(0, 1, 0.5)
y <- rbinom(n, 1, plogis(X %*% beta))

# 1. mode and covariance ------------------------------------------------
neg_log_post <- function(b, X, y) {
  eta <- drop(X %*% b)
  -sum(y * eta - log1p(exp(eta))) + sum(b^2) / (2 * PRIOR_SD^2)
}
fit <- fit_logistic(X, y)
opt <- optim(numeric(3), neg_log_post, X = X, y = y, method = "BFGS",
             hessian = TRUE, control = list(reltol = 1e-15))
ok("mode vs optim()", max(abs(fit$beta_hat - opt$par)), 1e-5)
ok("covariance vs optim() hessian",
   max(abs(fit$V - solve(opt$hessian))), 1e-4)

# 2. quadrature ---------------------------------------------------------
exact_prob <- function(m, s) {
  integrate(function(t) plogis(m + s * t) * dnorm(t), -Inf, Inf,
            rel.tol = 1e-13)$value
}
grid <- expand.grid(m = c(-30, -8, -1, 0, 1, 8, 30),
                    s = c(0.001, 0.1, 1, 4, 16, 32))
err <- mapply(function(m, s) abs(predictive_prob(m, s) - exact_prob(m, s)),
              grid$m, grid$s)
ok("quadrature vs integrate()", max(err), 1e-12)

# 3. LOO predictive -----------------------------------------------------
got <- loo_predictive_binary(X, y)$p
brute <- vapply(seq_len(n), function(i) {
  o <- optim(numeric(3), neg_log_post, X = X[-i, ], y = y[-i], method = "BFGS",
             hessian = TRUE, control = list(reltol = 1e-15))
  exact_prob(sum(X[i, ] * o$par),
             sqrt(drop(X[i, ] %*% solve(o$hessian) %*% X[i, ])))
}, numeric(1))
ok("LOO predictive vs brute force", max(abs(got - brute)), 1e-5)

# 4. pointwise scores ---------------------------------------------------
p <- got
s <- pointwise_scores_binary(y, list(p = p))
ok("elpd score", max(abs(s$elpd - log(ifelse(y == 1, p, 1 - p)))), 1e-12)
ok("sqe score", max(abs(s$sqe - (y - p)^2)), 1e-12)
ok("ae score", max(abs(s$ae - abs(y - p))), 1e-12)
# For a binary outcome the RPS equals the Brier score.
ok("rps equals the Brier score", max(abs(s$rps - (y - p)^2)), 1e-12)
ok("srps score",
   max(abs(s$srps - (-ifelse(y == 1, 1 - p, p) / (2 * p * (1 - p)) -
                       log(2 * p * (1 - p)) / 2))), 1e-12)
ok("acc score", max(abs(s$acc - as.numeric((p >= 0.5) == (y == 1)))), 1e-12)

# 5. bacc ---------------------------------------------------------------
acc_a <- rbinom(n, 1, 0.7)
acc_b <- rbinom(n, 1, 0.5)
got <- stratified_mean_se(acc_a - acc_b, y)
classes <- sort(unique(y))
d_by_class <- lapply(classes, function(cl) (acc_a - acc_b)[y == cl])
want_est <- mean(vapply(d_by_class, mean, numeric(1)))
want_se <- sqrt(sum(vapply(d_by_class, function(d) var(d) / length(d),
                           numeric(1)))) / length(classes)
ok("bacc estimate", abs(got[1] - want_est), 1e-12)
ok("bacc standard error", abs(got[2] - want_se), 1e-12)

# The guards of finite_se() and stratified_mean_se(). Each case returns NA,
# so that summarise.R drops the trial instead of counting an SE of zero.
is_na <- function(x) as.numeric(!all(is.na(x)))
same <- list(acc = c(1, 1, 0, 1))
ok("acc is NA when the models agree everywhere",
   is_na(loo_difference("acc", same, same, c(0, 0, 1, 1))[2]), 0)
ok("bacc SE is NA when the models agree everywhere",
   is_na(loo_difference("bacc", same, same, c(0, 0, 1, 1))[2]), 0)
other <- list(acc = c(0, 1, 0, 1))
ok("bacc is NA when a class is absent",
   is_na(loo_difference("bacc", same, other, c(0, 0, 0, 0))), 0)
ok("bacc is NA when a class holds one observation",
   is_na(loo_difference("bacc", same, other, c(0, 1, 1, 1))), 0)
ok("rmse is NA when the models agree everywhere",
   is_na(loo_difference("rmse", list(sqe = c(1, 2, 3)),
                        list(sqe = c(1, 2, 3)), c(0, 1, 1))[2]), 0)

# 6. one cell -----------------------------------------------------------
cell <- list(family = "binomial", n_obs = 40, beta_t = 0.5, out_dev = 5)
n_trial <- 20
set.seed(1)
train <- make_data_family(cell$family, n_trial, cell$n_obs, cell$n_obs,
                          make_beta(cell$beta_t), cell$out_dev)
test <- pool_sets(make_data_family(cell$family, 200, cell$n_obs, cell$n_obs,
                                   make_beta(cell$beta_t), cell$out_dev))
scores <- scores_for_family(SCORES_BINOMIAL, cell$family)
measures <- measures_for_family(names(MEASURE_SCORE), cell$family)
classes <- sort(unique(test$y))

rows <- lapply(seq_len(n_trial), function(t) {
  y_t <- train$y[, t]
  per_model <- lapply(c(A = "A", B = "B"), function(model) {
    cols <- model_columns(model, dim(train$X)[2])
    f <- fit_model_binary(train$X[, cols, t], y_t)
    test_i <- pointwise_scores_binary(test$y,
                                      predict_new_binary(f, test$X[, cols]),
                                      scores)
    list(loo = pointwise_scores_binary(y_t, f$loo, scores),
         test = lapply(test_i, mean),
         by_class = lapply(test_i, function(v) {
           vapply(classes, function(cl) mean(v[test$y == cl]), numeric(1))
         }))
  })
  est <- vapply(measures, loo_difference, numeric(2),
                loo_a = per_model$A$loo, loo_b = per_model$B$loo, y = y_t)
  tgt <- vapply(measures, target_difference, numeric(1),
                test_a = per_model$A$test, test_b = per_model$B$test,
                n_obs = cell$n_obs,
                y_test_var = mean((test$y - mean(test$y))^2),
                by_class_a = per_model$A$by_class,
                by_class_b = per_model$B$by_class)
  data.frame(trial = t, measure = measures, estimate = est[1, ],
             se = est[2, ], target = tgt, row.names = NULL)
})
res <- do.call(rbind, rows)

# Only acc and bacc may hold a dropped trial; every other measure is finite.
keep <- !(res$measure %in% c("acc", "bacc"))
ok("every measure but acc and bacc is finite",
   sum(!is.finite(as.matrix(res[keep, -(1:2)]))), 0)
ok("every SE is non-negative", sum(res$se < 0, na.rm = TRUE), 0)
ok("no SE equals zero", sum(res$se == 0, na.rm = TRUE), 0)
cat(sprintf("%-46s %9d\n", "dropped acc and bacc trials",
            sum(is.na(res$se[!keep]))))
# ic is -2 elpd on the same scale, so both must agree up to the factor.
wide <- reshape(res[res$measure %in% c("elpd", "ic", "mlpd"), ],
                direction = "wide", idvar = "trial", timevar = "measure")
ok("ic equals 2 elpd",
   max(abs(wide$estimate.ic - 2 * wide$estimate.elpd)), 1e-10)
ok("mlpd equals elpd / n",
   max(abs(wide$estimate.mlpd - wide$estimate.elpd / cell$n_obs)), 1e-10)
# For a binary outcome mse, brier and rps are the same measure.
wide2 <- reshape(res[res$measure %in% c("mse", "brier", "rps"), ],
                 direction = "wide", idvar = "trial", timevar = "measure")
ok("mse equals brier", max(abs(wide2$estimate.mse - wide2$estimate.brier)), 1e-12)
ok("mse equals rps", max(abs(wide2$estimate.mse - wide2$estimate.rps)), 1e-12)

cat("\nall checks passed\n")
print(aggregate(cbind(estimate, se, target) ~ measure, res, mean))

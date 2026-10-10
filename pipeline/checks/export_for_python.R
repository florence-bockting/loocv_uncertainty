source("pipeline/workflow/scripts/lib/dgp.R")
source("pipeline/workflow/scripts/lib/glm.R")
source("pipeline/workflow/scripts/lib/measures.R")

n_obs <- 16
beta_t <- 0.5
out_dev <- 20
n_trial <- 3
n_test <- 20
out <- "pipeline/checks/tmp"
dir.create(out, showWarnings = FALSE)

# 1. Data: training sets first, then test sets (same order as run_all.py)
set.seed(1)
beta <- make_beta_coeff(beta_t)
train <- make_data_gaussian(n_trial, n_obs, n_obs, beta, out_dev)
test <- make_data_gaussian(n_test, n_obs, n_obs, beta, out_dev)
test_pooled <- pool_sets(test)

# 2. R results for each training set t
loo_a <- loo_b <- matrix(NA, n_trial, n_obs)
target_a <- target_b <- numeric(n_trial)
for (t in 1:n_trial) {
  X <- train$X[, , t]
  y <- train$y[, t]
  loo_a[t, ] <- pointwise_scores(y, loo_predictive_gaussian(X[, -3], y), "elpd")$elpd
  loo_b[t, ] <- pointwise_scores(y, loo_predictive_gaussian(X, y), "elpd")$elpd
  target_a[t] <- n_obs * mean(pointwise_scores(test_pooled$y,
    predict_new_gaussian(fit_model_gaussian(X[, -3], y), test_pooled$X[, -3]), "elpd")$elpd)
  target_b[t] <- n_obs * mean(pointwise_scores(test_pooled$y,
    predict_new_gaussian(fit_model_gaussian(X, y), test_pooled$X), "elpd")$elpd)
}

# 3. Write the files (arrays in numpy order [T, n, d])
save <- function(x, name) writeLines(sprintf("%.17g", x), file.path(out, name))
save(aperm(train$X, c(2, 1, 3)), "X.txt")
save(train$y, "y.txt")
save(aperm(test$X, c(2, 1, 3)), "X_test.txt")
save(test$y, "y_test.txt")
save(t(loo_a), "loo_a.txt")
save(t(loo_b), "loo_b.txt")
save(target_a, "target_a.txt")
save(target_b, "target_b.txt")

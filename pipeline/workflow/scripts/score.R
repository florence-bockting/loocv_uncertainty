# Rule score: pointwise LOO scores and test-set mean scores of one model,
# for one chunk of trials.
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])
snakemake@source("lib/glm.R")
snakemake@source("lib/measures.R")

fit <- readRDS(snakemake@input[["fit"]])
train <- readRDS(snakemake@input[["train"]])
test <- readRDS(snakemake@input[["test"]])
cell <- snakemake@params[["cell"]]
# The scores that the measures of the family are computed from.
measures <- measures_for_family(unlist(snakemake@params[["measures"]]),
                                cell$family)
scores <- unique(unname(MEASURE_SCORE[measures]))
cols <- model_columns(snakemake@wildcards[["model"]], ncol(test$X))
X_test <- test$X[, cols, drop = FALSE]
binary <- cell$family == "binomial"

n_obs <- nrow(train$y)
k <- length(fit$trials)
loo <- sapply(scores, function(s) matrix(NA_real_, n_obs, k),
              simplify = FALSE)
test_mean <- sapply(scores, function(s) numeric(k), simplify = FALSE)
# bacc needs the test means within each class, not the pooled mean.
classes <- if (binary) sort(unique(test$y)) else NULL
test_by_class <- if (binary) {
  sapply(scores, function(s) matrix(NA_real_, length(classes), k),
         simplify = FALSE)
} else {
  NULL
}

start <- Sys.time()
for (i in seq_len(k)) {
  loo_i <- pointwise_scores_family(
    family = cell$family,
    y = train$y[, fit$trials[i]],
    pred = fit$loos[[i]],
    scores = scores
  )
  test_i <- pointwise_scores_family(
    family = cell$family,
    y = test$y,
    pred = predict_new_family(
      family = cell$family,
      fit = fit$fits[[i]],
      X_new = X_test
    ),
    scores = scores
  )

  for (s in scores) {
    loo[[s]][, i] <- loo_i[[s]]
    test_mean[[s]][i] <- mean(test_i[[s]])
    if (binary) {
      test_by_class[[s]][, i] <- vapply(
        classes,
        function(cl) mean(test_i[[s]][test$y == cl]),
        numeric(1)
      )
    }
  }
}
message(sprintf("%d trials in %.1f s", k,
                as.numeric(Sys.time() - start, units = "secs")))

saveRDS(list(trials = fit$trials, loo = loo, test = test_mean,
             test_by_class = test_by_class),
        snakemake@output[[1]])

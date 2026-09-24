# Rule compare: LOO estimate, SE and target of Ma - Mb for every measure,
# for one chunk of trials.
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])
snakemake@source("lib/measures.R")

a <- readRDS(snakemake@input[["a"]])
b <- readRDS(snakemake@input[["b"]])
train <- readRDS(snakemake@input[["train"]])
y_test_var <- readRDS(snakemake@input[["info"]])$y_test_var
cell <- snakemake@params[["cell"]]
measures <- measures_for_family(unlist(snakemake@params[["measures"]]),
                                cell$family)
stopifnot(identical(a$trials, b$trials))

rows <- lapply(seq_along(a$trials), function(i) {
  y <- train$y[, a$trials[i]]
  loo_a <- lapply(a$loo, function(m) m[, i])
  loo_b <- lapply(b$loo, function(m) m[, i])
  test_a <- lapply(a$test, function(v) v[i])
  test_b <- lapply(b$test, function(v) v[i])
  by_class_a <- lapply(a$test_by_class, function(m) m[, i])
  by_class_b <- lapply(b$test_by_class, function(m) m[, i])
  est <- vapply(measures, loo_difference, numeric(2),
                loo_a = loo_a, loo_b = loo_b, y = y)
  tgt <- vapply(measures, target_difference, numeric(1),
                test_a = test_a, test_b = test_b, n_obs = length(y),
                y_test_var = y_test_var, by_class_a = by_class_a,
                by_class_b = by_class_b)
  data.frame(trial = a$trials[i], measure = measures, estimate = est[1, ],
             se = est[2, ], target = tgt, row.names = NULL)
})
res <- do.call(rbind, rows)

res$cell <- snakemake@wildcards[["cell"]]
res$family <- cell$family
res$n_obs <- cell$n_obs
res$beta_t <- cell$beta_t
res$out_dev <- cell$out_dev
res$tau2 <- as.character(cell$tau2)
saveRDS(res, snakemake@output[[1]])

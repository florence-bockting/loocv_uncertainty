# Rule data: draws the training sets and the pooled test sets of one cell.
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])
snakemake@source("lib/dgp.R")

cell <- snakemake@params[["cell"]]
beta <- make_beta(cell$beta_t)
set.seed(snakemake@params[["seed"]])
# Draw order as in run_all.py: training sets first, then test sets.
train <- make_data_family(cell$family, snakemake@params[["n_trial"]],
                          cell$n_obs, snakemake@params[["n_obs_max"]],
                          beta, cell$out_dev)
# A binomial cell scores each test point by quadrature, which is far more
# costly than the closed-form normal score, so it uses fewer test sets.
# The Monte Carlo error of the target is then 1 / sqrt(n_test) of the LOO
# standard error, whatever n_obs is, because both shrink with n_obs.
# A poisson cell sums over a count grid on top of that quadrature, so it
# uses fewer test sets again.
n_test <- snakemake@params[["n_test"]]
n_test_family <- snakemake@params[[paste0("n_test_", cell$family)]]
if (!is.null(n_test_family)) n_test <- n_test_family
test <- pool_sets(make_data_family(cell$family, n_test, cell$n_obs,
                                   snakemake@params[["n_obs_max"]], beta,
                                   cell$out_dev))

saveRDS(train, snakemake@output[["train"]], compress = FALSE)
saveRDS(test, snakemake@output[["test"]], compress = FALSE)
saveRDS(list(y_test_var = mean((test$y - mean(test$y))^2)),
        snakemake@output[["info"]])

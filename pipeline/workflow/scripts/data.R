# Rule data: draws the training sets and the pooled test sets of one cell.
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])
snakemake@source("lib/dgp.R")

cell <- snakemake@params[["cell"]]
beta <- make_beta(cell$beta_t)
set.seed(snakemake@params[["seed"]])
# Draw order as in run_all.py: training sets first, then test sets.
train <- make_data(snakemake@params[["n_trial"]], cell$n_obs,
                   snakemake@params[["n_obs_max"]], beta, cell$out_dev)
test <- pool_sets(make_data(snakemake@params[["n_test"]], cell$n_obs,
                            snakemake@params[["n_obs_max"]], beta,
                            cell$out_dev))

saveRDS(train, snakemake@output[["train"]], compress = FALSE)
saveRDS(test, snakemake@output[["test"]], compress = FALSE)
saveRDS(list(y_test_var = mean((test$y - mean(test$y))^2)),
        snakemake@output[["info"]])

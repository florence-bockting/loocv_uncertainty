# Rule generate_data: draws the training sets and the pooled test sets of one cell.
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])
snakemake@source("lib/dgp.R")

cell <- snakemake@params[["cell"]]
beta <- make_beta_coeff(cell$beta_t)
set.seed(snakemake@params[["seed"]])
# Draw order as in run_all.py: training sets first, then test sets.
train <- make_data_family(
        family = cell$family,
        n_sets = snakemake@params[["n_sets_train"]],
        n_obs = cell$n_obs,
        n_obs_max = snakemake@params[["n_obs_max"]],
        beta = beta,
        out_dev = cell$out_dev,
        x_df = cell$x_df
)
# Each family has its own number of test sets, n_sets_test_<family>.
# A binomial cell scores each test point by quadrature, which is far more
# costly than the closed-form normal score, so it uses fewer test sets.
# The Monte Carlo error of the target is then 1 / sqrt(n_sets_test) of the
# LOO standard error, whatever n_obs is, because both shrink with n_obs.
# A poisson cell sums over a count grid on top of that quadrature, so it
# uses fewer test sets again.
n_sets_test <- snakemake@params[[paste0("n_sets_test_", cell$family)]]
test <- pool_sets(make_data_family(
        family = cell$family,
        n_sets = n_sets_test,
        n_obs = cell$n_obs,
        n_obs_max = snakemake@params[["n_obs_max"]],
        beta = beta,
        out_dev = cell$out_dev,
        x_df = cell$x_df
))

saveRDS(train, snakemake@output[["train"]], compress = FALSE)
saveRDS(test, snakemake@output[["test"]], compress = FALSE)
# needed for R2 computation
saveRDS(list(y_test_var = mean((test$y - mean(test$y))^2)),
        snakemake@output[["info"]])

# Rule fit_posterior_and_loo: fits one model to the training sets of one
# chunk of trials.
# Step 1 computes the full-data posterior, step 2 the LOO predictive.
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])
snakemake@source("lib/glm.R")

train <- readRDS(snakemake@input[["train"]])
trials <- seq(snakemake@params[["trials"]][[1]],
              snakemake@params[["trials"]][[2]])
cell <- snakemake@params[["cell"]]
cols <- model_columns(snakemake@wildcards[["model"]], dim(train$X)[2])

# Step 1: full-data posterior.
fits <- lapply(trials, function(t) {
  fit_model_family(
    family = cell$family,
    X = train$X[, cols, t],
    y = train$y[, t]
  )
})

# Step 2: LOO predictive, warm started at the full-data mode.
loos <- Map(function(t, fit) {
  loo_predictive_family(
    family = cell$family,
    X = train$X[, cols, t],
    y = train$y[, t],
    start = fit$beta_hat
  )
}, trials, fits)

saveRDS(list(trials = trials, fits = fits, loos = loos), snakemake@output[[1]])

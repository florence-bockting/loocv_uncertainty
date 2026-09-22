# Rule fit: fits one model to the training sets of one chunk of trials.
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])
snakemake@source("lib/predictive.R")

train <- readRDS(snakemake@input[["train"]])
trials <- seq(snakemake@params[["trials"]][[1]],
              snakemake@params[["trials"]][[2]])
tau2 <- parse_tau2(snakemake@params[["cell"]]$tau2)
cols <- model_columns(snakemake@wildcards[["model"]], dim(train$X)[2])

fits <- lapply(trials, function(t) {
  fit_model(train$X[, cols, t], train$y[, t], tau2)
})
saveRDS(list(trials = trials, fits = fits), snakemake@output[[1]])

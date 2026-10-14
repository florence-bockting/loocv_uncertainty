# Rule design_stats: a statistic of the training design for every trial
# of every gaussian cell, used by plot_joint_grouped.
#   mx2: mean(x2^2), drives the noise variance estimate of model A
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])

cells <- snakemake@params[["cells"]]
train_files <- setNames(unlist(snakemake@input[["train"]]), cells)

design <- do.call(rbind, lapply(cells, function(cell) {
  X <- readRDS(train_files[[cell]])$X
  trials <- seq_len(dim(X)[3])
  data.frame(
    cell = cell,
    trial = trials,
    mx2 = colMeans(X[, 3, ]^2)
  )
}))
saveRDS(design, snakemake@output[[1]])

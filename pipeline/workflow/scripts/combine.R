# Rule combine: stacks the trials of all cells and chunks into one table.
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])

res <- do.call(rbind, lapply(snakemake@input, readRDS))
rownames(res) <- NULL
saveRDS(res, snakemake@output[[1]])

# Rule compare_models: LOO estimate, SE and target of Ma - Mb for every measure,
# for every trial of every cell. Stacks all cells into one table.
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])
snakemake@source("lib/measures.R")

cells <- snakemake@params[["cells"]]
all_measures <- unlist(snakemake@params[["measures"]])

# The score files are in the order of the cells, then of the chunks
# (grid_files() in common.smk). train and info have one file per cell.
n_chunks <- vapply(cells, function(cell) as.integer(cell$n_chunks), integer(1))
chunk_cell <- factor(rep(names(cells), times = n_chunks), levels = names(cells))
a_files <- split(unlist(snakemake@input[["a"]]), chunk_cell)
b_files <- split(unlist(snakemake@input[["b"]]), chunk_cell)
train_files <- setNames(unlist(snakemake@input[["train"]]), names(cells))
info_files <- setNames(unlist(snakemake@input[["info"]]), names(cells))

compare_chunk <- function(a, b, train, y_test_var, measures) {
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
    data.frame(
      trial = a$trials[i],
      measure = measures,
      estimate = est[1, ],
      se = est[2, ],
      target = tgt,
      row.names = NULL
    )
  })
  do.call(rbind, rows)
}

compare_cell <- function(cell_name) {
  cell <- cells[[cell_name]]
  measures <- measures_for_family(measures = all_measures,
                                  family = cell$family)
  train <- readRDS(train_files[[cell_name]])
  y_test_var <- readRDS(info_files[[cell_name]])$y_test_var
  res <- do.call(rbind, Map(function(a_file, b_file) {
    compare_chunk(readRDS(a_file), readRDS(b_file), train, y_test_var,
                  measures)
  }, a_files[[cell_name]], b_files[[cell_name]]))
  res$cell <- cell_name
  res$family <- cell$family
  res$n_obs <- cell$n_obs
  res$beta_t <- cell$beta_t
  res$out_dev <- cell$out_dev
  res
}

res <- do.call(rbind, lapply(names(cells), compare_cell))
rownames(res) <- NULL
saveRDS(res, snakemake@output[[1]])

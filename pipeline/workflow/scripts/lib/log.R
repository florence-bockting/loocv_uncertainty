# Sends the messages, warnings and errors of a Snakemake R script to its
# log file.
start_log <- function(path) {
  con <- file(path, open = "wt")
  sink(con, type = "message")
  invisible(con)
}

# Converts the tau2 column of the cell table: "none" means unknown tau2.
parse_tau2 <- function(tau2) {
  if (identical(as.character(tau2), "none")) NULL else as.numeric(tau2)
}

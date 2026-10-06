# Sends the messages, warnings and errors of a Snakemake R script to its
# log file.
start_log <- function(path) {
  con <- file(path, open = "wt")
  sink(con, type = "message")
  invisible(con)
}

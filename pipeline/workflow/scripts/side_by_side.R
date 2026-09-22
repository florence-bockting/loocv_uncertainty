# Rule side_by_side: the paper figure (left) next to the rebuilt elpd
# figure (right), both rendered to PNG with pdftoppm.
snakemake@source("lib/log.R")
start_log(snakemake@log[[1]])

render <- function(pdf) {
  out <- tempfile()
  status <- system2("pdftoppm", c("-png", "-r", "100", "-singlefile",
                                  shQuote(pdf), shQuote(out)))
  stopifnot(status == 0)
  png::readPNG(paste0(out, ".png"))
}

imgs <- list(paper = render(snakemake@input[["paper"]]),
             rebuild = render(snakemake@input[["rebuild"]]))
titles <- c(paper = "Sivula et al. (2025)", rebuild = "R + Snakemake rebuild")

grDevices::pdf(snakemake@output[[1]], width = 16, height = 7)
par(mfrow = c(1, 2), mar = c(0, 0, 2, 0))
for (k in names(imgs)) {
  img <- imgs[[k]]
  asp <- dim(img)[1] / dim(img)[2]
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "",
       ylab = "", asp = 1)
  w <- min(1, 1 / asp)
  h <- w * asp
  rasterImage(img, (1 - w) / 2, (1 - h) / 2, (1 + w) / 2, (1 + h) / 2)
  title(main = titles[[k]])
}
dev.off()

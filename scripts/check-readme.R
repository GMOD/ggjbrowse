#!/usr/bin/env Rscript
# Runs the README's R blocks in order in one session, drawing each plot, so the
# README's "every example runs as written" stays true. The install block is
# skipped. Run from the repo root with the package installed.
blocks <- function(lines) {
  fence <- grep("^```", lines)
  starts <- fence[c(TRUE, FALSE)]
  ends <- fence[c(FALSE, TRUE)]
  r <- lines[starts] == "```r"
  Map(function(from, to) lines[seq_len(to - from - 1) + from], starts[r],
    ends[r])
}

readme <- blocks(readLines("README.md"))
env <- new.env()
setwd(tempdir())
grDevices::pdf(NULL)
for (code in readme) {
  if (any(grepl("install_github", code))) next
  message(code[1])
  for (expr in parse(text = code)) {
    value <- withVisible(eval(expr, env))
    if (value$visible && inherits(value$value, c("ggplot", "patchwork"))) {
      print(value$value)
    }
  }
}

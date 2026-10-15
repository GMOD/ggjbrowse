#!/usr/bin/env Rscript
# Runs the R blocks of the README and the guides in docs/ in order in one
# session, drawing each plot, so "every example runs as written" stays true.
# The install block is skipped. Run from the repo root with the package
# installed.
blocks <- function(lines) {
  fence <- grep("^```", lines)
  if (!length(fence)) return(list())
  starts <- fence[c(TRUE, FALSE)]
  ends <- fence[c(FALSE, TRUE)]
  r <- lines[starts] == "```r"
  Map(function(from, to) lines[seq_len(to - from - 1) + from], starts[r],
    ends[r])
}

pages <- c("README.md", file.path("docs",
  paste0(c("alignments", "grammar", "features", "wiggle", "variants", "stack",
    "synteny"), ".md")))
pages <- normalizePath(pages)
env <- new.env()
setwd(tempdir())
grDevices::pdf(NULL)
for (code in unlist(lapply(pages, function(p) blocks(readLines(p))),
  recursive = FALSE)) {
  if (any(grepl("install_github", code))) next
  message(code[1])
  for (expr in parse(text = code)) {
    value <- withVisible(eval(expr, env))
    if (value$visible && inherits(value$value, c("ggplot", "patchwork"))) {
      print(value$value)
    }
  }
}

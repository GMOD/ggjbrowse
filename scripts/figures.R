#!/usr/bin/env Rscript
# The README's quick-start figure, from the package's own test data: run after
# a change to what autoplot() draws, with the package installed.
# scripts/gallery.R draws the rest.
library(ggjbrowse)
library(ggplot2)

extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
gapdh <- autoplot(jbrowse_alignments(extdata("gapdh_direct_rna.bam"),
  "chr12:6,534,201-6,538,700"))
ggsave(file.path("man", "figures", "gapdh.png"), gapdh, width = 10,
  height = jbrowse_height(gapdh), dpi = 120)

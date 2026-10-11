#!/usr/bin/env Rscript
# The README's figures, from the package's own test data: run after a change
# to what autoplot() draws, with the package installed.
library(ggjbrowse)
library(ggplot2)

extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
out <- file.path("man", "figures")
save <- function(name, plot, height) {
  ggsave(file.path(out, name), plot, width = 10, height = height, dpi = 120)
}

spliced <- jbrowse_alignments(extdata("spliced.bam"), "ctgA", 400, 1300,
  reference = extdata("volvox.fa"), color = "strand")
save("spliced.png", autoplot(spliced), 6)

long <- jbrowse_alignments(extdata("volvox-long-reads-sv.bam"), "ctgA",
  1000, 2700, reference = extdata("volvox.fa"), color = "strand")
save("long-reads.png", autoplot(long), 5)

zoom <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA", 1000, 1120,
  reference = extdata("volvox.fa"), color = "mapq")
save("zoom.png", autoplot(zoom), 6)

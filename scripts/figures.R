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

meth <- jbrowse_alignments(extdata("methylation_clip.bam"), "20", 10000, 22000,
  reference = extdata("hg38_clip.fa"), base_color = "modifications",
  modifications = list(unmodified = "calls"))
save("methylation.png", autoplot(meth), 4)

genes <- jbrowse_features(extdata("volvox.sort.gff3.gz"), "ctgA", 1000, 10000,
  display = list(showOnlyGenes = TRUE))
reads <- jbrowse_alignments(extdata("spliced.bam"), "ctgA", 1000, 10000,
  reference = extdata("volvox.fa"), color = "strand")
save("genes-over-reads.png", autoplot(reads, features = genes), 7)

#!/usr/bin/env Rscript
# The README's figures, from the package's own test data: run after a change
# to what autoplot() draws, with the package installed.
library(ggjbrowse)
library(ggplot2)

extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
out <- file.path("man", "figures")
save <- function(name, plot, height = jbrowse_height(plot)) {
  ggsave(file.path(out, name), plot, width = 10, height = height, dpi = 120)
}

spliced <- jbrowse_alignments(extdata("spliced.bam"), "ctgA", 400, 1300,
  reference = extdata("volvox.fa"),
  display = list(color = list(field = "strand")))
save("spliced.png", autoplot(spliced))

long <- jbrowse_alignments(extdata("volvox-long-reads-sv.bam"), "ctgA",
  1000, 2700, reference = extdata("volvox.fa"),
  display = list(color = list(field = "strand")))
save("long-reads.png", autoplot(long))

zoom <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA", 1000, 1120,
  reference = extdata("volvox.fa"),
  display = list(color = list(field = "mapq")))
save("zoom.png", autoplot(zoom))

inv <- jbrowse_alignments(extdata("volvox-simple-inv-paired.bam"), "ctgA",
  18000, 32000,
  display = list(readConnections = "arc", facet = list(field = "strand")))
save("inversion.png", autoplot(inv))

meth <- jbrowse_alignments(extdata("methylation_clip.bam"), "20", 10000, 22000,
  reference = extdata("hg38_clip.fa"),
  display = list(baseColor = list(field = "modifications"),
    modifications = list(unmodified = "all")))
save("methylation.png", autoplot(meth))

genes <- jbrowse_features(extdata("volvox.sort.gff3.gz"), "ctgA", 1000, 10000,
  display = list(showOnlyGenes = TRUE))
reads <- jbrowse_alignments(extdata("spliced.bam"), "ctgA", 1000, 10000,
  reference = extdata("volvox.fa"),
  display = list(color = list(field = "strand")))
save("genes-over-reads.png", autoplot(reads, tracks = genes))

bw <- extdata("volvox_microarray_posneg.bw")
save("wiggle.png", patchwork::wrap_plots(ncol = 1,
  autoplot(jbrowse_wiggle(bw, "ctgA", 0, 50000, width = 100)),
  autoplot(jbrowse_wiggle(bw, "ctgA", 0, 50000, width = 100,
    display = list(mark = "line", interpolate = "linear"))),
  autoplot(jbrowse_wiggle(bw, "ctgA", 0, 50000,
    display = list(mark = "span")))), 6)

loci <- c("ctgA:1,001-4,000", "ctgA:6,001-9,000")
tracks <- list(
  jbrowse_features(extdata("volvox.sort.gff3.gz"), loci),
  jbrowse_wiggle(extdata("volvox_microarray_posneg.bw"), loci),
  jbrowse_alignments(extdata("volvox-sorted.bam"), loci,
    display = list(color = list(field = "strand"))))
save("stack.png", jbrowse_stack(tracks))

dogs <- jbrowse_variants(extdata("dog10k_cyp1a2_snvs.vcf.gz"), "chr30",
  38258000, 38265000)
save("variants.png", autoplot(dogs, layout = "columns"), 6)
save("variants-clustered.png",
  autoplot(jbrowse_cluster(dogs), layout = "columns"), 6)

paf <- extdata("volvox_inv_indels.paf")
synteny <- jbrowse_synteny(paf, "ctgA:1-49,186", "ctgA:1-50,001",
  view = list(color = list(field = "strand")))
target_genes <- jbrowse_features(extdata("volvox.sort.gff3.gz"), "ctgA:1-50,001",
  display = list(showOnlyGenes = TRUE))
save("synteny.png", autoplot(synteny, target = target_genes))
save("synteny-zoom.png", autoplot(jbrowse_synteny(paf, "ctgA:26,801-30,000",
  "ctgA:27,201-30,500", width = 1500)), 2.2)
ggsave(file.path(out, "dotplot.png"), autoplot(jbrowse_dotplot(paf,
  "ctgA:1-49,186", "ctgA:1-50,001", view = list(color = list(field = "strand")))),
  width = 6, height = 5.5, dpi = 120)

ggsave(file.path(out, "whole-genome.png"),
  autoplot(jbrowse_dotplot(extdata("peach_grape.paf.gz"))),
  width = 6, height = 5.5, dpi = 120)

multiway <- jbrowse_multiway(extdata("volvox_all_vs_all.paf"), "ctgA:1-50,001",
  "volvox", display = list(ribbonColor = list(field = "strand")))
save("multiway.png", autoplot(multiway), 3.2)

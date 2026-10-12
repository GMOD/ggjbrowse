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

gapdh <- "chr12:6,534,201-6,538,700"
rna <- jbrowse_alignments(extdata("gapdh_direct_rna.bam"), gapdh)
save("gapdh.png", autoplot(rna))
save("genes-over-reads.png", autoplot(rna,
  above = jbrowse_features(extdata("gapdh_gencode.gff3.gz"), gapdh)))

long <- jbrowse_alignments(extdata("colo829_tumor_chr3.bam"),
  "chr3:25,357,600-25,361,000",
  display = list(color = list(field = "strand")))
save("long-reads.png", autoplot(long))

deletion <- jbrowse_alignments(extdata("hg008t_cuzd1_deletion.bam"),
  "chr10:122,831,700-122,840,800",
  display = list(sort = list(type = "basePair", pos = 122836300)))
save("deletion.png", autoplot(deletion))

zoom <- jbrowse_alignments(extdata("hg008t_cuzd1_deletion.bam"),
  "chr10:122,835,300-122,835,420",
  display = list(sort = list(type = "basePair", pos = 122835344)))
save("zoom.png", autoplot(zoom))

inv <- jbrowse_alignments(extdata("hg00151_chr1_inversion.bam"),
  "chr1:197,786,900-197,789,700",
  display = list(readConnections = "arc", color = list(field = "strand"),
    facet = list(field = "split")))
save("inversion.png", autoplot(inv))

meth <- jbrowse_alignments(extdata("methylation_clip.bam"), "20", 10000, 22000,
  reference = extdata("hg38_clip.fa"),
  display = list(baseColor = list(field = "modifications"),
    modifications = list(unmodified = "all")))
save("methylation.png", autoplot(meth))

bw <- extdata("gapdh_phylop100way.bw")
save("wiggle.png", patchwork::wrap_plots(ncol = 1,
  autoplot(jbrowse_wiggle(bw, gapdh, width = 100)),
  autoplot(jbrowse_wiggle(bw, gapdh, width = 100,
    display = list(mark = "line", interpolate = "linear"))),
  autoplot(jbrowse_wiggle(bw, gapdh, display = list(mark = "span")))), 6)

cells <- c("gm12878", "k562", "hepg2", "huvec", "hmec", "nhek", "helas3",
  "a549", "mcf7", "hsmm", "hcm", "hct116")
dnase <- list(type = "MultiWiggleAdapter", subadapters = lapply(cells,
  function(cell) list(type = "BigWigAdapter", name = toupper(cell),
    bigWigLocation = list(uri = extdata(sprintf("hbb_dnase_%s.bw", cell))))))
hbb <- "chr11:5,270,001-5,320,000"
save("multi-wiggle.png", autoplot(jbrowse_wiggle(dnase, hbb)), 6)
clustered <- jbrowse_wiggle(dnase, hbb, display = list(mark = "span"))
clustered$tree <- hclust(dist(as.matrix(clustered)), "average")
save("multi-wiggle-clustered.png", autoplot(clustered), 4)

loci <- c("chr12:6,534,201-6,535,700", "chr12:6,537,201-6,538,700")
tracks <- list(
  jbrowse_features(extdata("gapdh_gencode.gff3.gz"), loci),
  jbrowse_wiggle(extdata("gapdh_phylop100way.bw"), loci),
  jbrowse_alignments(extdata("gapdh_direct_rna.bam"), loci,
    display = list(color = list(field = "strand"))))
save("stack.png", jbrowse_stack(tracks))

gff <- read.delim(gzfile(extdata("gapdh_gencode.gff3.gz")), header = FALSE,
  comment.char = "#")
exon <- gff[gff$V3 == "exon" &
  grepl("transcript_id=ENST00000229239.10;", gff$V9), ]
exon <- exon[order(exon$V4), ]
window <- 20
exons <- data.frame(chrom = exon$V1, start = exon$V4 - 1 - window,
  end = exon$V5 + window)
genes <- jbrowse_features(extdata("gapdh_gencode.gff3.gz"), exons)
reads <- jbrowse_alignments(extdata("gapdh_direct_rna.bam"), exons,
  display = list(color = list(field = "strand")))
save("collapsed-introns.png", jbrowse_stack(genes, reads))

people <- jbrowse_variants(extdata("1000g_dcst2.vcf.gz"),
  "1:155,000,001-155,050,000", samples = extdata("1000g_dcst2_samples.tsv"),
  display = list(rowColor = list(field = "superpopulation")))
save("variants.png", autoplot(people, layout = "columns"), 7)
people$tree <- hclust(dist(as.matrix(people)), "average")
save("variants-clustered.png", autoplot(people, layout = "columns"), 7)

chain <- extdata("hg002_chr8.chain")
inversion <- "5,250,001-14,250,000"
synteny <- jbrowse_synteny(chain, paste0("chr8_MATERNAL:", inversion),
  paste0("chr8_PATERNAL:", inversion),
  view = list(color = list(field = "strand")), swap = TRUE)
gene_track <- function(file, contig) jbrowse_features(extdata(file),
  paste0(contig, ":", inversion),
  display = list(displayMode = "compact", showLabels = FALSE,
    showOnlyGenes = TRUE))
synteny_plot <- autoplot(synteny,
  above = gene_track("hg002_mat_genes.gff3.gz", "chr8_MATERNAL"),
  below = gene_track("hg002_pat_genes.gff3.gz", "chr8_PATERNAL"))
save("synteny.png", synteny_plot)

amylase <- jbrowse_synteny(extdata("hprc_amylase.paf"),
  "GRCh38#0#chr1:103,520,894-103,832,637",
  "HG00133#1#CM090045.1:103,669,666-103,981,330",
  view = list(color = list(field = "strand")))
save("amylase.png", autoplot(amylase), 3)

save("synteny-zoom.png", autoplot(jbrowse_synteny(extdata("hprc_amylase.paf"),
  "GRCh38#0#chr1:103,619,501-103,622,500",
  "HG00133#1#CM090045.1:103,768,401-103,771,400", width = 1500)), 2.2)
ggsave(file.path(out, "dotplot.png"), autoplot(jbrowse_dotplot(chain,
  view = list(color = list(field = "strand")))),
  width = 7, height = 6, dpi = 120)

ggsave(file.path(out, "whole-genome.png"),
  autoplot(jbrowse_dotplot(extdata("peach_grape.paf.gz"))),
  width = 6, height = 5.5, dpi = 120)

liftover <- function(genome, assembly) list(type = "PairwiseIndexedPAFAdapter",
  pifGzLocation = list(uri = extdata(sprintf("hg38To%s.over.pif.gz", genome))),
  index = list(
    location = list(uri = extdata(sprintf("hg38To%s.over.pif.gz.csi", genome))),
    indexType = "CSI"),
  assemblyNames = c(assembly, "hg38"))
primates <- list(type = "MultiPairwiseSyntenyAdapter", adapters = list(
  liftover("PanTro6", "panTro6"), liftover("GorGor6", "gorGor6"),
  liftover("PonAbe3", "ponAbe3")))
multiway <- jbrowse_multiway(primates, "chr17:7,400,001-7,700,000", "hg38",
  lanes = c("panTro6", "gorGor6", "ponAbe3"),
  display = list(ribbonColor = list(field = "strand")))
save("multiway.png", autoplot(multiway), 4.5)

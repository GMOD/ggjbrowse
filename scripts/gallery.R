#!/usr/bin/env Rscript
# The README's gallery: each example drawn twice from the package's test data,
# by JBrowse's own renderer (`jb2export`, @jbrowse/img) and by ggjbrowse, side
# by side. Writes the images, the README's table and docs/gallery.md, which
# holds both recipes of each. Run from the repo root with the package installed
# and jbrowse-components built; name examples to redraw only those.
#   Rscript scripts/gallery.R [example ...]
# JBROWSE_COMPONENTS is the jbrowse-components checkout (default
# ~/src/jbrowse-components); JB2EXPORT replaces the command that runs
# jb2export. The hg38 and hg19 hubs and the primate chrom.sizes come over the
# network. Needs ImageMagick for the webp.
# Needs jb2export 5.0.0-beta or later: npm's latest, 4.3.0, cannot draw most of
# these. The clustered rows need jbrowse-components commit c50aaedcda, in no
# release yet.
suppressPackageStartupMessages({
  library(ggjbrowse)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a)) b else a
components <- Sys.getenv("JBROWSE_COMPONENTS", path.expand("~/src/jbrowse-components"))
jb2export <- strsplit(Sys.getenv("JB2EXPORT",
  paste("node", file.path(components, "products/jbrowse-img/esm/bin.js"))),
  " ")[[1]]
width <- 1200
dir <- "gallery"

# jb2export reads the checkout's copy of the data, and configs committed beside
# the images, so each command in docs/gallery.md runs from the repo root
inst <- function(f) file.path("inst", "extdata", f)
configs <- file.path(dir, "jb2export")
dir.create(configs, showWarnings = FALSE)
config_file <- function(name) file.path(configs, name)

# jb2export resolves a JSON config's localPath against the config's directory
local_file <- function(path) list(localPath = file.path("..", "..", path),
  locationType = "LocalPathLocation")

json_file <- function(x, name) {
  path <- config_file(name)
  writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE), path)
  path
}

# A chrom.sizes file of a pairwise file's contigs on one side, which jb2export
# draws an assembly from where no sequence is needed
chrom_sizes <- function(file, side, label, keep = NULL) {
  contigs <- jbrowse_contigs(file)
  contigs <- contigs[contigs$side == side, ]
  if (!is.null(keep)) contigs <- contigs[contigs$chrom %in% keep, ]
  path <- config_file(paste0(label, ".chrom.sizes"))
  write.table(contigs[c("chrom", "length")], path, quote = FALSE,
    sep = "\t", row.names = FALSE, col.names = FALSE)
  path
}

# A chrom.sizes file of one axis of a dotplot, its contigs in the order the
# plot lays them
axis_sizes <- function(layout, label) {
  path <- config_file(paste0(label, ".chrom.sizes"))
  write.table(data.frame(layout$chrom, layout$end - layout$start), path,
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  path
}

cells <- c("gm12878", "k562", "hepg2", "huvec", "hmec", "nhek", "helas3",
  "a549", "mcf7", "hsmm", "hcm", "hct116")

# Rows named for their cell line, as the multi-row track labels them
dnase_json <- function() json_file(lapply(cells, function(cell) list(
  type = "BigWigAdapter", name = toupper(cell),
  bigWigLocation = local_file(inst(sprintf("hbb_dnase_%s.bw", cell))))),
  "dnase.json")

fasta <- c("--fasta", inst("gapdh_hg38.fa.gz"))
hg38 <- c("--hub", "hg38")
grow <- c("heightMode:grow", "growMaxHeight=2000")
gapdh <- "chr12:6,534,201-6,538,700"
two_loci <- c("chr12:6,534,201-6,535,700", "chr12:6,537,201-6,538,700")
rna <- inst("gapdh_direct_rna.bam")
gencode <- inst("gapdh_gencode.gff3.gz")
phylop <- inst("gapdh_phylop100way.bw")
deletion <- inst("hg008t_cuzd1_deletion.bam")
inversion <- "5,250,001-14,250,000"
chain <- inst("hg002_chr8.chain")
amylase <- inst("hprc_amylase.paf")
peach_grape <- jbrowse_dotplot(inst("peach_grape.paf.gz"))
hg002_dot <- jbrowse_dotplot(chain, swap = TRUE)

chrom_sizes_assembly <- function(name, sizes) list(name = name,
  sequence = list(type = "ReferenceSequenceTrack",
    trackId = paste0(name, "-refseq"),
    adapter = list(type = "ChromSizesAdapter", chromSizesLocation = sizes)))

gene_track <- function(id, file, assembly) list(type = "FeatureTrack",
  trackId = id, name = id, assemblyNames = I(assembly),
  adapter = list(type = "Gff3TabixAdapter", gffGzLocation = local_file(file),
    index = list(location = local_file(paste0(file, ".tbi")),
      indexType = "TBI")),
  displays = list(list(type = "LinearBasicDisplay",
    displayId = paste0(id, "-display"), displayMode = "compact",
    showLabels = FALSE, showOnlyGenes = TRUE)))

synteny_genes_args <- function() {
  config <- json_file(list(
    assemblies = list(
      chrom_sizes_assembly("maternal",
        local_file(chrom_sizes(chain, "target", "maternal"))),
      chrom_sizes_assembly("paternal",
        local_file(chrom_sizes(chain, "query", "paternal")))),
    tracks = list(
      list(type = "SyntenyTrack", trackId = "hg002", name = "hg002",
        assemblyNames = c("maternal", "paternal"),
        adapter = list(type = "ChainAdapter",
          chainLocation = local_file(chain), targetAssembly = "maternal",
          queryAssembly = "paternal")),
      gene_track("maternal genes", inst("hg002_mat_genes.gff3.gz"),
        "maternal"),
      gene_track("paternal genes", inst("hg002_pat_genes.gff3.gz"),
        "paternal"))), "synteny_genes.config.json")
  spec <- json_file(list(views = list(list(type = "LinearSyntenyView",
    tracks = list(I("hg002")), color = list(field = "strand"),
    views = list(
      list(assembly = "maternal", loc = paste0("chr8_MATERNAL:", inversion),
        tracks = I("maternal genes")),
      list(assembly = "paternal", loc = paste0("chr8_PATERNAL:", inversion),
        tracks = I("paternal genes")))))), "synteny_genes.spec.json")
  c("--config", config, "--spec", spec)
}

multiway_args <- function(genomes) {
  assemblies <- lapply(c("hg38", genomes), function(genome)
    chrom_sizes_assembly(genome, list(uri = sprintf(
      "https://hgdownload.soe.ucsc.edu/goldenPath/%s/bigZips/%s.chrom.sizes",
      genome, genome))))
  liftover <- function(genome, assembly) {
    pif <- inst(sprintf("hg38To%s.over.pif.gz", genome))
    list(type = "PairwiseIndexedPAFAdapter", pifGzLocation = local_file(pif),
      index = list(location = local_file(paste0(pif, ".csi")),
        indexType = "CSI"),
      assemblyNames = c(assembly, "hg38"))
  }
  adapter <- list(type = "MultiPairwiseSyntenyAdapter", adapters = list(
    liftover("PanTro6", "panTro6"), liftover("GorGor6", "gorGor6"),
    liftover("PonAbe3", "ponAbe3")))
  config <- json_file(list(assemblies = assemblies, tracks = list(list(
    type = "SyntenyTrack", trackId = "primates", name = "liftOver to primates",
    assemblyNames = c("hg38", genomes), adapter = adapter,
    displays = list(list(type = "MultiWaySyntenyDisplay",
      displayId = "primates-display", height = 300,
      rows = list(domain = I(genomes)),
      ribbonColor = list(field = "strand")))))), "multiway.config.json")
  c("--config", config, "--assembly", "hg38", "--track", "primates")
}

population_tracks <- function() json_file(list(list(type = "VariantTrack",
  trackId = "people", name = "1000 Genomes", assemblyNames = I("hg19"),
  adapter = list(type = "VcfTabixAdapter",
    vcfGzLocation = local_file(inst("1000g_dcst2.vcf.gz")),
    index = list(location = local_file(inst("1000g_dcst2.vcf.gz.tbi"))),
    samplesTsvLocation = local_file(inst("1000g_dcst2_samples.tsv"))),
  displays = list(list(type = "LinearMultiSampleVariantDisplay",
    displayId = "people-display", height = 450, variantLayout = "columns",
    rowColor = "superpopulation")))), "people.tracks.json")

exon_locs <- function(transcript, window = 20) {
  gff <- read.delim(gzfile(gencode), header = FALSE, comment.char = "#")
  exon <- gff[gff$V3 == "exon" & grepl(paste0("transcript_id=", transcript, ";"),
    gff$V9, fixed = TRUE), ]
  exon <- exon[order(exon$V4), ]
  sprintf("%s:%d-%d", exon$V1, exon$V4 - window, exon$V5 + window)
}

# name -> `title`; `locus`, as the JBrowse location box writes it, for the
# caption; `jbrowse`, every jb2export argument; `r`, the ggjbrowse recipe,
# whose last value is the figure; `height`, its inches where the figure is not
# a track stack
examples <- list(
  genes_and_reads = list(
    title = "GENCODE genes over direct RNA reads", locus = gapdh,
    jbrowse = c(fasta, "--gffgz", gencode, "--bam", rna, "showLegend=true", grow,
      "--loc", gapdh),
    r = r"-(
reads <- jbrowse_alignments(e("gapdh_direct_rna.bam"),
  "chr12:6,534,201-6,538,700")
genes <- jbrowse_features(e("gapdh_gencode.gff3.gz"),
  "chr12:6,534,201-6,538,700")
autoplot(reads, above = genes)
)-"),
  long_reads_sv = list(
    title = "COLO829 der(3) split reads chained across its two chr3 breakpoints",
    locus = "chr3:25,358,500-25,359,800",
    jbrowse = c(hg38, "--bam", inst("colo829_tumor_chr3.bam"),
      "showLegend=true", grow,
      "arcs:down", "unit=chain", "color.field=strand", "filter.split=only",
      "--loc", "chr3:25,358,500-25,359,800"),
    r = r"-(
autoplot(jbrowse_alignments(e("colo829_tumor_chr3.bam"),
  "chr3:25,358,500-25,359,800",
  display = list(unit = "chain", readConnections = "arc",
    color = list(field = "strand"), filter = list(split = "only"))))
)-"),
  deletion = list(
    title = "Reads sorted by the base at the center, over a somatic deletion",
    locus = "chr10:122,831,700-122,840,800",
    jbrowse = c(hg38, "--bam", deletion, "showLegend=true", grow, "sort:base",
      "--loc", "chr10:122,831,700-122,840,800"),
    r = r"-(
autoplot(jbrowse_alignments(e("hg008t_cuzd1_deletion.bam"),
  "chr10:122,831,700-122,840,800",
  display = list(sort = list(type = "basePair", pos = 122836300))))
)-"),
  compact_reads = list(
    title = "Compact reads over the same deletion",
    locus = "chr10:122,831,700-122,840,800",
    jbrowse = c(hg38, "--bam", deletion, "showLegend=true", grow,
      "featureHeight:compact", "--loc", "chr10:122,831,700-122,840,800"),
    r = r"-(
autoplot(jbrowse_alignments(e("hg008t_cuzd1_deletion.bam"),
  "chr10:122,831,700-122,840,800", display = list(featureHeight = 3)))
)-"),
  base_zoom = list(
    title = "Base-level zoom at the deletion's left edge",
    locus = "chr10:122,835,300-122,835,420",
    jbrowse = c(hg38, "--bam", deletion, "showLegend=true", grow, "sort:base",
      "--loc", "chr10:122,835,300-122,835,420"),
    r = r"-(
autoplot(jbrowse_alignments(e("hg008t_cuzd1_deletion.bam"),
  "chr10:122,835,300-122,835,420",
  display = list(sort = list(type = "basePair", pos = 122835344))))
)-"),
  inversion = list(
    title = "Split reads chained across an inversion, with their arcs",
    locus = "chr1:197,786,900-197,789,700",
    jbrowse = c(hg38, "--bam", inst("hg00151_chr1_inversion.bam"),
      "showLegend=true", grow,
      "arcs:down", "unit=chain", "color.field=pairOrientation",
      "facet=split", "coverageHeight=80", "--loc",
      "chr1:197,786,900-197,789,700"),
    r = r"-(
autoplot(jbrowse_alignments(e("hg00151_chr1_inversion.bam"),
  "chr1:197,786,900-197,789,700",
  display = list(unit = "chain", readConnections = "arc",
    color = list(field = "pairOrientation"), facet = list(field = "split"))))
)-"),
  methylation = list(
    title = "CpG methylation from MM/ML tags", locus = "20:10,001-22,000",
    jbrowse = c("--fasta", inst("hg38_clip.fa"), "--bam",
      inst("methylation_clip.bam"), "baseColor:methylation",
      "showLegend=true", grow, "--loc", "20:10,001-22,000"),
    r = r"-(
autoplot(jbrowse_alignments(e("methylation_clip.bam"), "20:10,001-22,000",
  reference = e("hg38_clip.fa"),
  display = list(baseColor = list(field = "modifications"),
    modifications = list(unmodified = "all"))))
)-"),
  dnase_rows = list(
    title = "DNase in twelve cell lines over the beta-globin locus",
    locus = "chr11:5,270,001-5,320,000", height = 6,
    jbrowse = c(hg38, "--multiwig", dnase_json(), "showLegend=true",
      "name:DNase", "--loc", "chr11:5,270,001-5,320,000"),
    r = r"-(
cells <- c("gm12878", "k562", "hepg2", "huvec", "hmec", "nhek", "helas3",
  "a549", "mcf7", "hsmm", "hcm", "hct116")
dnase <- list(type = "MultiWiggleAdapter", subadapters = lapply(cells,
  function(cell) list(type = "BigWigAdapter", name = toupper(cell),
    bigWigLocation = list(uri = e(sprintf("hbb_dnase_%s.bw", cell))))))
autoplot(jbrowse_wiggle(dnase, "chr11:5,270,001-5,320,000"))
)-"),
  dnase_clustered = list(
    title = "DNase in twelve cell lines as density, clustered by signal",
    locus = "chr11:5,270,001-5,320,000", height = 4,
    jbrowse = c(hg38, "--multiwig", dnase_json(), "showLegend=true",
      "name:DNase", "mark=span", "runClustering=true",
      "--loc", "chr11:5,270,001-5,320,000"),
    r = r"-(
cells <- c("gm12878", "k562", "hepg2", "huvec", "hmec", "nhek", "helas3",
  "a549", "mcf7", "hsmm", "hcm", "hct116")
dnase <- list(type = "MultiWiggleAdapter", subadapters = lapply(cells,
  function(cell) list(type = "BigWigAdapter", name = toupper(cell),
    bigWigLocation = list(uri = e(sprintf("hbb_dnase_%s.bw", cell))))))
w <- jbrowse_wiggle(dnase, "chr11:5,270,001-5,320,000",
  display = list(mark = "span"))
w$tree <- hclust(dist(as.matrix(w)), "average")
autoplot(w)
)-"),
  conservation = list(
    title = "phyloP conservation across GAPDH as bars, a line and density",
    locus = gapdh, height = 6,
    jbrowse = c(fasta, "--bigwig", phylop, "showLegend=true",
      "--bigwig", phylop, "mark=line", "interpolate=linear",
      "--bigwig", phylop, "mark=span", "--loc", gapdh),
    r = r"-(
bw <- e("gapdh_phylop100way.bw")
gapdh <- "chr12:6,534,201-6,538,700"
patchwork::wrap_plots(ncol = 1,
  autoplot(jbrowse_wiggle(bw, gapdh, width = 1100)),
  autoplot(jbrowse_wiggle(bw, gapdh, width = 1100,
    display = list(mark = "line", interpolate = "linear"))),
  autoplot(jbrowse_wiggle(bw, gapdh, width = 1100,
    display = list(mark = "span"))))
)-"),
  two_regions = list(
    title = "Genes, conservation and strand-colored reads over two regions end to end",
    locus = paste(two_loci, collapse = " "),
    jbrowse = c(fasta, "--gffgz", gencode, "--bigwig", phylop, "showLegend=true",
      "--bam", rna, "showLegend=true", grow, "color.field=strand",
      "--loc", paste(two_loci, collapse = " ")),
    r = r"-(
loci <- c("chr12:6,534,201-6,535,700", "chr12:6,537,201-6,538,700")
genes <- jbrowse_features(e("gapdh_gencode.gff3.gz"), loci)
signal <- jbrowse_wiggle(e("gapdh_phylop100way.bw"), loci)
reads <- jbrowse_alignments(e("gapdh_direct_rna.bam"), loci,
  display = list(color = list(field = "strand")))
jbrowse_stack(genes, signal, reads)
)-"),
  variants_by_population = list(
    title = "100 people's genotypes across DCST2 (GRCh37), colored by superpopulation",
    locus = "1:155,000,001-155,050,000", height = 7,
    jbrowse = c("--hub", "hg19", "--tracks", population_tracks(), "--track",
      "people", "--loc", "1:155,000,001-155,050,000"),
    r = r"-(
v <- jbrowse_variants(e("1000g_dcst2.vcf.gz"), "1:155,000,001-155,050,000",
  samples = e("1000g_dcst2_samples.tsv"),
  display = list(rowColor = list(field = "superpopulation")))
autoplot(v, layout = "columns")
)-"),
  variants_clustered = list(
    title = "100 people's genotypes across DCST2 clustered, the tree beside the rows",
    locus = "1:155,000,001-155,050,000", height = 7,
    jbrowse = c("--hub", "hg19", "--tracks", population_tracks(), "--track",
      "people", "runClustering=true", "--loc", "1:155,000,001-155,050,000"),
    r = r"-(
v <- jbrowse_variants(e("1000g_dcst2.vcf.gz"), "1:155,000,001-155,050,000",
  samples = e("1000g_dcst2_samples.tsv"),
  display = list(rowColor = list(field = "superpopulation")))
v$tree <- hclust(dist(as.matrix(v)), "average")
autoplot(v, layout = "columns")
)-"),
  collapsed_introns = list(
    title = "GAPDH's nine exons side by side, with the junctions between them as arcs",
    locus = "ENST00000229239.10's exons, each padded by 20 bp", height = 7,
    jbrowse = c(fasta, "--gffgz", gencode, "--bam", rna, "sashimi:auto",
      "color.field=strand", "showLegend=true", grow, "--loc",
      paste(exon_locs("ENST00000229239.10"), collapse = " ")),
    r = r"-(
gff <- read.delim(gzfile(e("gapdh_gencode.gff3.gz")), header = FALSE,
  comment.char = "#")
exon <- gff[gff$V3 == "exon" &
  grepl("transcript_id=ENST00000229239.10;", gff$V9), ]
exon <- exon[order(exon$V4), ]
exons <- data.frame(chrom = exon$V1, start = exon$V4 - 1 - 20,
  end = exon$V5 + 20)
genes <- jbrowse_features(e("gapdh_gencode.gff3.gz"), exons)
reads <- jbrowse_alignments(e("gapdh_direct_rna.bam"), exons,
  display = list(color = list(field = "strand")))
jbrowse_stack(genes, reads)
)-"),
  inversion_synteny = list(
    title = "HG002's maternal chr8 over its paternal copy at 8p23.1",
    locus = sprintf("chr8_MATERNAL:%s chr8_PATERNAL:%s", inversion, inversion),
    jbrowse = c("synteny", "--chromSizes", chrom_sizes(chain, "target",
      "maternal"), paste0("loc:chr8_MATERNAL:", inversion), "--chain", chain,
      "--chromSizes", chrom_sizes(chain, "query", "paternal"),
      paste0("loc:chr8_PATERNAL:", inversion), "--colorBy", "strand"),
    r = r"-(
r <- "5,250,001-14,250,000"
autoplot(jbrowse_synteny(e("hg002_chr8.chain"), paste0("chr8_MATERNAL:", r),
  paste0("chr8_PATERNAL:", r),
  view = list(color = list(field = "strand")), swap = TRUE,
  labels = c(query = "paternal", target = "maternal")))
)-"),
  inversion_dotplot = list(
    title = "HG002's maternal chr8 against its paternal copy as a dotplot",
    locus = "chr8_MATERNAL chr8_PATERNAL", height = 6,
    jbrowse = c("dotplot", "--chromSizes", axis_sizes(hg002_dot$x, "maternal"),
      "--chain", chain, "--chromSizes", axis_sizes(hg002_dot$y, "paternal"),
      "--colorBy", "strand"),
    r = r"-(
autoplot(jbrowse_dotplot(e("hg002_chr8.chain"),
  view = list(color = list(field = "strand")), swap = TRUE,
  labels = c(query = "paternal", target = "maternal")))
)-"),
  whole_genome = list(
    title = "Peach against grape, every chromosome", locus = "whole genome",
    height = 5.5,
    jbrowse = c("dotplot", "--chromSizes", axis_sizes(peach_grape$x, "peach"),
      "--paf", inst("peach_grape.paf.gz"), "--chromSizes",
      axis_sizes(peach_grape$y, "grape")),
    r = r"-(
autoplot(jbrowse_dotplot(e("peach_grape.paf.gz"),
  labels = c(query = "peach", target = "grape")))
)-"),
  amylase_synteny = list(
    title = "GRCh38 over HG00133 at the amylase genes",
    locus = "GRCh38#0#chr1:103,520,894-103,832,637",
    jbrowse = c("synteny", "--chromSizes",
      chrom_sizes(amylase, "query", "GRCh38", "GRCh38#0#chr1"),
      "loc:GRCh38#0#chr1:103,520,894-103,832,637", "--paf", amylase,
      "--chromSizes", chrom_sizes(amylase, "target", "HG00133",
        "HG00133#1#CM090045.1"),
      "loc:HG00133#1#CM090045.1:103,669,666-103,981,330",
      "--colorBy", "strand"),
    r = r"-(
autoplot(jbrowse_synteny(e("hprc_amylase.paf"),
  "GRCh38#0#chr1:103,520,894-103,832,637",
  "HG00133#1#CM090045.1:103,669,666-103,981,330",
  view = list(color = list(field = "strand")),
  labels = c(query = "GRCh38", target = "HG00133")))
)-"),
  amylase_zoom = list(
    title = "3 kb at the amylase genes, the match block cut by an 80 bp deletion",
    locus = "GRCh38#0#chr1:103,619,501-103,622,500", height = 2.2,
    jbrowse = c("synteny", "--chromSizes",
      chrom_sizes(amylase, "query", "GRCh38", "GRCh38#0#chr1"),
      "loc:GRCh38#0#chr1:103,619,501-103,622,500", "--paf", amylase,
      "--chromSizes", chrom_sizes(amylase, "target", "HG00133",
        "HG00133#1#CM090045.1"),
      "loc:HG00133#1#CM090045.1:103,768,401-103,771,400",
      "--colorBy", "strand"),
    r = r"-(
autoplot(jbrowse_synteny(e("hprc_amylase.paf"),
  "GRCh38#0#chr1:103,619,501-103,622,500",
  "HG00133#1#CM090045.1:103,768,401-103,771,400",
  view = list(color = list(field = "strand")),
  labels = c(query = "GRCh38", target = "HG00133")))
)-"),
  synteny_genes = list(
    title = "HG002's maternal chr8 over its paternal copy, each with its genes",
    locus = sprintf("chr8_MATERNAL:%s chr8_PATERNAL:%s", inversion, inversion),
    jbrowse = synteny_genes_args(),
    r = r"-(
r <- "5,250,001-14,250,000"
s <- jbrowse_synteny(e("hg002_chr8.chain"), paste0("chr8_MATERNAL:", r),
  paste0("chr8_PATERNAL:", r),
  view = list(color = list(field = "strand")), swap = TRUE,
  labels = c(query = "paternal", target = "maternal"))
gene_track <- function(file, contig) jbrowse_features(e(file),
  paste0(contig, ":", r),
  display = list(displayMode = "compact", showLabels = FALSE,
    showOnlyGenes = TRUE))
autoplot(s, above = gene_track("hg002_mat_genes.gff3.gz", "chr8_MATERNAL"),
  below = gene_track("hg002_pat_genes.gff3.gz", "chr8_PATERNAL"))
)-"),
  multiway = list(
    title = "The TP53 neighbourhood on hg38 over chimp, gorilla and orangutan",
    locus = "chr17:7,400,001-7,700,000", height = 4.5,
    jbrowse = c(multiway_args(c("panTro6", "gorGor6", "ponAbe3")), "--loc",
      "chr17:7,400,001-7,700,000"),
    r = r"-(
liftover <- function(genome, assembly) list(type = "PairwiseIndexedPAFAdapter",
  pifGzLocation = list(uri = e(sprintf("hg38To%s.over.pif.gz", genome))),
  index = list(
    location = list(uri = e(sprintf("hg38To%s.over.pif.gz.csi", genome))),
    indexType = "CSI"),
  assemblyNames = c(assembly, "hg38"))
primates <- list(type = "MultiPairwiseSyntenyAdapter", adapters = list(
  liftover("PanTro6", "panTro6"), liftover("GorGor6", "gorGor6"),
  liftover("PonAbe3", "ponAbe3")))
autoplot(jbrowse_multiway(primates, "chr17:7,400,001-7,700,000", "hg38",
  lanes = c("panTro6", "gorGor6", "ponAbe3"),
  display = list(ribbonColor = list(field = "strand"))))
)-")
)
examples <- lapply(examples, function(x) {
  x$r <- trimws(x$r)
  x
})

png_of <- function(name, side) file.path(dir, sprintf("%s-%s.png", name, side))
webp_of <- function(name, side) sub("png$", "webp", png_of(name, side))

draw_jbrowse <- function(name, example) {
  out <- png_of(name, "jbrowse")
  status <- system2(jb2export[1], c(jb2export[-1],
    vapply(example$jbrowse, shQuote, ""), "--width", width, "--out", out),
    stdout = FALSE)
  if (status != 0 || !file.exists(out)) stop(name, ": jb2export failed")
}

draw_r <- function(name, example) {
  env <- new.env()
  env$e <- function(f) system.file("extdata", f, package = "ggjbrowse")
  figure <- eval(parse(text = example$r), env)
  height <- tryCatch(jbrowse_height(figure), error = function(e)
    example$height %||% 4)
  ggsave(png_of(name, "r"), figure, width = width / 120, height = height,
    dpi = 120)
}

to_webp <- function(name) {
  for (side in c("jbrowse", "r")) {
    system2("magick", c(png_of(name, side), "-quality", 85,
      webp_of(name, side)))
    file.remove(png_of(name, side))
  }
}

image_pair <- function(name, e, prefix) sprintf(
  "<tr><td width=\"50%%\"><img src=\"%s%s\" alt=\"%s in JBrowse\" width=\"100%%\"></td><td width=\"50%%\"><img src=\"%s%s\" alt=\"%s drawn by ggjbrowse\" width=\"100%%\"></td></tr>",
  prefix, webp_of(name, "jbrowse"), e$title, prefix, webp_of(name, "r"),
  e$title)

gallery_table <- function(examples) {
  row <- function(name, e) c(
    sprintf("<tr><td colspan=\"2\"><b><a href=\"docs/gallery.md#%s\">%s</a></b> <code>%s</code></td></tr>",
      name, e$title, e$locus),
    image_pair(name, e, ""))
  c("<!-- gallery -->", "<table>",
    "<tr><th width=\"50%\">JBrowse</th><th width=\"50%\">R</th></tr>",
    unlist(Map(row, names(examples), examples)), "</table>", "<!-- /gallery -->")
}

update_readme <- function(examples) {
  lines <- readLines("README.md")
  from <- grep("^<!-- gallery -->$", lines)
  to <- grep("^<!-- /gallery -->$", lines)
  stopifnot(length(from) == 1, length(to) == 1)
  writeLines(c(lines[seq_len(from - 1)], gallery_table(examples),
    lines[-seq_len(to)]), "README.md")
}

# Quoted where a shell would read it otherwise, as zsh reads `#` as a glob
shell_word <- function(x) ifelse(grepl("^[A-Za-z0-9_./:=,+-]+$", x), x,
  shQuote(x))

jb2export_command <- function(name, args) {
  words <- shell_word(c(args, "--width", width, "--out", paste0(name, ".png")))
  lines <- vapply(split(words, cumsum(startsWith(words, "--") &
    words != "--out" & seq_along(words) > 1)), paste, "", collapse = " ")
  paste0(c("jb2export ", rep("  ", length(lines) - 1)), lines,
    c(rep(" \\", length(lines) - 1), ""))
}

recipe <- function(name, e) c(
  sprintf("<a id=\"%s\"></a>", name), "",
  paste("##", e$title), "",
  sprintf("`%s`", e$locus), "",
  "<table>", image_pair(name, e, "../"), "</table>", "",
  "```r", e$r, "```", "",
  "```sh", jb2export_command(name, e$jbrowse), "```", "")

write_recipes <- function(examples) writeLines(c(
  "# Gallery recipes", "",
  "Each figure in the README's gallery, drawn from the package's own data by",
  "ggjbrowse and by JBrowse.", "",
  "- **ggjbrowse:** each R recipe runs as written after",
  "",
  "  ```r",
  "  library(ggjbrowse)",
  "  e <- function(f) system.file(\"extdata\", f, package = \"ggjbrowse\")",
  "  ```",
  "",
  "- **JBrowse:** each shell recipe runs `jb2export` from `@jbrowse/img`,",
  "  JBrowse's own renderer, from the root of a checkout of this repo. It needs",
  "  5.0.0-beta or later, and the clustered rows need jbrowse-components",
  "  commit `c50aaedcda`, in no release yet.",
  "- `scripts/gallery.R` draws both and writes this page, so edit an example",
  "  there.",
  "",
  unlist(Map(recipe, names(examples), examples))),
  file.path("docs", "gallery.md"))

wanted <- commandArgs(trailingOnly = TRUE)
todo <- if (length(wanted)) examples[wanted] else examples
for (name in names(todo)) {
  message(name)
  draw_jbrowse(name, todo[[name]])
  draw_r(name, todo[[name]])
  to_webp(name)
}
update_readme(examples)
write_recipes(examples)

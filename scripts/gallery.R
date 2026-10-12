#!/usr/bin/env Rscript
# The README's gallery: each example drawn twice from the package's test data,
# by JBrowse's own renderer (`jb2export`, @jbrowse/img) and by ggjbrowse, side
# by side. Run with the package installed and jbrowse-components built; name
# examples to redraw only those.
#   Rscript scripts/gallery.R [example ...]
# JBROWSE_COMPONENTS is the jbrowse-components checkout (default
# ~/src/jbrowse-components); JB2EXPORT replaces the command that runs
# jb2export. The hg38 and hg19 hubs and the clips' reference sequences come
# over the network. Needs ImageMagick for the webp.
suppressPackageStartupMessages({
  library(ggjbrowse)
  library(ggplot2)
})

extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
`%||%` <- function(a, b) if (is.null(a)) b else a
components <- Sys.getenv("JBROWSE_COMPONENTS", path.expand("~/src/jbrowse-components"))
jb2export <- strsplit(Sys.getenv("JB2EXPORT",
  paste("node", file.path(components, "products/jbrowse-img/esm/bin.js"))),
  " ")[[1]]
width <- 1200
dir <- "gallery"
work <- tempfile("gallery")
dir.create(work)

# A chrom.sizes file of a pairwise file's contigs on one side, which jb2export
# draws an assembly from where no sequence is needed
chrom_sizes <- function(file, side, label, keep = NULL) {
  contigs <- jbrowse_contigs(file)
  contigs <- contigs[contigs$side == side, ]
  if (!is.null(keep)) contigs <- contigs[contigs$chrom %in% keep, ]
  path <- file.path(work, paste0(label, ".chrom.sizes"))
  write.table(contigs[c("chrom", "length")], path, quote = FALSE,
    sep = "\t", row.names = FALSE, col.names = FALSE)
  path
}

# A chrom.sizes file of one axis of a dotplot, its contigs in the order the
# plot lays them
axis_sizes <- function(layout, label) {
  path <- file.path(work, paste0(label, ".chrom.sizes"))
  write.table(data.frame(layout$chrom, layout$end - layout$start), path,
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  path
}

# Rows named for their cell line, as the multi-row track labels them
rows_json <- function(files, names) {
  path <- tempfile(fileext = ".json", tmpdir = work)
  writeLines(jsonlite::toJSON(unname(Map(function(f, n) list(type = "BigWigAdapter",
    name = n, bigWigLocation = list(localPath = f,
      locationType = "LocalPathLocation")), files, names)), auto_unbox = TRUE,
    pretty = TRUE), path)
  path
}

fasta <- c("--fasta", extdata("gapdh_hg38.fa.gz"))
hg38 <- c("--hub", "hg38")
gapdh <- "chr12:6,534,201-6,538,700"
two_loci <- c("chr12:6,534,201-6,535,700", "chr12:6,537,201-6,538,700")
rna <- extdata("gapdh_direct_rna.bam")
gencode <- extdata("gapdh_gencode.gff3.gz")
phylop <- extdata("gapdh_phylop100way.bw")
deletion <- extdata("hg008t_cuzd1_deletion.bam")
cells <- c("gm12878", "k562", "hepg2", "huvec", "hmec", "nhek", "helas3",
  "a549", "mcf7", "hsmm", "hcm", "hct116")
dnase <- extdata(sprintf("hbb_dnase_%s.bw", cells))
inversion <- "5,250,001-14,250,000"
chain <- extdata("hg002_chr8.chain")
amylase <- extdata("hprc_amylase.paf")
peach_grape <- jbrowse_dotplot(extdata("peach_grape.paf.gz"))

# name -> `title`; `locus`, as the JBrowse location box writes it, for the
# caption; `jbrowse`, every jb2export argument; `r`, the ggjbrowse figure;
# `height`, its inches where the figure is not a track stack
examples <- list(
  genes_and_reads = list(
    title = "GENCODE genes over direct RNA reads", locus = gapdh,
    jbrowse = c(fasta, "--gffgz", gencode, "--bam", rna, "--loc", gapdh),
    r = function() autoplot(jbrowse_alignments(rna, gapdh),
      above = jbrowse_features(gencode, gapdh))),
  sashimi = list(
    title = "Splice junctions over direct RNA reads", locus = gapdh,
    jbrowse = c(fasta, "--bam", rna, "--loc", gapdh),
    r = function() autoplot(jbrowse_alignments(rna, gapdh))),
  strand_sections = list(
    title = "Reads split into strand sections", locus = gapdh,
    jbrowse = c(fasta, "--bam", rna, "facet=strand", "--loc", gapdh),
    r = function() autoplot(jbrowse_alignments(rna, gapdh,
      display = list(facet = list(field = "strand"))))),
  long_reads_sv = list(
    title = "Long reads across two breakpoints of the COLO829 der(3)",
    locus = "chr3:25,357,600-25,361,000",
    jbrowse = c(hg38, "--bam", extdata("colo829_tumor_chr3.bam"),
      "color.field=strand", "--loc", "chr3:25,357,600-25,361,000"),
    r = function() autoplot(jbrowse_alignments(
      extdata("colo829_tumor_chr3.bam"), "chr3:25,357,600-25,361,000",
      display = list(color = list(field = "strand"))))),
  deletion = list(
    title = "Reads sorted by the base at the center, over a somatic deletion",
    locus = "chr10:122,831,700-122,840,800",
    jbrowse = c(hg38, "--bam", deletion, "sort:base", "--loc",
      "chr10:122,831,700-122,840,800"),
    r = function() autoplot(jbrowse_alignments(deletion,
      "chr10:122,831,700-122,840,800",
      display = list(sort = list(type = "basePair", pos = 122836300))))),
  compact_reads = list(
    title = "Compact reads over the same deletion",
    locus = "chr10:122,831,700-122,840,800",
    jbrowse = c(hg38, "--bam", deletion, "featureHeight:compact", "--loc",
      "chr10:122,831,700-122,840,800"),
    r = function() autoplot(jbrowse_alignments(deletion,
      "chr10:122,831,700-122,840,800", display = list(featureHeight = 3)))),
  base_zoom = list(
    title = "Base-level zoom at the deletion's left edge",
    locus = "chr10:122,835,300-122,835,420",
    jbrowse = c(hg38, "--bam", deletion, "sort:base", "--loc",
      "chr10:122,835,300-122,835,420"),
    r = function() autoplot(jbrowse_alignments(deletion,
      "chr10:122,835,300-122,835,420",
      display = list(sort = list(type = "basePair", pos = 122835344))))),
  inversion = list(
    title = "Split reads chained across an inversion, with their arcs",
    locus = "chr1:197,786,900-197,789,700",
    jbrowse = c(hg38, "--bam", extdata("hg00151_chr1_inversion.bam"),
      "arcs:down", "unit=chain", "color.field=pairOrientation",
      "facet=split", "coverageHeight=80", "--loc",
      "chr1:197,786,900-197,789,700"),
    r = function() autoplot(jbrowse_alignments(
      extdata("hg00151_chr1_inversion.bam"), "chr1:197,786,900-197,789,700",
      display = list(unit = "chain", readConnections = "arc",
        color = list(field = "pairOrientation"),
        facet = list(field = "split"))))),
  methylation = list(
    title = "CpG methylation from MM/ML tags", locus = "20:10,001-22,000",
    jbrowse = c("--fasta", extdata("hg38_clip.fa"), "--bam",
      extdata("methylation_clip.bam"), "baseColor:methylation",
      "showLegend=true", "--loc", "20:10,001-22,000"),
    r = function() autoplot(jbrowse_alignments(extdata("methylation_clip.bam"),
      "20:10,001-22,000", reference = extdata("hg38_clip.fa"),
      display = list(baseColor = list(field = "modifications"),
        modifications = list(unmodified = "all"))))),
  dnase_rows = list(
    title = "DNase in twelve cell lines over the beta-globin locus",
    locus = "chr11:5,270,001-5,320,000", height = 6,
    jbrowse = c(hg38, "--multiwig", rows_json(dnase, toupper(cells)),
      "name:DNase", "--loc", "chr11:5,270,001-5,320,000"),
    r = function() autoplot(jbrowse_wiggle(list(type = "MultiWiggleAdapter",
      subadapters = lapply(seq_along(cells), function(i) list(
        type = "BigWigAdapter", name = toupper(cells[i]),
        bigWigLocation = list(uri = dnase[i])))), "chr11:5,270,001-5,320,000"))),
  conservation = list(
    title = "phyloP conservation across GAPDH", locus = gapdh, height = 3,
    jbrowse = c(fasta, "--bigwig", phylop, "--loc", gapdh),
    r = function() autoplot(jbrowse_wiggle(phylop, gapdh, width = 100))),
  two_regions = list(
    title = "Genes, conservation and reads over two regions end to end",
    locus = paste(two_loci, collapse = " "),
    jbrowse = c(fasta, "--gffgz", gencode, "--bigwig", phylop, "--bam", rna,
      "--loc", paste(two_loci, collapse = " ")),
    r = function() jbrowse_stack(jbrowse_features(gencode, two_loci),
      jbrowse_wiggle(phylop, two_loci), jbrowse_alignments(rna, two_loci))),
  genotype_matrix = list(
    title = "100 people's genotypes across DCST2 (GRCh37)",
    locus = "1:155,000,001-155,050,000", height = 7,
    jbrowse = c("--hub", "hg19", "--vcfgz", extdata("1000g_dcst2.vcf.gz"),
      "display:multivariantmatrix", "height:450", "--loc",
      "1:155,000,001-155,050,000"),
    r = function() autoplot(jbrowse_variants(extdata("1000g_dcst2.vcf.gz"),
      "1:155,000,001-155,050,000",
      samples = extdata("1000g_dcst2_samples.tsv")), layout = "columns")),
  inversion_synteny = list(
    title = "HG002's maternal chr8 over its paternal copy at 8p23.1",
    locus = sprintf("chr8_MATERNAL:%s chr8_PATERNAL:%s", inversion, inversion),
    jbrowse = c("synteny", "--chromSizes", chrom_sizes(chain, "target",
      "maternal"), paste0("loc:chr8_MATERNAL:", inversion), "--chain", chain,
      "--chromSizes", chrom_sizes(chain, "query", "paternal"),
      paste0("loc:chr8_PATERNAL:", inversion), "--colorBy", "strand"),
    r = function() autoplot(jbrowse_synteny(chain,
      paste0("chr8_PATERNAL:", inversion), paste0("chr8_MATERNAL:", inversion),
      view = list(color = list(field = "strand")), swap = TRUE,
      labels = c(query = "paternal", target = "maternal")))),
  whole_genome = list(
    title = "Peach against grape, every chromosome", locus = "whole genome",
    height = 5.5,
    jbrowse = c("dotplot", "--chromSizes", axis_sizes(peach_grape$x, "peach"),
      "--paf", extdata("peach_grape.paf.gz"), "--chromSizes",
      axis_sizes(peach_grape$y, "grape")),
    r = function() autoplot(jbrowse_dotplot(extdata("peach_grape.paf.gz"),
      labels = c(query = "peach", target = "grape")))),
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
    r = function() autoplot(jbrowse_synteny(amylase,
      "GRCh38#0#chr1:103,520,894-103,832,637",
      "HG00133#1#CM090045.1:103,669,666-103,981,330",
      view = list(color = list(field = "strand")),
      labels = c(query = "GRCh38", target = "HG00133"))))
)

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
  figure <- example$r()
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

gallery_table <- function(examples) {
  row <- function(name, e) c(
    sprintf("<tr><td colspan=\"2\"><b>%s</b> <code>%s</code></td></tr>",
      e$title, e$locus),
    sprintf("<tr><td width=\"50%%\"><img src=\"gallery/%s-jbrowse.webp\" alt=\"%s in JBrowse\" width=\"100%%\"></td><td width=\"50%%\"><img src=\"gallery/%s-r.webp\" alt=\"%s drawn by ggjbrowse\" width=\"100%%\"></td></tr>",
      name, e$title, name, e$title))
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

wanted <- commandArgs(trailingOnly = TRUE)
todo <- if (length(wanted)) examples[wanted] else examples
for (name in names(todo)) {
  message(name)
  draw_jbrowse(name, todo[[name]])
  draw_r(name, todo[[name]])
  to_webp(name)
}
update_readme(examples)

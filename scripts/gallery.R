#!/usr/bin/env Rscript
# The README's gallery: each example drawn twice from the package's test data,
# by JBrowse's own renderer (`jb2export`, @jbrowse/img) and by ggjbrowse, side
# by side. Writes the images, the README's table and docs/gallery.md, which
# holds both recipes of each. Run from the repo root with the package installed.
#   Rscript scripts/gallery.R [example ...]  redraw all, or the examples named
#   Rscript scripts/gallery.R --docs         rewrite the table and docs only
#   Rscript scripts/gallery.R --check        run each R recipe, then as --docs
# Drawing needs jbrowse-components built, at JBROWSE_COMPONENTS (default
# ~/src/jbrowse-components), or JB2EXPORT naming the command that runs
# jb2export; `jb2export_version` below says which. The hg38 and hg19 hubs and
# the primate chrom.sizes come over the network. Needs ImageMagick for the webp.
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
jb2export_version <- c(
  "It needs 5.0.0-beta or later, since 4.3.0 cannot draw most of these, and",
  "the clustered rows and the alignments' band settings, such as",
  "`coverage.height=`, need jbrowse-components commit `586faa2bd3`, in no",
  "release yet.")

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
  contigs <- contigs[contigs$assembly == side, ]
  if (!is.null(keep)) contigs <- contigs[contigs$chrom %in% keep, ]
  path <- config_file(paste0(label, ".chrom.sizes"))
  write.table(contigs[c("chrom", "end")], path, quote = FALSE,
    sep = "\t", row.names = FALSE, col.names = FALSE)
  path
}

# A chrom.sizes file of one axis of a dotplot, its contigs in the order the
# plot lays them
axis_sizes <- function(result, axis, label) {
  layout <- result$regions[result$regions$axis == axis, ]
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

methyldackel_json <- function() json_file(
  lapply(c("CpG", "CHG", "CHH"), function(context) list(
    type = "BigWigAdapter", name = context, bigWigLocation = local_file(
      inst(sprintf("arabidopsis_methyldackel_%s.bw", context))))),
  "methyldackel.json")

fasta <- c("--fasta", inst("gapdh_hg38.fa.gz"))
hg38 <- c("--hub", "hg38")
grow <- c("heightMode:grow", "growMaxHeight=2000")
gapdh <- "chr12:6,534,201-6,538,700"
two_loci <- c("chr12:6,534,201-6,535,700", "chr12:6,537,201-6,538,700")
aa_loci <- c("chr12:6,534,811-6,534,880", "chr12:6,536,476-6,536,545")
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
    showLabels = FALSE, filter = list("geneLike"))))

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
    tracks = list(I("hg002")), encoding = list(color = list(field = "strand")),
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
      rows = list(sort = I(genomes)),
      ribbons = list(encoding = list(color = list(field = "strand")))))))),
    "multiway.config.json")
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

# name -> `title`; `caption`, what to see in it; `locus`, as the JBrowse
# location box writes it; `jbrowse`, every jb2export argument; `r`, the ggjbrowse recipe,
# whose last value is the figure; `height`, its inches where the figure is not
# a track stack
examples <- list(
  genes_and_reads = list(
    caption = "Junction arcs join each splice over the coverage band, and the reads split at every intron.",
    title = "GENCODE genes over direct RNA reads", locus = gapdh,
    jbrowse = c(fasta, "--gffgz", gencode, "--bam", rna, "legend.show=true", grow,
      "--loc", gapdh),
    r = r"-(
reads <- jbrowse_alignments(e("gapdh_direct_rna.bam"),
  "chr12:6,534,201-6,538,700")
genes <- jbrowse_features(e("gapdh_gencode.gff3.gz"),
  "chr12:6,534,201-6,538,700")
jbrowse_stack(genes, reads)
)-"),
  long_reads_sv = list(
    caption = "Coverage halves at the first breakpoint, and each read's segments share a row, overlapping in dark grey and colored by strand past it, joined by duplication-type arcs.",
    title = "COLO829 der(3) split reads chained across its two chr3 breakpoints",
    locus = "chr3:25,358,500-25,359,800",
    jbrowse = c(hg38, "--bam", inst("colo829_tumor_chr3.bam"),
      "legend.show=true", grow,
      "arcs:down", "unit=chain", "encoding.color.field=strand", "filter.split=only",
      "--loc", "chr3:25,358,500-25,359,800"),
    r = r"-(
autoplot(jbrowse_alignments(e("colo829_tumor_chr3.bam"),
  "chr3:25,358,500-25,359,800",
  display = list(unit = "chain", readConnections = list(show = TRUE), encoding = list(color = list(field = "strand")), filter = list(list(field = "split", equal = TRUE)))))
)-"),
  deletion = list(
    caption = "Coverage halves across the deletion, and the reads sorted to the top span it with dark deletion bars, each labelled with its length.",
    title = "Reads sorted by the base at the center, over a somatic deletion",
    locus = "chr10:122,831,700-122,840,800",
    jbrowse = c(hg38, "--bam", deletion, "legend.show=true", grow, "sort:base",
      "--loc", "chr10:122,831,700-122,840,800"),
    r = r"-(
autoplot(jbrowse_alignments(e("hg008t_cuzd1_deletion.bam"),
  "chr10:122,831,700-122,840,800",
  display = list(sort = list(field = "base", pos = 122836300))))
)-"),
  compact_reads = list(
    caption = "Sorted by the base at the center, at 3 pixels a row the pileup takes far less height, the deletion's dark bars banded through it.",
    title = "Compact reads over the same deletion",
    locus = "chr10:122,831,700-122,840,800",
    jbrowse = c(hg38, "--bam", deletion, "legend.show=true", grow, "sort:base",
      "featureHeight:compact", "--loc", "chr10:122,831,700-122,840,800"),
    r = r"-(
autoplot(jbrowse_alignments(e("hg008t_cuzd1_deletion.bam"),
  "chr10:122,831,700-122,840,800",
  display = list(pileup = list(mark = list(height = 3)),
    sort = list(field = "base", pos = 122836300))))
)-"),
  base_zoom = list(
    caption = "At 120 bp the pileup writes out each mismatch's base, and every read sorted to the top carries a G a few bases before its deletion.",
    title = "Base-level zoom at the deletion's left edge",
    locus = "chr10:122,835,300-122,835,420",
    jbrowse = c(hg38, "--bam", deletion, "legend.show=true", grow, "sort:base",
      "--loc", "chr10:122,835,300-122,835,420"),
    r = r"-(
autoplot(jbrowse_alignments(e("hg008t_cuzd1_deletion.bam"),
  "chr10:122,835,300-122,835,420",
  display = list(sort = list(field = "base", pos = 122835344))))
)-"),
  inversion = list(
    caption = "The split reads, in their own section above the unsplit ones, are colored inverted between the two breakpoints, where inversion-type arcs join their segments.",
    title = "Split reads chained across an inversion, with their arcs",
    locus = "chr1:197,786,900-197,789,700",
    jbrowse = c(hg38, "--bam", inst("hg00151_chr1_inversion.bam"),
      "legend.show=true", grow,
      "arcs:down", "unit=chain", "encoding.color.field=pairOrientation",
      "facet=split", "coverage.height=80", "--loc",
      "chr1:197,786,900-197,789,700"),
    r = r"-(
autoplot(jbrowse_alignments(e("hg00151_chr1_inversion.bam"),
  "chr1:197,786,900-197,789,700",
  display = list(unit = "chain", readConnections = list(show = TRUE), encoding = list(color = list(field = "pairOrientation")), facet = list(field = "split"),
    coverage = list(height = 80))))
)-"),
  methylation = list(
    caption = "Each read's 5mC calls in red where methylated and blue where not.",
    title = "CpG methylation from MM/ML tags", locus = "20:10,001-22,000",
    jbrowse = c("--fasta", inst("hg38_clip.fa"), "--bam",
      inst("methylation_clip.bam"), "baseColor:methylation",
      "legend.show=true", grow, "--loc", "20:10,001-22,000"),
    r = r"-(
autoplot(jbrowse_alignments(e("methylation_clip.bam"), "20:10,001-22,000",
  assembly = e("hg38_clip.fa"),
  display = list(pileup = list(bases = list(show = TRUE,
    encoding = list(color = "modifications"), unmodified = list(show = TRUE))))))
)-"),
  bisulfite = list(
    caption = "RefSeq genes and RepeatMasker over MethylDackel's three context bigWigs and each read's methylated CHH cytosines in red: CHH is dense over the LTR/Copia element on the right and rare over the gene on the left, which carries CpG methylation only.",
    title = "CHH methylation from bisulfite reads",
    locus = "NC_003070.9:4,398,001-4,412,000",
    jbrowse = c("--fasta", inst("arabidopsis_tair10.fa.gz"), "--bigbed",
      inst("arabidopsis_refseq.bb"), "name:RefSeq genes", "--bigbed",
      inst("arabidopsis_rmsk.bb"), "name:RepeatMasker", "--multiwig",
      methyldackel_json(), "name:MethylDackel", "legend.show=true", "--bam",
      inst("arabidopsis_wgbs.bam"), "baseColor:bisulfite",
      '{"pileup":{"bases":{"filter":[{"field":"context","oneOf":["CHH"]}]}}}',
      "legend.show=true", grow, "--loc",
      "NC_003070.9:4,398,001-4,412,000"),
    r = r"-(
locus <- "NC_003070.9:4,398,001-4,412,000"
genes <- jbrowse_features(e("arabidopsis_refseq.bb"), locus)
repeats <- jbrowse_features(e("arabidopsis_rmsk.bb"), locus)
methyldackel <- list(type = "MultiWiggleAdapter",
  subadapters = lapply(c("CpG", "CHG", "CHH"), function(context)
    list(type = "BigWigAdapter", name = context, bigWigLocation = list(
      uri = e(sprintf("arabidopsis_methyldackel_%s.bw", context))))))
reads <- jbrowse_alignments(e("arabidopsis_wgbs.bam"), locus,
  assembly = e("arabidopsis_tair10.fa.gz"),
  display = list(pileup = list(bases = list(show = TRUE,
    encoding = list(color = "bisulfite"),
    filter = list(list(field = "context", oneOf = "CHH"))))))
jbrowse_stack(`RefSeq genes` = genes, RepeatMasker = repeats,
  MethylDackel = jbrowse_wiggle(methyldackel, locus), `CHH reads` = reads)
)-"),
  dnase_rows = list(
    caption = "Most lines share one peak near 5,281,000, and K562 alone shows several more to either side.",
    title = "DNase in twelve cell lines over the beta-globin locus",
    locus = "chr11:5,270,001-5,320,000", height = 6,
    jbrowse = c(hg38, "--multiwig", dnase_json(), "legend.show=true",
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
    caption = "K562 splits off alone from the other eleven lines.",
    title = "DNase in twelve cell lines as density, clustered by signal",
    locus = "chr11:5,270,001-5,320,000", height = 4,
    jbrowse = c(hg38, "--multiwig", dnase_json(), "legend.show=true",
      "name:DNase", "mark=rect", "runClustering=true",
      "--loc", "chr11:5,270,001-5,320,000"),
    r = r"-(
cells <- c("gm12878", "k562", "hepg2", "huvec", "hmec", "nhek", "helas3",
  "a549", "mcf7", "hsmm", "hcm", "hct116")
dnase <- list(type = "MultiWiggleAdapter", subadapters = lapply(cells,
  function(cell) list(type = "BigWigAdapter", name = toupper(cell),
    bigWigLocation = list(uri = e(sprintf("hbb_dnase_%s.bw", cell))))))
w <- jbrowse_wiggle(dnase, "chr11:5,270,001-5,320,000",
  display = list(mark = "rect"))
w$tree <- hclust(dist(as.matrix(w)), "average")
autoplot(w)
)-"),
  conservation = list(
    caption = "The exons stand out as runs of high scores in bars, in a mean line over its min-max band, and in a density strip.",
    title = "phyloP conservation across GAPDH as bars, a line and density",
    locus = gapdh, height = 6,
    jbrowse = c(fasta, "--bigwig", phylop, "legend.show=true",
      "--bigwig", phylop, "mark.type=line", "mark.interpolate=linear",
      "--bigwig", phylop, "mark=rect", "--loc", gapdh),
    r = r"-(
bw <- e("gapdh_phylop100way.bw")
gapdh <- "chr12:6,534,201-6,538,700"
patchwork::wrap_plots(ncol = 1,
  autoplot(jbrowse_wiggle(bw, gapdh, width = 1100)),
  autoplot(jbrowse_wiggle(bw, gapdh, width = 1100,
    display = list(mark = list(type = "line", interpolate = "linear")))),
  autoplot(jbrowse_wiggle(bw, gapdh, width = 1100,
    display = list(mark = "rect"))))
)-"),
  two_regions = list(
    caption = "A grey divider separates the two loci, and each tick gives a position in its own region.",
    title = "Genes, conservation and strand-colored reads over two regions end to end",
    locus = paste(two_loci, collapse = " "),
    jbrowse = c(fasta, "--gffgz", gencode, "--bigwig", phylop, "legend.show=true",
      "--bam", rna, "legend.show=true", grow, "encoding.color.field=strand",
      "--loc", paste(two_loci, collapse = " ")),
    r = r"-(
loci <- c("chr12:6,534,201-6,535,700", "chr12:6,537,201-6,538,700")
genes <- jbrowse_features(e("gapdh_gencode.gff3.gz"), loci)
signal <- jbrowse_wiggle(e("gapdh_phylop100way.bw"), loci)
reads <- jbrowse_alignments(e("gapdh_direct_rna.bam"), loci,
  display = list(encoding = list(color = list(field = "strand"))))
jbrowse_stack(genes, signal, reads)
)-"),
  amino_acids = list(
    caption = "GAPDH's coding transcripts open MGKVKVGVNG, the last codon split by the intron and lettered red on each side, the frame running on in the next exon.",
    title = "Amino acids over GAPDH's first coding exons, read in frame across the intron",
    locus = paste(aa_loci, collapse = " "),
    jbrowse = c(fasta, "--gffgz", gencode, grow,
      "--loc", paste(aa_loci, collapse = " ")),
    r = r"-(
loci <- c("chr12:6,534,811-6,534,880", "chr12:6,536,476-6,536,545")
genes <- jbrowse_features(e("gapdh_gencode.gff3.gz"), loci,
  assembly = e("gapdh_hg38.fa.gz"), width = 1200)
jbrowse_stack(genes)
)-"),
  variants_by_population = list(
    caption = "20 people from each superpopulation genotyped at 249 variants, a column per record, the East Asians homozygous for the alt allele at many of the first 40.",
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
    caption = "Most East Asians gather near the top, the Africans and Europeans further down.",
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
  variants_banded = list(
    caption = "Each superpopulation clusters apart in its own band, and the records crowding a pixel fade by how many share it.",
    title = "100 people's genotypes across DCST2 in a band per superpopulation, each band clustered",
    locus = "1:155,000,001-155,050,000", height = 7,
    jbrowse = c("--hub", "hg19", "--tracks", population_tracks(), "--track",
      "people", "variantLayout=genomic", "facet=superpopulation",
      "runClustering=true", "--loc", "1:155,000,001-155,050,000"),
    r = r"-(
v <- jbrowse_variants(e("1000g_dcst2.vcf.gz"), "1:155,000,001-155,050,000",
  samples = e("1000g_dcst2_samples.tsv"),
  display = list(rowColor = list(field = "superpopulation"),
    facet = "superpopulation"))
m <- as.matrix(v)
v$tree <- lapply(split(rownames(m), v$rows$band), function(rows)
  hclust(dist(m[rows, ]), "average"))
autoplot(v)
)-"),
  collapsed_introns = list(
    caption = "An arc over the coverage joins each pair of exons, and the direct RNA reads run across.",
    title = "GAPDH's nine exons side by side, with the junctions between them as arcs",
    locus = "ENST00000229239.10's exons, each padded by 20 bp", height = 7,
    jbrowse = c(fasta, "--gffgz", gencode, "--bam", rna, "sashimi:auto",
      "encoding.color.field=strand", "legend.show=true", grow, "--loc",
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
  display = list(encoding = list(color = list(field = "strand"))))
jbrowse_stack(genes, reads)
)-"),
  inversion_synteny = list(
    caption = "Forward blocks in tomato flank the 4 Mb inversion, whose reverse ribbon crosses in blue.",
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
  display = list(encoding = list(color = list(field = "strand"))), swap = TRUE,
  labels = c(query = "paternal", target = "maternal")))
)-"),
  inversion_dotplot = list(
    caption = "One forward diagonal runs the length of chr8, broken only by the inversion's short reverse segment near 10 Mb.",
    title = "HG002's maternal chr8 against its paternal copy as a dotplot",
    locus = "chr8_MATERNAL chr8_PATERNAL", height = 6,
    jbrowse = c("dotplot", "--chromSizes", axis_sizes(hg002_dot, "x", "maternal"),
      "--chain", chain, "--chromSizes", axis_sizes(hg002_dot, "y", "paternal"),
      "--colorBy", "strand"),
    r = r"-(
autoplot(jbrowse_dotplot(e("hg002_chr8.chain"),
  display = list(encoding = list(color = list(field = "strand"))), swap = TRUE,
  labels = c(query = "paternal", target = "maternal")))
)-"),
  whole_genome = list(
    caption = "Short diagonal runs scatter over most pairs of peach and grape chromosomes.",
    title = "Peach against grape, every chromosome", locus = "whole genome",
    height = 5.5,
    jbrowse = c("dotplot", "--chromSizes", axis_sizes(peach_grape, "x", "peach"),
      "--paf", inst("peach_grape.paf.gz"), "--chromSizes",
      axis_sizes(peach_grape, "y", "grape")),
    r = r"-(
autoplot(jbrowse_dotplot(e("peach_grape.paf.gz"),
  labels = c(query = "peach", target = "grape")))
)-"),
  amylase_synteny = list(
    caption = "Forward blocks flank a reversed copy in the middle, with alignments of the other copies crossing over it.",
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
  display = list(encoding = list(color = list(field = "strand"))),
  labels = c(query = "GRCh38", target = "HG00133")))
)-"),
  amylase_zoom = list(
    caption = "The deletion shows as a purple wedge narrowing to a point on HG00133's row.",
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
  display = list(encoding = list(color = list(field = "strand"))),
  labels = c(query = "GRCh38", target = "HG00133")))
)-"),
  synteny_genes = list(
    caption = "Forward blocks in tomato, the 4 Mb inversion crossing in blue, and the genes above and below.",
    title = "HG002's maternal chr8 over its paternal copy, each with its genes",
    locus = sprintf("chr8_MATERNAL:%s chr8_PATERNAL:%s", inversion, inversion),
    jbrowse = synteny_genes_args(),
    r = r"-(
r <- "5,250,001-14,250,000"
s <- jbrowse_synteny(e("hg002_chr8.chain"), paste0("chr8_MATERNAL:", r),
  paste0("chr8_PATERNAL:", r),
  display = list(encoding = list(color = list(field = "strand"))), swap = TRUE,
  labels = c(query = "paternal", target = "maternal"))
gene_track <- function(file, contig) jbrowse_features(e(file),
  paste0(contig, ":", r),
  display = list(displayMode = "compact", showLabels = FALSE,
    filter = "geneLike"))
autoplot(s, above = gene_track("hg002_mat_genes.gff3.gz", "chr8_MATERNAL"),
  below = gene_track("hg002_pat_genes.gff3.gz", "chr8_PATERNAL"))
)-"),
  multiway = list(
    caption = "Each lane runs in its own genome's coordinates, every liftOver ribbon in hg38's orientation, with gaps where a genome lacks sequence the one above has.",
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
  display = list(ribbons = list(encoding = list(color = list(field = "strand"))))))
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

recipe_figure <- function(example) {
  env <- new.env()
  env$e <- function(f) system.file("extdata", f, package = "ggjbrowse")
  eval(parse(text = example$r), env)
}

check_r <- function(name, example) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  print(recipe_figure(example))
}

draw_r <- function(name, example) {
  figure <- recipe_figure(example)
  height <- tryCatch(jbrowse_height(figure), error = function(e)
    example$height %||% 4)
  ggsave(png_of(name, "r"), figure, width = width / 120, height = height,
    dpi = 120)
}

to_webp <- function(name, sides = c("jbrowse", "r")) {
  for (side in sides) {
    system2("magick", c(png_of(name, side), "-quality", 85,
      webp_of(name, side)))
    file.remove(png_of(name, side))
  }
}

image_pair <- function(name, e, prefix) sprintf(
  "<tr><td width=\"50%%\"><img src=\"%s%s\" alt=\"%s in JBrowse. %s\" width=\"100%%\"></td><td width=\"50%%\"><img src=\"%s%s\" alt=\"%s drawn by ggjbrowse. %s\" width=\"100%%\"></td></tr>",
  prefix, webp_of(name, "jbrowse"), e$title, e$caption, prefix,
  webp_of(name, "r"), e$title, e$caption)

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
  e$caption, "",
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
  "  JBrowse's own renderer, from the root of a checkout of this repo.",
  paste0("  ", jb2export_version),
  "- `scripts/gallery.R` draws both and writes this page, so edit an example",
  "  there.",
  "",
  unlist(Map(recipe, names(examples), examples))),
  file.path("docs", "gallery.md"))

wanted <- commandArgs(trailingOnly = TRUE)
if (identical(wanted, "--check")) {
  for (name in names(examples)) {
    message(name)
    check_r(name, examples[[name]])
  }
} else if (!identical(wanted, "--docs")) {
  r_only <- "--r" %in% wanted
  wanted <- setdiff(wanted, "--r")
  todo <- if (length(wanted)) examples[wanted] else examples
  for (name in names(todo)) {
    message(name)
    if (!r_only) draw_jbrowse(name, todo[[name]])
    draw_r(name, todo[[name]])
    to_webp(name, if (r_only) "r" else c("jbrowse", "r"))
  }
}
update_readme(examples)
write_recipes(examples)

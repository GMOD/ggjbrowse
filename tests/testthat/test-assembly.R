extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
aliases <- tempfile(fileext = ".txt")
writeLines(c("# ucsc\tensembl\trefseq", "chr12\t12\tNC_000012.12"), aliases)
hg38 <- jbrowse_assembly(extdata("gapdh_hg38.fa.gz"), aliases, "hg38")
locus <- "chr12:6,534,201-6,538,700"

test_that("an assembly reads each file under its own name and the tables under the asked one", {
  expect_error(jbrowse_alignments(extdata("gapdh_direct_rna.bam"),
    "12:6,534,201-6,538,700"), "12 is not a contig of this file",
    class = "ggjbrowse_missing_contig")
  reads <- jbrowse_alignments(extdata("gapdh_direct_rna.bam"), locus)
  by_alias <- jbrowse_alignments(extdata("gapdh_direct_rna.bam"),
    "12:6,534,201-6,538,700", assembly = hg38)
  expect_equal(by_alias$reads[c("start", "end", "row")],
    reads$reads[c("start", "end", "row")])
  expect_equal(unique(by_alias$reads$chrom), "12")
  genes <- jbrowse_features(extdata("gapdh_gencode.gff3.gz"),
    "NC_000012.12:6,534,201-6,538,700", assembly = hg38)
  expect_gt(nrow(genes$features), 0)
  expect_equal(unique(genes$features$chrom), "NC_000012.12")
  signal <- jbrowse_wiggle(extdata("gapdh_phylop100way.bw"),
    "12:6,534,201-6,538,700", assembly = hg38)
  expect_gt(nrow(signal$marks), 0)
  expect_no_warning(jbrowse_stack(
    jbrowse_features(extdata("gapdh_gencode.gff3.gz"),
      "12:6,534,201-6,538,700", assembly = hg38),
    signal, by_alias))
})

test_that("a CRAM decodes against the assembly's sequence, and the sequence band reads it", {
  cram <- jbrowse_alignments(extdata("gapdh_direct_rna.cram"),
    "NC_000012.12:6,534,201-6,534,500", assembly = hg38)
  bam <- jbrowse_alignments(extdata("gapdh_direct_rna.bam"),
    "chr12:6,534,201-6,534,500")
  expect_equal(cram$reads$start, bam$reads$start)
  expect_equal(cram$sequence, jbrowse_sequence(extdata("gapdh_hg38.fa.gz"),
    "chr12:6,534,201-6,534,500"))
  expect_equal(jbrowse_sequence(hg38, "12:6,534,301-6,534,320"),
    jbrowse_sequence(extdata("gapdh_hg38.fa.gz"), "chr12:6,534,301-6,534,320"))
  expect_s3_class(ggplot2::ggplot_build(autoplot(cram, "sequence")),
    "ggplot_built")
})

test_that("an assembly config as a JBrowse config writes it serves as one", {
  config <- list(name = "hg38", uri = extdata("gapdh_hg38.fa.gz"),
    refNameAliases = list(uri = aliases))
  expect_equal(jbrowse_contigs(config),
    data.frame(assembly = "hg38", chrom = "chr12", start = 0, end = 6570000))
  expect_equal(jbrowse_contigs(extdata("gapdh_hg38.fa.gz"))$chrom, "chr12")
  bam <- jbrowse_contigs(extdata("volvox-sorted.bam"))
  expect_equal(bam, data.frame(assembly = NA_character_, chrom = "ctgA",
    start = 0, end = 50001, has_data = TRUE))
  expect_equal(jbrowse_contigs(extdata("gapdh_gencode.gff3.gz"),
    assembly = config)$end, 6570000)
  expect_equal(nrow(jbrowse_features(extdata("volvox.sort.gff3.gz"),
    jbrowse_contigs(extdata("volvox.sort.gff3.gz"))[2, ])$regions), 1)
  expect_equal(nrow(jbrowse_alignments(extdata("gapdh_direct_rna.bam"),
    "12:6,534,201-6,538,700", assembly = config)$reads), 73)
})

test_that("an assembly's locations read as a track adapter's do", {
  local <- function(f) list(localPath = extdata(f),
    locationType = "LocalPathLocation")
  config <- list(name = "hg38", sequence = list(
    type = "ReferenceSequenceTrack", trackId = "hg38-seq",
    adapter = list(type = "BgzipFastaAdapter",
      fastaLocation = local("gapdh_hg38.fa.gz"),
      faiLocation = local("gapdh_hg38.fa.gz.fai"),
      gziLocation = local("gapdh_hg38.fa.gz.gzi"))))
  expect_equal(jbrowse_contigs(config)$end, 6570000)
  location <- assembly_config(config)$sequence$adapter$fastaLocation
  expect_match(location$uri, "^file://.*\\?v=")
})

test_that("an alias file rewritten in place reads anew", {
  rewritten <- tempfile(fileext = ".txt")
  writeLines("chr12\t12", rewritten)
  assembly <- jbrowse_assembly(extdata("gapdh_hg38.fa.gz"), rewritten)
  bam <- extdata("gapdh_direct_rna.bam")
  expect_gt(nrow(jbrowse_alignments(bam, "12:6,534,201-6,534,700",
    assembly = assembly)$reads), 0)
  writeLines("chr12\ttwelve", rewritten)
  Sys.setFileTime(rewritten, Sys.time() + 5)
  expect_no_warning(twelve <- jbrowse_alignments(bam,
    "twelve:6,534,201-6,534,700", assembly = assembly))
  expect_gt(nrow(twelve$reads), 0)
})

test_that("refNameAliases written as the alias file's path alone serves", {
  aliases_file <- tempfile(fileext = ".txt")
  writeLines("chr12\t12", aliases_file)
  config <- list(name = "hg38", uri = extdata("gapdh_hg38.fa.gz"),
    refNameAliases = aliases_file)
  expect_gt(nrow(jbrowse_alignments(extdata("gapdh_direct_rna.bam"),
    "12:6,534,201-6,534,700", assembly = config)$reads), 0)
})

test_that("an inline alias list of one alias stays a list", {
  config <- list(name = "hg38", uri = extdata("gapdh_hg38.fa.gz"),
    refNameAliases = list(adapter = list(type = "FromConfigAdapter",
      features = list(list(refName = "chr12", uniqueId = "chr12",
        aliases = "12")))))
  expect_equal(nrow(jbrowse_alignments(extdata("gapdh_direct_rna.bam"),
    "12:6,534,201-6,538,700", assembly = config)$reads), 73)
})

test_that("a calls VCF with no GT reads as features with its FILTER, INFO and FORMAT", {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.3",
    "##FILTER=<ID=Too_low_VAF,Description=\"Low VAF\">",
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">",
    "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">",
    "##FORMAT=<ID=DP,Number=1,Type=Float,Description=\"Tumor depth\">",
    "##FORMAT=<ID=AD,Number=1,Type=Float,Description=\"Tumor alt depth\">",
    paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
      "FORMAT", "TUMOR", sep = "\t"),
    paste("12", "6534300", "d_1", "C", "<DEL>", ".", "PASS",
      "SVTYPE=DEL;SVLEN=-400", "DP:AD", "40:12", sep = "\t"),
    paste("12", "6535000", "i_2", "G", "<INS>", ".", "Too_low_VAF",
      "SVTYPE=INS", "DP:AD", "35:2", sep = "\t")), vcf)
  expect_equal(nrow(jbrowse_variants(vcf, "12:6,534,201-6,538,700")$variants),
    0)
  calls <- jbrowse_features(vcf, locus, assembly = hg38)
  f <- calls$features
  expect_equal(f$name, c("d_1", "i_2"))
  expect_equal(f$filter, c("PASS", "Too_low_VAF"))
  expect_equal(f$alt, c("<DEL>", "<INS>"))
  expect_equal(f$SVTYPE, c("DEL", "INS"))
  expect_equal(f$SVLEN, c(-400L, NA))
  expect_equal(unique(f$chrom), "chr12")
  expect_equal(calls$samples$feature, 1:2)
  expect_equal(calls$samples$sample, c("TUMOR", "TUMOR"))
  expect_equal(calls$samples$AD, c(12, 2))
  expect_equal(calls$fields$format$key, c("DP", "AD"))
  expect_equal(calls$fields$info$column, c("SVTYPE", "SVLEN"))
})

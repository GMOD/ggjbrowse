extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
gff <- extdata("volvox.sort.gff3.gz")
bam <- extdata("volvox-sorted.bam")

test_that("a locus reads as JBrowse's location box writes it, 1-based", {
  by_locus <- jbrowse_features(gff, "ctgA:1,001-9,000")
  by_bounds <- jbrowse_features(gff, "ctgA", 1000, 9000)
  expect_equal(by_locus$features, by_bounds$features)
  expect_equal(jbrowse_features(gff, "ctgA:1001..9000")$features,
    by_bounds$features)
  expect_error(jbrowse_features(gff, "ctgA"), "ctgA:1,001-9,000")
})

test_that("several regions keep their genomic positions and draw end to end", {
  loci <- c("ctgA:1,001-4,000", "ctgA:6,001-9,000")
  both <- jbrowse_alignments(bam, loci)
  alone <- jbrowse_alignments(bam, loci[2])
  r <- both$regions
  second <- both$reads[both$reads$.region == 2, ]
  expect_equal(second$start, alone$reads$start)
  expect_equal(unique(second$chrom), "ctgA")
  m <- both$mismatches
  expect_true(all(m$position >= both$reads$start[m$read] &
    m$position < both$reads$end[m$read]))
  expect_true(all(m$.region == both$reads$.region[m$read]))
  drawn <- on_axis(both)$reads
  expect_equal(min(drawn$start[drawn$.region == 2]),
    max(min(alone$reads$start), 6000) - r$start[2] + r$offset[2])
  expect_equal(ggplot2::fortify(both)$start, drawn$start)
  expect_equal(c(both$region$start, both$region$end),
    c(min(r$cum_start), max(r$cum_end)))
})

test_that("a contig the file lacks warns only where no region is on one it has", {
  expect_warning(jbrowse_features(gff, "chr1:1-5,000"),
    "chr1 is not a contig of this file, whose contigs are ctgA, ctgB")
  expect_no_warning(jbrowse_features(gff, c("ctgA:1-5,000", "chr1:1-5,000")))
})

test_that("a CRAM read without its reference says to pass one", {
  expect_error(jbrowse_alignments(extdata("volvox-sorted.cram"), "ctgA:1-100"),
    "pass that FASTA or 2bit as `reference =`")
})

test_that("jbrowse_stack stacks each track's bands as the display shows them", {
  locus <- "ctgA:1,001-9,000"
  genes <- jbrowse_features(gff, locus)
  reads <- jbrowse_alignments(bam, locus)
  expect_length(jbrowse_stack(genes, reads), 3)
  hidden <- jbrowse_alignments(bam, locus,
    display = list(showCoverage = FALSE))
  expect_length(jbrowse_stack(list(genes, hidden)), 2)
})

test_that("an axis of several regions is labelled in each region's own positions", {
  w <- jbrowse_wiggle(extdata("volvox_microarray_posneg.bw"),
    c("ctgA:1,001-4,000", "ctgA:6,001-9,000"))
  b <- ggplot2::ggplot_build(ggplot2::autoplot(w))
  expect_equal(b$layout$panel_params[[1]]$x$get_labels(),
    c("1,000", "2,000", "3,000", "7,000", "8,000", "9,000"))
})

test_that("jbrowse_stack refuses tracks over different regions", {
  genes <- jbrowse_features(gff, "ctgB:1-5,000")
  reads <- jbrowse_alignments(bam, "ctgA:1,001-9,000")
  expect_error(jbrowse_stack(genes, reads),
    "track 2 reads ctgA:1,001-9,000 and track 1 ctgB:1-5,000")
})

test_that("jbrowse_stack draws a ggplot of one's own on the tracks' axis", {
  genes <- jbrowse_features(gff, "ctgA:1,001-9,000")
  own <- ggplot2::ggplot(data.frame(x = c(2000, 5000), y = 1:2),
    ggplot2::aes(x, y)) + ggplot2::geom_point()
  fig <- jbrowse_stack(genes, own)
  expect_length(fig, 2)
  b <- ggplot2::ggplot_build(fig[[2]])
  expect_equal(b$layout$panel_params[[1]]$x.range, c(1000, 9000))
})

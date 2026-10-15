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

test_that("loci print with commas however their starts differ in size", {
  expect_equal(loci_of_layout(data.frame(chrom = "1", start = c(0, 155e6),
    end = c(1000, 155050000))), c("1:1-1,000", "1:155,000,001-155,050,000"))
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
    "pass that FASTA or 2bit as `assembly =`")
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

test_that("print names the loci read, however many regions", {
  loci <- c("ctgA:1,001-4,000", "ctgA:6,001-9,000")
  expect_output(print(jbrowse_alignments(bam, loci)),
    "ctgA:1,001-4,000, ctgA:6,001-9,000")
  expect_output(print(jbrowse_features(gff, loci)),
    "ctgA:1,001-4,000, ctgA:6,001-9,000")
  expect_output(print(jbrowse_alignments(bam, loci[1])), "ctgA:1,001-4,000")
})

test_that("each band of a stack marks the seam between regions with a bar", {
  loci <- c("ctgA:1,001-4,000", "ctgA:6,001-9,000")
  genes <- jbrowse_features(gff, loci)
  reads <- jbrowse_alignments(bam, loci)
  r <- reads$regions
  for (p in jbrowse_stack(genes, reads)$patches$plots) {
    b <- ggplot2::ggplot_build(p)
    bars <- Filter(function(d) !is.null(d$xintercept), b$data)
    expect_length(bars, 1)
    expect_true(bars[[1]]$xintercept > r$cum_end[1] &&
      bars[[1]]$xintercept < r$offset[2])
  }
  one <- ggplot2::ggplot_build(ggplot2::autoplot(jbrowse_features(gff,
    loci[1])))
  expect_false(any(vapply(one$data, function(d) !is.null(d$xintercept), NA)))
})

test_that("a plain ggplot of several regions draws them, and one's own layers, end to end", {
  loci <- c("ctgA:1,001-4,000", "ctgA:6,001-9,000")
  a <- jbrowse_alignments(bam, loci)
  p <- ggplot2::ggplot(a) + geom_jbrowse_coverage() +
    ggplot2::geom_point(ggplot2::aes(position, 1), data = a$mismatches)
  b <- ggplot2::ggplot_build(p)
  expect_equal(b$layout$panel_params[[1]]$x$get_labels(),
    c("1,000", "2,000", "3,000", "7,000", "8,000", "9,000"))
  r <- a$regions
  expect_lte(max(b$data[[2]]$x), max(r$cum_end))
  m <- a$mismatches
  second <- m$.region == 2 & m$position >= r$start[2] & m$position < r$end[2]
  expect_equal(sort(b$data[[2]]$x[b$data[[2]]$x > r$cum_end[1]]),
    sort(a$mismatches$position[second] - r$start[2] + r$offset[2]))
  own <- ggplot2::ggplot(a) + geom_jbrowse_coverage() +
    ggplot2::scale_x_continuous(name = "mine")
  expect_equal(ggplot2::ggplot_build(own)$plot$scales$get_scales("x")$name,
    "mine")
  expect_s3_class(scale_x_jbrowse(a)[[1]], "ScaleContinuousPosition")
  named <- scale_x_jbrowse(a, name = "Position")[[1]]
  expect_equal(named$name, "Position")
  expect_equal(named$breaks, scale_x_jbrowse(a)[[1]]$breaks)
})

test_that("jbrowse_stack takes a contig's other spelling and warns at another contig", {
  a <- jbrowse_features(gff, "ctgA:1,001-3,000")
  b <- jbrowse_features(gff, "ctgB:1,001-3,000")
  expect_warning(jbrowse_stack(a, b), "track 2 reads ctgB where track 1 reads ctgA")
  prefixed <- a
  prefixed$region$chrom <- "chrctgA"
  expect_no_warning(jbrowse_stack(a, prefixed))
  loci <- c("ctgA:1,001-2,000", "ctgA:3,001-4,000")
  two <- jbrowse_features(gff, loci)
  mixed <- jbrowse_features(gff, c("ctgA:1,001-2,000", "ctgB:3,001-4,000"))
  expect_error(jbrowse_stack(two, mixed), "same contigs")
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

test_that("a result's table handed to a geom or stacked in a plot of one's own draws on the axis", {
  loci <- c("ctgA:1,001-4,000", "ctgA:6,001-9,000")
  a <- jbrowse_alignments(bam, loci)
  axis <- c(a$region$start, a$region$end)
  minus <- subset(a$reads, strand == -1)
  p <- ggplot2::autoplot(a, "pileup") + geom_jbrowse_reads(data = minus)
  drawn <- ggplot2::layer_data(p, length(p$layers))
  expect_equal(nrow(drawn), nrow(minus))
  expect_true(all(drawn$xmin >= axis[1] & drawn$xmax <= axis[2]))

  own <- ggplot2::ggplot(a$coverage, ggplot2::aes(start, depth)) +
    ggplot2::geom_step() + ggplot2::coord_cartesian(ylim = c(0, 5))
  fig <- jbrowse_stack(a, own)
  band <- fig[[length(fig)]]
  built <- ggplot2::ggplot_build(band)
  expect_true(all(built$data[[1]]$x >= axis[1] & built$data[[1]]$x <= axis[2]))
  expect_equal(built$layout$panel_params[[1]]$y.range, c(0, 5))
  expect_null(own$data$.on_axis)
})

test_that("an alignments track refuses tracks above it over other regions", {
  expect_error(ggplot2::autoplot(jbrowse_alignments(bam, "ctgA:1,001-2,000"),
    above = jbrowse_features(gff, "ctgA:20,001-30,000")), "same regions")
})

test_that("a split read across two regions arcs over the gap, kept once", {
  r <- jbrowse_alignments(
    system.file("extdata", "volvox-long-reads-sv.bam", package = "ggjbrowse"),
    c("ctgA:1,800-3,000", "ctgA:33,900-35,100"),
    display = list(readConnections = "arc"))
  drawn <- on_axis(r)$arcs
  axis <- r$regions
  feet <- cbind(pmin(drawn$start, drawn$end), pmax(drawn$start, drawn$end))
  expect_true(all(feet[, 1] < axis$cum_end[1] & feet[, 2] > axis$cum_start[2]))
  expect_equal(nrow(drawn), nrow(unique(r$arcs[c("start", "end")])))
  expect_lt(nrow(drawn), nrow(r$arcs))
})

test_that("a line breaks at the seam between regions", {
  w <- jbrowse_wiggle(extdata("volvox_microarray_posneg.bw"),
    c("ctgA:1-10000", "ctgA:20001-30000"), display = list(mark = "line"))
  regions_per_run <- tapply(w$marks$.region, w$marks$run,
    function(r) length(unique(r)))
  expect_true(all(regions_per_run == 1))
})

test_that("a splice junction across two regions arcs over the seam, kept once", {
  r <- jbrowse_alignments(extdata("gapdh_direct_rna.bam"),
    c("chr12:6,534,451-6,534,700", "chr12:6,534,751-6,534,900"))
  drawn <- on_axis(r)$junctions
  axis <- r$regions
  across <- drawn$start < axis$cum_end[1] & drawn$end > axis$cum_start[2]
  expect_true(any(across))
  expect_false(any(duplicated(drawn[c("start", "end", "strand", "count")])))
  expect_true(all(drawn$start >= 0 & drawn$end <= axis$cum_end[2]))
})

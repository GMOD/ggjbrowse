extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
gff <- extdata("volvox.sort.gff3.gz")
bam <- extdata("volvox-sorted.bam")

test_that("a locus reads as JBrowse's location box writes it, 1-based", {
  by_locus <- jbrowse_features(gff, "ctgA:1,001-9,000")
  by_frame <- jbrowse_features(gff,
    data.frame(chrom = "ctgA", start = 1000, end = 9000))
  expect_equal(by_locus$features, by_frame$features)
  expect_equal(jbrowse_features(gff, "ctgA:1001..9000")$features,
    by_frame$features)
  expect_equal(by_locus$regions[c("chrom", "start", "end")],
    data.frame(chrom = "ctgA", start = 1000, end = 9000))
})

test_that("a contig named alone reads whole, at the length the file states", {
  whole <- jbrowse_alignments(bam, "ctgA")
  expect_equal(whole$regions$end, 50001)
  expect_equal(nrow(whole$reads), nrow(jbrowse_alignments(bam,
    "ctgA:1-50,001")$reads))
  expect_error(jbrowse_alignments(bam, "ctgZ"),
    class = "ggjbrowse_bad_region")
  expect_error(jbrowse_alignments(bam, "ctgZ"), "whose contigs are ctgA")
})

test_that("a contig whose length no file states reads to its records' end", {
  contigs <- jbrowse_contigs(gff)
  expect_true(all(is.na(contigs$end)))
  whole <- jbrowse_features(gff, "ctgA")
  expect_equal(whole$regions$end, max(whole$features$end))
  file <- jbrowse_features(gff)
  expect_equal(file$regions$end[file$regions$chrom == "ctgA"],
    whole$regions$end)
  expect_lt(max(file$regions$cum_end), 1e6)
  framed <- jbrowse_features(gff, contigs[1, ])
  expect_equal(framed$regions$end, whole$regions$end)
})

test_that("a reversed or 0-based locus is refused, an empty one read", {
  for (bad in c("ctgA:0-100", "ctgA:2,000-1,000", "ctgA:500-498")) {
    expect_error(regions_of(bad), class = "ggjbrowse_bad_region")
  }
  expect_error(regions_of("ctgA:2,000-1,000"), "ctgA:2,000-1,000 is reversed")
  expect_error(regions_of(data.frame(chrom = character(0), start = numeric(0),
    end = numeric(0))), "no rows", class = "ggjbrowse_bad_region")
  expect_equal(regions_of("ctgA:500-500")$end, 500)
  expect_equal(regions_of(data.frame(chrom = "ctgA", start = 5, end = 5))$end,
    5)
})

test_that("a frame's NA end is the contig's end where the file states it", {
  framed <- jbrowse_alignments(bam, data.frame(chrom = "ctgA", start = 0,
    end = NA))
  expect_equal(framed$regions$end, 50001)
})

test_that("a GRanges reads as its 1-based ranges, a 1-based frame refused", {
  skip_if_not_installed("GenomicRanges")
  gr <- GenomicRanges::GRanges("ctgA",
    IRanges::IRanges(c(1001, 6001), c(4000, 9000)))
  by_ranges <- jbrowse_features(gff, gr)
  expect_equal(by_ranges$regions$start, c(1000, 6000))
  expect_equal(by_ranges$features,
    jbrowse_features(gff, c("ctgA:1,001-4,000", "ctgA:6,001-9,000"))$features)
  expect_error(jbrowse_features(gff, as.data.frame(gr)),
    class = "ggjbrowse_bad_region")
  one_based <- data.frame(chrom = "ctgA", start = 1001, end = 4000,
    width = 3000)
  expect_error(jbrowse_features(gff, one_based), "0-based and half-open")
  expect_no_error(jbrowse_features(gff, transform(one_based, start = 1000)))
  skip_if_not_installed("tibble")
  expect_no_warning(jbrowse_features(gff,
    tibble::tibble(chrom = "ctgA", start = 1000, end = 4000)))
})

test_that("loci print with commas however their starts differ in size", {
  expect_equal(loci_of_layout(data.frame(chrom = "1", start = c(0, 155e6, 3e9),
    end = c(1000, 155050000, 3000000100))), c("1:1-1,000",
    "1:155,000,001-155,050,000", "1:3,000,000,001-3,000,000,100"))
})

test_that("several regions keep their genomic positions and draw end to end", {
  loci <- c("ctgA:1,001-4,000", "ctgA:6,001-9,000")
  both <- jbrowse_alignments(bam, loci)
  alone <- jbrowse_alignments(bam, loci[2])
  r <- both$regions
  second <- both$reads[both$reads$region == 2, ]
  expect_equal(second$start, alone$reads$start)
  expect_equal(unique(second$chrom), "ctgA")
  m <- both$mismatches
  read <- of(m, both$reads, "read")
  expect_false(anyNA(read$read))
  expect_true(all(m$position >= read$start & m$position < read$end))
  drawn <- on_axis(both)$reads
  expect_equal(min(drawn$start[drawn$region == 2]),
    max(min(alone$reads$start), 6000) - r$start[2] + r$offset[2])
  expect_equal(ggplot2::fortify(both)$start, drawn$start)
  expect_equal(xlim_of(both),
    c(min(r$cum_start), max(r$cum_end)))
})

test_that("a contig the file lacks errors where no region is on one it has", {
  expect_error(jbrowse_features(gff, "chr1:1-5,000"),
    "chr1 is not a contig of this file, whose contigs are ctgA, ctgB",
    class = "ggjbrowse_missing_contig")
  expect_warning(partial <- jbrowse_features(gff,
    c("ctgA:1-5,000", "chr1:1-5,000")), "chr1 is not a contig")
  expect_equal(partial$regions$missing[1], "")
  expect_match(partial$regions$missing[2], "chr1 is not a contig")
  expect_match(partial$warnings, "chr1 is not a contig")
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
    display = list(coverage = list(show = FALSE)))
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
  second <- m$region == 2 & m$position >= r$start[2] & m$position < r$end[2]
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
  prefixed$regions$chrom <- "chrctgA"
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
  axis <- xlim_of(a)
  minus <- subset(a$mismatches, base == "A")
  p <- ggplot2::autoplot(a, "pileup") + geom_jbrowse_mismatches(data = minus)
  drawn <- ggplot2::layer_data(p, length(p$layers))
  expect_lte(nrow(drawn), nrow(minus))
  expect_gt(nrow(drawn), 0)
  expect_true(all(drawn$xmin >= axis[1] & drawn$xmax <= axis[2]))

  own <- ggplot2::ggplot(a$coverage, ggplot2::aes(start, depth)) +
    ggplot2::geom_step() + ggplot2::coord_cartesian(ylim = c(0, 5))
  fig <- jbrowse_stack(a, own)
  band <- fig[[length(fig)]]
  built <- ggplot2::ggplot_build(band)
  expect_true(all(built$data[[1]]$x >= axis[1] & built$data[[1]]$x <= axis[2]))
  expect_equal(built$layout$panel_params[[1]]$y.range, c(0, 5))
  expect_null(attr(own$data, "on_axis"))
})

test_that("an alignments track refuses tracks stacked with it over other regions", {
  expect_error(jbrowse_stack(jbrowse_features(gff, "ctgA:20,001-30,000"),
    jbrowse_alignments(bam, "ctgA:1,001-2,000")), "same regions")
})

test_that("a frame of one's own lands in the region its chrom and positions fall in", {
  loci <- c("ctgA:1,001-4,000", "ctgA:6,001-9,000")
  genes <- jbrowse_features(gff, loci)
  r <- genes$regions
  peaks <- data.frame(chrom = c("ctgA", "ctgA", "ctgA", "ctgB"),
    start = c(1500, 6500, 4500, 1500), end = c(1600, 6600, 4600, 1600),
    y = 1:4)
  p <- ggplot2::ggplot(genes) + ggplot2::geom_rect(ggplot2::aes(
    xmin = start, xmax = end, ymin = y, ymax = y + 1), data = peaks)
  drawn <- ggplot2::layer_data(p, 1)
  expect_equal(drawn$ymin, 1:2)
  expect_equal(drawn$xmin, c(1500, 6500 - r$start[2] + r$offset[2]))
  sites <- data.frame(chrom = "ctgA", position = c(2000, 7000), y = 1)
  own <- ggplot2::ggplot(sites, ggplot2::aes(position, y)) +
    ggplot2::geom_point()
  band <- jbrowse_stack(genes, own)[[2]]
  expect_equal(ggplot2::layer_data(band, 1)$x,
    c(2000, 7000 - r$start[2] + r$offset[2]))
  skip_if_not_installed("tibble")
  marks <- tibble::tibble(chrom = "ctgA", start = c(2000, 7000), y = 1)
  expect_no_warning(drawn <- ggplot2::layer_data(ggplot2::ggplot(genes) +
    ggplot2::geom_point(ggplot2::aes(start, y), data = marks), 1))
  expect_equal(drawn$x, c(2000, 7000 - r$start[2] + r$offset[2]))
  one <- jbrowse_features(gff, "ctgA:1,001-9,000")
  expect_equal(ggplot2::layer_data(ggplot2::ggplot(one) +
    ggplot2::geom_rect(ggplot2::aes(xmin = start, xmax = end, ymin = y,
      ymax = y + 1), data = peaks), 1)$ymin, 1:3)
})

test_that("a result's $ matches names exactly, so $region is not $regions", {
  g <- jbrowse_features(gff, "ctgA:1-1,000")
  expect_null(g$region)
  expect_s3_class(g$regions, "data.frame")
})

test_that("an arc between two contigs rises half as high as its feet lie apart", {
  axis <- region_layout(data.frame(chrom = c("chr3", "chr10"),
    start = c(1000, 5000), end = c(2000, 6000)))
  arcs <- data.frame(start = c(1500, 1200), end = c(5500, 1800),
    y = c(0, 300), startRegion = c(0L, 0L), endRegion = c(1L, 0L),
    region = c(1L, 1L))
  placed <- arcs_on_axis(arcs, axis)
  expect_equal(placed$y[1], (placed$end[1] - placed$start[1]) / 2)
  expect_equal(placed$y[2], 300)
})

test_that("a split read across two regions arcs over the gap, listed once", {
  r <- jbrowse_alignments(
    system.file("extdata", "volvox-long-reads-sv.bam", package = "ggjbrowse"),
    c("ctgA:1,800-3,000", "ctgA:33,900-35,100"),
    display = list(readConnections = list(show = TRUE)))
  drawn <- on_axis(r)$arcs
  axis <- r$regions
  feet <- cbind(pmin(drawn$start, drawn$end), pmax(drawn$start, drawn$end))
  expect_true(all(feet[, 1] < axis$cum_end[1] & feet[, 2] > axis$cum_start[2]))
  expect_equal(nrow(drawn), nrow(r$arcs))
  expect_equal(nrow(r$arcs), nrow(unique(r$arcs[c("start", "end")])))
})

test_that("a line breaks at the seam between regions", {
  w <- jbrowse_wiggle(extdata("volvox_microarray_posneg.bw"),
    c("ctgA:1-10000", "ctgA:20001-30000"), display = list(mark = "line"))
  regions_per_run <- tapply(w$marks$region, w$marks$run,
    function(r) length(unique(r)))
  expect_true(all(regions_per_run == 1))
})

test_that("a splice junction across two regions arcs over the seam, listed once", {
  r <- jbrowse_alignments(extdata("gapdh_direct_rna.bam"),
    c("chr12:6,534,451-6,534,700", "chr12:6,534,751-6,534,900"))
  drawn <- on_axis(r)$junctions
  axis <- r$regions
  across <- drawn$start < axis$cum_end[1] & drawn$end > axis$cum_start[2]
  expect_true(any(across))
  expect_false(any(duplicated(drawn[c("start", "end", "strand", "count")])))
  expect_true(all(drawn$start >= 0 & drawn$end <= axis$cum_end[2]))
})

test_that("a read crossing a seam keeps one row in both regions", {
  a <- jbrowse_alignments(bam, c("ctgA:1,001-3,000", "ctgA:3,001-5,000"))
  r <- a$reads
  shared <- intersect(r$read[r$region == 1], r$read[r$region == 2])
  expect_gt(length(shared), 10)
  first <- r[r$region == 1 & r$read %in% shared, ]
  second <- r[r$region == 2 & r$read %in% shared, ]
  expect_equal(second$row[match(first$read, second$read)], first$row)
  expect_equal(second$id[match(first$read, second$read)], first$id)
})

test_that("a facet value only the second region holds stacks its section", {
  a <- jbrowse_alignments(bam, c("ctgA:100-110", "ctgA:6000-6100"),
    display = list(facet = list(field = "strand")))
  expect_setequal(as.character(a$reads$section), a$sections$label)
  expect_equal(nrow(a$sections), 2)
  expect_equal(unique(as.character(a$reads$section[a$reads$region == 1])),
    "Forward strand")
  expect_length(alignment_bands(a)$plots, 4)
})

test_that("regions on several contigs the file lacks warn once", {
  expect_warning(r <- jbrowse_features(gff,
    c("ctgA:1-5,000", "chr1:1-100", "chr2:1-100")),
    "2 contigs are not the file's and read empty: chr1, chr2")
  expect_length(r$warnings, 1)
  expect_equal(r$regions$missing != "", c(FALSE, TRUE, TRUE))
})

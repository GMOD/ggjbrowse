extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
paired <- extdata("volvox-simple-inv-paired.bam")
read <- function(...) jbrowse_alignments(paired, "ctgA:18,001-32,000", ...)

test_that("a facet splits the reads into sections, each packed on its own", {
  plain <- read()
  expect_equal(plain$sections$label, "")
  s <- read(display = list(facet = list(field = "strand")))
  expect_equal(nrow(s$sections), 2)
  expect_equal(levels(s$reads$section), s$sections$label)
  expect_equal(nrow(s$reads), nrow(plain$reads))
  expect_equal(nlevels(droplevels(s$reads$section)), 2)
  for (r in split(s$reads, s$reads$section)) {
    expect_equal(length(unique(r$strand)), 1)
    expect_equal(min(r$row), 0)
  }
  expect_equal(s$mismatches$section, s$reads$section[s$mismatches$read])
  flipped <- read(display = list(facet = list(field = "strand",
    domain = rev(s$sections$key))))
  expect_equal(flipped$sections$label, rev(s$sections$label))
  p <- ggplot2::ggplot_build(ggplot2::autoplot(s, "pileup"))
  expect_equal(nrow(p$layout$layout), 2)
  heights <- jbrowse_stack(s)$patches$layout$heights
  expect_equal(heights[c(1, 3)], c(45, 45))
  rows <- vapply(split(s$reads$row, s$reads$section), max, 0) + 1
  expect_equal(heights[c(2, 4)], unname(rows) * 8)
})

test_that("read connections come back as the display's band draws them", {
  expect_equal(nrow(read()$arcs), 0)
  a <- read(display = list(readConnections = "arc"))
  expect_gt(nrow(a$arcs), 100)
  expect_setequal(unique(a$arcs$shape), "arc")
  expect_true(all(c("pairLL", "pairRR") %in% a$arcs$category))
  expect_setequal(unique(a$arcs$color), unique(a$arc_legend$color))
  expect_equal(attr(a$arc_legend, "title"), "Arc colors")
  abnormal <- read(display = list(readConnections = "arc",
    showOrdinaryPairs = FALSE))
  expect_false("normalInsert" %in% abnormal$arcs$category)
  domes <- ggplot2::ggplot_build(ggplot2::autoplot(a, "arcs"))$data[[2]]
  expect_equal(length(unique(domes$group)), nrow(a$arcs))

  cloud <- read(display = list(readConnections = "cloud",
    facet = list(field = "strand")))
  expect_true(all(cloud$arcs$shape %in% c("flat", "flatSplit")))
  expect_equal(cloud$arc_span_max, max(cloud$arcs$span))
  expect_equal(attr(cloud$arc_legend, "title"), "Read cloud colors")
  expect_setequal(as.character(unique(cloud$arcs$section)),
    cloud$sections$label)
  track <- ggplot2::autoplot(cloud)
  expect_s3_class(track, "patchwork")
  expect_equal(length(track$patches$plots) + 1,
    3 * nrow(cloud$sections))
})

test_that("an arc rises as high as it is wide, the band reaching its pixels at the figure's scale", {
  a <- read(display = list(readConnections = "arc"))
  span <- 32000 - 18000
  expect_equal(arc_top(a), a$settings$readConnectionsHeight / (1500 / span))
  expect_lt(median(a$arcs$y), arc_top(a))
  expect_gt(max(a$arcs$y), arc_top(a))
  reach <- function(width) {
    ggplot2::autoplot(a, "arcs", width = width)$coordinates$limits$y[2]
  }
  expect_equal(reach(3000), arc_top(a) / 2)
})

test_that("a chain puts a read's mates and split segments on one row, a hairline across it", {
  plain <- read()
  chained <- read(display = list(unit = "chain"))
  expect_equal(nrow(plain$chains), 0)
  expect_gt(nrow(chained$chains), 10)
  expect_equal(chained$reads$start, plain$reads$start)
  rows <- tapply(chained$reads$row, chained$reads$name, function(r) length(unique(r)))
  expect_true(all(rows == 1))
  expect_true(all(chained$chains$end > chained$chains$start))
  built <- ggplot2::ggplot_build(ggplot2::autoplot(chained, "pileup"))$data
  lines <- Filter(function(d) all(c("xend", "yend") %in% names(d)) &&
    nrow(d) == nrow(chained$chains), built)
  expect_length(lines, 1)
  expect_equal(-lines[[1]]$yend, chained$chains$row + 0.5)
})

test_that("a chain over an inversion frames the split segments' strands against each other", {
  chained <- jbrowse_alignments(extdata("hg00151_chr1_inversion.bam"),
    "chr1:197,786,900-197,789,700", display = list(unit = "chain",
      color = list(field = "pairOrientation")))
  expect_gt(length(unique(chained$reads$color)), 2)
  expect_true(any(grepl("same strand", chained$legend$label)))
  expect_true(any(grepl("inverted", chained$legend$label)))
})

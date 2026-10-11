extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
paired <- extdata("volvox-simple-inv-paired.bam")
read <- function(...) jbrowse_alignments(paired, "ctgA", 18000, 32000, ...)

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
  expect_equal(length(track$patches$plots) + 1, 3)
})

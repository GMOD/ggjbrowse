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
  expect_equal(on_axis(s)$mismatches$section,
    of(s$mismatches, s$reads, "read")$section)
  flipped <- read(display = list(facet = list(field = "strand",
    sort = rev(s$sections$key))))
  expect_equal(flipped$sections$label, rev(s$sections$label))
  p <- ggplot2::ggplot_build(ggplot2::autoplot(s, "pileup"))
  expect_equal(nrow(p$layout$layout), 2)
  heights <- jbrowse_stack(s)$patches$layout$heights
  expect_equal(heights[c(1, 3)], c(45, 45))
  rows <- vapply(split(s$reads$row, s$reads$section), max, 0) + 1
  expect_equal(heights[c(2, 4)], unname(rows) * 9)
})

test_that("read connections come back as the display's band draws them", {
  expect_equal(nrow(read()$arcs), 0)
  a <- read(display = list(readConnections = list(show = TRUE)))
  expect_gt(nrow(a$arcs), 100)
  expect_setequal(unique(a$arcs$shape), "arc")
  expect_true(all(c("LL", "RR") %in% a$arcs$level))
  expect_setequal(unique(a$arcs$color), unique(a$keys$arcs$color))
  expect_equal(attr(a$keys$arcs, "title"), "Arc colors")
  abnormal <- read(display = list(readConnections = list(show = TRUE,
    ordinaryPairs = list(show = FALSE))))
  expect_true("normal" %in% a$arcs$level)
  expect_false("normal" %in% abnormal$arcs$level)
  domes <- ggplot2::ggplot_build(ggplot2::autoplot(a, "arcs"))$data[[2]]
  expect_equal(length(unique(domes$group)), nrow(a$arcs))

  cloud <- read(display = list(readConnections = list(show = TRUE, mark = "rule"),
    facet = list(field = "strand")))
  expect_true(all(cloud$arcs$shape %in% c("flat", "flatSplit")))
  expect_equal(cloud$arc_span_max, max(cloud$arcs$span))
  expect_equal(attr(cloud$keys$arcs, "title"), "Read cloud colors")
  expect_setequal(as.character(unique(cloud$arcs$section)),
    cloud$sections$label)
  track <- ggplot2::autoplot(cloud)
  expect_s3_class(track, "patchwork")
  expect_equal(length(track$patches$plots) + 1,
    3 * nrow(cloud$sections))
})

test_that("an arc rises as high as it is wide, the band reaching its pixels at the figure's scale", {
  a <- read(display = list(readConnections = list(show = TRUE)))
  span <- 32000 - 18000
  s <- a$arc_scale
  expect_equal(s$reachPx, s$heightPx - s$marginPx)
  expect_equal(arc_top(a), s$reachPx / (1500 / span))
  expect_lt(median(a$arcs$y), arc_top(a))
  expect_gt(max(a$arcs$y), arc_top(a))
  band <- function(width) {
    ggplot2::autoplot(a, "arcs", width = width)$coordinates$limits$y[2]
  }
  expect_equal(band(3000), s$heightPx / (3000 / span))
  expect_equal(arc_top(a, 3000), arc_top(a) / 2)
})

test_that("a chain puts a read's mates and split segments on one row, joined by connectors", {
  plain <- read()
  chained <- read(display = list(unit = "chain"))
  expect_equal(nrow(plain$connectors), 0)
  expect_gt(nrow(chained$connectors), 10)
  expect_true(all(plain$reads$name %in% chained$reads$name))
  rows <- tapply(chained$reads$row, chained$reads$name, function(r) length(unique(r)))
  expect_true(all(rows == 1))
  c <- chained$connectors
  expect_equal(c$end_row, c$row)
  name_of <- function(read) chained$reads$name[match(read, chained$reads$read)]
  expect_equal(name_of(c$mate), name_of(c$read))
  expect_true(all(c$kind %in% c("straight", "curve", "loop")))
})

test_that("each connector draws the curve the engine places, end to end on its rows", {
  chained <- read(display = list(unit = "chain"))
  c <- chained$connectors
  dips <- which(c$kind != "straight" & c$dip > 0)
  expect_gt(length(dips), 0)
  paths <- connector_paths(on_axis(chained), 1200)
  by_path <- split(paths, paths$path)[seq_len(nrow(c))]
  expect_equal(vapply(by_path, function(p) p$x[1], 0), c$start,
    ignore_attr = TRUE)
  expect_equal(vapply(by_path, function(p) p$x[nrow(p)], 0), c$end,
    ignore_attr = TRUE)
  expect_equal(vapply(by_path, function(p) p$y[1], 0), c$row + 0.5,
    ignore_attr = TRUE)
  depth <- function(paths, i) {
    p <- paths[paths$path == i, ]
    max(p$y) - (max(c$row[i], c$end_row[i]) + 0.5)
  }
  wide <- vapply(dips, function(i) depth(paths, i), 0)
  expect_true(all(wide > 0))
  narrow <- connector_paths(on_axis(chained), 60)
  expect_true(all(vapply(dips, function(i) depth(narrow, i), 0) < wide))
  built <- ggplot2::ggplot_build(ggplot2::autoplot(chained, "pileup"))$data
  drawn <- Filter(function(d) "group" %in% names(d) &&
    nrow(d) == nrow(paths) && "linetype" %in% names(d), built)
  expect_length(drawn, 1)
})

test_that("a chain over an inversion frames the split segments' strands against each other", {
  chained <- jbrowse_alignments(extdata("hg00151_chr1_inversion.bam"),
    "chr1:197,786,900-197,789,700", display = list(unit = "chain", encoding = list(color = list(field = "insertSizeAndOrientation"))))
  expect_gt(length(unique(chained$reads$color)), 2)
  expect_true(any(grepl("same strand", chained$keys$reads$label)))
  expect_true(any(grepl("inverted", chained$keys$reads$label)))
  expect_true(all(c("unpairedSplitInverted", "unpairedSplitSameStrand") %in%
    chained$reads$level))
})

test_that("each level paints one color, which the read key holds", {
  for (field in c("strand", "mapq", "insertSizeAndOrientation")) {
    a <- read(display = list(encoding = list(color = list(field = field))))
    colors <- tapply(a$reads$color, a$reads$level, unique)
    expect_true(all(lengths(colors) == 1))
    expect_true(all(unlist(colors) %in% a$keys$reads$color))
  }
})

test_that("insertSizeAndOrientation colors by level, in the band normalInsertSize pins", {
  color <- function(...) list(color = list(field = "insertSizeAndOrientation",
    ...))
  wide <- read(display = list(encoding = color(),
    normalInsertSize = list(min = 1, max = 1e7)))
  expect_true("LL" %in% wide$reads$level)
  expect_false(any(c("short", "long") %in% wide$reads$level))
  narrow <- read(display = list(encoding = color(),
    normalInsertSize = list(min = 1, max = 2)))
  expect_true("long" %in% narrow$reads$level)
  teal <- read(display = list(encoding = color(scale = list(domain = "LL",
    range = "#008080"))))
  expect_equal(unique(teal$reads$color[teal$reads$level == "LL"]), "#008080")
  expect_false("#008080" %in% teal$reads$color[teal$reads$level != "LL"])
})

test_that("the arcs and the read cloud hang from the band's top under placement down, and rise under up", {
  y_transform <- function(x) {
    ggplot2::ggplot_build(ggplot2::autoplot(x, "arcs"))$layout$
      panel_scales_y[[1]]$get_transformation()$name
  }
  arcs <- read(display = list(readConnections = list(show = TRUE)))
  expect_equal(arcs$display$readConnections$placement, "down")
  expect_match(y_transform(arcs), "reverse")
  cloud <- read(display = list(readConnections = list(show = TRUE,
    mark = "rule")))
  expect_match(y_transform(cloud), "log-10.*reverse")
  up <- read(display = list(readConnections = list(show = TRUE,
    mark = "rule", placement = "up")))
  expect_match(y_transform(up), "log-10")
  expect_no_match(y_transform(up), "reverse")
  expect_equal(y_transform(read(display = list(readConnections =
    list(show = TRUE, placement = "up")))), "identity")
})

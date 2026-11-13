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
  expect_true(all(c$kind %in% c("straight", "curve", "bow", "loop")))
})

test_that("a chain's connector over a hidden segment bows over its row, where a flat curve runs along it", {
  a <- jbrowse_alignments(extdata("volvox-simple-inv-bow.bam"),
    c("ctgA:15,000-19,950", "ctgA:30,010-33,000"),
    display = list(unit = "chain"))
  c <- a$connectors
  bows <- which(c$kind == "bow")
  expect_length(bows, 2)
  expect_equal(c$dip[bows], c(0, 0))
  expect_equal(c$end_row[bows], c$row[bows])
  flat <- a
  flat$connectors$kind[bows] <- "curve"
  bowed <- connector_paths(on_axis(a), NULL)
  run <- connector_paths(on_axis(flat), NULL)
  for (i in bows) {
    expect_lt(min(bowed$y[bowed$path == i]), c$row[i] + 0.5)
    expect_equal(range(run$y[run$path == i]), rep(c$row[i] + 0.5, 2))
  }
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
  zoomed_out <- connector_paths(on_axis(chained), 60)
  narrow <- vapply(dips, function(i) depth(zoomed_out, i), 0)
  expect_true(all(narrow <= wide))
  expect_true(any(narrow < wide))
  built <- ggplot2::ggplot_build(ggplot2::autoplot(chained, "pileup"))$data
  drawn <- Filter(function(d) "group" %in% names(d) &&
    nrow(d) == nrow(paths) && "linetype" %in% names(d), built)
  expect_length(drawn, 1)
})

polylines <- function(g) if (inherits(g, "polyline")) list(g) else
  unlist(lapply(g$children, polylines), recursive = FALSE)

drawn_strokes <- function(plot) {
  lines <- polylines(ggplot2::layer_grob(plot)[[1]])
  unique(do.call(rbind, lapply(lines, function(l) data.frame(
    px = l$gp$lwd / 96 * PX_PER_INCH,
    alpha = grDevices::col2rgb(l$gp$col, alpha = TRUE)[4, ] / 255,
    lty = l$gp$lty))))
}

test_that("each connector strokes the engine's width in px at its opacity", {
  chained <- read(display = list(unit = "chain"))
  n <- nrow(chained$connectors)
  chained$connectors$width <- seq(1, 4, length.out = n)
  chained$connectors$opacity <- seq(0.2, 1, length.out = n)
  drawn <- drawn_strokes(ggplot2::ggplot(chained) + geom_jbrowse_connectors())
  drawn <- drawn[order(drawn$px), ]
  expect_equal(drawn$px, chained$connectors$width)
  expect_equal(drawn$alpha, chained$connectors$opacity, tolerance = 1 / 255)
})

test_that("every connector draws solid, over a hidden segment too", {
  a <- jbrowse_alignments(extdata("volvox-simple-inv-bow.bam"),
    c("ctgA:15,000-19,950", "ctgA:30,010-33,000"),
    display = list(unit = "chain"))
  expect_gt(sum(a$connectors$kind == "bow"), 0)
  drawn <- drawn_strokes(ggplot2::ggplot(a) + geom_jbrowse_connectors())
  expect_true(all(drawn$lty %in% c("1", "solid")))
})

test_that("every connector's curve stays within the pileup's rows, where its y axis ends", {
  chained <- read(display = list(unit = "chain"))
  rows <- max(chained$reads$row) + 1
  pitch <- chained$read_pitch$featureHeight + chained$read_pitch$featureSpacing
  for (width in list(NULL, 1200)) {
    expect_lte(max(connector_paths(on_axis(chained), width)$y), rows)
  }
  plot <- ggplot2::autoplot(chained, "pileup")
  expect_equal(plot$coordinates$limits$y, c(0, rows))
  expect_equal(alignment_px(on_axis(chained), alignment_options(chained))[["pileup"]],
    max(30, rows * pitch))
})

test_that("a faceted pileup draws each section's connectors in its own panel, its y axis on its rows", {
  a <- read(display = list(facet = list(field = "pairOrientation"),
    pileup = list(connectors = list(show = TRUE))))
  plot <- ggplot2::autoplot(a, "pileup")
  layer <- which(vapply(plot$layers,
    function(l) inherits(l$geom, "GeomPath"), TRUE))
  built <- ggplot2::ggplot_build(plot)
  panels <- built$layout$layout
  drawn <- built$data[[layer]]
  paths <- connector_paths(on_axis(a), NULL)
  expect_equal(as.vector(table(panels$section[drawn$PANEL])),
    as.vector(table(paths$section)))
  bare <- plot
  bare$layers <- plot$layers[-layer]
  ranges <- function(b) lapply(b$layout$panel_params, function(p) p$y.range)
  expect_equal(ranges(built), ranges(ggplot2::ggplot_build(bare)))
})

test_that("a section's curve dipping past its last row draws over the sections below it, clipped at its panel's sides", {
  a <- read(display = list(facet = list(field = "pairOrientation"),
    pileup = list(connectors = list(show = TRUE))))
  rr <- a$sections$label[2]
  paths <- connector_paths(on_axis(a), NULL)
  expect_gt(max(paths$y[paths$section == rr]),
    max(a$reads$row[a$reads$section == rr]) + 1)
  grDevices::pdf(NULL, width = 9, height = 5)
  on.exit(grDevices::dev.off())
  grid::grid.newpage()
  grid::grid.draw(ggplot2::ggplotGrob(ggplot2::autoplot(a, "pileup")))
  grid::grid.force()
  spill <- grid::getGrob(grid::grid.get("panel-1-2", grep = TRUE),
    "connectors", grep = TRUE)
  y <- as.numeric(spill$children[[1]]$y)
  expect_true(spill$vp$clip)
  expect_equal(spill$vp$xscale, c(0, 1))
  expect_lt(min(y), 0)
  expect_gte(min(y), spill$vp$yscale[1])
  grid::grid.newpage()
  rr_only <- section_of(on_axis(a), rr)
  grid::grid.draw(ggplot2::ggplotGrob(ggplot2::ggplot(rr_only) +
    geom_jbrowse_connectors() + ggplot2::scale_y_reverse() +
    ggplot2::coord_cartesian(ylim = c(0, 2))))
  grid::grid.force()
  spill <- grid::grid.get("connectors", grep = TRUE)
  in_window <- with(rr_only$connectors, sum(pmax(row, end_row) < 2))
  expect_gt(in_window, 0)
  expect_lt(in_window, nrow(rr_only$connectors))
  expect_length(unique(spill$children[[1]]$id), in_window)
})

test_that("a track draws every section's connector curves over every section's band, as the display paints them last", {
  a <- read(display = list(facet = list(field = "pairOrientation"),
    pileup = list(connectors = list(show = TRUE))))
  grDevices::pdf(NULL, width = 9, height = 7)
  on.exit(grDevices::dev.off())
  print(ggplot2::autoplot(a))
  grid::grid.force()
  drawn <- grid::grid.ls(print = FALSE)
  curves <- grep("^connectors", drawn$name)
  bands <- which(startsWith(drawn$name, "panel") & drawn$gDepth == 1)
  expect_length(curves, nrow(a$sections))
  expect_gt(min(curves), max(bands))
  dips <- vapply(drawn$name[curves],
    function(name) grid::grid.get(name)$vp$yscale[1], 0)
  expect_lt(min(dips), 0)
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

test_that("a pileup keys the connectors it draws, and the read-connection band keys them once", {
  colo <- extdata("colo829_tumor_chr3.bam")
  chain <- function(...) jbrowse_alignments(colo, "chr3:25,358,500-25,359,800",
    display = list(unit = "chain", facet = list(field = "split"), ...))
  hidden <- "Split read (segment not in view)"

  a <- chain()
  expect_equal(a$keys$arcs$label, hidden)
  expect_equal(attr(a$keys$arcs, "title"), "Connector colors")
  key <- drawn_guides(autoplot(a, "pileup"))[["Connector colors"]]
  expect_equal(key$label, hidden)
  expect_setequal(key$color, unique(a$connectors$color))
  expect_equal(sum(drawn_labels(autoplot(a)) == hidden), 1)

  plain <- jbrowse_alignments(colo, "chr3:25,358,500-25,359,800")
  expect_equal(nrow(plain$connectors), 0)
  expect_false("Connector colors" %in%
    names(drawn_guides(autoplot(plain, "pileup"))))
  expect_false("Connector colors" %in% drawn_labels(autoplot(plain)))

  banded <- chain(readConnections = list(show = TRUE))
  expect_gt(nrow(banded$arcs), 0)
  rows <- drawn_labels(autoplot(banded))
  expect_equal(sum(rows == hidden), 1)
  expect_equal(sum(rows == "Arc colors"), 1)
  expect_false("Connector colors" %in% rows)
  expect_equal(sum(rows == "Inter-chromosomal"), 1)
})

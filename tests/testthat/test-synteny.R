paf <- system.file("extdata", "volvox_inv_indels.paf", package = "ggjbrowse")
gff <- system.file("extdata", "volvox.sort.gff3.gz", package = "ggjbrowse")
query <- "ctgA:1-49,186"
target <- "ctgA:1-50,001"
s <- jbrowse_synteny(paf, query, target)

test_that("a PAF reads as the synteny band lays out and colors it", {
  expect_equal(nrow(s$alignments), 7)
  expect_equal(sum(s$alignments$strand == -1), 2)
  expect_equal(s$assemblies, c("query", "target"))
  expect_equal(s$display$opacity, 0.25)
  expect_true(all(c("match", "insertion") %in% s$ribbons$kind))
  match <- s$ribbons[s$ribbons$kind == "match", ]
  expect_equal(nrow(match), 7)
  reversed <- s$alignments$strand[match$alignment] == -1
  expect_equal(match$query_start > match$query_end, reversed)
  expect_true(all(match$target_start < match$target_end))
  expect_equal(unique(as.character(match$color)), "#ff0000")
  expect_equal(s$keys$ribbons$label, c("match", "insertion", "deletion"))
  expect_true("deletion" %in% s$ribbons$kind)
  expect_s3_class(s$ribbons$color, "jbrowse_color")
  expect_equal(nrow(s$offscreen), 0)
})

test_that("the view's settings color by strand, set the opacity and pick the CIGAR mode", {
  strand <- jbrowse_synteny(paf, query, target,
    display = list(color = list(field = "strand"), opacity = 0.6,
      cigarMode = "off"))
  expect_equal(unique(strand$ribbons$kind), "match")
  expect_equal(length(unique(strand$ribbons$color)), 2)
  expect_equal(strand$display$opacity, 0.6)
  expect_equal(strand$keys$ribbons$label[1:2], c("forward", "reverse"))
  grey <- jbrowse_synteny(paf, query, target, display = list(color = "grey"))
  expect_equal(unique(as.character(grey$alignments$color)), "#808080")
})

test_that("the band draws every ribbon between the rows, with both rulers", {
  p <- autoplot(s)
  built <- ggplot2::ggplot_build(p)
  polygons <- built$data[[1]]
  expect_equal(length(unique(polygons$group)), sum(!s$ribbons$thin))
  expect_true(all(polygons$y %in% c(0, 1)))
  expect_true(all(polygons$x >= 0 & polygons$x <= s$width))
  scale <- built$layout$panel_scales_x[[1]]
  expect_match(scale$name, "bp")
  expect_s3_class(scale$secondary.axis, "AxisSecondary")
  curved <- ggplot2::ggplot_build(autoplot(s, curves = TRUE))$data[[1]]
  expect_gt(nrow(curved), nrow(polygons))
})

test_that("tracks over each row stack around the band", {
  genes <- jbrowse_features(gff, target, display = list(showOnlyGenes = TRUE))
  fig <- autoplot(s, below = genes)
  expect_equal(length(fig$patches$plots) + 1, 2)
  expect_equal(fig$patches$layout$heights[1], s$display$height)
  expect_gt(jbrowse_height(fig), 1)
  expect_error(autoplot(s, above = genes), "top row.s regions")
})

test_that("a PAF reads as the dotplot lays out and colors it", {
  d <- jbrowse_dotplot(paf, query, target)
  expect_equal(nrow(d$alignments), 7)
  expect_gt(nrow(d$segments), 7)
  expect_equal(unique(as.character(d$segments$color)), "#000000")
  expect_equal(d$display$lineWidth, 2.5)
  expect_equal(d$display$opacity, 1)
  expect_equal(d$skipped$count, 0)
  built <- ggplot2::ggplot_build(autoplot(d))
  expect_equal(nrow(built$data[[3]]), nrow(d$segments))
  plain <- jbrowse_dotplot(paf, query, target,
    display = list(drawCigar = FALSE, color = list(field = "strand")))
  expect_equal(plain$segments$op, rep("match", 7))
  expect_equal(length(unique(plain$segments$color)), 2)
})

test_that("two regions on an axis lay end to end, each ribbon on its own", {
  halves <- c("ctgA:1-25,000", "ctgA:25,001-50,001")
  two <- jbrowse_synteny(paf, query, halves)
  expect_setequal(unique(two$ribbons$target_region), 2:3)
  axis <- on_axis(two)$ribbons
  expect_true(all(axis$target_px_start >= 0 & axis$target_px_start <= two$width))
  d <- jbrowse_dotplot(paf, query, halves)
  moved <- on_axis(d)$segments
  second <- moved$y_region == 3
  expect_true(any(second))
  expect_true(all(moved$y_start[second] >= d$regions$offset[3]))
  expect_equal(length(layout_seams(on_axis_of(d$regions, "y"))), 1)
})

test_that("one regions frame lays out every axis, and the tables index its rows", {
  halves <- c("ctgA:1-25,000", "ctgA:25,001-50,001")
  two <- jbrowse_synteny(paf, query, halves)
  expect_named(two$regions, c("axis", "assembly", "chrom", "start", "end",
    "offset", "cum_start", "cum_end", "missing"))
  expect_equal(two$regions$axis, c("query", "target", "target"))
  expect_equal(two$regions$assembly, c("query", "target", "target"))
  expect_true(all(two$regions$axis[two$ribbons$query_region] == "query"))
  expect_true(all(two$regions$axis[two$ribbons$target_region] == "target"))
  inside <- function(region, start, end) {
    r <- two$regions[region, ]
    pmin(start, end) >= r$start & pmax(start, end) <= r$end
  }
  expect_true(all(inside(two$ribbons$target_region, two$ribbons$target_start,
    two$ribbons$target_end)))
  expect_false(any(c("query", "target") %in% names(two)))
  d <- jbrowse_dotplot(paf, query, halves)
  expect_equal(d$regions$axis, c("x", "y", "y"))
  expect_false(any(c("x", "y") %in% names(d)))
  expect_true(all(d$regions$axis[d$segments$y_region] == "y"))
})

test_that("drawing places tables on the axis once, marked by an attribute", {
  placed <- on_axis(s)
  expect_true(isTRUE(attr(placed$ribbons, "on_axis")))
  expect_null(placed$ribbons$.on_axis)
  expect_identical(on_axis(placed)$ribbons, placed$ribbons)
  d <- jbrowse_dotplot(paf, query, c("ctgA:1-25,000", "ctgA:25,001-50,001"))
  once <- on_axis(d)
  expect_null(once$segments$.on_axis)
  expect_equal(on_axis(once)$segments, once$segments)
})

test_that("the display slot holds the settings each view drew with", {
  expect_equal(s$display, list(minAlignmentLength = 0, cigarMode = "full",
    drawLocationMarkers = FALSE, fadeThinAlignmentsMode = "auto",
    opacity = 0.25, height = 100))
  strand <- jbrowse_synteny(paf, query, target,
    display = list(color = list(field = "strand"), cigarMode = "off"))
  expect_equal(strand$display$color, list(field = "strand"))
  expect_equal(strand$display$cigarMode, "off")
  expect_null(strand[["field"]])
  expect_null(strand$alpha)
  expect_null(strand$height)
  ref <- jbrowse_synteny(paf, query, target,
    display = list(color = list(field = "reference")))
  expect_equal(ref$display$color$field, "query")
  faded <- jbrowse_synteny(paf, query, target,
    display = list(opacity = list(field = "identity")))
  expect_equal(faded$display$opacity, list(field = "identity"))
  d <- jbrowse_dotplot(paf, query, target, display = list(color = "grey"))
  expect_equal(d$height, 1500)
  expect_equal(jbrowse_dotplot(paf, query, target, width = 900)$height, 900)
  expect_equal(d$display, list(minAlignmentLength = 0, minIdentity = 0,
    drawCigar = TRUE, color = "grey", opacity = 1, lineWidth = 2.5))
  expect_null(d$line_width)
  expect_error(jbrowse_synteny(paf, query, target, view = list()),
    "unused argument")
})

test_that("a comparative result holds the slots every result holds", {
  results <- list(s, jbrowse_dotplot(paf, query, target),
    jbrowse_multiway(system.file("extdata", "volvox_all_vs_all.paf",
      package = "ggjbrowse"), "ctgA:1-50,001", "volvox"))
  for (x in results) {
    expect_true(all(c("regions", "width", "display", "keys", "ramps",
      "y_scale", "fields", "warnings") %in% names(x)))
    expect_equal(x$width, 1500)
    expect_null(x$y_scale)
    expect_equal(x$fields, list())
    expect_equal(x$warnings, character(0))
    expect_false(any(c("legend", "field", "alpha", "line_width", "curves",
      "query", "target", "x", "y", "anchor") %in% names(x)))
  }
})

test_that("warnings the engine raises land in warnings", {
  real <- call_engine
  local_mocked_bindings(call_engine = function(table, args) {
    if (table == "synteny") warning("the engine says so", call. = FALSE)
    real(table, args)
  })
  expect_warning(w <- jbrowse_synteny(paf, query, target), "engine says so")
  expect_equal(w$warnings, "the engine says so")
})

test_that("an unknown extension and a missing contig say so", {
  expect_error(jbrowse_synteny("aln.tsv", query, target), "expected .paf")
  expect_error(jbrowse_synteny(paf, "ctgY:1-1,000", "ctgZ:1-1,000"),
    "ctgY is not a contig of the file's query side",
    class = "ggjbrowse_missing_contig")
  expect_error(jbrowse_dotplot(paf, "ctgY:1-1,000", "ctgZ:1-1,000"),
    class = "ggjbrowse_missing_contig")
  expect_warning(empty <- jbrowse_synteny(paf, query, "ctgZ:1-1,000"),
    "ctgZ is not a contig of the file's target side")
  expect_equal(nrow(empty$ribbons), 0)
  expect_warning(beside <- jbrowse_synteny(paf, query,
    c(target, "ctgZ:1-1,000")), "ctgZ is not a contig")
  expect_equal(beside$regions$missing[1:2], c("", ""))
  expect_match(beside$regions$missing[3], "ctgZ is not a contig")
  expect_match(beside$warnings, "ctgZ is not a contig")
  expect_equal(nrow(beside$alignments), 7)
  expect_error(jbrowse_synteny(paf, "ctgZ", target),
    class = "ggjbrowse_bad_region")
})

test_that("each axis takes every region form a linear reader takes", {
  contigs <- jbrowse_contigs(paf)
  whole <- jbrowse_synteny(paf, contigs, contigs)
  expect_equal(whole$regions[c("assembly", "chrom", "start", "end")],
    contigs[c("assembly", "chrom", "start", "end")], ignore_attr = TRUE)
  expect_equal(jbrowse_synteny(paf, whole$regions, whole$regions)$regions,
    whole$regions)
  expect_equal(jbrowse_synteny(paf)$regions, whole$regions)
  expect_equal(jbrowse_synteny(paf, "ctgA", factor("ctgA"))$regions,
    whole$regions)
  frame <- data.frame(chrom = "ctgA", start = 0, end = 49186)
  expect_equal(jbrowse_synteny(paf, frame, tibble::tibble(chrom = "ctgA",
    start = 0, end = 50001))$regions, whole$regions)
  expect_error(jbrowse_synteny(paf, data.frame(assembly = "grape",
    chrom = "ctgA", start = 0, end = 10)), "none of them this axis's query",
    class = "ggjbrowse_bad_region")
  expect_error(jbrowse_synteny(paf, data.frame(chrom = "ctgA", start = 1,
    end = 10, width = 10)), class = "ggjbrowse_bad_region")
  skip_if_not_installed("GenomicRanges")
  g <- GenomicRanges::GRanges("ctgA", IRanges::IRanges(1, 49186))
  expect_equal(jbrowse_synteny(paf, g, target)$regions$end[1], 49186)
})

test_that("a setting the view does not take raises ggjbrowse_bad_display", {
  expect_error(jbrowse_synteny(paf, query, target,
    display = list(cigarMode = "some")), "display\\$cigarMode",
    class = "ggjbrowse_bad_display")
  expect_error(jbrowse_synteny(paf, query, target,
    display = list(fadeThinAlignmentsMode = TRUE)),
    class = "ggjbrowse_bad_display")
})

test_that("identity keys as a bar from 0% to 100%, and mapq as chips", {
  id <- jbrowse_synteny(paf, query, target,
    display = list(color = list(field = "identity")))
  expect_length(id$ramps, 1)
  expect_equal(attr(id$ramps[[1]], "title"), "Identity")
  expect_equal(attr(id$ramps[[1]], "labels"), c("0%", "100%"))
  g <- ggplot2::ggplotGrob(autoplot(id))
  expect_false(inherits(g$grobs[[which(g$layout$name == "guide-box-right")]],
    "zeroGrob"))
  d <- jbrowse_dotplot(paf, query, target,
    display = list(color = list(field = "identity")))
  expect_equal(d$ramps, id$ramps)
  mapq <- jbrowse_synteny(paf, query, target,
    display = list(color = list(field = "mapq")))
  expect_length(mapq$ramps, 0)
  expect_gt(nrow(mapq$keys$ribbons), 1)
})

swap <- system.file("extdata", "volvox_contig_swap.paf", package = "ggjbrowse")

test_that("the off-screen strip draws as the band paints it and keys its contigs", {
  skip_if(swap == "")
  expect_warning(o <- jbrowse_synteny(swap, "ctgA:1-50,001", "ctgC:1-1,000"),
    "ctgC is not a contig")
  expect_equal(o$offscreen$mate_chrom, "ctgB")
  expect_true(o$offscreen$drawn)
  expect_equal(o$offscreen_contigs$label, "ctgB · 6.08Kbp in 1 alignments")
  expect_equal(o$offscreen_strip$marks$height, 6)
  expect_equal(o$offscreen_strip$labels$label, "ctgB")
  t <- suppressWarnings(jbrowse_synteny(swap, "ctgA:1-50,001", "ctgC:1-1,000",
    display = list(color = list(field = "target"))))
  expect_false(identical(t$offscreen$color, o$offscreen$color))
  built <- ggplot2::ggplot_build(autoplot(t))
  rects <- built$data[[3]]
  expect_equal(rects$ymax, 1)
  expect_equal(rects$ymin, 1 - 6 / t$display$height)
  expect_equal(as.character(rects$fill), as.character(t$offscreen$color))
  expect_false(is.null(offscreen_guides(t)))
  expect_null(offscreen_guides(s))
})

test_that("a file's contigs list with their lengths, and an axis with no regions shows them all", {
  contigs <- jbrowse_contigs(paf)
  expect_named(contigs, c("assembly", "chrom", "start", "end"))
  expect_equal(contigs$assembly, c("query", "target"))
  expect_equal(contigs$chrom, c("ctgA", "ctgA"))
  expect_equal(contigs$end, c(49186, 50001))
  whole <- jbrowse_synteny(paf)
  expect_equal(whole$regions$end, c(49186, 50001))
  expect_equal(nrow(whole$alignments), 7)
  grape <- system.file("extdata", "peach_grape.paf.gz", package = "ggjbrowse")
  d <- jbrowse_dotplot(grape)
  h <- on_axis_of(d$regions, "x")
  v <- on_axis_of(d$regions, "y")
  expect_equal(nrow(h), 13)
  expect_equal(h$chrom[1:3], c("Pp01", "Pp02", "Pp03"))
  expect_equal(v$chrom[c(1, 2, 10)], c("chr1", "chr2", "chr10"))
  expect_equal(nrow(v), 24)
  expect_gt(nrow(d$segments), 100)
  built <- ggplot2::ggplot_build(autoplot(d))
  labels <- built$layout$panel_params[[1]]$x$get_labels()
  expect_true(all(labels %in% h$chrom))
  named <- jbrowse_dotplot(grape, "Pp01", "chr1")
  expect_equal(named$regions$end[1], jbrowse_contigs(grape)$end[
    jbrowse_contigs(grape)$chrom == "Pp01"])
  expect_equal(named$regions$axis, c("x", "y"))
})

test_that("labels name each axis's assembly in its title", {
  grape <- system.file("extdata", "peach_grape.paf.gz", package = "ggjbrowse")
  named <- c(query = "peach", target = "grape")
  d <- jbrowse_dotplot(grape, labels = named)
  expect_equal(d$assemblies, c("peach", "grape"))
  expect_equal(jbrowse_dotplot(grape, labels = c("peach", "grape"))$assemblies,
    c("peach", "grape"))
  title <- function(p) ggplot2::ggplot_build(p)$layout$resolve_label(
    ggplot2::ggplot_build(p)$layout$panel_scales_x[[1]], p$labels)$primary
  expect_match(title(autoplot(d)), "^peach, 13 regions")
  swapped <- jbrowse_dotplot(grape, swap = TRUE, labels = named)
  expect_equal(swapped$assemblies, c("grape", "peach"))
  expect_equal(jbrowse_synteny(paf, labels = c(target = "volvox"))$assemblies,
    c("query", "volvox"))
  expect_error(jbrowse_dotplot(grape, labels = "peach"), "named vector")
})

test_that("swap puts the file's target on the top row and the x axis", {
  swapped <- jbrowse_synteny(paf, swap = TRUE)
  expect_equal(swapped$assemblies, c("target", "query"))
  expect_equal(on_axis_of(swapped$regions, "query")$end, 50001)
  expect_equal(on_axis_of(swapped$regions, "query")$assembly, "target")
  expect_equal(nrow(swapped$alignments), 7)
  expect_error(jbrowse_synteny(list(type = "PAFAdapter",
    pafLocation = list(uri = paf), assemblyNames = "one"), swap = TRUE),
    "two assembly names")
  d <- jbrowse_dotplot(paf, swap = TRUE)
  expect_equal(on_axis_of(d$regions, "x")$end, 50001)
})

test_that("comparative axes lay their regions end to end with no gap, as the browser does", {
  halves <- c("ctgA:1-25,000", "ctgA:25,001-50,001")
  two <- jbrowse_synteny(paf, query, halves)
  expect_equal(two$regions$cum_start[3], two$regions$cum_end[2])
  genes <- jbrowse_features(gff, halves)
  fig <- autoplot(two, below = genes)
  expect_equal(length(fig$patches$plots) + 1, 2)
})

test_that("a binned key lists every bin, present in the data or not", {
  mapq <- jbrowse_synteny(paf, query, target,
    display = list(color = list(field = "mapq")))
  expect_equal(nrow(mapq$keys$ribbons), 4)
  expect_equal(length(unique(mapq$alignments$color)), 1)
  built <- ggplot2::ggplot_build(autoplot(mapq))
  fill <- built$plot$scales$get_scales("fill")
  expect_setequal(fill$get_breaks(), as.character(mapq$keys$ribbons$color))
  expect_true(all(as.character(mapq$keys$ribbons$color) %in% fill$get_limits()))
  reads <- jbrowse_alignments(system.file("extdata", "spliced.bam",
    package = "ggjbrowse"), "ctgA:401-1,300")
  fill <- ggplot2::ggplot_build(autoplot(reads, "pileup"))$plot$scales$get_scales("fill")
  expect_true(all(fill$get_limits() %in% as.character(reads$reads$color)))
})

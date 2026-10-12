paf <- system.file("extdata", "volvox_inv_indels.paf", package = "ggjbrowse")
gff <- system.file("extdata", "volvox.sort.gff3.gz", package = "ggjbrowse")
query <- "ctgA:1-49,186"
target <- "ctgA:1-50,001"
s <- jbrowse_synteny(paf, query, target)

test_that("a PAF reads as the synteny band lays out and colors it", {
  expect_equal(nrow(s$alignments), 7)
  expect_equal(sum(s$alignments$strand == -1), 2)
  expect_equal(s$assemblies, c("query", "target"))
  expect_equal(s$alpha, 0.25)
  expect_true(all(c("match", "insertion") %in% s$ribbons$kind))
  match <- s$ribbons[s$ribbons$kind == "match", ]
  expect_equal(nrow(match), 7)
  reversed <- s$alignments$strand[match$alignment] == -1
  expect_equal(match$query_start > match$query_end, reversed)
  expect_true(all(match$target_start < match$target_end))
  expect_equal(unique(as.character(match$color)), "#ff0000")
  expect_equal(s$legend$label, c("match", "insertion"))
  expect_s3_class(s$ribbons$color, "jbrowse_color")
  expect_equal(nrow(s$offscreen), 0)
})

test_that("the view's settings color by strand, set the opacity and pick the CIGAR mode", {
  strand <- jbrowse_synteny(paf, query, target,
    view = list(color = list(field = "strand"), opacity = 0.6,
      cigarMode = "off"))
  expect_equal(unique(strand$ribbons$kind), "match")
  expect_equal(length(unique(strand$ribbons$color)), 2)
  expect_equal(strand$alpha, 0.6)
  expect_equal(strand$legend$label[1:2], c("forward", "reverse"))
  grey <- jbrowse_synteny(paf, query, target, view = list(color = "grey"))
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
  fig <- autoplot(s, target = genes)
  expect_equal(length(fig$patches$plots) + 1, 2)
  expect_equal(fig$patches$layout$heights[1], s$height)
  expect_gt(jbrowse_height(fig), 1)
  expect_error(autoplot(s, query = genes), "query row's regions")
})

test_that("a PAF reads as the dotplot lays out and colors it", {
  d <- jbrowse_dotplot(paf, query, target)
  expect_equal(nrow(d$alignments), 7)
  expect_gt(nrow(d$segments), 7)
  expect_equal(unique(as.character(d$segments$color)), "#000000")
  expect_equal(d$line_width, 2.5)
  expect_equal(d$alpha, 1)
  expect_equal(d$skipped$count, 0)
  built <- ggplot2::ggplot_build(autoplot(d))
  expect_equal(nrow(built$data[[3]]), nrow(d$segments))
  plain <- jbrowse_dotplot(paf, query, target,
    view = list(drawCigar = FALSE, color = list(field = "strand")))
  expect_equal(plain$segments$op, rep("match", 7))
  expect_equal(length(unique(plain$segments$color)), 2)
})

test_that("two regions on an axis lay end to end, each ribbon on its own", {
  halves <- c("ctgA:1-25,000", "ctgA:25,001-50,001")
  two <- jbrowse_synteny(paf, query, halves)
  expect_setequal(unique(two$ribbons$target_region), 1:2)
  axis <- on_axis(two)$ribbons
  expect_true(all(axis$target_px_start >= 0 & axis$target_px_start <= two$width))
  d <- jbrowse_dotplot(paf, query, halves)
  moved <- on_axis(d)$segments
  second <- moved$y_region == 2
  expect_true(all(moved$y_start[second] > d$y$offset[2] - d$y$start[2]))
  expect_equal(length(layout_seams(d$y)), 1)
})

test_that("an unknown extension and a missing contig say so", {
  expect_error(jbrowse_synteny("aln.tsv", query, target), "expected .paf")
  expect_equal(nrow(jbrowse_synteny(paf, query, "ctgZ:1-1,000")$alignments), 0)
})

test_that("identity keys as a bar from 0% to 100%, and mapq as chips", {
  id <- jbrowse_synteny(paf, query, target,
    view = list(color = list(field = "identity")))
  expect_length(id$ramps, 1)
  expect_equal(attr(id$ramps[[1]], "title"), "Identity")
  expect_equal(attr(id$ramps[[1]], "labels"), c("0%", "100%"))
  g <- ggplot2::ggplotGrob(autoplot(id))
  expect_false(inherits(g$grobs[[which(g$layout$name == "guide-box-right")]],
    "zeroGrob"))
  d <- jbrowse_dotplot(paf, query, target,
    view = list(color = list(field = "identity")))
  expect_equal(d$ramps, id$ramps)
  mapq <- jbrowse_synteny(paf, query, target,
    view = list(color = list(field = "mapq")))
  expect_length(mapq$ramps, 0)
  expect_gt(nrow(mapq$legend), 1)
})

swap <- system.file("extdata", "volvox_contig_swap.paf", package = "ggjbrowse")

test_that("the off-screen strip draws as the band paints it and keys its contigs", {
  skip_if(swap == "")
  o <- jbrowse_synteny(swap, "ctgA:1-50,001", "ctgC:1-1,000")
  expect_equal(o$offscreen$mate_chrom, "ctgB")
  expect_true(o$offscreen$drawn)
  expect_equal(o$offscreen_contigs$label, "ctgB · 6.08Kbp in 1 alignments")
  expect_equal(o$offscreen_strip$marks$height, 6)
  expect_equal(o$offscreen_strip$labels$label, "ctgB")
  t <- jbrowse_synteny(swap, "ctgA:1-50,001", "ctgC:1-1,000",
    view = list(color = list(field = "target")))
  expect_false(identical(t$offscreen$color, o$offscreen$color))
  built <- ggplot2::ggplot_build(autoplot(t))
  rects <- built$data[[3]]
  expect_equal(rects$ymax, 1)
  expect_equal(rects$ymin, 1 - 6 / t$height)
  expect_equal(as.character(rects$fill), as.character(t$offscreen$color))
  expect_false(is.null(offscreen_guides(t)))
  expect_null(offscreen_guides(s))
})

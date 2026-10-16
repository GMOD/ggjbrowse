paf <- system.file("extdata", "volvox_all_vs_all.paf", package = "ggjbrowse")
gff <- system.file("extdata", "volvox.sort.gff3.gz", package = "ggjbrowse")
m <- jbrowse_multiway(paf, "ctgA:1-50,001", "volvox")

test_that("an all-vs-all PAF reads as the multi-way display stacks it", {
  expect_equal(nrow(m$lanes), 3)
  expect_equal(m$lanes$assembly[1], "volvox")
  expect_null(m$assemblies)
  expect_setequal(m$lanes$assembly[-1], c("volvox_ins", "volvox_del"))
  expect_equal(m$lanes$lane, 1:3)
  expect_true(all(m$ribbons$upper_lane < m$ribbons$lower_lane))
  expect_setequal(unique(m$ribbons$gutter), 1:2)
  expect_true(all(m$ribbons$alignment %in% seq_len(nrow(m$alignments))))
  expect_type(m$ribbons$color, "character")
  expect_true(all(m$segments$x2 >= m$segments$x1))
  expect_true(all(m$segments$x1 >= 0 & m$segments$x2 <= m$width))
  expect_match(m$lanes$locus[1], "^ctgA:1-50,001$")
  expect_equal(m$regions$axis, m$lanes$assembly)
  expect_equal(m$regions$assembly, m$lanes$assembly)
  expect_null(m$anchor)
  expect_equal(m$regions$axis[m$segments$region],
    m$lanes$assembly[m$segments$lane])
  expect_equal(m$regions$axis[m$ribbons$upper_region],
    m$lanes$assembly[m$ribbons$upper_lane])
  expect_equal(m$regions$axis[m$ribbons$lower_region],
    m$lanes$assembly[m$ribbons$lower_lane])
  expect_false(anyNA(m$blocks$region))
  expect_error(jbrowse_multiway(paf, "ctgA:1-50,001"), "anchor lane's assembly")
})

test_that("an anchor of two regions lays them out as a track over them does", {
  two <- jbrowse_multiway(paf, c("ctgA:1-20,000", "ctgA:25,001-50,001"),
    "volvox")
  anchor <- on_axis_of(two$regions, "volvox")
  expect_equal(nrow(anchor), 2)
  genes <- jbrowse_features(gff, c("ctgA:1-20,000", "ctgA:25,001-50,001"))
  expect_equal(anchor$offset, genes$regions$offset)
  expect_equal(two$segments$region[two$segments$lane == 1], 1:2)
  expect_equal(two$lanes$locus[1], "ctgA:1-20,000 ctgA:25,001-50,001")
  fig <- autoplot(two, above = genes)
  expect_equal(length(fig$patches$plots) + 1, 2)
})

test_that("the anchor takes loci, frames and GRanges, and refuses whole contigs", {
  frame <- jbrowse_multiway(paf, data.frame(chrom = "ctgA", start = 0,
    end = 50001), "volvox")
  expect_equal(frame$regions, m$regions)
  expect_error(jbrowse_multiway(paf, NULL, "volvox"),
    class = "ggjbrowse_bad_region")
  expect_error(jbrowse_multiway(paf, "ctgA", "volvox"), "no extent",
    class = "ggjbrowse_bad_region")
  expect_error(jbrowse_multiway(paf, "ctgA:1-50,001", "volvox",
    lanes = c(volvox_del = "ctgA:x")), class = "ggjbrowse_bad_region")
})

test_that("assigning alignments drops their ribbons", {
  kept <- m
  kept$alignments <- m$alignments[1, ]
  expect_true(all(kept$ribbons$alignment == m$alignments$alignment[1]))
  expect_equal(nrow(kept$ribbons), sum(m$ribbons$alignment == 1))
  expect_equal(nrow(kept$segments), nrow(m$segments))
})

test_that("lanes draw in the order given, and a named lane pins its contig", {
  ordered <- jbrowse_multiway(paf, "ctgA:1-50,001", "volvox",
    lanes = c("volvox_del", "volvox_ins"))
  expect_equal(ordered$lanes$assembly, c("volvox", "volvox_del", "volvox_ins"))
  one <- jbrowse_multiway(paf, "ctgA:1-50,001", "volvox", lanes = "volvox_ins")
  expect_equal(nrow(one$lanes), 2)
  expect_equal(unique(one$ribbons$gutter), 1)
  pinned <- jbrowse_multiway(paf, "ctgA:1-50,001", "volvox",
    lanes = c(volvox_del = "ctgA"))
  expect_equal(pinned$lanes$chrom[2], "ctgA")
})

test_that("the display's settings color the ribbons", {
  strand <- jbrowse_multiway(paf, "ctgA:1-50,001", "volvox",
    display = list(ribbonColor = list(field = "strand")))
  expect_equal(strand$keys$ribbons$label[1:2],
    c("Same orientation as lane above", "Inverted vs lane above"))
  expect_equal(strand$display$ribbonColor$field, "strand")
  grey <- jbrowse_multiway(paf, "ctgA:1-50,001", "volvox",
    display = list(ribbonColor = "grey"))
  expect_equal(unique(as.character(grey$ribbons$color)), "#808080")
  expect_equal(grey$display[c("hideUnlabelled", "bridgeSkippedLanes",
    "splitStrands", "showGeneLabels", "inlineLaneNames", "drawCurves",
    "height")], list(hideUnlabelled = FALSE, bridgeSkippedLanes = TRUE,
    splitStrands = TRUE, showGeneLabels = TRUE, inlineLaneNames = FALSE,
    drawCurves = FALSE, height = 240L))
  expect_equal(grey$display$ribbonColor$value, "grey")
  ordered <- jbrowse_multiway(paf, "ctgA:1-50,001", "volvox",
    lanes = c("volvox_del", "volvox_ins"))
  expect_equal(unlist(ordered$display$rows$kept), c("volvox_del", "volvox_ins"))
  expect_null(grey[["field"]])
  expect_null(grey$curves)
})

test_that("the stack draws every lane and ribbon, with tracks over the anchor and bottom lanes", {
  p <- autoplot(m)
  built <- ggplot2::ggplot_build(p)
  expect_s3_class(built, "ggplot_built")
  polygons <- ggplot2::layer_data(p, 1)
  expect_equal(length(unique(polygons$group)), sum(!m$ribbons$thin))
  genes <- jbrowse_features(gff, m$lanes$locus[1],
    display = list(showOnlyGenes = TRUE))
  fig <- autoplot(m, above = genes)
  expect_equal(length(fig$patches$plots) + 1, 2)
  expect_s3_class(fig$theme$axis.text.x, "element_blank")
  expect_s3_class(fig$patches$plots[[1]]$theme$axis.text.x, "element_blank")
  expect_gt(jbrowse_height(fig), 1)
  curved <- ggplot2::layer_data(autoplot(m, curves = TRUE), 1)
  expect_gt(nrow(curved), nrow(polygons))
})

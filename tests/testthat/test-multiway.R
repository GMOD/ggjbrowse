paf <- system.file("extdata", "volvox_all_vs_all.paf", package = "ggjbrowse")
gff <- system.file("extdata", "volvox.sort.gff3.gz", package = "ggjbrowse")
m <- jbrowse_multiway(paf, "ctgA:1-50,001", "volvox")

test_that("an all-vs-all PAF reads as the multi-way display stacks it", {
  expect_equal(nrow(m$lanes), 3)
  expect_equal(m$lanes$assembly[1], "volvox")
  expect_setequal(m$lanes$assembly[-1], c("volvox_ins", "volvox_del"))
  expect_equal(m$lanes$lane, 1:3)
  expect_true(all(m$ribbons$upper_lane < m$ribbons$lower_lane))
  expect_setequal(unique(m$ribbons$gutter), 1:2)
  expect_true(all(m$ribbons$alignment %in% seq_len(nrow(m$alignments))))
  expect_s3_class(m$ribbons$color, "jbrowse_color")
  expect_true(all(m$segments$x2 >= m$segments$x1))
  expect_true(all(m$segments$x1 >= 0 & m$segments$x2 <= m$width))
  expect_match(m$lanes$locus[1], "^ctgA:1-50,001$")
  expect_error(jbrowse_multiway(paf, "ctgA:1-50,001"), "anchor lane's assembly")
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
  expect_equal(strand$legend$label[1:2],
    c("Same orientation as lane above", "Inverted vs lane above"))
  grey <- jbrowse_multiway(paf, "ctgA:1-50,001", "volvox",
    display = list(ribbonColor = "grey"))
  expect_equal(unique(as.character(grey$ribbons$color)), "#808080")
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
  expect_gt(jbrowse_height(fig), 1)
  curved <- ggplot2::layer_data(autoplot(m, curves = TRUE), 1)
  expect_gt(nrow(curved), nrow(polygons))
})

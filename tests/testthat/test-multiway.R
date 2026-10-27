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
  mixed <- jbrowse_multiway(paf, "ctgA:1-50,001", "volvox",
    lanes = c(volvox_del = "ctgA", "volvox_ins"))
  expect_equal(mixed$lanes$assembly, c("volvox", "volvox_del", "volvox_ins"))
  expect_equal(lane_choice(c(volvox_del = "ctgA", "volvox_ins"))$pins,
    list(list(lane = "volvox_del", refName = "ctgA")))
  expect_equal(lane_choice(list(a = "ctgA:11-20", b = NULL))$pins,
    list(list(lane = "a", refName = "ctgA", at = 15)))
})

test_that("the display's settings color the ribbons", {
  strand <- jbrowse_multiway(paf, "ctgA:1-50,001", "volvox",
    display = list(ribbons = list(encoding = list(color = list(field = "strand")))))
  expect_equal(strand$keys$ribbons$label[1:2],
    c("Same orientation as lane above", "Inverted vs lane above"))
  expect_equal(strand$display$ribbons$encoding$color$field, "strand")
  grey <- jbrowse_multiway(paf, "ctgA:1-50,001", "volvox",
    display = list(ribbons = list(encoding = list(color = list(value = "grey")))))
  expect_equal(unique(as.character(grey$ribbons$color)), "#808080")
  expect_equal(grey$display$labels$lanes$placement, "above")
  expect_equal(grey$display$height, 240L)
  expect_true(grey$display$labels$name$show)
  expect_equal(grey$display$ribbons$mark$interpolate, "linear")
  expect_equal(grey$display$ribbons$encoding$color$value, "grey")
  faded <- jbrowse_multiway(paf, "ctgA:1-50,001", "volvox",
    display = list(ribbons = list(encoding = list(opacity = list(value = 0.8)))))
  expect_equal(unique(faded$ribbons$alpha), 0.8)
  inverted <- jbrowse_multiway(paf, "ctgA:1-50,001", "volvox",
    display = list(ribbons = list(filter = list(list(field = "strand", equal = -1)))))
  expect_equal(nrow(inverted$ribbons), 0)
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
    display = list(filter = "geneLike"))
  fig <- autoplot(m, above = genes)
  expect_equal(length(fig$patches$plots) + 1, 2)
  expect_s3_class(fig$theme$axis.text.x, "element_blank")
  expect_s3_class(fig$patches$plots[[1]]$theme$axis.text.x, "element_blank")
  expect_equal(fig$patches$plots[[1]]$coordinates$limits$y, c(0, genes$height))
  expect_gt(jbrowse_height(fig), 1)
  bottom <- nrow(m$lanes)
  expect_gt(min(m$segments$x1[m$segments$lane == bottom]), 1)
  under <- jbrowse_features(gff, m$lanes$locus[bottom],
    display = list(filter = "geneLike"))
  band <- lane_band(under, m, bottom, "below")$plots[[1]]
  xlim <- ggplot2::ggplot_build(band)$layout$panel_params[[1]]$x.range
  expect_equal(xlim, c(0, m$width), tolerance = 1e-6)
  curved <- ggplot2::layer_data(autoplot(m, curves = TRUE), 1)
  expect_gt(nrow(curved), nrow(polygons))
  bumped <- m
  bumped$display$ribbons$mark$interpolate <- "bump"
  expect_equal(nrow(ggplot2::layer_data(autoplot(bumped), 1)), nrow(curved))
})

# The PAF is reference-anchored: HG01361 aligns to GRCh38 alone, so a lane
# under HG00133 has no records, and a read of it fails.
test_that("an inversion draws as its records state it, its ribbons crossed", {
  inv <- jbrowse_multiway(
    system.file("extdata", "hprc_amylase.paf", package = "ggjbrowse"),
    "chr1:103,500,001-103,900,000", "GRCh38", lanes = "HG00133",
    display = list(ribbons = list(encoding = list(color = list(field = "strand")))))
  expect_true(any(inv$blocks$strand == -1))
  expect_false(any(inv$lanes$flipped))
  expect_true(all(inv$segments$x2 >= inv$segments$x1))
  r <- inv$ribbons[inv$ribbons$kind == "match", ]
  crossed <- sign(r$upper_x2 - r$upper_x1) != sign(r$lower_x2 - r$lower_x1)
  inverted <- inv$keys$ribbons$color[
    inv$keys$ribbons$label == "Inverted vs lane above"]
  expect_true(any(crossed))
  expect_equal(crossed, r$color == inverted)
  expect_s3_class(ggplot2::ggplot_build(autoplot(inv)), "ggplot_built")
})

test_that("an anchor contig or kept lane the file lacks is refused, not read empty", {
  expect_error(jbrowse_multiway(paf, "ctgZ:1-50,001", "volvox"),
    "ctgZ is not a contig", class = "ggjbrowse_missing_contig")
  expect_error(jbrowse_multiway(paf, "ctgA:1-50,001", "volvox",
    lanes = "volvox_typo"), "no lane for")
})

extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
gff <- extdata("volvox.sort.gff3.gz")
g <- jbrowse_features(gff, "ctgA", 0, 25000)

test_that("a GFF3 region reads as JBrowse's feature display lays it out", {
  f <- g$features
  expect_gt(nrow(f), 50)
  for (r in split(f, f$top)) {
    r <- r[order(r$start), ]
    expect_true(all(utils::tail(r$start, -1) >= utils::head(r$end, -1)))
  }
  r <- g$rects
  expect_true(all(r$start >= f$start[r$feature] & r$end <= f$end[r$feature]))
  expect_true(all(r$top >= f$top[r$feature] &
    r$top + r$height <= f$bottom[r$feature] + 1e-3))
  eden <- which(f$name == "EDEN")
  expect_equal(g$subfeatures$name[g$subfeatures$feature == eden],
    c("EDEN.1", "EDEN.2", "EDEN.3"))
  expect_true("EDEN" %in% g$labels$text)
  expect_match(r$color, "^#[0-9a-f]{6}$")
  expect_gt(nrow(g$chevrons), 0)
  expect_gt(nrow(g$arrows), 0)
})

test_that("the display's settings reach the tables", {
  bare <- jbrowse_features(gff, "ctgA", 0, 25000,
    display = list(showLabels = "none", displayDirectionalChevrons = FALSE))
  expect_equal(nrow(bare$labels), 0)
  expect_equal(nrow(bare$chevrons), 0)
  expect_lt(bare$height, g$height)
  genes <- jbrowse_features(gff, "ctgA", 0, 25000,
    display = list(showOnlyGenes = TRUE))
  expect_true(all(genes$features$gene))
  typed <- jbrowse_features(gff, "ctgA", 0, 25000,
    display = list(color = list(field = "type")))
  purple <- jbrowse_features(gff, "ctgA", 0, 25000,
    display = list(color = "purple"))
  expect_equal(unique(as.character(purple$rects$color)), "#800080")
  expect_gt(nrow(typed$legend), 1)
  expect_true(all(typed$legend$color %in% typed$rects$color))
  expect_equal(attr(typed$legend, "title"), "type")
})

test_that("a wider layout fits more names", {
  narrow <- jbrowse_features(gff, "ctgA", 0, 25000, width = 300)
  expect_lt(nrow(narrow$labels), nrow(g$labels))
})

test_that("a BigBed reads through the same display", {
  b <- jbrowse_features(extdata("volvox.bb"), "ctgA", 0, 25000)
  expect_gt(nrow(b$features), 0)
  expect_gt(nrow(b$rects), 0)
})

test_that("a track's adapter config reads with its own settings", {
  bb <- list(type = "BigBedAdapter",
    bigBedLocation = list(localPath = extdata("volvox.bb")),
    aggregateField = "geneName")
  b <- jbrowse_features(bb, "ctgA", 0, 25000)
  eden <- which(b$features$name == "EDEN")
  expect_equal(b$features$type[eden], "gene")
  expect_setequal(b$subfeatures$name[b$subfeatures$feature == eden],
    c("EDEN.1", "EDEN.2", "EDEN.3"))
})

test_that("autoplot draws every box, line and label", {
  p <- ggplot2::ggplot_build(ggplot2::autoplot(g))
  rows <- vapply(p$data, nrow, 1L)
  expect_equal(rows[1], nrow(g$lines))
  expect_equal(rows[3], nrow(g$rects))
  expect_equal(rows[length(rows)], nrow(g$labels))
  expect_setequal(unique(p$data[[3]]$fill), unique(g$rects$color))
})

test_that("labels take the display's pixel size at the panel's scale", {
  p <- ggplot2::autoplot(g)
  text <- ggplot2::layer_grob(p, i = length(p$layers))[[1]]
  expect_s3_class(text, "ggjbrowse_display_text")
  font_at <- function(w, h) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    grid::pushViewport(grid::viewport(width = grid::unit(w, "mm"),
      height = grid::unit(h, "mm")))
    unique(grid::makeContent(text)$children[[1]]$gp$fontsize)
  }
  pt <- g$glyph$fontSize * ggplot2::.pt
  expect_equal(font_at(100, 1000), pt * 100 / 800)
  expect_equal(font_at(1000, 10), pt * 10 / g$height)
})

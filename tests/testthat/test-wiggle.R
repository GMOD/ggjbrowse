bw <- system.file("extdata", "volvox_microarray_posneg.bw",
  package = "ggjbrowse")
w <- jbrowse_wiggle(bw, "ctgA", 0, 50000)

test_that("a BigWig region reads as the wiggle display scales and colors it", {
  expect_gt(nrow(w$bins), 100)
  expect_equal(nrow(w$marks), nrow(w$bins))
  expect_equal(w$rendering, "xyplot")
  expect_lte(w$domain[1], min(w$bins$score))
  expect_gte(w$domain[2], max(w$bins$score))
  sides <- split(w$marks$color, w$marks$score < 0)
  expect_equal(lengths(lapply(sides, unique)), c(`FALSE` = 1, `TRUE` = 1))
  expect_false(identical(sides[[1]][1], sides[[2]][1]))
  pinned <- jbrowse_wiggle(bw, "ctgA", 0, 50000,
    display = list(scales = list(y = list(domainMin = -50, domainMax = 300))))
  expect_equal(pinned$domain, c(-50, 300))
})

test_that("a zoomed-out min-max extent draws the max, mean and min bands on each side", {
  wide <- jbrowse_wiggle(bw, "ctgA", 0, 50000, width = 100)
  expect_equal(length(unique(wide$marks$layer)), 6)
  expect_true(all(wide$bins$min <= wide$bins$score &
    wide$bins$score <= wide$bins$max))
  line <- jbrowse_wiggle(bw, "ctgA", 0, 50000, width = 100,
    display = list(mark = "line"))
  band <- line$marks[line$marks$band, ]
  expect_gte(nrow(band), nrow(line$bins))
  expect_true(all(band$low <= band$score))
  expect_true(all(band$opacity < 1))
  layers <- ggplot2::ggplot_build(ggplot2::autoplot(line))$data
  expect_equal(nrow(layers[[1]]), nrow(band))
})

test_that("each plot type draws every mark", {
  for (mark in c("bar", "point", "line", "span")) {
    x <- jbrowse_wiggle(bw, "ctgA", 0, 50000, display = list(mark = mark,
      extent = "none"))
    built <- ggplot2::ggplot_build(ggplot2::autoplot(x))
    expect_equal(sum(vapply(built$data[2:6], nrow, 1L)),
      nrow(x$marks) * if (mark == "line") 2 else 1)
  }
})

test_that("an empty region draws, whatever the plot type", {
  for (mark in c("bar", "line", "span")) {
    empty <- jbrowse_wiggle(bw, "ctgB", 0, 1000, display = list(mark = mark))
    expect_equal(nrow(empty$marks), 0)
    expect_s3_class(ggplot2::ggplot_build(ggplot2::autoplot(empty)),
      "ggplot_built")
  }
})

test_that("a stepped line joins its bins, and a rule is labelled in its color", {
  x <- jbrowse_wiggle(bw, "ctgA", 0, 50000, display = list(mark = "line",
    extent = "none", scales = list(y = list(rules = list(
      list(value = 500, label = "gain"), list(value = 0))))))
  expect_equal(x$rules$label, c("gain", ""))
  expect_equal(as.character(x$rules$color), c("#787878", "#787878"))
  built <- ggplot2::ggplot_build(ggplot2::autoplot(x))
  path <- built$data[[5]]
  expect_equal(length(unique(path$group)), 1)
  expect_equal(path$x[2], path$x[3])
  labels <- built$data[[8]]
  expect_equal(labels$label, "gain")
})

test_that("a BigWig stacks over reads as one of the track's tracks", {
  bam <- system.file("extdata", "volvox-sorted.bam", package = "ggjbrowse")
  reads <- jbrowse_alignments(bam, "ctgA", 1000, 6000)
  track <- ggplot2::autoplot(reads,
    tracks = jbrowse_wiggle(bw, "ctgA", 1000, 6000))
  expect_equal(length(track$patches$plots) + 1, 3)
})

test_that("the score axis takes the display's log or symlog transform", {
  axis <- function(type) {
    w <- jbrowse_wiggle(bw, "ctgA", 0, 50000,
      display = list(scales = list(y = list(type = type))))
    scale_y_jbrowse(w)$get_transformation()$name
  }
  expect_equal(axis("log"), "log-10")
  expect_equal(axis("symlog"), "symlog")
  expect_equal(axis("linear"), "identity")
})

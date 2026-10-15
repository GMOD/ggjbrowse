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

test_that("density keys its ramp as a bar beside the plot, and bars key none", {
  d <- jbrowse_wiggle(bw, "ctgA", 0, 50000, display = list(mark = "span"))
  expect_length(d$ramps, 1)
  ramp <- d$ramps[[1]]
  expect_equal(attr(ramp, "title"), "Score")
  expect_equal(attr(ramp, "labels"), as.character(d$domain))
  expect_equal(range(ramp$offset), c(0, 1))
  expect_length(w$ramps, 0)
  keyed <- function(x) {
    g <- ggplot2::ggplotGrob(ggplot2::autoplot(x))
    !inherits(g$grobs[[which(g$layout$name == "guide-box-right")]], "zeroGrob")
  }
  expect_true(keyed(d))
  expect_false(keyed(w))
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
    drawn <- sum(vapply(built$data[2:6], nrow, 1L))
    if (mark == "line") expect_gte(drawn, nrow(x$marks) * 2)
    else expect_equal(drawn, nrow(x$marks))
  }
})

test_that("a line changes color where it crosses the origin, not at the next point", {
  d <- data.frame(x = c(0, 10, 20), score = c(2, -2, -4), run = 0,
    color = c("blue", "red", "red"))
  cut <- cut_at_origin(d, 0)
  expect_equal(cut$x, c(0, 5, 10, 20))
  expect_equal(cut$score, c(2, 0, -2, -4))
  expect_equal(cut$color, c("blue", "red", "red", "red"))
  expect_identical(cut_at_origin(d[2:3, ], 0), d[2:3, ])
  apart <- d
  apart$run <- c(0, 1, 1)
  expect_identical(cut_at_origin(apart, 0), apart)
})

test_that("a stepped line's riser across the origin is cut there too", {
  line <- jbrowse_wiggle(bw, "ctgA", 0, 50000, width = 100,
    display = list(mark = "line"))
  steps <- ggplot2::ggplot_build(ggplot2::autoplot(line))$data[[5]]
  at <- which(steps$y == 0)
  at <- at[at > 1 & at < nrow(steps) & steps$group[at - 1] == steps$group[at + 1]]
  expect_gt(length(at), 0)
  expect_true(all(steps$y[at - 1] * steps$y[at + 1] < 0))
  expect_equal(steps$colour[at], steps$colour[at + 1])
})

test_that("an empty region draws, whatever the plot type", {
  for (mark in c("bar", "line", "span")) {
    expect_warning(
      empty <- jbrowse_wiggle(bw, "ctgB", 0, 1000, display = list(mark = mark)),
      "ctgB is not a contig of this file")
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
    above = jbrowse_wiggle(bw, "ctgA", 1000, 6000))
  expect_equal(length(track$patches$plots) + 1, 3)
})

test_that("the score axis takes the display's log or symlog transform", {
  axis <- function(type) {
    w <- jbrowse_wiggle(bw, "ctgA", 0, 50000,
      display = list(scales = list(y = list(type = type))))
    scale_y_jbrowse(w)$get_transformation()$name
  }
  expect_equal(axis("log"), "log-10")
  w <- jbrowse_wiggle(bw, "ctgA", 0, 50000,
    display = list(scales = list(y = list(type = "symlog"))))
  expect_equal(scale_y_jbrowse(w, name = "Score")$name, "Score")
  expect_equal(scale_y_jbrowse(w, name = "Score")$get_transformation()$name,
    "symlog")
  expect_equal(axis("symlog"), "symlog")
  expect_equal(axis("linear"), "identity")
})

test_that("a row layout packs its sources into the display's height, a slice each with no gap", {
  bw <- normalizePath(system.file("extdata", "volvox_microarray_posneg.bw",
    package = "ggjbrowse"))
  multi <- function(n) list(type = "MultiWiggleAdapter",
    subadapters = lapply(seq_len(n), function(i) list(type = "BigWigAdapter",
      name = sprintf("s%03d", i), bigWigLocation = list(uri = bw))))
  w <- jbrowse_wiggle(multi(120), "ctgA", 0, 50000,
    display = list(mark = "span"))
  expect_true(w$row_layout)
  expect_equal(jbrowse_stack(w)$patches$layout$heights, w$height)
  p <- autoplot(w)
  expect_null(p$facet$params$rows)
  tiles <- ggplot2::layer_data(p, 3)
  expect_setequal(unique(tiles$ymin), 0:119)
  expect_equal(tiles$ymax - tiles$ymin, rep(1, nrow(tiles)))
  expect_s3_class(p$theme$axis.text.y, "element_blank")
  two <- autoplot(jbrowse_wiggle(multi(2), "ctgA", 0, 50000))
  expect_false(inherits(two$theme$axis.text.y, "element_blank"))
})

test_that("several regions share one score axis and color ramp", {
  loci <- c("ctgA:1-10000", "ctgA:20001-30000")
  d <- jbrowse_wiggle(bw, loci, display = list(mark = "span"))
  m <- d$marks
  colors <- tapply(m$color, m$score, function(c) length(unique(c)))
  expect_true(all(colors == 1))
  expect_true(any(duplicated(m$score[!duplicated(m[c("score", ".region")])])))
  expect_equal(as.numeric(attr(d$ramps[[1]], "labels")), d$domain)
  alone <- jbrowse_wiggle(bw, loci[1], display = list(mark = "span"))
  expect_equal(d$marks[d$marks$.region == 1, c("start", "end", "score")],
    alone$marks[c("start", "end", "score")])
  expect_false(isTRUE(all.equal(d$domain, alone$domain)))
})

test_that("as.matrix() bins each source into the display's pixel columns", {
  sources <- list(type = "MultiWiggleAdapter", subadapters = lapply(1:3,
    function(i) list(type = "BigWigAdapter", name = paste0("rep", i),
      bigWigLocation = list(uri = normalizePath(bw)))))
  m3 <- jbrowse_wiggle(sources, "ctgA:1-50,000", display = list(mark = "span"))
  m <- as.matrix(m3)
  expect_equal(dim(m), c(3, 800))
  expect_equal(rownames(m), c("rep1", "rep2", "rep3"))
  b <- m3$bins[m3$bins$source == "rep1", ]
  covering <- function(k) {
    from <- pmax(0, trunc(b$start / m3$bp_per_px))
    to <- pmax(trunc(b$end / m3$bp_per_px), from + 1)
    b$score[from < k & to >= k]
  }
  expect_equal(unname(m[1, c(1, 2, 400)]),
    vapply(c(1, 2, 400), function(k) mean(covering(k)), 0))
  two <- jbrowse_wiggle(sources, c("ctgA:1-10,000", "ctgA:20,001-30,000"),
    display = list(mark = "span"))
  expect_equal(colnames(as.matrix(two))[c(1, 401)], c("ctgA:1", "ctgA:20001"))

  m3$tree <- stats::hclust(stats::dist(m + c(0, 5, 1)), "average")
  drawn <- on_axis(m3)
  leaves <- m3$tree$labels[m3$tree$order]
  expect_equal(drawn$sources$name, leaves)
  expect_equal(unique(drawn$marks$row[drawn$marks$source == leaves[1]]), 0)
  p <- autoplot(m3)
  trees <- Filter(function(l) inherits(l$geom, "GeomTree"), p$layers)
  expect_length(trees, 1)
  expect_no_error(ggplot2::ggplotGrob(p))
  m3$tree <- stats::hclust(stats::dist(m[1:2, ]))
  expect_error(autoplot(m3), "label it by x\\$sources\\$name")

  m3$tree <- NULL
  m3$sources <- m3$sources[2:3, ]
  expect_setequal(unique(on_axis(m3)$marks$source), c("rep2", "rep3"))
  expect_equal(rownames(as.matrix(m3)), c("rep2", "rep3"))
})

test_that("a mark narrower than a pixel draws a pixel wide", {
  sparse <- w
  sparse$marks <- w$marks[1, ]
  sparse$marks$end <- sparse$marks$start + 1
  drawn <- ggplot2::layer_data(ggplot2::ggplot(sparse) +
    geom_jbrowse_wiggle(), 2)
  expect_equal(drawn$xmax - drawn$xmin, w$bp_per_px)
  wide <- ggplot2::layer_data(ggplot2::ggplot(w) + geom_jbrowse_wiggle(), 2)
  expect_equal(wide$xmax - wide$xmin, w$marks$end - w$marks$start)
})

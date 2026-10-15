extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
people <- extdata("1000g_dcst2.vcf.gz")
dcst2 <- "1:155,000,001-155,050,000"
read_people <- function(display = NULL) {
  jbrowse_variants(people, dcst2, samples = extdata("1000g_dcst2_samples.tsv"),
    display = display)
}

layer_of <- function(p, geom) {
  i <- which(vapply(p$layers, function(l) inherits(l$geom, geom), NA))
  ggplot2::layer_data(p, i[1])
}

test_that("the rows take a saved arrangement: its order, labels and focus", {
  plain <- read_people()
  last <- utils::tail(plain$rows$name, 2)
  v <- read_people(list(rows = list(domain = rev(last),
    labels = stats::setNames(list("Proband"), last[2]))))
  expect_equal(v$rows$name[1:2], rev(last))
  expect_equal(v$rows$label[1], "Proband")
  expect_equal(v$rows$label[-1], v$rows$name[-1])
  p <- ggplot2::autoplot(v, layout = "columns")
  labels <- ggplot2::ggplot_build(p)$data
  expect_equal(labels[[length(labels)]]$label[1], "Proband")
  focused <- read_people(list(rows = list(kept = I(last))))
  expect_setequal(focused$rows$name, last)
  expect_setequal(as.character(unique(focused$cells$row)), last)
})

test_that("a facet stacks the rows in bands that hold whatever order they take", {
  v <- read_people(list(facet = list(field = "superpopulation",
    domain = I("European"))))
  expect_equal(v$bands$label[1], "European")
  expect_equal(rle(v$rows$band)$values, v$bands$key)
  expect_true(all(v$rows$band[1:20] == "European"))
  shuffled <- v
  shuffled$rows <- v$rows[rev(seq_len(nrow(v$rows))), ]
  drawn <- on_axis(shuffled)$rows
  expect_equal(rle(drawn$band)$values, v$bands$key)
  expect_equal(drawn$name[1:20], rev(v$rows$name[1:20]))
  p <- ggplot2::autoplot(v)
  bands <- layer_of(p, "GeomBandLabels")
  expect_equal(bands$label, v$bands$label)
  expect_no_error(ggplot2::ggplotGrob(p))
  expect_null(row_bands(read_people()$rows, read_people()$bands))
})

test_that("under bands a forest draws a clade per band, and a whole tree none", {
  v <- read_people(list(facet = "superpopulation"))
  m <- as.matrix(v)
  v$tree <- lapply(split(rownames(m), v$rows$band), function(r)
    stats::hclust(stats::dist(m[r, ]), "average"))
  p <- ggplot2::autoplot(v, layout = "columns")
  expect_equal(nrow(layer_of(p, "GeomTree")), 3 * (nrow(v$rows) - 5))
  drawn <- on_axis(v)$rows
  expect_equal(rle(drawn$band)$lengths, rep(20, 5))
  expect_no_error(ggplot2::ggplotGrob(p))
  whole <- v
  whole$tree <- stats::hclust(stats::dist(m), "average")
  expect_equal(nrow(layer_of(ggplot2::autoplot(whole), "GeomTree")), 0)
  part <- read_people()
  part$tree <- list(stats::hclust(stats::dist(m[5:10, ])))
  expect_equal(on_axis(part)$rows$name[1:6],
    part$tree[[1]]$labels[part$tree[[1]]$order])
})

test_that("zoomed out, records sharing a pixel fade as the display fades them", {
  v <- read_people()
  alpha <- function(...) {
    p <- ggplot2::ggplot(v) + geom_jbrowse_genotypes(...)
    ggplot2::layer_data(p)$alpha
  }
  far <- alpha(width = 300)
  expect_lt(min(far), 1)
  expect_gt(min(far), 0)
  expect_gt(mean(alpha(width = 1e5)), mean(far))
  expect_true(all(alpha(width = 300, fade = FALSE) == 1))
  expect_true(all(alpha(width = 300, layout = "columns") == 1))
  v$variants <- v$variants[c(1, 1), ]
  v$variants$end[1] <- v$variants$start[1] + 30000
  expect_equal(density_alpha(v$variants, 50)[1], 1)
})

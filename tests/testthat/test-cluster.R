extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
vcf <- extdata("dog10k_cyp1a2_snvs.vcf.gz")
locus <- "chr30:38,258,001-38,265,000"

clustered <- function() {
  v <- jbrowse_variants(vcf, locus)
  v$tree <- stats::hclust(stats::dist(stats::xtabs(dosage ~ row + variant,
    v$cells)), "average")
  v
}

test_that("the cells' dosages make the matrix the README clusters", {
  v <- jbrowse_variants(vcf, locus)
  m <- stats::xtabs(dosage ~ row + variant, v$cells)
  expect_equal(dim(m), c(nrow(v$rows), length(unique(v$cells$variant))))
  het <- v$cells[v$cells$genotype == "0/1", ][1, ]
  expect_equal(m[as.character(het$row), het$variant], 0.5)
  hom <- v$cells[v$cells$genotype == "1/1", ][1, ]
  expect_equal(m[as.character(hom$row), hom$variant], 1)
})

test_that("a tree of one's own draws beside the rows, each leaf at its row", {
  v <- clustered()
  n <- nrow(v$rows)
  drawn <- on_axis(v)
  segments <- tree_segments(stats::as.hclust(v$tree), drawn$rows$name)
  expect_equal(nrow(segments), 3 * (n - 1))
  leaves <- segments[segments$h == 0 & segments$y == segments$yend &
    segments$y == round(segments$y), ]
  expect_setequal(leaves$y, seq_len(n))
  expect_equal(max(segments$hend), 1)
  p <- ggplot2::autoplot(v, layout = "columns")
  trees <- Filter(function(l) inherits(l$geom, "GeomTree"), p$layers)
  expect_length(trees, 1)
  expect_length(Filter(function(l) inherits(l$geom, "GeomTree"),
    ggplot2::autoplot(jbrowse_variants(vcf, locus))$layers), 0)
  expect_no_error(ggplot2::ggplotGrob(p))
  own <- ggplot2::ggplot(v) + geom_jbrowse_genotypes(layout = "columns") +
    geom_jbrowse_tree(ggplot2::aes(colour = I("red")), span = 20) +
    ggplot2::coord_cartesian(clip = "off")
  expect_equal(nrow(ggplot2::layer_data(own, 2)), 3 * (n - 1))
  expect_no_error(ggplot2::ggplotGrob(own))
  stacked <- jbrowse_stack(jbrowse_features(extdata("dog10k_cyp1a2_snvs.vcf.gz"),
    locus), v)
  expect_no_error(patchwork::patchworkGrob(stacked))
})

test_that("rows draw in x$rows' order, or a tree's leaf order where x$tree holds one", {
  v <- jbrowse_variants(vcf, locus)
  y_of <- function(x) {
    p <- ggplot2::ggplot(x) + geom_jbrowse_genotypes()
    d <- p$layers[[1]]$layer_data(p$data)
    stats::setNames(d$y, d$sample)
  }
  reversed <- v
  reversed$rows <- v$rows[rev(seq_len(nrow(v$rows))), ]
  y <- y_of(reversed)
  expect_true(all(y[names(y) == reversed$rows$name[1]] == 1))
  expect_true(all(y[names(y) == v$rows$name[1]] == nrow(v$rows)))
  labels <- ggplot2::ggplot_build(ggplot2::autoplot(reversed))$data
  expect_equal(labels[[length(labels)]]$label, reversed$rows$name)
  kept <- v
  kept$rows <- v$rows[1:5, ]
  expect_setequal(names(y_of(kept)), v$rows$name[1:5])
  own <- v
  m <- stats::xtabs(dosage ~ row + variant, v$cells)
  own$tree <- stats::hclust(stats::dist(m), "complete")
  expect_equal(on_axis(own)$rows$name,
    own$tree$labels[own$tree$order])
  expect_equal(levels(on_axis(own)$cells$row), on_axis(own)$rows$name)
  own$tree <- stats::as.dendrogram(own$tree)
  expect_no_error(ggplot2::ggplotGrob(ggplot2::autoplot(own)))
  wrong <- v
  wrong$tree <- stats::hclust(stats::dist(m[1:5, ]))
  expect_error(ggplot2::autoplot(wrong), "label it by x\\$rows\\$name")
})

test_that("a matrix with no rows draws an empty panel", {
  v <- jbrowse_variants(vcf, locus)
  v$rows <- v$rows[0, ]
  expect_no_error(ggplot2::ggplotGrob(ggplot2::autoplot(v)))
})

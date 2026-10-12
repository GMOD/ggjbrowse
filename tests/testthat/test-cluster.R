extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
vcf <- extdata("dog10k_cyp1a2_snvs.vcf.gz")
locus <- "chr30:38,258,001-38,265,000"

clustered <- function() {
  v <- jbrowse_variants(vcf, locus)
  v$tree <- stats::hclust(stats::dist(as.matrix(v)), "average")
  v
}

test_that("as.matrix() scores each ALT as the display's clustering does", {
  v <- jbrowse_variants(vcf, locus)
  m <- as.matrix(v)
  alts <- lengths(strsplit(v$variants$alt, ","))
  expect_equal(dim(m), c(nrow(v$rows), sum(alts)))
  expect_equal(rownames(m), as.character(v$rows$name))
  at <- function(cell, column = as.character(cell$variant)) {
    unname(m[as.character(cell$row), column])
  }
  biallelic <- v$cells[alts[v$cells$variant] == 1, ]
  expect_equal(at(biallelic[biallelic$genotype == "0/1", ][1, ]), 1)
  expect_equal(at(biallelic[biallelic$genotype == "1/1", ][1, ]), 2)
  second <- v$cells[v$cells$genotype == "0/2", ][1, ]
  expect_equal(at(second, paste0(second$variant, c(".1", ".2"))), c(0, 1))
  nocall <- v$cells[v$cells$genotype == "./.", ][1, ]
  expect_true(all(is.na(at(nocall, grep(paste0("^", nocall$variant, "(\\.|$)"),
    colnames(m))))))
  drawn <- paste(v$cells$row, v$cells$variant)
  skipped <- which(!paste(v$rows$name[1], seq_along(alts)) %in% drawn &
    alts == 1)[1]
  expect_equal(unname(m[1, as.character(skipped)]), 0)
})

test_that("as.matrix() of haplotype rows is each haplotype's ALT", {
  h <- jbrowse_variants(vcf, locus, display = list(unit = "haplotype"))
  m <- as.matrix(h)
  expect_equal(dim(m), c(nrow(h$rows), nrow(h$variants)))
  cells <- h$cells[!is.na(h$cells$dosage), ]
  expect_equal(unname(m[cbind(as.character(cells$row),
    as.character(cells$variant))]), cells$dosage)
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
  m <- as.matrix(v)
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

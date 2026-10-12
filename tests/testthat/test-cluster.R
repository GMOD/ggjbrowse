extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
vcf <- extdata("dog10k_cyp1a2_snvs.vcf.gz")
locus <- "chr30:38,258,001-38,265,000"

test_that("a genotype counts each ALT's dosage, reference where no cell is", {
  v <- jbrowse_variants(vcf, locus)
  m <- genotype_matrix(v)
  alts <- lengths(strsplit(v$variants$alt, ","))
  expect_equal(dim(m), c(nrow(v$rows), sum(alts)))
  expect_false(anyNA(m))
  first <- cumsum(c(0, alts))
  at <- function(genotype) {
    cell <- v$cells[v$cells$genotype == genotype, ][1, ]
    unname(m[as.character(cell$row),
      first[cell$variant] + seq_len(alts[cell$variant])])
  }
  expect_equal(at("0/1")[1], 1)
  expect_equal(at("1|1")[1], 2)
  expect_equal(at("0/2"), c(0, 1))
  expect_equal(at("2/2"), c(0, 2))
  no_call <- v$cells[v$cells$genotype == "./.", ][1, ]
  column <- first[no_call$variant] + 1
  others <- setdiff(v$rows$name, as.character(no_call$row))
  expect_equal(unname(m[as.character(no_call$row), column]),
    mean(m[others, column]))
  hv <- jbrowse_variants(vcf, locus, display = list(unit = "haplotype"))
  h <- genotype_matrix(hv)
  phased <- hv$cells[hv$cells$genotype == "0|1", ][1, ]
  column <- first[phased$variant] + 1
  expect_equal(unname(h[paste(phased$sample, c("HP0", "HP1")), column]),
    c(0, 2))
})

test_that("clustering orders the rows as JBrowse's R script does, cells following", {
  v <- jbrowse_variants(vcf, locus)
  c1 <- jbrowse_cluster(v)
  expected <- stats::hclust(stats::dist(genotype_matrix(v)), "average")
  expect_equal(c1$rows$name, v$rows$name[expected$order])
  expect_equal(c1$tree$merge, expected$merge)
  expect_equal(levels(c1$cells$row), c1$rows$name)
  expect_equal(as.character(c1$cells$row), as.character(v$cells$row))
  expect_equal(jbrowse_cluster(v, "complete")$tree$method, "complete")
  expect_error(jbrowse_cluster(jbrowse_wiggle(
    extdata("volvox_microarray_posneg.bw"), "ctgA:1-5,000")),
    "jbrowse_variants")
})

test_that("a clustered matrix draws its tree beside the rows, each leaf at its row", {
  c1 <- jbrowse_cluster(jbrowse_variants(vcf, locus))
  n <- nrow(c1$rows)
  segments <- tree_segments(c1$tree)
  expect_equal(nrow(segments), 3 * (n - 1))
  leaves <- segments[segments$h == 0 & segments$y == segments$yend &
    segments$y == round(segments$y), ]
  expect_setequal(leaves$y, seq_len(n))
  expect_equal(max(segments$hend), 1)
  p <- ggplot2::autoplot(c1, layout = "columns")
  trees <- Filter(function(l) inherits(l$geom, "GeomTree"), p$layers)
  expect_length(trees, 1)
  expect_length(Filter(function(l) inherits(l$geom, "GeomTree"),
    ggplot2::autoplot(jbrowse_variants(vcf, locus))$layers), 0)
  expect_no_error(ggplot2::ggplotGrob(p))
})

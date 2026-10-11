extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
bam <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA", 1000, 3000)
library(ggplot2)

test_that("a display color column takes an identity scale, any other column its own", {
  default <- ggplot_build(ggplot(bam) + geom_jbrowse_reads())
  expect_s3_class(default$plot$scales$get_scales("fill"),
    "ScaleDiscreteIdentity")
  expect_setequal(default$data[[1]]$fill, unique(as.character(bam$reads$color)))
  mapped <- ggplot_build(ggplot(bam) + geom_jbrowse_reads(aes(fill = mapq)))
  expect_s3_class(mapped$plot$scales$get_scales("fill"), "ScaleContinuous")
  expect_gt(length(unique(mapped$data[[1]]$fill)), 1)
})

test_that("a geom reads its table from the plot's result or from its own data", {
  from_plot <- layer_data(ggplot(bam) + geom_jbrowse_mismatches(), 1)
  from_data <- layer_data(ggplot() + geom_jbrowse_mismatches(data = bam), 1)
  expect_equal(nrow(from_plot), nrow(bam$mismatches))
  expect_equal(from_plot$x, from_data$x)
  expect_error(ggplot_build(ggplot(bam$reads) + geom_jbrowse_reads()),
    "needs ggplot")
})

test_that("the display's keys are scales, each listing what was painted", {
  p <- ggplot(bam) + geom_jbrowse_mismatches() + scale_colour_jbrowse(bam, "bases")
  scale <- ggplot_build(p)$plot$scales$get_scales("colour")
  expect_equal(scale$name, "Base")
  expect_error(scale_fill_jbrowse(bam, "nope"), "keys reads, bases")
  s <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA", 1000, 3000,
    color = "strand")
  grDevices::pdf(NULL)
  key <- ggplotGrob(ggplot(s) + geom_jbrowse_reads() + scale_fill_jbrowse(s))
  grDevices::dev.off()
  expect_true(any(grepl("guide-box", key$layout$name)))
})

test_that("a feature track's boxes take any feature column", {
  g <- jbrowse_features(extdata("volvox.sort.gff3.gz"), "ctgA", 0, 25000)
  built <- ggplot_build(ggplot(g) + geom_jbrowse_glyphs(aes(fill = type)))
  boxes <- built$data[[3]]
  expect_equal(nrow(boxes), nrow(g$rects))
  expect_equal(length(unique(boxes$fill)), length(unique(g$rects$type)))
})

test_that("a multi-sample VCF reads as the display paints its genotypes", {
  v <- jbrowse_variants(extdata("dog10k_cyp1a2_snvs.vcf.gz"), "chr30",
    38258000, 38265000)
  expect_equal(nrow(v$variants), 490)
  expect_equal(nrow(v$rows), 39)
  expect_equal(levels(v$cells$row), v$rows$name)
  expect_false("0/0" %in% v$cells$genotype)
  expect_equal(v$cells$start, v$variants$start[v$cells$variant])
  hom <- unique(v$cells$color[v$cells$genotype == "1/1"])
  expect_equal(as.character(hom), as.character(v$variants$color[1]))
  for (layout in c("genomic", "columns")) {
    built <- ggplot_build(autoplot(v, layout = layout))
    expect_equal(nrow(built$data[[1]]), nrow(v$cells))
  }
  haplotypes <- jbrowse_variants(extdata("dog10k_cyp1a2_snvs.vcf.gz"),
    "chr30", 38258000, 38265000, display = list(unit = "haplotype"))
  expect_equal(nrow(haplotypes$rows), 78)
})

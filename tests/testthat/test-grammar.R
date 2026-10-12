extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
bam <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA", 1000, 3000)
library(ggplot2)

test_that("autoplot refuses an argument it does not take", {
  genes <- jbrowse_features(extdata("volvox.sort.gff3.gz"), "ctgA:1,001-3,000")
  expect_error(autoplot(genes, tracks = bam), "tracks")
  expect_error(autoplot(bam, trakcs = genes), "trakcs")
  expect_error(autoplot(jbrowse_wiggle(extdata("volvox_microarray_posneg.bw"),
    "ctgA:1,001-3,000"), layout = "columns"), "layout")
})

test_that("a display color column takes an identity scale, any other column its own", {
  default <- ggplot_build(ggplot(bam) + geom_jbrowse_reads())
  expect_s3_class(default$plot$scales$get_scales("fill"),
    "ScaleDiscreteIdentity")
  expect_setequal(default$data[[1]]$fill, unique(as.character(bam$reads$color)))
  mapped <- ggplot_build(ggplot(bam) + geom_jbrowse_reads(aes(fill = mapq)))
  expect_s3_class(mapped$plot$scales$get_scales("fill"), "ScaleContinuous")
  expect_gt(length(unique(mapped$data[[1]]$fill)), 1)
  recolored <- bam
  recolored$reads$color <- ifelse(bam$reads$strand == 1, "#ff0000", "navy")
  expect_false(inherits(recolored$reads$color, "jbrowse_color"))
  b <- ggplot_build(ggplot(recolored) + geom_jbrowse_reads())
  expect_s3_class(b$plot$scales$get_scales("fill"), "ScaleDiscreteIdentity")
  expect_setequal(b$data[[1]]$fill, c("#ff0000", "navy"))
  labelled <- bam
  labelled$reads$color <- ifelse(bam$reads$strand == 1, "fwd", "rev")
  b <- ggplot_build(ggplot(labelled) + geom_jbrowse_reads())
  expect_false(inherits(b$plot$scales$get_scales("fill"),
    "ScaleDiscreteIdentity"))
})

test_that("a geom reads its table from the plot's result or from its own data", {
  from_plot <- layer_data(ggplot(bam) + geom_jbrowse_mismatches(), 1)
  from_data <- layer_data(ggplot() + geom_jbrowse_mismatches(data = bam), 1)
  expect_equal(nrow(from_plot), nrow(bam$mismatches))
  expect_equal(from_plot$x, from_data$x)
  expect_error(ggplot_build(ggplot(bam$reads) + geom_jbrowse_reads()),
    "needs ggplot")
})

test_that("a data frame as a geom's data stands in for the table it draws", {
  s <- jbrowse_alignments(extdata("spliced.bam"), "ctgA", 400, 1300,
    reference = extdata("volvox.fa"))
  blocks <- ggjbrowse:::read_blocks(s)
  forward <- s$reads[s$reads$strand == 1, ]
  drawn <- layer_data(ggplot(s) + geom_jbrowse_reads(data = forward))
  expect_equal(nrow(drawn), sum(blocks$read %in% forward$read))
  skips <- layer_data(ggplot(s) + geom_jbrowse_skips(data = s$gaps))
  expect_equal(nrow(skips), sum(s$gaps$type == "skip"))
  common <- s$mismatches[s$mismatches$frequency > 0.5, ]
  expect_equal(nrow(layer_data(ggplot(s) +
    geom_jbrowse_mismatches(data = common))), nrow(common))
  expect_equal(nrow(layer_data(ggplot(s) +
    geom_jbrowse_junctions(data = s$junctions[1, ]))), 48)
  expect_error(geom_jbrowse_letters(data = s$mismatches), "several")
  expect_error(ggplot_build(ggplot() + geom_jbrowse_reads(data = forward)),
    "replaces one of the plot's result's tables")
})

test_that("the display's keys are scales, each listing what was painted", {
  p <- ggplot(bam) + geom_jbrowse_mismatches() + scale_colour_jbrowse(bam, "bases")
  scale <- ggplot_build(p)$plot$scales$get_scales("colour")
  expect_equal(scale$name, "Base")
  expect_error(scale_fill_jbrowse(bam, "nope"), "keys reads, bases")
  s <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA", 1000, 3000,
    display = list(color = list(field = "strand",
      range = c("#FF0000", "#0000FF"))))
  expect_setequal(as.character(unique(s$reads$color)), c("#ff0000", "#0000ff"))
  key <- get_guide_data(ggplot(s) + geom_jbrowse_reads() +
    scale_fill_jbrowse(s), "fill")
  expect_equal(nrow(key), 2)
  expect_setequal(key$fill, c("#ff0000", "#0000ff"))
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
    expect_equal(as.character(built$data[[1]]$fill),
      as.character(v$legend$color[1]))
    expect_equal(nrow(built$data[[2]]), nrow(v$cells))
  }
  p <- autoplot(v, layout = "columns")
  expect_equal(ggplot_build(p)$data[[3]]$label, v$rows$name)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  find_labels <- function(g) {
    if (inherits(g, "ggjbrowse_row_labels")) return(g)
    for (child in c(g$grobs, g$children)) {
      found <- find_labels(child)
      if (!is.null(found)) return(found)
    }
    NULL
  }
  grob <- find_labels(ggplotGrob(p))
  drawn_pt <- function(row_pt) {
    grid::pushViewport(grid::viewport(height = grid::unit(39 * row_pt, "pt")))
    on.exit(grid::popViewport())
    drawn <- grid::makeContent(grob)$children
    if (length(drawn)) drawn[[1]]$gp$fontsize else 0
  }
  expect_equal(drawn_pt(20), 8)
  expect_equal(drawn_pt(5), 0.85 * 5)
  expect_equal(drawn_pt(2), 0)
  haplotypes <- jbrowse_variants(extdata("dog10k_cyp1a2_snvs.vcf.gz"),
    "chr30", 38258000, 38265000, display = list(unit = "haplotype"))
  expect_equal(nrow(haplotypes$rows), 78)
})

test_that("a key row names every label its color stands for", {
  paired <- jbrowse_alignments(extdata("volvox-simple-inv-paired.bam"), "ctgA",
    18000, 32000, display = list(readConnections = "arc"))
  legend <- paired$arc_legend
  expect_true(anyDuplicated(legend$color) > 0)
  key <- get_guide_data(ggplot(paired) + geom_jbrowse_arcs() +
    scale_colour_jbrowse(paired, "arcs"), "colour")
  expect_equal(nrow(key), length(unique(legend$color)))
  for (label in legend$label) {
    expect_true(any(grepl(label, key$.label, fixed = TRUE)))
  }
  expect_equal(as.character(jbrowse_color("#FF00AA")), "#ff00aa")
})

test_that("the multi-layer geoms take any column of their tables", {
  paired <- jbrowse_alignments(extdata("volvox-simple-inv-paired.bam"), "ctgA",
    18000, 32000, display = list(readConnections = "arc"))
  domes <- layer_data(ggplot(paired) + geom_jbrowse_arcs(aes(alpha = span)), 1)
  expect_equal(length(unique(domes$group)), nrow(paired$arcs))
  expect_gt(length(unique(domes$alpha)), 1)
  g <- jbrowse_features(extdata("volvox.sort.gff3.gz"), "ctgA", 0, 25000)
  built <- ggplot_build(ggplot(g) + geom_jbrowse_glyphs(aes(alpha = end - start)))
  expect_equal(nrow(built$data[[2]]), 2 * nrow(g$chevrons))
  labelled <- jbrowse_features(extdata("volvox.sort.gff3.gz"), "ctgA", 0, 25000,
    width = 5000, display = list(subfeatureLabels = "below"))
  expect_true(all(labelled$labels$feature >= 1))
  sub <- labelled$labels[!is.na(labelled$labels$subfeature), ]
  expect_gt(nrow(sub), 0)
  expect_equal(labelled$subfeatures$name[sub$subfeature], sub$text)
  expect_equal(labelled$subfeatures$feature[sub$subfeature], sub$feature)
})

test_that("the junction band keys the strands its arcs are drawn in", {
  spliced <- jbrowse_alignments(extdata("spliced.bam"), "ctgA", 400, 1300)
  key <- get_guide_data(autoplot(spliced, "junctions", min_junction_reads = 1),
    "colour")
  expect_equal(key$.label, "Forward")
})

test_that("rows take their samples file's colors under rowColor, keyed by the attribute", {
  v <- jbrowse_variants(extdata("dog10k_amy2b_svs.vcf.gz"),
    "chr6:47,340,001-47,420,000",
    samples = extdata("dog10k_amy2b_samples.tsv"),
    display = list(rowColor = list(field = "group")))
  expect_equal(attr(v$keys$rows, "title"), "Group")
  expect_equal(v$keys$rows$label, c("Breed dog", "Village dog", "Gray wolf"))
  expect_setequal(as.character(unique(v$rows$color)),
    as.character(v$keys$rows$color))
  built <- ggplot_build(autoplot(v))
  bars <- built$data[[3]]
  expect_equal(nrow(bars), nrow(v$rows))
  plain <- jbrowse_variants(extdata("dog10k_amy2b_svs.vcf.gz"),
    "chr6:47,340,001-47,420,000")
  expect_true(all(plain$rows$color == ""))
  expect_false(any(vapply(autoplot(plain)$layers, function(l)
    inherits(l$geom, "GeomSegment"), NA)))
})

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

test_that("an identity key lists every domain color, drawn or not, as the browser does", {
  keyed <- jbrowse_features(gff, "ctgA", 0, 25000, display = list(color = list(
    scale = "identity", domain = c("goldenrod", "purple"),
    labels = c("Drawn", "Not drawn"))))
  expect_false("#800080" %in% keyed$rects$color)
  expect_equal(keyed$legend$label, c("Drawn", "Not drawn"))
  guide <- ggplot2::get_guide_data(ggplot2::autoplot(keyed), "fill")
  expect_equal(guide$.label, c("Drawn", "Not drawn"))
})

test_that("a fitted track gives up what the browser gives up to fit it", {
  fit <- jbrowse_features(gff, "ctgA", 0, 25000,
    display = list(heightMode = "fit", height = 100))
  eden <- which(fit$features$name == "EDEN")
  expect_lte(max(fit$features$bottom), 100)
  expect_false("EDEN" %in% fit$labels$text)
  expect_equal(sum(fit$subfeatures$feature == eden), 1)
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

test_that("a VCF reads as a variant track's feature display draws it", {
  v <- jbrowse_features(extdata("dog10k_cyp1a2_snvs.vcf.gz"),
    "chr30:38,258,001-38,265,000")
  expect_gt(nrow(v$features), 100)
  expect_gt(nrow(v$rects), 0)
})

test_that("a feature narrower than a pixel draws a display pixel wide", {
  v <- jbrowse_features(extdata("dog10k_cyp1a2_snvs.vcf.gz"),
    "chr30:38,258,001-38,265,000", width = 500)
  expect_true(all(v$rects$end - v$rects$start < v$bp_per_px))
  boxes <- ggplot2::layer_data(ggplot2::ggplot(v) + geom_jbrowse_glyphs(), 3)
  expect_equal(nrow(boxes), nrow(v$rects))
  expect_true(all(boxes$xmax - boxes$xmin >= v$bp_per_px - 1e-9))
})

test_that("a gene's exons with the introns collapsed keep each feature on one row", {
  exons <- c("ctgA:1,001-1,600", "ctgA:2,901-4,000", "ctgA:4,901-5,600",
    "ctgA:6,901-9,100")
  g <- jbrowse_features(gff, exons)
  expect_setequal(unique(g$features$.region), 1:4)
  tops <- tapply(g$features$top, g$features$id, function(t) length(unique(t)))
  expect_true(all(tops == 1))
  eden <- g$features[g$features$name == "EDEN", ]
  expect_equal(eden$.region, 1:4)
  for (name in c("rects", "lines", "labels", "subfeatures")) {
    d <- g[[name]]
    expect_equal(d$.region, g$features$.region[d$feature])
  }
  subs <- g$labels[!is.na(g$labels$subfeature), ]
  expect_equal(g$subfeatures$feature[subs$subfeature], subs$feature)
  apart <- lapply(exons, function(l) jbrowse_features(gff, l))
  expect_gt(length(unique(vapply(apart, function(a)
    a$features$top[a$features$name == "EDEN"], 0))), 1)
})

test_that("a figure at jbrowse_width() draws labels the same size at any view width", {
  gff <- extdata("gapdh_gencode.gff3.gz")
  label_pt <- function(px) {
    fig <- jbrowse_stack(jbrowse_features(gff, "chr12:6,534,201-6,538,700",
      width = px))
    grDevices::pdf(NULL, width = jbrowse_width(fig),
      height = jbrowse_height(fig))
    on.exit(grDevices::dev.off())
    grid::grid.draw(patchwork::patchworkGrob(fig))
    grid::grid.force()
    first_label(grid::grid.grab())
  }
  first_label <- function(g) {
    if (inherits(g, "ggjbrowse_display_text")) {
      return(grid::makeContent(g)$children[[1]]$gp$fontsize[1])
    }
    for (child in if (inherits(g, "gTree")) g$children) {
      pt <- first_label(child)
      if (!is.null(pt)) return(pt)
    }
  }
  expect_equal(label_pt(2500), label_pt(800), tolerance = 0.05)
  expect_equal(label_pt(800), 11 / 150 * 72.27, tolerance = 0.05)
})

test_that("the display's facet packs a section per value, each on its own rows", {
  s <- jbrowse_features(gff, "ctgA", 0, 25000,
    display = list(facet = "strand"))
  expect_gt(nrow(s$sections), 1)
  expect_equal(nrow(s$features), nrow(g$features))
  expect_equal(levels(s$features$section), s$sections$label)
  chip <- s$glyph$sectionLabelHeight
  expect_gt(chip, 0)
  for (i in seq_len(nrow(s$sections))) {
    f <- s$features[s$features$section == s$sections$label[i], ]
    expect_true(all(f$top >= s$sections$top[i] + chip))
    expect_true(all(f$bottom <= s$sections$top[i] + s$sections$height[i]))
    expect_length(unique(f$strand), 1)
  }
  for (table in c("rects", "lines", "chevrons", "arrows", "labels")) {
    d <- s[[table]]
    expect_equal(d$section, s$features$section[d$feature])
  }
  keys <- s$sections$key
  hidden <- jbrowse_features(gff, "ctgA", 0, 25000,
    display = list(facet = list(field = "strand", hidden = keys[1])))
  expect_false(keys[1] %in% hidden$sections$key)
  expect_lt(nrow(hidden$features), nrow(s$features))
  reordered <- jbrowse_features(gff, "ctgA", 0, 25000,
    display = list(facet = list(field = "strand", domain = rev(keys))))
  expect_equal(reordered$sections$key, rev(keys))
  expect_equal(g$sections$label, "")
})

test_that("autoplot draws each section in its own panel over its rows", {
  s <- jbrowse_features(gff, "ctgA", 0, 25000,
    display = list(facet = "strand"))
  built <- ggplot2::ggplot_build(ggplot2::autoplot(s))
  expect_equal(nrow(built$layout$layout), nrow(s$sections))
  chip <- s$glyph$sectionLabelHeight
  ranges <- lapply(built$layout$panel_params, function(q)
    sort(-q$y$continuous_range))
  expect_equal(ranges, lapply(seq_len(nrow(s$sections)), function(i)
    c(s$sections$top[i] + chip, s$sections$top[i] + s$sections$height[i])))
  labels <- built$data[[length(built$data)]]
  expect_equal(labels$px_tall,
    (s$sections$height - chip)[as.integer(s$labels$section)])
})

test_that("a zoomed-in CDS carries its residues, in frame across an intron", {
  gapdh <- extdata("gapdh_gencode.gff3.gz")
  loci <- c("chr12:6,534,811-6,534,880", "chr12:6,536,476-6,536,545")
  reference <- extdata("gapdh_hg38.fa.gz")
  near <- jbrowse_features(gapdh, loci, reference = reference, width = 1200)
  aa <- near$amino_acids
  first <- aa[aa$.region == 1 & aa$top == aa$top[1], ]
  expect_equal(paste(first$amino_acid, collapse = ""), "MGKVKVGVNG")
  expect_equal(first$protein_index, 0:9)
  expect_equal(first$text, first$amino_acid)
  expect_equal(first$color == "#ff0000", first$stop_or_non_triplet)
  expect_equal(first$end[10] - first$start[10], 2)
  expect_true(all(aa$y > aa$top & aa$y < aa$top + aa$height))
  expect_true(any(near$rects$color != near$rects$color[1]))

  unlettered <- jbrowse_features(gapdh, loci, reference = reference,
    width = 400)
  expect_equal(nrow(unlettered$amino_acids), nrow(aa))
  expect_true(all(unlettered$amino_acids$text == ""))
  expect_equal(nrow(jbrowse_features(gapdh, loci, width = 1200)$amino_acids),
    0)
  expect_equal(nrow(jbrowse_features(gapdh, loci, reference = reference,
    width = 1200, amino_acids = FALSE)$amino_acids), 0)
})

test_that("the residues draw as haloed monospace letters", {
  near <- jbrowse_features(extdata("gapdh_gencode.gff3.gz"),
    "chr12:6,534,811-6,534,880", reference = extdata("gapdh_hg38.fa.gz"),
    width = 1200)
  p <- ggplot2::ggplot(near) + geom_jbrowse_amino_acids()
  d <- ggplot2::layer_data(p)
  expect_equal(nrow(d), sum(nzchar(near$amino_acids$text)))
  expect_equal(unique(d$family), "mono")
  text <- ggplot2::layer_grob(p)[[1]]
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  drawn <- grid::makeContent(text)$children
  expect_length(drawn, 9)
  expect_equal(unique(drawn[[1]]$gp$col), near$glyph$aminoAcidHalo)
  expect_setequal(tolower(drawn[[9]]$gp$col),
    unique(near$amino_acids$color))
})

test_that("a track's name in jbrowse_stack titles its first band", {
  bam <- system.file("extdata", "volvox-sorted.bam", package = "ggjbrowse")
  gff <- system.file("extdata", "volvox.sort.gff3.gz", package = "ggjbrowse")
  genes <- jbrowse_features(gff, "ctgA", 1000, 3000)
  reads <- jbrowse_alignments(bam, "ctgA", 1000, 3000)
  titles <- function(fig) vapply(c(fig$patches$plots, list(fig)),
    function(p) p$labels$title %||% "", "")
  expect_equal(titles(jbrowse_stack(genes, reads)), rep("", 3))
  expect_equal(titles(jbrowse_stack(Genes = genes, reads)),
    c("Genes", "", ""))
  expect_equal(titles(jbrowse_stack(genes, list(a = reads, b = reads))),
    c("", "a", "", "b", ""))
})

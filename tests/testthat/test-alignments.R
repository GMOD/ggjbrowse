extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
bam <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA", 1000, 6000)

test_that("a BAM region reads as JBrowse lays it out", {
  expect_equal(nrow(bam$reads), 982)
  for (r in split(bam$reads, bam$reads$row)) {
    r <- r[order(r$start), ]
    expect_true(all(utils::tail(r$start, -1) >= utils::head(r$end, -1)))
  }
  m <- bam$mismatches
  expect_true(all(m$position >= bam$reads$start[m$read] &
    m$position < bam$reads$end[m$read]))
  expect_equal(m$row, bam$reads$row[m$read])
  expect_setequal(unique(m$base), c("A", "C", "G", "T"))
  expect_true(all(m$frequency >= 0 & m$frequency <= 1))
  expect_gt(max(bam$coverage$depth), 20)
})

test_that("a CRAM reads as the same tables as its BAM", {
  cram <- jbrowse_alignments(extdata("volvox-sorted.cram"), "ctgA", 1000, 6000,
    reference = extdata("volvox.fa"))
  key <- function(x) sort(paste(x$reads$name, x$reads$start, x$reads$end,
    x$reads$row))
  expect_equal(key(cram), key(bam))
  mm <- function(x) sort(paste(x$mismatches$position, x$mismatches$base,
    x$reads$name[x$mismatches$read]))
  expect_equal(mm(cram), mm(bam))
})

test_that("insertions and clips sit on their reads, in the track's layers", {
  l <- jbrowse_alignments(extdata("volvox-long-reads-sv.bam"), "ctgA", 1000, 20000)
  expect_gt(nrow(l$insertions), 0)
  expect_true(all(l$insertions$position >= l$reads$start[l$insertions$read] &
    l$insertions$position <= l$reads$end[l$insertions$read]))
  expect_equal(nchar(l$insertions$sequence), l$insertions$length)
  expect_setequal(unique(l$clips$type), c("soft", "hard"))
  expect_equal(l$clips$row, l$reads$row[l$clips$read])
  p <- ggplot2::ggplot_build(ggplot2::autoplot(l, "pileup"))
  layers <- vapply(p$data, nrow, 1L)
  expect_equal(layers[8], nrow(l$insertions))
  expect_equal(layers[9], nrow(l$clips))
  # sequencing-error insertions on 20 long reads fade toward the floor
  faded <- p$data[[7]]$alpha
  expect_true(any(faded < 1))
  expect_true(all(faded >= (1200 / 19000)^2 - 1e-9))
  opaque <- ggplot2::ggplot_build(ggplot2::autoplot(l, "pileup",
    fade_low_frequency = FALSE))$data[[7]]$alpha
  expect_true(all(opaque == 1))
})

test_that("the coverage band stacks each position's alleles as the display does", {
  a <- bam$allele_coverage
  expect_gt(nrow(a), 0)
  expect_true(all(a$position %in% bam$mismatches$position))
  expect_true(all(a$from >= 0 & a$to > a$from & a$to <= 1 + 1e-6))
  cov <- ggplot2::ggplot_build(ggplot2::autoplot(bam, "coverage"))
  alleles <- cov$data[[2]]
  expect_equal(nrow(alleles), nrow(a))
  expect_true(all(alleles$y >= 0 & alleles$yend > alleles$y))
  depth <- bam$coverage$depth[findInterval(alleles$x, bam$coverage$start)]
  expect_true(all(alleles$yend <= depth + 1e-6))
  strict <- ggplot2::ggplot_build(ggplot2::autoplot(bam, "coverage",
    min_allele_frequency = 0.5))$data[[2]]
  expect_lt(nrow(strict), nrow(alleles))
  expect_true(all((strict$yend - strict$y) /
    bam$coverage$depth[findInterval(strict$x, bam$coverage$start)] >= 0.5))
})

sv <- jbrowse_alignments(extdata("volvox-long-reads-sv.bam"), "ctgA", 1000,
  20000)

test_that("the interbase tables count the reads' insertions and clips", {
  ib <- sv$interbase_coverage
  expect_setequal(unique(ib$type), c("insertion", "softclip", "hardclip"))
  events <- c(paste(sv$insertions$position, "insertion"),
    paste(sv$clips$position, ifelse(sv$clips$type == "soft", "softclip",
      "hardclip")))
  expect_equal(as.vector(table(events)[paste(ib$position, ib$type)]),
    ib$to - ib$from)
  expect_equal(sv$indicators$type, "softclip")
  expect_equal(sv$keys$interbase$label,
    c("Insertion", "Soft clip", "Hard clip"))
})

test_that("the coverage band draws on the display's depth axis, inset as the display insets it", {
  scale <- sv$coverage_scale
  built <- ggplot2::ggplot_build(ggplot2::autoplot(sv, "coverage"))
  y <- built$layout$panel_params[[1]]$y
  expect_equal(y$breaks[!is.na(y$breaks)], scale$ticks)
  inset <- sv$glyph$coverageInset / (sv$settings$coverageHeight -
    2 * sv$glyph$coverageInset)
  span <- diff(scale$domain)
  expect_equal(y$continuous_range,
    scale$domain + c(-1, 1) * inset * span)

  hidden <- jbrowse_alignments(extdata("volvox-long-reads-sv.bam"), "ctgA",
    1000, 20000, display = list(showInterbaseIndicators = FALSE))
  expect_length(ggplot2::ggplot_build(ggplot2::autoplot(hidden,
    "coverage"))$data, 3)
})

# The panel the band draws in, 1 inch tall, and where the grob's `y` lands in
# it, in points from the bottom
grob_y_pt <- function(grob, y) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  grid::pushViewport(grid::viewport(height = grid::unit(1, "in")))
  grid::convertY(y, "pt", valueOnly = TRUE)
}

test_that("a stack hangs from the triangles' tip, as deep a share of the depth area as the domain's top is twice", {
  b <- ggplot2::ggplot_build(ggplot2::ggplot(sv) + geom_jbrowse_interbase() +
    geom_jbrowse_coverage())
  d <- b$data[[1]]
  ib <- sv$interbase_coverage
  top <- 2 * sv$coverage_scale$domain[2]
  expect_equal(d$hang, ib$from / top)
  expect_equal(d$hang_end, ib$to / top)
  grob <- ggplot2::layer_grob(ggplot2::ggplot(sv) + geom_jbrowse_interbase() +
    geom_jbrowse_coverage(), 1)[[1]]
  # in a 1 inch panel: 72.27 points, a display px 1/150 of an inch
  inch <- 72.27
  pt <- function(px) px / 150 * inch
  deepest <- which.max(ib$to)
  expect_equal(grob_y_pt(grob, grob$y)[deepest],
    inch - pt(4) - ib$from[deepest] / top * (inch - 2 * pt(5)))
  expect_equal(grob_y_pt(grob, grob$height)[deepest],
    (ib$to - ib$from)[deepest] / top * (inch - 2 * pt(5)))
})

test_that("a log depth axis leaves the hang linear", {
  log <- jbrowse_alignments(extdata("volvox-long-reads-sv.bam"), "ctgA", 1000,
    20000, display = list(scales = list(y = list(type = "log"))))
  expect_equal(log$coverage_scale$domain[1], 1)
  built <- ggplot2::ggplot_build(ggplot2::autoplot(log, "coverage"))
  hang <- built$data[[4]]
  expect_equal(hang$hang, log$interbase_coverage$from /
    (2 * log$coverage_scale$domain[2]))
  # an allele keeps its share of the bar's drawn height from the bar's base,
  # one read: on a log axis its top is at depth^to
  a <- ggplot2::layer_data(ggplot2::ggplot(log) + geom_jbrowse_alleles())
  alleles <- log$allele_coverage
  expect_gt(nrow(alleles), 0)
  expect_equal(a$ymax, alleles$depth^alleles$to)
  expect_equal(a$ymin, alleles$depth^alleles$from)
})

test_that("every section shares the display's one depth axis", {
  s <- jbrowse_alignments(extdata("volvox-long-reads-sv.bam"), "ctgA", 1000,
    20000, display = list(facet = "strand"))
  expect_gt(length(unique(s$interbase_coverage$section)), 1)
  built <- ggplot2::ggplot_build(ggplot2::autoplot(s, "coverage"))
  ranges <- lapply(built$layout$panel_params, function(p) p$y$continuous_range)
  expect_true(all(vapply(ranges, identical, TRUE, ranges[[1]])))
  expect_gte(s$coverage_scale$domain[2], max(s$coverage$depth))
})

test_that("several regions' indicators land on the axis that lays them end to end", {
  m <- jbrowse_alignments(extdata("volvox-long-reads-sv.bam"),
    c("ctgA:2600-2800", "ctgA:2690-2710"))
  expect_equal(m$indicators$position, c(2700, 2700))
  drawn <- ggplot2::layer_data(ggplot2::ggplot(m) + geom_jbrowse_indicators())
  expect_equal(sort(drawn$x),
    sort(2700 - m$regions$start + m$regions$offset))
})

test_that("the display's sort object lays the rows out", {
  site <- bam$mismatches$position[1]
  s <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA", 1000, 6000,
    display = list(sort = list(type = "basePair", pos = site)))
  covering <- which(s$reads$start <= site & s$reads$end > site)
  mismatched <- unique(s$mismatches$read[s$mismatches$position == site])
  expect_gt(length(mismatched), 0)
  expect_lt(max(s$reads$row[mismatched]),
    min(s$reads$row[setdiff(covering, mismatched)]))
  l <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA", 1000, 6000,
    display = list(sort = list(type = "length")))
  expect_false(identical(l$reads$row, bam$reads$row))
})

test_that("splice junctions arc over the coverage, counted and stranded", {
  s <- jbrowse_alignments(extdata("spliced.bam"), "ctgA", 400, 1300,
    reference = extdata("volvox.fa"))
  expect_equal(nrow(s$junctions), 2)
  expect_equal(sum(s$junctions$count), sum(s$gaps$type == "skip"))
  expect_true(all(s$junctions$strand %in% c(1, -1, 0)))
  arcs <- ggplot2::ggplot_build(ggplot2::autoplot(s, "junctions",
    min_reads = 1))
  expect_equal(length(unique(arcs$data[[1]]$group)), 2)
  expect_equal(nrow(arcs$data[[2]]), 2)
  track <- ggplot2::autoplot(s)
  expect_s3_class(track, "patchwork")
  expect_s3_class(ggplot2::autoplot(bam), "patchwork")
})

test_that("base modifications paint under the layer, keyed and tinted as the display paints them", {
  bam <- extdata("methylation_clip.bam")
  plain <- jbrowse_alignments(bam, "20", 10000, 20000)
  expect_equal(nrow(plain$modifications), 0)
  m <- jbrowse_alignments(bam, "20", 10000, 20000,
    reference = extdata("hg38_clip.fa"),
    display = list(baseColor = list(field = "modifications")))
  expect_gt(nrow(m$modifications), 100)
  expect_setequal(unique(m$modifications$name), "5mC")
  expect_equal(m$modifications$row, m$reads$row[m$modifications$read])
  expect_true(all(m$modifications$opacity >= 0.1 & m$modifications$opacity <= 1))
  expect_true("5mC" %in% m$legend$label)
  expect_setequal(unique(m$reads$category), c("modFwd", "modRev"))
  expect_gt(nrow(m$modification_coverage), 0)
  band <- ggplot2::ggplot_build(ggplot2::autoplot(m, "coverage"))$data[[3]]
  expect_equal(nrow(band), nrow(m$modification_coverage))
  p <- ggplot2::ggplot_build(ggplot2::autoplot(m, "pileup"))
  expect_equal(nrow(p$data[[6]]), nrow(m$modifications))
  expect_true(all(p$data[[6]]$fill %in% m$legend$color))
  calls <- jbrowse_alignments(bam, "20", 10000, 20000,
    display = list(baseColor = list(field = "modifications"),
      modifications = list(unmodified = "calls")))
  expect_gt(sum(calls$modifications$unmodified), 0)
  every <- jbrowse_alignments(bam, "20", 10000, 20000,
    display = list(baseColor = list(field = "modifications"),
      modifications = list(unmodified = "all")))
  expect_gt(sum(every$modifications$unmodified),
    10 * sum(calls$modifications$unmodified))
})

test_that("modification coloring greys the mismatched bases and drops their key", {
  bam <- extdata("methylation_clip.bam")
  fa <- extdata("hg38_clip.fa")
  plain <- jbrowse_alignments(bam, "20", 10000, 12000, reference = fa)
  expect_gt(length(unique(plain$mismatches$color)), 1)
  expect_equal(nrow(plain$keys$bases), 5)
  m <- jbrowse_alignments(bam, "20", 10000, 12000, reference = fa,
    display = list(baseColor = list(field = "modifications")))
  expect_gt(nrow(m$mismatches), 0)
  expect_setequal(as.character(m$mismatches$color), "#888")
  expect_setequal(as.character(m$allele_coverage$color), "#888")
  expect_equal(nrow(m$keys$bases), 0)
})

test_that("bisulfite reads paint their methylated cytosines by context", {
  bam <- extdata("arabidopsis_wgbs.bam")
  fa <- extdata("arabidopsis_tair10.fa.gz")
  locus <- "NC_003070.9:4,406,001-4,408,000"
  bisulfite <- function(modifications = NULL) jbrowse_alignments(bam, locus,
    reference = fa, display = list(baseColor = list(field = "bisulfite"),
      modifications = modifications))
  cg <- bisulfite()
  expect_gt(nrow(cg$modifications), 100)
  expect_false(any(cg$modifications$unmodified))
  expect_setequal(cg$modifications$probability, 1)
  expect_equal(cg$modifications$row, cg$reads$row[cg$modifications$read])
  expect_true("5mC methylated" %in% cg$legend$label)
  expect_gt(nrow(cg$modification_coverage), 0)
  expect_setequal(as.character(cg$mismatches$color), "#888")
  expect_equal(nrow(cg$keys$bases), 0)
  seq <- jbrowse_sequence(fa, "NC_003070.9", 4406000, 4408000)
  base_at <- function(p, offset) toupper(substring(seq, p - 4406000 + 1 +
    offset, p - 4406000 + 1 + offset))
  forward <- cg$modifications[base_at(cg$modifications$position, 0) == "C", ]
  expect_gt(nrow(forward), 0)
  expect_setequal(base_at(forward$position, 1), "G")
  chh <- bisulfite(list(cytosineContext = "CHH"))
  expect_false(any(chh$modifications$position %in% cg$modifications$position))
  calls <- bisulfite(list(unmodified = "calls"))
  expect_gt(sum(calls$modifications$unmodified), 0)
  expect_equal(sum(!calls$modifications$unmodified), nrow(cg$modifications))
  expect_true("Unmethylated" %in% calls$legend$label)
  p <- ggplot2::ggplot_build(ggplot2::autoplot(cg, "pileup"))
  expect_true(any(vapply(p$data, nrow, 0) == nrow(cg$modifications)))
  expect_error(jbrowse_alignments(bam, locus,
    display = list(baseColor = list(field = "bisulfite"))), "reference")
})

test_that("base quality paints each aligned base on its read, keyed as a bar", {
  bam <- extdata("methylation_clip.bam")
  q <- jbrowse_alignments(bam, "20", 10000, 10300,
    display = list(baseColor = list(field = "baseQuality")))
  expect_gt(nrow(q$qualities), 1000)
  expect_equal(q$qualities$row, q$reads$row[q$qualities$read])
  colors <- tapply(as.character(q$qualities$color), q$qualities$score,
    function(x) length(unique(x)))
  expect_true(all(colors == 1))
  expect_equal(vapply(q$ramps, attr, "", "title"), "Base quality")
  expect_equal(attr(q$ramps[[1]], "labels")[1], "0")
  p <- ggplot2::ggplot_build(ggplot2::autoplot(q, "pileup"))
  expect_equal(nrow(p$data[[3]]), nrow(q$qualities))
  expect_equal(nrow(jbrowse_alignments(bam, "20", 10000, 10300)$qualities), 0)
})

test_that("a bgzipped FASTA reads as a reference, and MD tags stand in for one", {
  fa <- extdata("gapdh_hg38.fa.gz")
  expect_equal(jbrowse_sequence(fa, "chr12", 6499990, 6500000),
    strrep("N", 10))
  expect_match(jbrowse_sequence(fa, "chr12", 6534500, 6534600),
    "^[ACGTacgt]{100}$")
  bam <- extdata("gapdh_direct_rna.bam")
  locus <- "chr12:6,534,201-6,538,700"
  tagged <- jbrowse_alignments(bam, locus)
  expect_gt(nrow(tagged$mismatches), 1000)
  expect_true(all(tagged$junctions$motif == ""))
  read <- jbrowse_alignments(bam, locus, reference = fa)
  expect_equal(read$mismatches, tagged$mismatches)
  common <- read$junctions[read$junctions$count >= 10, ]
  expect_gt(nrow(common), 5)
  expect_true(all(common$motif == "GT-AG (canonical)"))
})

test_that("the reference reads as its bases, drawn as a band only when asked", {
  fa <- extdata("volvox.fa")
  seq <- jbrowse_sequence(fa, "ctgA", 1000, 1020)
  expect_match(seq, "^[ACGTNacgtn]{20}$")
  z <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA", 1000, 1100,
    reference = fa)
  expect_equal(nchar(jbrowse_sequence(fa, "ctgA", 1000, 1100)), 100)
  s <- ggplot2::ggplot_build(ggplot2::autoplot(z, "sequence"))
  expect_equal(nrow(s$data[[2]]), 100)
  track <- ggplot2::autoplot(z)
  expect_s3_class(track, "patchwork")
  expect_equal(length(track$patches$plots) + 1, 2)
})

test_that("a spliced read splits into its blocks at its skips", {
  s <- jbrowse_alignments(extdata("spliced.bam"), "ctgA", 400, 1300)
  skips <- s$gaps[s$gaps$type == "skip", ]
  expect_gt(nrow(skips), 0)
  blocks <- ggjbrowse:::read_blocks(s)
  covered <- tapply(blocks$end - blocks$start, blocks$read, sum)
  skipped <- tapply(skips$end - skips$start, skips$read, sum)
  span <- s$reads$end - s$reads$start
  i <- seq_len(nrow(s$reads))
  gap <- ifelse(is.na(skipped[as.character(i)]), 0, skipped[as.character(i)])
  expect_equal(as.vector(covered[as.character(i)] + gap), span)
  kept <- s
  kept$reads <- s$reads[s$reads$read %% 2 == 0, ]
  expect_equal(ggjbrowse:::read_blocks(kept),
    blocks[blocks$read %% 2 == 0, ], ignore_attr = TRUE)
  drawn <- ggplot2::layer_data(ggplot2::ggplot(kept) + geom_jbrowse_reads())
  expect_equal(nrow(drawn), sum(blocks$read %% 2 == 0))
})

test_that("the pileup and coverage draw in the display's colors", {
  colors <- jbrowse_colors()
  expect_type(colors$bases, "character")
  expect_setequal(names(colors$bases), c("A", "C", "G", "T", "N"))
  expect_equal(unique(bam$reads$color), tolower(colors$read))
  expect_equal(bam$legend$label, "Reads")
  p <- ggplot2::ggplot_build(ggplot2::autoplot(bam))
  expect_equal(unique(p$data[[1]]$fill), tolower(colors$read))
  expect_equal(nrow(p$data[[7]]), nrow(bam$mismatches))
  cov <- ggplot2::ggplot_build(ggplot2::autoplot(bam, "coverage"))
  expect_equal(nrow(cov$data[[1]]), nrow(bam$coverage))
})

test_that("reads fill as the display's color object paints them, with its key", {
  s <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA", 1000, 6000,
    display = list(color = list(field = "strand")))
  expect_equal(attr(s$legend, "title"), "Read colors")
  expect_equal(s$legend$label, c("Forward strand", "Reverse strand"))
  expect_equal(unique(s$reads$color[s$reads$strand == 1]), s$legend$color[1])
  expect_equal(unique(s$reads$color[s$reads$strand == -1]), s$legend$color[2])
  expect_setequal(unique(s$reads$category), c("fwdStrand", "revStrand"))
  p <- ggplot2::ggplot_build(ggplot2::autoplot(s))
  expect_setequal(unique(p$data[[1]]$fill), s$legend$color)

  own <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA", 1000, 6000,
    display = list(color = list(field = "strand",
      range = c("#ff0000", "#0000ff"), title = "Mine")))
  expect_equal(unique(own$reads$color[own$reads$strand == 1]), "#ff0000")
  expect_equal(attr(own$legend, "title"), "Mine")

  m <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA", 1000, 6000,
    display = list(color = list(field = "mapq")))
  expect_equal(attr(m$legend, "title"), "Mapping quality")
  expect_true(all(m$reads$color %in% m$legend$color))
  expect_equal(unique(m$reads$category), "mapq")
})

test_that("the display's filter object keeps the reads it keeps", {
  r <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA", 1000, 6000,
    display = list(filter = list(flagInclude = 16)))
  expect_true(all(r$reads$strand == -1))
  expect_gt(nrow(r$reads), 0)
  expect_lt(nrow(r$reads), nrow(bam$reads))
})

test_that("svSignal only and exclude split the reads between them", {
  locus <- "chr3:25,357,600-25,361,000"
  bam <- extdata("colo829_tumor_chr3.bam")
  read <- function(level) jbrowse_alignments(bam, locus,
    display = list(filter = list(svSignal = level)))
  all <- jbrowse_alignments(bam, locus)
  only <- read("only")
  none <- read("exclude")
  expect_gt(nrow(only$reads), 0)
  expect_gt(nrow(none$reads), 0)
  expect_equal(nrow(only$reads) + nrow(none$reads), nrow(all$reads))
})

test_that("a remote BAM reads over HTTP", {
  skip_on_cran()
  skip_if_offline("jbrowse.org")
  u <- "https://jbrowse.org/code/jb2/main/test_data/volvox/volvox-sorted.bam"
  r <- jbrowse_alignments(u, "ctgA", 1000, 2000)
  expect_gt(nrow(r$reads), 100)
})

test_that("a file overwritten in place is read again", {
  dir <- tempfile()
  dir.create(dir)
  path <- file.path(dir, "x.bam")
  put <- function(f) {
    file.copy(extdata(f), path, overwrite = TRUE)
    file.copy(extdata(paste0(f, ".bai")), paste0(path, ".bai"), overwrite = TRUE)
  }
  put("volvox-sorted.bam")
  first <- jbrowse_alignments(path, "ctgA", 400, 1300)
  Sys.sleep(1.1)
  put("spliced.bam")
  second <- jbrowse_alignments(path, "ctgA", 400, 1300)
  fresh <- jbrowse_alignments(extdata("spliced.bam"), "ctgA", 400, 1300)
  expect_false(nrow(first$reads) == nrow(second$reads))
  expect_equal(second$reads$name, fresh$reads$name)
})

test_that("a path with %, ? and # in it reads", {
  odd <- if (.Platform$OS.type == "windows") "dir 50%41 a#c" else "dir 50%41 a?b#c"
  dir <- file.path(tempfile(), odd)
  dir.create(dir, recursive = TRUE)
  path <- file.path(dir, "x.bam")
  file.copy(extdata("volvox-sorted.bam"), path)
  file.copy(extdata("volvox-sorted.bam.bai"), paste0(path, ".bai"))
  expect_equal(nrow(jbrowse_alignments(path, "ctgA", 1000, 6000)$reads), 982)
})

test_that("a region over the byte limit is refused with its size", {
  expect_error(
    jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA", 0, 50000,
      byte_limit = 1000),
    "over the 1000-byte limit")
})

test_that("reads draw at the display's read height, and the pileup band is as tall as its rows", {
  bam <- extdata("volvox-sorted.bam")
  read_at <- function(h) jbrowse_alignments(bam, "ctgA:1,001-9,000",
    display = list(featureHeight = h))
  normal <- read_at(7)
  compact <- read_at(3)
  expect_equal(compact$settings$featureHeight, 3)
  body <- function(x) {
    d <- ggplot2::layer_data(ggplot2::ggplot(x) + geom_jbrowse_reads(), 1)
    unique(round(d$ymax - d$ymin, 6))
  }
  expect_equal(body(normal), 7 / 8)
  expect_equal(body(compact), 1)
  rows <- max(normal$reads$row) + 1
  tall <- jbrowse_stack(normal)
  short <- jbrowse_stack(compact)
  expect_equal(tall$patches$layout$heights, c(45, rows * 8))
  expect_equal(short$patches$layout$heights, c(45, rows * 3))
  expect_lt(jbrowse_height(short), jbrowse_height(tall))
})

test_that("a short pileup's figure is still tall enough for its keys", {
  r <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA:1,001-1,400",
    display = list(featureHeight = 1, color = list(field = "strand")))
  fig <- jbrowse_stack(r)
  bands <- sum(fig$patches$layout$heights) / 150 + 0.86
  expect_lt(bands, 2)
  expect_gt(jbrowse_height(fig), 2.2)
  expect_equal(jbrowse_height(fig & ggplot2::theme_classic()),
    jbrowse_height(fig), tolerance = 0.2)
})

test_that("a base modification call draws at least a display pixel wide", {
  m <- jbrowse_alignments(extdata("methylation_clip.bam"), "20", 10000, 22000,
    reference = extdata("hg38_clip.fa"),
    display = list(baseColor = list(field = "modifications")))
  d <- ggplot2::layer_data(ggplot2::ggplot(m) + geom_jbrowse_modifications(), 1)
  expect_true(all(d$xmax - d$xmin >= 12000 / 1500 - 1e-9))
})

test_that("a chain's overlapping segments come back as overlap spans", {
  colo <- extdata("colo829_tumor_chr3.bam")
  locus <- "chr3:25,358,500-25,359,800"
  chained <- jbrowse_alignments(colo, locus,
    display = list(unit = "chain", filter = list(split = "only")))
  expect_gt(nrow(chained$overlaps), 0)
  expect_true(all(chained$overlaps$end > chained$overlaps$start))
  expect_equal(nrow(jbrowse_alignments(colo, locus)$overlaps), 0)
})

test_that("a stacked pileup draws at the width its stack's tracks were read at", {
  locus <- "chr12:6,534,201-6,538,700"
  reads <- jbrowse_alignments(extdata("gapdh_direct_rna.bam"), locus)
  genes <- jbrowse_features(extdata("gapdh_gencode.gff3.gz"), locus,
    width = 2500)
  built <- function(p) lapply(ggplot2::ggplot_build(p)$data, function(d)
    d[setdiff(names(d), c("PANEL", "group"))])
  pileup <- built(jbrowse_stack(genes, reads))
  expect_equal(pileup, built(autoplot(reads, "pileup", width = 2500)))
  expect_false(identical(pileup, built(autoplot(reads, "pileup"))))
})

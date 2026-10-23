extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
bam <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA:1,001-6,000")

test_that("a BAM region reads as JBrowse lays it out", {
  expect_equal(nrow(bam$reads), 982)
  for (r in split(bam$reads, bam$reads$row)) {
    r <- r[order(r$start), ]
    expect_true(all(utils::tail(r$start, -1) >= utils::head(r$end, -1)))
  }
  m <- bam$mismatches
  read <- of(m, bam$reads, "read")
  expect_true(all(m$position >= read$start & m$position < read$end))
  expect_null(m$row)
  expect_equal(on_axis(bam)$mismatches$row, read$row)
  expect_setequal(unique(m$base), c("A", "C", "G", "T"))
  expect_true(all(m$frequency >= 0 & m$frequency <= 1))
  expect_gt(max(bam$coverage$depth), 20)
})

test_that("a CRAM reads as the same tables as its BAM", {
  cram <- jbrowse_alignments(extdata("volvox-sorted.cram"), "ctgA:1,001-6,000",
    assembly = extdata("volvox.fa"))
  key <- function(x) sort(paste(x$reads$name, x$reads$start, x$reads$end,
    x$reads$row))
  expect_equal(key(cram), key(bam))
  mm <- function(x) sort(paste(x$mismatches$position, x$mismatches$base,
    of(x$mismatches, x$reads, "read")$name))
  expect_equal(mm(cram), mm(bam))
})

test_that("insertions and clips sit on their reads, in the track's layers", {
  l <- jbrowse_alignments(extdata("volvox-long-reads-sv.bam"), "ctgA:1,001-20,000")
  expect_gt(nrow(l$insertions), 0)
  read <- of(l$insertions, l$reads, "read")
  expect_true(all(l$insertions$position >= read$start &
    l$insertions$position <= read$end))
  expect_equal(nchar(l$insertions$sequence), l$insertions$length)
  expect_setequal(unique(l$clips$type), c("soft", "hard"))
  expect_equal(on_axis(l)$clips$row, of(l$clips, l$reads, "read")$row)
  p <- ggplot2::ggplot_build(ggplot2::autoplot(l, "pileup"))
  layers <- vapply(p$data, nrow, 1L)
  expect_equal(layers[8], nrow(l$insertions))
  expect_equal(layers[9], nrow(l$clips))
  faded <- p$data[[7]]$alpha
  expect_equal(length(faded), nrow(l$mismatches))
  expect_true(any(faded < 1))
  expect_true(all(faded >= 1200 / 19000 - 1e-9))
  l$mismatches$fade <- 1
  opaque <- ggplot2::ggplot_build(ggplot2::autoplot(l, "pileup"))$data[[7]]$alpha
  expect_true(all(opaque == 1))
})

test_that("insertions fade and widen as the display's insertion shader draws them", {
  l <- jbrowse_alignments(extdata("volvox-long-reads-sv.bam"), "ctgA:1,001-20,000")
  px <- 1500 / 19000
  ins <- ggplot2::layer_data(ggplot2::autoplot(l, "pileup"), 8)
  short <- l$insertions$length < 10
  expect_true(all(ins$alpha[short] <= px * l$insertions$length[short] + 1e-9))
  expect_equal(ins$alpha[!short], pmin(1, l$insertions$length[!short] * px))
  bp <- diff(xlim_of(l))
  widths <- function(x) mark_alpha(x, bp, "insertions",
    data.frame(length = c(1, 9, 12, 30, 200), fade = 1))$widthPx
  expect_equal(widths(l), c(1, 1, 4, 22, 28))
  l$read_pitch$featureHeight <- 4
  expect_equal(widths(l)[5], 5)
  expect_equal(mark_alpha(l, bp, "deletions",
    data.frame(length = c(0.001, 0.5, 3), fade = 1))$alpha, c(0, 0.5, 1))
  expect_equal(mark_alpha(l, bp, "clips", data.frame(fade = numeric()))$alpha,
    numeric())
})

test_that("mismatches fade by base quality where the display does", {
  q <- jbrowse_alignments(extdata("gapdh_direct_rna.bam"),
    "chr12:6,534,201-6,538,700",
    display = list(pileup = list(mismatches = list(
      encoding = list(opacity = "baseQuality")))))
  q$mismatches$fade <- 1
  drawn <- ggplot2::layer_data(ggplot2::ggplot(q) + geom_jbrowse_mismatches())
  expect_equal(drawn$alpha, pmin(1, q$mismatches$quality / 50))
  expect_lt(min(drawn$alpha), 1)
})

test_that("mismatches$show off drops the mismatches and indels and keeps the skips", {
  off <- jbrowse_alignments(extdata("gapdh_direct_rna.bam"),
    "chr12:6,534,201-6,538,700",
    display = list(pileup = list(mismatches = list(show = FALSE))))
  expect_equal(nrow(off$mismatches), 0)
  expect_equal(nrow(off$insertions), 0)
  expect_gt(nrow(off$gaps), 0)
  expect_equal(unique(off$gaps$type), "skip")
  expect_gt(nrow(off$clips), 0)
})

test_that("the coverage band stacks each position's alleles as the display does", {
  a <- bam$allele_coverage
  expect_gt(nrow(a), 0)
  expect_true(all(a$position %in% bam$mismatches$position))
  expect_true(all(a$from >= 0 & a$to > a$from & a$to <= 1 + 1e-6))
  cov <- ggplot2::ggplot_build(ggplot2::autoplot(bam, "coverage"))
  alleles <- cov$data[[2]]
  expect_equal(nrow(alleles), nrow(a))
  expect_true(all(alleles$ymin >= 0 & alleles$ymax > alleles$ymin))
  depth <- bam$coverage$depth[findInterval(alleles$xmin, bam$coverage$start)]
  expect_true(all(alleles$ymax <= depth + 1e-6))
  strict <- bam
  strict$allele_coverage$shown <- with(strict$allele_coverage, to - from >= 0.5)
  strict <- ggplot2::ggplot_build(ggplot2::autoplot(strict, "coverage"))$data[[2]]
  expect_lt(nrow(strict), nrow(alleles))
  expect_true(all((strict$ymax - strict$ymin) /
    bam$coverage$depth[findInterval(strict$xmin, bam$coverage$start)] >= 0.5))
})

sv <- jbrowse_alignments(extdata("volvox-long-reads-sv.bam"), "ctgA:1,001-20,000")

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
  scale <- sv$y_scale
  built <- ggplot2::ggplot_build(ggplot2::autoplot(sv, "coverage"))
  y <- built$layout$panel_params[[1]]$y
  expect_equal(y$breaks[!is.na(y$breaks)], scale$ticks)
  inset <- sv$glyph$coverageInset / (sv$display$coverage$height -
    2 * sv$glyph$coverageInset)
  span <- diff(scale$domain)
  expect_equal(y$continuous_range,
    scale$domain + c(-1, 1) * inset * span)

  hidden <- jbrowse_alignments(extdata("volvox-long-reads-sv.bam"), "ctgA:1,001-20,000",
     display = list(coverage = list(interbase = list(show = FALSE))))
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
  top <- 2 * sv$y_scale$domain[2]
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
  log <- jbrowse_alignments(extdata("volvox-long-reads-sv.bam"), "ctgA:1,001-20,000",
     display = list(scales = list(y = list(type = "log"))))
  expect_equal(log$y_scale$domain[1], 1)
  built <- ggplot2::ggplot_build(ggplot2::autoplot(log, "coverage"))
  hang <- built$data[[4]]
  expect_equal(hang$hang, log$interbase_coverage$from /
    (2 * log$y_scale$domain[2]))
  # an allele keeps its share of the bar's drawn height from the bar's base,
  # one read: on a log axis its top is at depth^to
  a <- ggplot2::layer_data(ggplot2::ggplot(log) + geom_jbrowse_alleles())
  alleles <- log$allele_coverage
  expect_gt(nrow(alleles), 0)
  expect_equal(a$ymax, alleles$depth^alleles$to)
  expect_equal(a$ymin, alleles$depth^alleles$from)
})

test_that("every section shares the display's one depth axis", {
  s <- jbrowse_alignments(extdata("volvox-long-reads-sv.bam"), "ctgA:1,001-20,000",
     display = list(facet = "strand"))
  expect_gt(length(unique(s$interbase_coverage$section)), 1)
  built <- ggplot2::ggplot_build(ggplot2::autoplot(s, "coverage"))
  ranges <- lapply(built$layout$panel_params, function(p) p$y$continuous_range)
  expect_true(all(vapply(ranges, identical, TRUE, ranges[[1]])))
  expect_gte(s$y_scale$domain[2], max(s$coverage$depth))
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
  s <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA:1,001-6,000",
    display = list(sort = list(field = "base", pos = site)))
  covering <- which(s$reads$start <= site & s$reads$end > site)
  mismatched <- unique(s$mismatches$read[s$mismatches$position == site])
  expect_gt(length(mismatched), 0)
  expect_lt(max(s$reads$row[mismatched]),
    min(s$reads$row[setdiff(covering, mismatched)]))
  l <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA:1,001-6,000",
    display = list(sort = list(field = "length")))
  expect_false(identical(l$reads$row, bam$reads$row))
})

test_that("splice junctions arc over the coverage, counted and stranded", {
  s <- jbrowse_alignments(extdata("spliced.bam"), "ctgA:401-1,300",
    assembly = extdata("volvox.fa"))
  expect_equal(nrow(s$junctions), 2)
  expect_equal(sum(s$junctions$count), sum(s$gaps$type == "skip"))
  expect_true(all(s$junctions$strand %in% c(1, -1, 0)))
  s$junctions$shown <- s$junctions$count >= 1
  arcs <- ggplot2::ggplot_build(ggplot2::autoplot(s, "junctions"))
  expect_equal(length(unique(arcs$data[[1]]$group)), 2)
  expect_length(arcs$data, 1)
  counted <- jbrowse_alignments(extdata("spliced.bam"), "ctgA:401-1,300",
    assembly = extdata("volvox.fa"),
    display = list(sashimi = list(labels = list(show = TRUE))))
  counted$junctions$shown <- TRUE
  counts <- ggplot2::layer_data(ggplot2::autoplot(counted, "junctions"), 2)
  expect_equal(nrow(counts), 2)
  track <- ggplot2::autoplot(s)
  expect_s3_class(track, "patchwork")
  expect_s3_class(ggplot2::autoplot(bam), "patchwork")
})

test_that("base modifications paint under the layer, keyed and tinted as the display paints them", {
  bam <- extdata("methylation_clip.bam")
  plain <- jbrowse_alignments(bam, "20:10,001-20,000")
  expect_equal(nrow(plain$modifications), 0)
  m <- jbrowse_alignments(bam, "20:10,001-20,000",
    assembly = extdata("hg38_clip.fa"),
    display = list(pileup = list(bases = list(show = TRUE, encoding = list(color = "modifications")))))
  expect_gt(nrow(m$modifications), 100)
  expect_setequal(unique(m$modifications$name), "5mC")
  expect_equal(on_axis(m)$modifications$row,
    of(m$modifications, m$reads, "read")$row)
  expect_true(all(m$modifications$opacity >= 0.1 & m$modifications$opacity <= 1))
  expect_true("5mC" %in% m$keys$reads$label)
  expect_setequal(unique(m$reads$category), c("modFwd", "modRev"))
  expect_gt(nrow(m$modification_coverage), 0)
  band <- ggplot2::ggplot_build(ggplot2::autoplot(m, "coverage"))$data[[3]]
  expect_equal(nrow(band), nrow(m$modification_coverage))
  p <- ggplot2::ggplot_build(ggplot2::autoplot(m, "pileup"))
  expect_equal(nrow(p$data[[6]]), nrow(m$modifications))
  expect_true(all(p$data[[6]]$fill %in% m$keys$reads$color))
  calls <- jbrowse_alignments(bam, "20:10,001-20,000",
    display = list(pileup = list(bases = list(show = TRUE, encoding = list(color = "modifications"))),
      modifications = list(unmodified = "calls")))
  expect_gt(sum(calls$modifications$unmodified), 0)
  every <- jbrowse_alignments(bam, "20:10,001-20,000",
    display = list(pileup = list(bases = list(show = TRUE, encoding = list(color = "modifications"))),
      modifications = list(unmodified = "all")))
  expect_gt(sum(every$modifications$unmodified),
    10 * sum(calls$modifications$unmodified))
})

test_that("modification coloring greys the mismatched bases and drops their key", {
  bam <- extdata("methylation_clip.bam")
  fa <- extdata("hg38_clip.fa")
  plain <- jbrowse_alignments(bam, "20:10,001-12,000", assembly = fa)
  expect_gt(length(unique(plain$mismatches$color)), 1)
  expect_setequal(as.character(plain$keys$bases$color),
    unique(as.character(c(plain$mismatches$color,
      plain$allele_coverage$color))))
  m <- jbrowse_alignments(bam, "20:10,001-12,000", assembly = fa,
    display = list(pileup = list(bases = list(show = TRUE, encoding = list(color = "modifications")))))
  expect_gt(nrow(m$mismatches), 0)
  expect_setequal(as.character(m$mismatches$color), "#888888")
  expect_setequal(as.character(m$allele_coverage$color), "#888888")
  expect_equal(nrow(m$keys$bases), 0)
})

test_that("bisulfite reads paint their methylated cytosines by context", {
  bam <- extdata("arabidopsis_wgbs.bam")
  fa <- extdata("arabidopsis_tair10.fa.gz")
  locus <- "NC_003070.9:4,406,001-4,408,000"
  bisulfite <- function(modifications = NULL) jbrowse_alignments(bam, locus,
    assembly = fa, display = list(pileup = list(bases = list(show = TRUE, encoding = list(color = "bisulfite"))),
      modifications = modifications))
  cg <- bisulfite()
  expect_gt(nrow(cg$modifications), 100)
  expect_false(any(cg$modifications$unmodified))
  expect_setequal(cg$modifications$probability, 1)
  expect_equal(on_axis(cg)$modifications$row,
    of(cg$modifications, cg$reads, "read")$row)
  expect_true("5mC methylated" %in% cg$keys$reads$label)
  expect_gt(nrow(cg$modification_coverage), 0)
  expect_setequal(as.character(cg$mismatches$color), "#888888")
  expect_equal(nrow(cg$keys$bases), 0)
  seq <- jbrowse_sequence(fa, "NC_003070.9:4,406,001-4,408,000")
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
  expect_true("Unmethylated" %in% calls$keys$reads$label)
  p <- ggplot2::ggplot_build(ggplot2::autoplot(cg, "pileup"))
  expect_true(any(vapply(p$data, nrow, 0) == nrow(cg$modifications)))
  expect_error(jbrowse_alignments(bam, locus,
    display = list(pileup = list(bases = list(show = TRUE, encoding = list(color = "bisulfite"))))), "reference")
})

test_that("base quality paints each aligned base on its read, keyed as a bar", {
  bam <- extdata("methylation_clip.bam")
  q <- jbrowse_alignments(bam, "20:10,001-10,300",
    display = list(pileup = list(bases = list(show = TRUE, encoding = list(color = "baseQuality")))))
  expect_gt(nrow(q$qualities), 1000)
  expect_equal(on_axis(q)$qualities$row,
    of(q$qualities, q$reads, "read")$row)
  colors <- tapply(as.character(q$qualities$color), q$qualities$score,
    function(x) length(unique(x)))
  expect_true(all(colors == 1))
  expect_equal(vapply(q$ramps, attr, "", "title"), "Base quality")
  expect_equal(attr(q$ramps[[1]], "labels")[1], "0")
  p <- ggplot2::ggplot_build(ggplot2::autoplot(q, "pileup"))
  expect_equal(nrow(p$data[[3]]), nrow(q$qualities))
  expect_equal(nrow(jbrowse_alignments(bam, "20:10,001-10,300")$qualities), 0)
})

test_that("a bgzipped FASTA reads as a reference, and MD tags stand in for one", {
  fa <- extdata("gapdh_hg38.fa.gz")
  expect_equal(jbrowse_sequence(fa, "chr12:6,499,991-6,500,000"),
    strrep("N", 10))
  expect_match(jbrowse_sequence(fa, "chr12:6,534,501-6,534,600"),
    "^[ACGTacgt]{100}$")
  bam <- extdata("gapdh_direct_rna.bam")
  locus <- "chr12:6,534,201-6,538,700"
  tagged <- jbrowse_alignments(bam, locus)
  expect_gt(nrow(tagged$mismatches), 1000)
  expect_true(all(tagged$junctions$motif == ""))
  read <- jbrowse_alignments(bam, locus, assembly = fa)
  expect_equal(read$mismatches, tagged$mismatches)
  common <- read$junctions[read$junctions$count >= 10, ]
  expect_gt(nrow(common), 5)
  expect_true(all(common$motif == "GT-AG (canonical)"))
})

test_that("the reference reads as its bases, drawn as a band only when asked", {
  fa <- extdata("volvox.fa")
  seq <- jbrowse_sequence(fa, "ctgA:1,001-1,020")
  expect_match(seq, "^[ACGTNacgtn]{20}$")
  z <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA:1,001-1,100",
    assembly = fa)
  expect_equal(nchar(jbrowse_sequence(fa, "ctgA:1,001-1,100")), 100)
  expect_equal(nchar(jbrowse_sequence(fa, c("ctgB", "ctgA:1-10"))),
    c(6079, 10))
  s <- ggplot2::ggplot_build(ggplot2::autoplot(z, "sequence"))
  expect_equal(nrow(s$data[[2]]), 100)
  track <- ggplot2::autoplot(z)
  expect_s3_class(track, "patchwork")
  expect_equal(length(track$patches$plots) + 1, 2)
})

test_that("the sequence band reads each region's own bases at its offset", {
  fa <- extdata("volvox.fa")
  two <- jbrowse_alignments(extdata("volvox-sorted.bam"),
    c("ctgA:1,001-1,040", "ctgA:3,001-3,040"), assembly = fa)
  band <- ggplot2::ggplot_build(ggplot2::autoplot(two, "sequence"))$data[[2]]
  second <- band[band$x > two$regions$offset[2], ]
  second <- second[order(second$x), ]
  expect_equal(paste(second$label, collapse = ""),
    toupper(jbrowse_sequence(fa, "ctgA:3,001-3,040")))
  expect_equal(min(second$x), two$regions$offset[2] + 0.5)
  past_end <- jbrowse_alignments(extdata("volvox-sorted.bam"),
    "ctgA:49,901-50,100", assembly = fa)
  tail_band <- ggplot2::ggplot_build(ggplot2::autoplot(past_end,
    "sequence"))$data[[2]]
  expect_equal(nrow(tail_band), 101)
})

test_that("a region the assembly holds too few bases of warns", {
  fa <- extdata("volvox.fa")
  expect_warning(s <- jbrowse_sequence(fa, "ctgB:6,070-6,090"),
    "holds 10 of the 21 bases of ctgB:6,070-6,090")
  expect_equal(nchar(s), 10)
  expect_warning(other <- jbrowse_alignments(extdata("volvox-sorted.bam"),
    "ctgA:1,001-1,100", assembly = extdata("gapdh_hg38.fa.gz")),
    "holds no bases of ctgA:1,001-1,100")
  band <- ggplot2::ggplot_build(ggplot2::autoplot(other, "sequence"))
  expect_equal(nrow(band$data[[2]]), 0)
})

test_that("a spliced read splits into its blocks at its skips", {
  s <- jbrowse_alignments(extdata("spliced.bam"), "ctgA:401-1,300")
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
  expect_setequal(names(base_colors()), c("A", "C", "G", "T", "N"))
  expect_equal(bam$keys$reads$label, "Reads")
  expect_equal(unique(as.character(bam$reads$color)), bam$keys$reads$color)
  p <- ggplot2::ggplot_build(ggplot2::autoplot(bam))
  expect_equal(unique(p$data[[1]]$fill), bam$keys$reads$color)
  expect_equal(nrow(p$data[[7]]), nrow(bam$mismatches))
  cov <- ggplot2::ggplot_build(ggplot2::autoplot(bam, "coverage"))
  expect_equal(nrow(cov$data[[1]]), nrow(bam$coverage))
})

test_that("every mark takes the engine's color, and every key the browser's wording", {
  a <- jbrowse_alignments(extdata("gapdh_direct_rna.bam"),
    "chr12:6,534,201-6,538,700")
  for (table in c("mismatches", "gaps", "insertions", "clips", "coverage",
    "allele_coverage", "interbase_coverage", "indicators", "junctions")) {
    expect_type(a[[table]]$color, "character")
    expect_false(anyNA(a[[table]]$color))
  }
  expect_false(identical(unique(a$gaps$color[a$gaps$type == "skip"]),
    unique(a$gaps$color[a$gaps$type == "deletion"])))
  expect_equal(attr(a$keys$junctions, "title"), "Sashimi arcs")
  expect_equal(a$keys$junctions$label, "Forward transcript strand")
  expect_equal(attr(a$keys$interbase, "title"), "Insertions and clips")
  expect_equal(attr(a$keys$bases, "title"), "Base")
  expect_setequal(as.character(a$keys$bases$color),
    unique(as.character(c(a$mismatches$color, a$allele_coverage$color))))
})

test_that("reads fill as the display's color object paints them, with its key", {
  s <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA:1,001-6,000",
    display = list(encoding = list(color = list(field = "strand"))))
  expect_equal(attr(s$keys$reads, "title"), "Read colors")
  expect_equal(s$keys$reads$label, c("Forward strand", "Reverse strand"))
  expect_equal(unique(s$reads$color[s$reads$strand == 1]), s$keys$reads$color[1])
  expect_equal(unique(s$reads$color[s$reads$strand == -1]), s$keys$reads$color[2])
  expect_setequal(unique(s$reads$category), c("fwdStrand", "revStrand"))
  p <- ggplot2::ggplot_build(ggplot2::autoplot(s))
  expect_setequal(unique(p$data[[1]]$fill), s$keys$reads$color)

  own <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA:1,001-6,000",
    display = list(encoding = list(color = list(field = "strand",
      scale = list(range = c("#ff0000", "#0000ff")),
      legend = list(title = "Mine")))))
  expect_equal(unique(own$reads$color[own$reads$strand == 1]), "#ff0000")
  expect_equal(attr(own$keys$reads, "title"), "Mine")

  m <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA:1,001-6,000",
    display = list(encoding = list(color = list(field = "mapq"))))
  expect_equal(attr(m$keys$reads, "title"), "Mapping quality")
  expect_true(all(m$reads$color %in% m$keys$reads$color))
  expect_equal(unique(m$reads$category), "mapq")
})

test_that("the display's filter list keeps the reads it keeps", {
  r <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA:1,001-6,000",
    display = list(filter = list(list(field = "flags.REVERSE", equal = TRUE))))
  expect_true(all(r$reads$strand == -1))
  expect_gt(nrow(r$reads), 0)
  expect_lt(nrow(r$reads), nrow(bam$reads))
})

test_that("svSignal only and exclude split the reads between them", {
  locus <- "chr3:25,357,600-25,361,000"
  bam <- extdata("colo829_tumor_chr3.bam")
  read <- function(equal) jbrowse_alignments(bam, locus,
    display = list(filter = list(list(field = "svSignal", equal = equal))))
  all <- jbrowse_alignments(bam, locus)
  only <- read(TRUE)
  none <- read(FALSE)
  expect_gt(nrow(only$reads), 0)
  expect_gt(nrow(none$reads), 0)
  expect_equal(nrow(only$reads) + nrow(none$reads), nrow(all$reads))
})

test_that("a remote BAM reads over HTTP", {
  skip_on_cran()
  skip_if_offline("jbrowse.org")
  u <- "https://jbrowse.org/code/jb2/main/test_data/volvox/volvox-sorted.bam"
  r <- jbrowse_alignments(u, "ctgA:1,001-2,000")
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
  first <- jbrowse_alignments(path, "ctgA:401-1,300")
  Sys.sleep(1.1)
  put("spliced.bam")
  second <- jbrowse_alignments(path, "ctgA:401-1,300")
  fresh <- jbrowse_alignments(extdata("spliced.bam"), "ctgA:401-1,300")
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
  expect_equal(nrow(jbrowse_alignments(path, "ctgA:1,001-6,000")$reads), 982)
})

test_that("a region over a limit is refused naming R's argument", {
  bam <- extdata("volvox-sorted.bam")
  expect_error(jbrowse_alignments(bam, "ctgA:1-50,000", byte_limit = 1000),
    "ctgA:1-50,000 needs [0-9,]+ bytes of alignments, past `byte_limit` \\(1,000 bytes\\)",
    class = "ggjbrowse_too_large")
  expect_error(jbrowse_alignments(bam, c("ctgA:1-20,000", "ctgA:30,001-40,000"),
    span_limit = 1000), paste0("ctgA:1-20,000, ctgA:30,001-40,000 span 30,000 ",
      "bp of alignments in total, past `span_limit` \\(1,000 bp\\)"),
    class = "ggjbrowse_too_large")
  limit <- tryCatch(jbrowse_features(extdata("volvox.sort.gff3.gz"),
    "ctgA:1-50,000", byte_limit = 10), error = identity)
  expect_s3_class(limit, "ggjbrowse_too_large")
  expect_equal(limit$limit, "byte_limit")
})

test_that("a display slot the schema refuses is named", {
  expect_error(jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA:1-400",
    display = list(unit = "bogus")), "`display\\$unit` is \"bogus\"",
    class = "ggjbrowse_bad_display")
  expect_error(jbrowse_features(extdata("volvox.sort.gff3.gz"), "ctgA:1-400",
    display = list(displayMode = "huge")), "`display\\$displayMode`",
    class = "ggjbrowse_bad_display")
})

test_that("a whole file reads the contigs its index shows reads on, to the assembly's lengths", {
  cram <- extdata("gapdh_direct_rna.cram")
  fa <- extdata("gapdh_hg38.fa.gz")
  listed <- file_contigs(alignments_adapter(cram, NULL), assembly_config(fa))
  listed <- listed[!listed$has_data %in% FALSE, ]
  expect_equal(listed$chrom, "chr12")
  expect_equal(listed$end, 6570000)
  whole <- jbrowse_alignments(cram, assembly = fa)
  expect_equal(max(whole$regions$end), 6570000)
  expect_equal(nrow(whole$reads), nrow(jbrowse_alignments(cram,
    "chr12:1-6,570,000", assembly = fa)$reads))
})

test_that("a whole file spanning past span_limit is refused before any read", {
  expect_error(jbrowse_alignments(extdata("gapdh_direct_rna.bam")),
    "the whole file spans 3,217,346,917 bp over 3366 contigs, past `span_limit`")
})

test_that("regions spanning past span_limit in all are refused", {
  loci <- c("ctgA:1-3,000", "ctgA:10,001-13,000")
  expect_error(jbrowse_alignments(extdata("volvox-sorted.bam"), loci,
    span_limit = 5000), "6,000 bp of alignments in total, past `span_limit`")
  expect_no_error(jbrowse_alignments(extdata("volvox-sorted.bam"), loci,
    span_limit = 6000))
})

test_that("reads draw at the display's read height, and the pileup band is as tall as its rows", {
  bam <- extdata("volvox-sorted.bam")
  read_at <- function(h) jbrowse_alignments(bam, "ctgA:1,001-9,000",
    display = list(pileup = list(readHeight = h)))
  normal <- read_at(7)
  compact <- read_at(3)
  expect_equal(compact$display$pileup$readHeight, 3)
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
    display = list(pileup = list(readHeight = 1), encoding = list(color = list(field = "strand"))))
  fig <- jbrowse_stack(r)
  bands <- sum(fig$patches$layout$heights) / 150 + 0.86
  expect_lt(bands, 2)
  expect_gt(jbrowse_height(fig), 2.2)
  expect_equal(jbrowse_height(fig & ggplot2::theme_classic()),
    jbrowse_height(fig), tolerance = 0.2)
})

test_that("a base modification call draws at least a display pixel wide", {
  m <- jbrowse_alignments(extdata("methylation_clip.bam"), "20:10,001-22,000",
    assembly = extdata("hg38_clip.fa"),
    display = list(pileup = list(bases = list(show = TRUE, encoding = list(color = "modifications")))))
  d <- ggplot2::layer_data(ggplot2::ggplot(m) + geom_jbrowse_modifications(), 1)
  expect_true(all(d$xmax - d$xmin >= 12000 / 1500 - 1e-9))
})

test_that("a chain's overlapping segments come back as overlap spans", {
  colo <- extdata("colo829_tumor_chr3.bam")
  locus <- "chr3:25,358,500-25,359,800"
  chained <- jbrowse_alignments(colo, locus,
    display = list(unit = "chain", filter = list(list(field = "split", equal = TRUE))))
  expect_gt(nrow(chained$overlaps), 0)
  expect_true(all(chained$overlaps$end > chained$overlaps$start))
  key <- chained$keys$reads
  expect_equal(key$label[nrow(key)], "Pair/chain reads overlap here")
  expect_equal(key$color[nrow(key)], chained$overlaps$color[1])
  expect_equal(nrow(jbrowse_alignments(colo, locus)$overlaps), 0)
  thin <- jbrowse_alignments(colo, locus, display = list(unit = "chain",
    filter = list(list(field = "split", equal = TRUE)), pileup = list(readHeight = 2)))
  expect_equal(nrow(thin$overlaps), 0)
  expect_false("Pair/chain reads overlap here" %in% thin$keys$reads$label)
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

test_that("an empty display slot reads as the display's defaults", {
  bam <- extdata("arabidopsis_wgbs.bam")
  fa <- extdata("arabidopsis_tair10.fa.gz")
  locus <- "NC_003070.9:4,406,001-4,407,000"
  read <- function(modifications) jbrowse_alignments(bam, locus,
    assembly = fa, display = list(pileup = list(bases = list(show = TRUE, encoding = list(color = "bisulfite"))),
      modifications = modifications))
  expect_equal(read(list())$modifications, read(NULL)$modifications)
})

test_that("the sequence band asks for the reference it draws", {
  a <- jbrowse_alignments(extdata("arabidopsis_wgbs.bam"),
    "NC_003070.9:4,406,001-4,406,300")
  expect_error(autoplot(a, "sequence"), "assembly =")
})

test_that("every section keys the track's read colors, insertions and clips", {
  a <- jbrowse_alignments(extdata("hg00151_chr1_inversion.bam"),
    "chr1:197,786,900-197,789,700", display = list(unit = "chain", encoding = list(color = list(field = "pairOrientation")), facet = list(field = "split")))
  flat <- function(p) c(if (inherits(p, "patchwork"))
    unlist(lapply(p$patches$plots, flat), recursive = FALSE), list(p))
  keyed <- lapply(flat(autoplot(a)), function(p)
    unlist(lapply(drawn_guides(p), `[[`, "color")))
  for (key in a$keys[c("reads", "interbase")]) {
    sections <- Filter(function(k) any(k %in% key$color), keyed)
    expect_length(sections, 2)
    for (k in sections) expect_true(all(key$color %in% k))
  }
})

test_that("a stack draws a key two tracks share once, with every entry either holds", {
  locus <- "chr1:197,786,900-197,789,700"
  tumor <- jbrowse_alignments(extdata("hg00151_chr1_inversion.bam"), locus)
  normal <- jbrowse_alignments(extdata("hg00151_chr1_inversion.bam"), locus,
    display = list(filter = list(
      list(field = "flags.SUPPLEMENTARY", equal = FALSE))))
  expect_true("Hard clip" %in% tumor$keys$interbase$label)
  expect_false("Hard clip" %in% normal$keys$interbase$label)
  collected <- function(fig) {
    grob <- patchwork::patchworkGrob(fig)
    text <- character(0)
    walk <- function(g) {
      if (inherits(g, "text") || inherits(g, "titleGrob")) {
        text <<- c(text, g$label)
      }
      for (child in c(g$children, g$grobs)) walk(child)
    }
    for (box in grob$grobs[startsWith(grob$layout$name, "guide-box")]) {
      walk(box)
    }
    text
  }
  for (fig in list(jbrowse_stack(tumor, normal), jbrowse_stack(normal, tumor))) {
    text <- collected(fig)
    expect_equal(sum(text == "Insertions and clips"), 1)
    expect_equal(sum(text == "Hard clip"), 1)
    expect_equal(sum(text == "Base"), 1)
  }

  locus <- "chr12:6,534,201-6,538,700"
  forward <- jbrowse_alignments(extdata("gapdh_direct_rna.bam"), locus,
    display = list(pileup = list(show = FALSE)))
  expect_equal(unique(forward$junctions$strand), 1)
  reverse <- forward
  reverse$junctions$strand <- -1
  reverse$junctions$color <- "#8f8fd8"
  reverse$keys$junctions <- key_frame("Reverse transcript strand", "#8f8fd8",
    "Sashimi arcs")
  text <- collected(jbrowse_stack(forward, reverse))
  expect_equal(sum(text == "Sashimi arcs"), 1)
  expect_equal(sum(text %in% c("Forward transcript strand",
    "Reverse transcript strand")), 2)
})

test_that("the junction band draws the junctions the display shows, at its heights and strokes", {
  a <- jbrowse_alignments(extdata("gapdh_direct_rna.bam"),
    "chr12:6,534,201-6,538,700")
  j <- a$junctions
  expect_gt(sum(!j$shown), 0)
  drawn <- ggplot2::layer_data(ggplot2::ggplot(a) + geom_jbrowse_junctions())
  expect_equal(length(unique(drawn$group)), sum(j$shown))
  shown <- j[j$shown, ]
  expect_equal(tapply(drawn$y, drawn$group, max), shown$arcHeight,
    tolerance = 1e-3, ignore_attr = TRUE)
  expect_equal(tapply(drawn$linewidth, drawn$group, unique),
    shown$strokeWidth / PX_PER_INCH * 72.27 / ggplot2::.pt,
    ignore_attr = TRUE)
  a$junctions$shown <- TRUE
  every <- ggplot2::layer_data(ggplot2::ggplot(a) + geom_jbrowse_junctions())
  expect_equal(length(unique(every$group)), nrow(j))
})

test_that("the coverage band stacks the alleles the display shows", {
  a <- jbrowse_alignments(extdata("gapdh_direct_rna.bam"),
    "chr12:6,534,201-6,538,700",
    display = list(coverage = list(minAlleleFraction = 0.2)))
  expect_gt(sum(!a$allele_coverage$shown), 0)
  drawn <- ggplot2::layer_data(ggplot2::ggplot(a) + geom_jbrowse_alleles())
  expect_equal(nrow(drawn), sum(a$allele_coverage$shown))
  a$allele_coverage$shown <- TRUE
  every <- ggplot2::layer_data(ggplot2::ggplot(a) + geom_jbrowse_alleles())
  expect_equal(nrow(every), nrow(a$allele_coverage))
})

test_that("a result carries the display's resolved config and read pitch", {
  expect_equal(bam$display$type, "LinearAlignmentsDisplay")
  floor <- Filter(function(p) identical(p$field, "junction.reads"),
    bam$display$filter)
  expect_equal(floor[[1]]$gte, 2)
  expect_false(bam$display$readConnections$show)
  expect_equal(bam$read_pitch, list(featureHeight = 7L, featureSpacing = 1L))
  compact <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA:1,001-2,000",
     display = list(pileup = list(readHeight = 3)))
  expect_equal(compact$read_pitch$featureSpacing, 0)
  expect_equal(read_inset(compact), 0)
  off <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA:1,001-2,000",
     display = list(pileup = list(fadeLowFrequency = FALSE)))
  expect_false(fades(off))
  expect_true(fades(bam))
})

test_that("the junction band hides non-canonical junctions where the display does", {
  read <- function(...) jbrowse_alignments(extdata("gapdh_direct_rna.bam"),
    "chr12:6,534,201-6,538,700", assembly = extdata("gapdh_hg38.fa.gz"),
    display = list(filter = list(list(field = "junction.reads", gte = 1), ...)))
  arcs <- function(x) length(unique(ggplot2::layer_data(
    ggplot2::ggplot(x) + geom_jbrowse_junctions(), 1)$group))
  shown <- read()
  odd <- shown$junctions$motif == "non-canonical" & shown$junctions$shown
  expect_gt(sum(odd), 0)
  canonical <- function(equal) list(field = "junction.canonical", equal = equal)
  expect_equal(arcs(read(canonical(TRUE))), arcs(shown) - sum(odd))
  expect_equal(arcs(read(canonical(FALSE))), sum(odd))
})

test_that("coverage stops at the region's ends, where the reads do not", {
  expect_lt(min(bam$reads$start), 1000)
  expect_gt(max(bam$reads$end), 6000)
  expect_equal(range(c(bam$coverage$start, bam$coverage$end)), c(1000, 6000))
  expect_equal(nrow(bam$coverage), 5000)
})

test_that("a band the display hides reads as an empty table", {
  hidden <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA:1,001-6,000",
     display = list(coverage = list(show = FALSE)))
  expect_equal(nrow(hidden$coverage), 0)
  expect_equal(nrow(hidden$allele_coverage), 0)
  expect_equal(hidden$reads$start, bam$reads$start)
})

test_that("a read's id is its record's, where mates share a name", {
  paired <- jbrowse_alignments(extdata("volvox-simple-inv-paired.bam"),
    "ctgA:18,001-32,000")$reads
  expect_lt(length(unique(paired$name)), nrow(paired))
  expect_false(anyDuplicated(paired$id) > 0)
  later <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA:3,001-6,000")$reads
  shared <- match(later$id, bam$reads$id)
  expect_gt(sum(!is.na(shared)), 100)
  expect_equal(later$start[!is.na(shared)], bam$reads$start[na.omit(shared)])
  expect_equal(later$name[!is.na(shared)], bam$reads$name[na.omit(shared)])
})

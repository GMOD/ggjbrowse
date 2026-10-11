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
})

test_that("the pileup and coverage draw in the display's colors", {
  colors <- jbrowse_colors()
  expect_type(colors$bases, "character")
  expect_setequal(names(colors$bases), c("A", "C", "G", "T", "N"))
  expect_equal(unique(bam$reads$color), tolower(colors$read))
  expect_equal(bam$legend$label, "Reads")
  p <- ggplot2::ggplot_build(ggplot2::autoplot(bam))
  expect_equal(unique(p$data[[1]]$fill), tolower(colors$read))
  expect_equal(nrow(p$data[[4]]), nrow(bam$mismatches))
  cov <- ggplot2::ggplot_build(ggplot2::autoplot(bam, "coverage"))
  expect_equal(nrow(cov$data[[1]]), nrow(bam$coverage))
})

test_that("reads fill as the display's color object paints them, with its key", {
  s <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA", 1000, 6000,
    color = "strand")
  expect_equal(attr(s$legend, "title"), "Read colors")
  expect_equal(s$legend$label, c("Forward strand", "Reverse strand"))
  expect_equal(unique(s$reads$color[s$reads$strand == 1]), s$legend$color[1])
  expect_equal(unique(s$reads$color[s$reads$strand == -1]), s$legend$color[2])
  expect_setequal(unique(s$reads$category), c("fwdStrand", "revStrand"))
  p <- ggplot2::ggplot_build(ggplot2::autoplot(s))
  expect_setequal(unique(p$data[[1]]$fill), s$legend$color)

  own <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA", 1000, 6000,
    color = list(field = "strand", range = c("#ff0000", "#0000ff"),
      title = "Mine"))
  expect_equal(unique(own$reads$color[own$reads$strand == 1]), "#ff0000")
  expect_equal(attr(own$legend, "title"), "Mine")

  m <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA", 1000, 6000,
    color = list(field = "mapq"))
  expect_equal(attr(m$legend, "title"), "Mapping quality")
  expect_true(all(m$reads$color %in% m$legend$color))
  expect_equal(unique(m$reads$category), "mapq")
})

test_that("the display's filter object keeps the reads it keeps", {
  r <- jbrowse_alignments(extdata("volvox-sorted.bam"), "ctgA", 1000, 6000,
    filter = list(flagInclude = 16))
  expect_true(all(r$reads$strand == -1))
  expect_gt(nrow(r$reads), 0)
  expect_lt(nrow(r$reads), nrow(bam$reads))
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

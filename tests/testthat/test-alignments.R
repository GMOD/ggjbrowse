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
  p <- ggplot2::ggplot_build(ggplot2::autoplot(bam, color_by = "strand"))
  expect_setequal(unique(toupper(p$data[[1]]$fill)),
    toupper(c(colors$forward, colors$reverse)))
  expect_equal(nrow(p$data[[4]]), nrow(bam$mismatches))
  cov <- ggplot2::ggplot_build(ggplot2::autoplot(bam, "coverage"))
  expect_equal(nrow(cov$data[[1]]), nrow(bam$coverage))
})

test_that("a remote BAM reads over HTTP", {
  skip_on_cran()
  skip_if_offline("jbrowse.org")
  u <- "https://jbrowse.org/code/jb2/main/test_data/volvox/volvox-sorted.bam"
  r <- jbrowse_alignments(u, "ctgA", 1000, 2000)
  expect_gt(nrow(r$reads), 100)
})

bam <- system.file("extdata", "volvox-sorted.bam", package = "ggjbrowse")

test_that("a handler that catches an engine warning leaves the engine working", {
  caught <- tryCatch(
    jbrowse_alignments(bam, "ctgA:1,001-9,000", display = list(colour = "red")),
    warning = function(w) conditionMessage(w))
  expect_match(caught, "does not declare colour")
  expect_gt(nrow(jbrowse_alignments(bam, "ctgA:1,001-9,000")$reads), 0)
})

test_that("an adapter guess ignores a URL's query", {
  expect_equal(alignments_adapter("https://host/x.cram?token=abc", NULL)$type,
    "CramAdapter")
  expect_equal(alignments_adapter("https://host/x.bam?token=abc", NULL)$type,
    "BamAdapter")
  expect_equal(sequence_adapter("https://host/x.2bit?x=1")$type,
    "TwoBitAdapter")
  expect_equal(sequence_adapter("https://host/x.fa.gz?x=1")$type,
    "BgzipFastaAdapter")
})

test_that("a nested column reads as a list column, a missing row NA", {
  wire <- jsonlite::fromJSON(paste0('{"t": {"n": [1, 2, 3], "x": {"$list": ',
    '{"lengths": [2, -1, 1], "values": [4, 5, 6]}}}}'), simplifyVector = TRUE,
    simplifyDataFrame = FALSE)
  d <- as_frame(from_wire(wire)$t)
  expect_equal(nrow(d), 3)
  expect_identical(d$x, list(c(4L, 5L), NA_integer_, 6L))
  same <- from_wire(list(x = list("$list" = list(lengths = c(2L, 2L),
    values = c("a", "b", "c", "d")))))
  expect_identical(same$x, list(c("a", "b"), c("c", "d")))
})

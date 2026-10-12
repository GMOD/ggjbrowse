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

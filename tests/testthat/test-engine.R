bam <- system.file("extdata", "volvox-sorted.bam", package = "ggjbrowse")

test_that("a handler that catches an engine warning leaves the engine working", {
  caught <- tryCatch(
    jbrowse_alignments(bam, "ctgA:1,001-9,000", display = list(colour = "red")),
    warning = function(w) conditionMessage(w))
  expect_match(caught, "does not declare colour")
  expect_gt(nrow(jbrowse_alignments(bam, "ctgA:1,001-9,000")$reads), 0)
})

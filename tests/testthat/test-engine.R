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

test_that("an index beside a URL takes the URL's query after its extension", {
  bam <- alignments_adapter("https://host/x.bam?token=abc", NULL)
  expect_equal(bam$index$location$uri, "https://host/x.bam.bai?token=abc")
  cram <- alignments_adapter("https://host/x.cram?token=abc", NULL)
  expect_equal(cram$craiLocation$uri, "https://host/x.cram.crai?token=abc")
  fasta <- sequence_adapter("https://host/x.fa.gz?x=1")
  expect_equal(c(fasta$faiLocation$uri, fasta$gziLocation$uri),
    c("https://host/x.fa.gz.fai?x=1", "https://host/x.fa.gz.gzi?x=1"))
  vcf <- features_adapter("https://host/x.vcf.gz?x=1",
    "https://host/x.vcf.gz.csi?x=1")
  expect_equal(vcf$index$indexType, "CSI")
})

test_that("a nested column reads as a list column, a missing row NA", {
  wire <- jsonlite::fromJSON(paste0('{"t": {"n": [1, 2, 3], "x": {"$list": ',
    '{"lengths": [2, -1, 1], "values": [4, 5, 6]}}}}'), simplifyVector = TRUE,
    simplifyDataFrame = FALSE)
  d <- as_frame(from_wire(wire, NULL)$t)
  expect_equal(nrow(d), 3)
  expect_identical(d$x, list(c(4L, 5L), NA_integer_, 6L))
  same <- from_wire(list(x = list("$list" = list(lengths = c(2L, 2L),
    values = c("a", "b", "c", "d")))), NULL)
  expect_identical(same$x, list(c("a", "b"), c("c", "d")))
})

test_that("numbers read off the engine's bytes, NA where it wrote none", {
  bytes <- c(
    writeBin(c(7L, NA, -3L), raw(), size = 4, endian = "little"),
    writeBin(c(0.5, NA), raw(), size = 8, endian = "little"),
    writeBin(c(1L, 0L, NA), raw(), size = 4, endian = "little"),
    writeBin(c(1L, -1L, 0L), raw(), size = 4, endian = "little"))
  con <- rawConnection(bytes)
  on.exit(close(con))
  wire <- jsonlite::fromJSON(paste0('{"i": {"$i32": [0, 3]}, ',
    '"f": {"$f64": [12, 2]}, "l": {"$lgl": [28, 3]}, ',
    '"d": {"$dict": {"levels": ["red", "blue"], "codes": {"$i32": [40, 3]}}}}'),
    simplifyVector = TRUE, simplifyDataFrame = FALSE)
  x <- from_wire(wire, con)
  expect_identical(x$i, c(7L, NA, -3L))
  expect_identical(x$f, c(0.5, NA))
  expect_identical(x$l, c(TRUE, FALSE, NA))
  expect_identical(x$d, c("blue", NA, "red"))
})

test_that("a column keeps the type it had when the tables crossed as JSON", {
  bam <- jbrowse_alignments(system.file("extdata", "volvox-sorted.bam",
    package = "ggjbrowse"), "ctgA", 1000, 2000)
  expect_type(bam$reads$start, "integer")
  expect_type(bam$reads$strand, "integer")
  expect_type(bam$mismatches$frequency, "double")
  expect_type(bam$reads$category, "character")
  expect_type(bam$reads$name, "character")
})

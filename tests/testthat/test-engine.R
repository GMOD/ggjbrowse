extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
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
  expect_equal(beside("https://host/x.bam#part", ".bai"),
    "https://host/x.bam.bai#part")
  expect_equal(features_adapter("https://host/x.vcf.gz#part", NULL)$type,
    "VcfTabixAdapter")
  expect_equal(bare_name("/data/a?b.vcf.gz"), "/data/a?b.vcf.gz")
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
    package = "ggjbrowse"), "ctgA:1,001-2,000")
  expect_type(bam$reads$start, "integer")
  expect_type(bam$reads$strand, "integer")
  expect_type(bam$mismatches$frequency, "double")
  expect_type(bam$reads$category, "character")
  expect_type(bam$reads$name, "character")
})

test_that("the main process is no fork, on every platform", {
  expect_false(forked())
})

test_that("a forked worker refuses to start the engine, and inherits a started one", {
  skip_on_os("windows")
  bam <- system.file("extdata", "volvox-sorted.bam", package = "ggjbrowse")
  read <- function() nrow(jbrowse_alignments(bam, "ctgA:1-5,000")$reads)
  in_fork <- function(f) parallel::mccollect(parallel::mcparallel(f()))[[1]]
  started <- read()
  expect_equal(in_fork(read), started)
  expect_equal(in_fork(function() {
    jbrowse_reset()
    read()
  }), started)
  expect_equal(in_fork(function() tryCatch({
    engine(drop = TRUE)
    read()
  }, ggjbrowse_forked = function(e) "refused")), "refused")
})

test_that("an engine error reads once, without JavaScript's prefix", {
  vcf <- system.file("extdata", "1000g_dcst2.vcf.gz", package = "ggjbrowse")
  e <- tryCatch(jbrowse_variants(vcf, "1:1-1000",
    display = list(unit = "allele")), error = conditionMessage)
  expect_match(e, "^`display\\$unit` is \"allele\"")
})

test_that("the engine's warnings warn once a call", {
  bam <- system.file("extdata", "volvox-sorted.bam", package = "ggjbrowse")
  w <- character(0)
  withCallingHandlers(jbrowse_alignments(bam,
    c("ctgA:1-1,000", "ctgA:2,001-3,000"),
    display = list(sort = list(type = "basePair"))),
    warning = function(c) {
      w <<- c(w, conditionMessage(c))
      invokeRestart("muffleWarning")
    })
  expect_equal(w, "the sort by basePair names no pos, so it is dropped")
})

# Copies of one BAM under new names, each a file the engine has not read
bam_copies <- function(n) {
  dir <- tempfile()
  dir.create(dir)
  paths <- file.path(dir, sprintf("copy%03d.bam", seq_len(n)))
  file.copy(bam, paths)
  file.copy(paste0(bam, ".bai"), paste0(paths, ".bai"))
  paths
}

# The URLs the engine reads whole while `code` runs, each as often as read
whole_reads <- function(code) {
  reads <- character(0)
  real <- read_range
  local_mocked_bindings(read_range = function(url, start, end) {
    if (start == 0 && end < 0) reads <<- c(reads, url)
    real(url, start, end)
  }, .env = asNamespace("ggjbrowse"))
  force(code)
  table(sub("\\?.*$", "", basename(reads)))
}

test_that("a read through an assembly opens each file once", {
  fa <- extdata("gapdh_hg38.fa.gz")
  for (file in extdata(c("gapdh_direct_rna.bam", "gapdh_direct_rna.cram"))) {
    jbrowse_reset()
    reads <- whole_reads(jbrowse_alignments(file, "chr12:6,534,201-6,538,700",
      assembly = fa))
    expect_gt(length(reads), 1)
    expect_true(all(reads == 1))
    expect_equal(jbrowse_cache()$adapters, 2)
  }
})

test_that("the engine keeps decoded records and chunks under the limit, and the files read last open", {
  limit <- jbrowse_cache()$limit
  on.exit(jbrowse_cache(limit = limit))
  jbrowse_reset()
  jbrowse_cache(limit = 4 * 2^20)
  paths <- bam_copies(12)
  for (path in paths) {
    jbrowse_alignments(path, "ctgA:1-2,000")
    held <- jbrowse_cache()
    expect_lte(held$decoded, 4 * 2^20)
    expect_lte(held$chunks, 4 * 2^20)
  }
  expect_lt(held$adapters, length(paths))
  indexes <- function(path) whole_reads(jbrowse_alignments(path,
    "ctgA:1-2,000"))[paste0(basename(path), ".bai")]
  expect_true(is.na(indexes(paths[11])))
  expect_equal(unname(indexes(paths[1])), 1)
  expect_error(jbrowse_cache(limit = -1), "number of bytes")
})

test_that("a lower limit keeps the files the last call read open", {
  limit <- jbrowse_cache()$limit
  on.exit(jbrowse_cache(limit = limit))
  jbrowse_reset()
  jbrowse_alignments(bam, "ctgA:1-2,000")
  expect_equal(jbrowse_cache(limit = 0)$adapters, 1)
})

test_that("a limit of 0 keeps no chunks, and a higher one keeps them", {
  limit <- jbrowse_cache()$limit
  on.exit(jbrowse_cache(limit = limit))
  jbrowse_reset()
  jbrowse_cache(limit = 0)
  jbrowse_alignments(bam, "ctgA:1-2,000")
  expect_equal(jbrowse_cache()$chunks, 0)
  jbrowse_reset()
  jbrowse_cache(limit = 2^30)
  jbrowse_alignments(bam, "ctgA:1-2,000")
  expect_gt(jbrowse_cache()$chunks, 0)
})

test_that("a limit of Inf is none, and a restarted engine takes it again", {
  limit <- jbrowse_cache()$limit
  on.exit(jbrowse_cache(limit = limit))
  expect_equal(jbrowse_cache(limit = Inf)$limit, Inf)
  engine(drop = TRUE)
  expect_equal(jbrowse_cache()$limit, Inf)
})

test_that("an engine that cannot empty its caches restarts, and says so", {
  jbrowse_alignments(bam, "ctgA:1-2,000")
  local_mocked_bindings(in_engine = function(...) stop("unreachable"))
  expect_warning(jbrowse_reset(),
    "could not empty its caches, so it restarts: unreachable")
  expect_null(engine(start = FALSE))
})

test_that("jbrowse_reset() empties the caches and keeps the limit", {
  limit <- jbrowse_cache()$limit
  on.exit(jbrowse_cache(limit = limit))
  jbrowse_cache(limit = 2^30)
  jbrowse_alignments(bam, "ctgA:1-5,000")
  expect_gt(jbrowse_cache()$decoded, 0)
  jbrowse_reset()
  expect_equal(jbrowse_cache(), list(limit = 2^30, decoded = 0, chunks = 0,
    parsed = 0, adapters = 0, assemblies = 0))
  expect_gt(nrow(jbrowse_alignments(bam, "ctgA:1-5,000")$reads), 0)
})

test_that("memory stays bounded over many files, and a reset frees it for reuse", {
  skip_on_cran()
  skip_if_not(file.exists("/proc/self/status"))
  rss <- function() {
    s <- readLines("/proc/self/status")
    as.numeric(sub("\\D+(\\d+).*", "\\1", grep("^VmRSS", s, value = TRUE))) /
      1024
  }
  limit <- jbrowse_cache()$limit
  on.exit(jbrowse_cache(limit = limit))
  jbrowse_cache(limit = 32 * 2^20)
  paths <- bam_copies(160)
  pass <- function() {
    for (i in seq_along(paths)) {
      jbrowse_alignments(paths[i], "ctgA:1-2,000")
      if (i == 20) warm <<- rss()
    }
    rss()
  }
  warm <- NA
  first <- pass()
  expect_lt(first - warm, 128)
  jbrowse_reset()
  expect_lt(pass() - first, 64)
})

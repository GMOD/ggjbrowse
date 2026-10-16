extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
bam <- extdata("volvox-sorted.bam")
loci <- c("ctgA:1,001-3,000", "ctgA:2,001-5,000")

test_that("every linear result holds the same non-table slots", {
  results <- list(
    jbrowse_alignments(bam, loci),
    jbrowse_features(extdata("volvox.sort.gff3.gz"), loci),
    jbrowse_variants(extdata("1000g_dcst2.vcf.gz"), "1:155,000,001-155,010,000"),
    jbrowse_wiggle(extdata("volvox_microarray_posneg.bw"), loci))
  for (x in results) {
    expect_true(all(c("regions", "width", "display", "keys", "ramps",
      "y_scale", "fields", "warnings") %in% names(x)))
    expect_named(x$regions, c("chrom", "start", "end", "offset", "cum_start",
      "cum_end", "missing"))
    expect_equal(x$width, 1500)
    expect_false(any(c("region", "missing_contig", "legend", "arc_legend",
      "settings", "coverage_scale", "domain") %in% names(x)))
    for (name in names(x)) {
      d <- x[[name]]
      if (is.data.frame(d) && name != "regions") expect_null(d$.region)
    }
  }
})

test_that("a read keeps one key in every copy, and its marks name it by that key", {
  a <- jbrowse_alignments(bam, loci)
  copies <- a$reads[a$reads$read %in% a$reads$read[duplicated(a$reads$read)], ]
  expect_gt(nrow(copies), 0)
  expect_true(all(tapply(copies$id, copies$read, function(i)
    length(unique(i))) == 1))
  expect_null(a$mismatches$row)
  expect_equal(a$mismatches$chrom, rep("ctgA", nrow(a$mismatches)))
  expect_null(jbrowse_variants(extdata("1000g_dcst2.vcf.gz"),
    "1:155,000,001-155,010,000")$cells$chrom)
  expect_false(anyNA(of(a$mismatches, a$reads, "read")$read))
})

test_that("assigning reads drops the marks of the reads it leaves out", {
  a <- jbrowse_alignments(bam, loci)
  kept <- a
  kept$reads <- subset(a$reads, strand == 1)
  expect_lt(nrow(kept$mismatches), nrow(a$mismatches))
  expect_true(all(paste(kept$mismatches$read, kept$mismatches$region) %in%
    paste(kept$reads$read, kept$reads$region)))
  forward <- of(a$mismatches, a$reads, "read")$strand == 1
  expect_equal(nrow(kept$mismatches), sum(forward))
  expect_equal(nrow(kept$coverage), nrow(a$coverage))
  expect_equal(attr(kept, "filtered"), "reads")

  first <- a
  first$reads <- a$reads[a$reads$region == 1, ]
  expect_setequal(unique(first$mismatches$region), 1)
  expect_equal(nrow(first$mismatches), sum(a$mismatches$region == 1))
})

test_that("a summary table drawn after its reads were filtered warns once", {
  rlang::local_options(rlib_warning_verbosity = "verbose")
  a <- jbrowse_alignments(bam, loci[1])
  expect_no_warning(ggplot2::ggplot_build(ggplot2::autoplot(a, "coverage")))
  a$reads <- subset(a$reads, strand == 1)
  warned <- character(0)
  withCallingHandlers(ggplot2::ggplot_build(ggplot2::autoplot(a, "coverage")),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("display = list(filter", warned, fixed = TRUE)))
  expect_no_warning(ggplot2::ggplot_build(ggplot2::autoplot(a, "pileup")))
})

test_that("reassigning a read's row moves its marks", {
  a <- jbrowse_alignments(bam, loci[1])
  flipped <- a
  flipped$reads$row <- max(a$reads$row) - a$reads$row
  expect_equal(nrow(flipped$mismatches), nrow(a$mismatches))
  drawn <- ggplot2::layer_data(ggplot2::ggplot(flipped) +
    geom_jbrowse_mismatches())
  expect_equal(sort(unique(drawn$ymin)), sort(unique(max(a$reads$row) -
    on_axis(a)$mismatches$row + read_inset(a))))
})

test_that("a dplyr tibble assigns and draws as a data frame does", {
  skip_if_not_installed("dplyr")
  a <- jbrowse_alignments(bam, loci[1])
  kept <- a
  kept$reads <- dplyr::filter(tibble::as_tibble(a$reads), strand == 1)
  expect_type(kept$reads$color, "character")
  base <- a
  base$reads <- subset(a$reads, strand == 1)
  expect_equal(nrow(kept$mismatches), nrow(base$mismatches))
  expect_s3_class(jbrowse_stack(kept), "patchwork")
  expect_error(geom_jbrowse_reads(data = kept$reads),
    class = "ggjbrowse_parent_data")
})

test_that("assigning variants drops their cells, and as.matrix() follows", {
  v <- jbrowse_variants(extdata("1000g_dcst2.vcf.gz"),
    "1:155,000,001-155,050,000")
  first <- v
  first$variants <- v$variants[1:10, ]
  expect_setequal(unique(first$cells$variant), v$variants$variant[1:10])
  alts <- lengths(strsplit(first$variants$alt, ","))
  m <- as.matrix(first)
  expect_equal(dim(m), c(nrow(v$rows), sum(alts)))
  expect_equal(m, as.matrix(v)[, colnames(m)])
  fewer <- v
  fewer$rows <- v$rows[1:5, ]
  expect_setequal(as.character(unique(fewer$cells$row)), v$rows$name[1:5])
})

test_that("assigning features drops their transcripts and marks", {
  g <- jbrowse_features(extdata("volvox.sort.gff3.gz"), loci)
  eden <- g
  eden$features <- g$features[g$features$name == "EDEN", ]
  expect_setequal(unique(eden$rects$feature), eden$features$feature)
  expect_setequal(unique(eden$subfeatures$feature), eden$features$feature)
  expect_true(all(is.na(eden$labels$subfeature) |
    eden$labels$subfeature %in% eden$subfeatures$subfeature))
  one <- eden
  one$subfeatures <- eden$subfeatures[1, ]
  expect_true(all(is.na(one$rects$subfeature) |
    one$rects$subfeature == one$subfeatures$subfeature))
})

test_that("assigning sources drops the marks of the sources left out", {
  sources <- list(type = "MultiWiggleAdapter", subadapters = lapply(1:2,
    function(i) list(type = "BigWigAdapter", name = paste0("rep", i),
      bigWigLocation = list(uri = normalizePath(
        extdata("volvox_microarray_posneg.bw"))))))
  w <- jbrowse_wiggle(sources, "ctgA:1-10,000")
  w$sources <- w$sources[w$sources$name == "rep2", ]
  expect_setequal(unique(w$marks$source), "rep2")
  expect_setequal(unique(w$bins$source), "rep2")
})

test_that("a saved result redraws after its files are gone", {
  dir <- withr::local_tempdir()
  copy <- function(f) {
    file.copy(extdata(f), dir)
    file.path(dir, f)
  }
  bam <- copy("volvox-sorted.bam")
  copy("volvox-sorted.bam.bai")
  fa <- copy("volvox.fa")
  copy("volvox.fa.fai")
  vcf <- copy("1000g_dcst2.vcf.gz")
  copy("1000g_dcst2.vcf.gz.tbi")
  a <- jbrowse_alignments(bam, c("ctgA:1,001-1,100", "ctgA:2,001-2,050"),
    assembly = fa)
  expect_equal(nchar(a$sequence), c(100, 50))
  expect_equal(a$sequence[1], jbrowse_sequence(fa, "ctgA:1,001-1,100"))
  v <- jbrowse_variants(vcf, "1:155,000,001-155,050,000", width = 300)
  expect_true(all(v$variants$opacity > 0 & v$variants$opacity <= 1))
  expect_lt(min(v$variants$opacity), 1)
  saved <- file.path(dir, "saved.rds")
  saveRDS(list(a = a, v = v), saved)
  unlink(c(bam, fa, vcf, paste0(c(bam, fa, vcf), c(".bai", ".fai", ".tbi"))))
  jbrowse_reset()
  back <- readRDS(saved)
  bases <- ggplot2::layer_data(ggplot2::autoplot(back$a, "sequence"), 1)
  expect_equal(nrow(bases), 150)
  expect_s3_class(ggplot2::ggplotGrob(ggplot2::autoplot(back$a)), "gtable")
  cells <- ggplot2::layer_data(ggplot2::autoplot(back$v), 2)
  expect_setequal(unique(cells$alpha), unique(of(on_axis(v)$cells,
    v$variants, "variant")$opacity))
})

test_that("autoplot and jbrowse_stack draw at a width, warning where narrower than a layout", {
  gff <- extdata("volvox.sort.gff3.gz")
  genes <- jbrowse_features(gff, loci[1], width = 1200)
  reads <- jbrowse_alignments(bam, loci[1], width = 900)
  expect_equal(jbrowse_width(jbrowse_stack(genes, reads)),
    jbrowse_width(jbrowse_stack(genes, reads, width = 1200)))
  expect_gt(jbrowse_width(jbrowse_stack(genes, reads, width = 2400)),
    jbrowse_width(jbrowse_stack(genes, reads)) + 7)
  expect_warning(jbrowse_stack(genes, reads, width = 600),
    "narrower than the 1200 px it was laid out for")
  expect_no_warning(jbrowse_stack(reads, width = 600))
  expect_warning(ggplot2::autoplot(genes, width = 600), "1200 px")
  wide <- ggplot2::autoplot(genes, width = 3000)
  expect_equal(attr(wide$data, "jbrowse")$width, 3000)
  expect_equal(genes$width, 1200)
})

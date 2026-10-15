extdata <- function(f) system.file("extdata", f, package = "ggjbrowse")
people <- extdata("1000g_dcst2.vcf.gz")
dcst2 <- "1:155,000,001-155,050,000"
read_people <- function(display = NULL) {
  jbrowse_variants(people, dcst2, samples = extdata("1000g_dcst2_samples.tsv"),
    display = display)
}

layer_of <- function(p, geom) {
  i <- which(vapply(p$layers, function(l) inherits(l$geom, geom), NA))
  ggplot2::layer_data(p, i[1])
}

test_that("the rows take a saved arrangement: its order, labels and focus", {
  plain <- read_people()
  last <- utils::tail(plain$rows$name, 2)
  v <- read_people(list(rows = list(domain = rev(last),
    labels = stats::setNames(list("Proband"), last[2]))))
  expect_equal(v$rows$name[1:2], rev(last))
  expect_equal(v$rows$label[1], "Proband")
  expect_equal(v$rows$label[-1], v$rows$name[-1])
  p <- ggplot2::autoplot(v, layout = "columns")
  labels <- ggplot2::ggplot_build(p)$data
  expect_equal(labels[[length(labels)]]$label[1], "Proband")
  focused <- read_people(list(rows = list(kept = I(last))))
  expect_setequal(focused$rows$name, last)
  expect_setequal(as.character(unique(focused$cells$row)), last)
})

test_that("a facet stacks the rows in bands that hold whatever order they take", {
  v <- read_people(list(facet = list(field = "superpopulation",
    domain = I("European"))))
  expect_equal(v$bands$label[1], "European")
  expect_equal(rle(v$rows$band)$values, v$bands$key)
  expect_true(all(v$rows$band[1:20] == "European"))
  shuffled <- v
  shuffled$rows <- v$rows[rev(seq_len(nrow(v$rows))), ]
  drawn <- on_axis(shuffled)$rows
  expect_equal(rle(drawn$band)$values, v$bands$key)
  expect_equal(drawn$name[1:20], rev(v$rows$name[1:20]))
  p <- ggplot2::autoplot(v)
  bands <- layer_of(p, "GeomBandLabels")
  expect_equal(bands$label, v$bands$label)
  expect_no_error(ggplot2::ggplotGrob(p))
  expect_null(row_bands(read_people()$rows, read_people()$bands))
})

test_that("under bands a forest draws a clade per band, and a whole tree none", {
  v <- read_people(list(facet = "superpopulation"))
  m <- as.matrix(v)
  v$tree <- lapply(split(rownames(m), v$rows$band), function(r)
    stats::hclust(stats::dist(m[r, ]), "average"))
  p <- ggplot2::autoplot(v, layout = "columns")
  expect_equal(nrow(layer_of(p, "GeomTree")), 3 * (nrow(v$rows) - 5))
  drawn <- on_axis(v)$rows
  expect_equal(rle(drawn$band)$lengths, rep(20, 5))
  expect_no_error(ggplot2::ggplotGrob(p))
  whole <- v
  whole$tree <- stats::hclust(stats::dist(m), "average")
  expect_equal(nrow(layer_of(ggplot2::autoplot(whole), "GeomTree")), 0)
  part <- read_people()
  part$tree <- list(stats::hclust(stats::dist(m[5:10, ])))
  expect_equal(on_axis(part)$rows$name[1:6],
    part$tree[[1]]$labels[part$tree[[1]]$order])
})

test_that("zoomed out, records sharing a pixel fade as the display fades them", {
  v <- read_people()
  alpha <- function(...) {
    p <- ggplot2::ggplot(v) + geom_jbrowse_genotypes(...)
    ggplot2::layer_data(p)$alpha
  }
  far <- alpha(width = 300)
  expect_lt(min(far), 1)
  expect_gt(min(far), 0)
  expect_gt(mean(alpha(width = 1e5)), mean(far))
  expect_true(all(alpha(width = 300, fade = FALSE) == 1))
  expect_true(all(alpha(width = 300, layout = "columns") == 1))
  v$variants <- v$variants[c(1, 1), ]
  v$variants$end[1] <- v$variants$start[1] + 30000
  expect_equal(density_alpha(v$variants, 50)[1], 1)
})

dogs <- extdata("dog10k_cyp1a2_snvs.vcf.gz")
cyp1a2 <- "chr30:38,258,001-38,265,000"

test_that("each record carries its FILTER and INFO, typed by the header", {
  v <- jbrowse_variants(dogs, cyp1a2)$variants
  site <- v[v$start == 38258773, ]
  expect_equal(site$alt, "C,A")
  expect_equal(site$filter, "VQSRTrancheSNP99.00to99.90")
  expect_identical(site$AC[[1]], c(33L, 1L))
  expect_identical(site$AN, 78L)
  expect_identical(site$culprit, "MQ")
  expect_type(v$AF, "list")
  expect_type(v$POSITIVE_TRAIN_SITE, "logical")
  expect_equal(sum(v$POSITIVE_TRAIN_SITE), 5)
  expect_equal(sum(v$filter == "PASS"), 22)
  expect_true(all(is.na(unlist(v$RAW_MQandDP))))
  expect_equal(sum(is.na(v$BaseQRankSum)), 4)
})

test_that("the INFO columns agree with bcftools query", {
  skip_if(Sys.which("bcftools") == "")
  # R's LD_LIBRARY_PATH can shadow the htslib bcftools was linked against
  out <- suppressWarnings(system2("bcftools", c("query", "-f",
    shQuote("%POS\\t%FILTER\\t%INFO/AC\\t%INFO/AN\\t%INFO/MQ\\t%INFO/POSITIVE_TRAIN_SITE\\n"),
    dogs), stdout = TRUE, stderr = FALSE, env = "LD_LIBRARY_PATH="))
  skip_if(!is.null(attr(out, "status")), "bcftools does not run")
  query <- utils::read.table(text = out, sep = "\t", na.strings = ".",
    col.names = c("pos", "filter", "AC", "AN", "MQ", "flag"),
    colClasses = c("integer", "character", "character", "integer",
      "numeric", "character"))
  v <- jbrowse_variants(dogs, cyp1a2)$variants
  expect_equal(v$start + 1, query$pos)
  expect_equal(v$filter, query$filter)
  expect_equal(vapply(v$AC, paste, "", collapse = ","), query$AC)
  expect_equal(v$AN, query$AN)
  expect_equal(v$MQ, query$MQ)
  expect_equal(v$POSITIVE_TRAIN_SITE, !is.na(query$flag))
})

test_that("an INFO key named like a record column takes an info_ prefix", {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2",
    '##INFO=<ID=type,Number=1,Type=String,Description="Named like a column">',
    '##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
      "FORMAT", "s1", sep = "\t"),
    paste("ctgA", "100", ".", "A", "G", "50", ".", "type=snv;DB;EXTRA=x,y",
      "GT", "0/1", sep = "\t"),
    paste("ctgA", "200", ".", "C", "T", "50", "q10;s50", "type=other", "GT",
      "1/1", sep = "\t")), vcf)
  v <- jbrowse_variants(vcf, "ctgA:1-1,000")
  expect_equal(v$info_fields$key, c("type", "DB", "EXTRA"))
  expect_equal(v$info_fields$column, c("info_type", "DB", "EXTRA"))
  expect_equal(v$variants$type, c("SNV", "SNV"))
  expect_equal(v$variants$info_type, c("snv", "other"))
  expect_equal(v$variants$filter, c(NA, "q10;s50"))
  expect_equal(v$variants$DB, c(TRUE, FALSE))
  expect_equal(v$variants$EXTRA, list(c("x", "y"), NA_character_))
  apart <- jbrowse_variants(vcf, data.frame(chrom = "ctgA", start = c(0, 150),
    end = c(150, 1000)))
  expect_equal(apart$variants$EXTRA, v$variants$EXTRA)
  expect_equal(apart$info_fields$key, c("type", "DB", "EXTRA"))
})

test_that("an INFO column keeps its type with no records or no values", {
  expect_equal(info_column(list(), "A", "Integer"), list())
  expect_identical(info_column(c(1L, NA), "1", "Integer"), c(1L, NA))
  expect_identical(info_column(list(c(1L, 2L), NA), "A", "Integer"),
    list(c(1L, 2L), NA_integer_))
  expect_identical(info_column(c(NA, NA), "R", "Float"),
    list(NA_real_, NA_real_))
})

test_that("a list column stays a list whatever lengths its rows share", {
  v <- jbrowse_variants(dogs, "chr30:38,259,001-38,259,200")$variants
  expect_true(all(lengths(v$AC) == 1))
  expect_type(v$AC, "list")
  expect_identical(v$AC[[1]], v$AC[[1]][1])
})

test_that("two regions bind their INFO columns, each record once", {
  halves <- data.frame(chrom = "chr30", start = c(38258000, 38262000),
    end = c(38262000, 38265000))
  both <- jbrowse_variants(dogs, halves)
  whole <- jbrowse_variants(dogs, cyp1a2)
  expect_equal(both$variants$AC, whole$variants$AC)
  expect_equal(both$variants$filter, whole$variants$filter)
  expect_equal(both$info_fields, whole$info_fields)
})

test_that("a first region with no records leaves every row and the second region's cells", {
  vcf <- system.file("extdata", "1000g_dcst2.vcf.gz", package = "ggjbrowse")
  both <- jbrowse_variants(vcf, c("1:1,001-2,000", "1:155,000,001-155,050,000"))
  alone <- jbrowse_variants(vcf, "1:155,000,001-155,050,000")
  expect_equal(both$rows$name, alone$rows$name)
  expect_equal(nrow(both$cells), nrow(alone$cells))
  expect_true(all(both$variants$region == 2))
  drawn <- ggplot2::layer_data(ggplot2::autoplot(both), 2)
  expect_equal(nrow(drawn), nrow(alone$cells))
})

test_that("a VCF with no region reads whole", {
  whole <- jbrowse_variants(dogs)
  expect_equal(whole$variants$start,
    jbrowse_variants(dogs, "chr30:1-40,000,000")$variants$start)
  expect_gt(nrow(whole$cells), 1000)
})

test_that("the result keeps the display's resolved unit, row height and config", {
  v <- jbrowse_variants(people, dcst2)
  expect_equal(v$display$type, "LinearMultiSampleVariantDisplay")
  expect_equal(v$unit, "sample")
  expect_equal(v$row_height, 0)
  phased <- jbrowse_variants(people, dcst2,
    display = list(unit = "haplotype", rowHeight = 4))
  expect_equal(phased$unit, "haplotype")
  expect_equal(phased$row_height, 4)
})

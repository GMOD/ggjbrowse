alignments_adapter <- function(file, index) {
  if (grepl("\\.cram$", file, ignore.case = TRUE)) {
    list(type = "CramAdapter", cramLocation = location(file),
      craiLocation = location(index %||% paste0(file, ".crai")))
  } else {
    index <- index %||% paste0(file, ".bai")
    list(type = "BamAdapter", bamLocation = location(file),
      index = list(location = location(index),
        indexType = if (grepl("\\.csi$", index)) "CSI" else "BAI"))
  }
}

sequence_adapter <- function(reference) {
  if (grepl("\\.2bit$", reference, ignore.case = TRUE)) {
    list(type = "TwoBitAdapter", twoBitLocation = location(reference))
  } else {
    list(type = "IndexedFastaAdapter", fastaLocation = location(reference),
      faiLocation = location(paste0(reference, ".fai")))
  }
}

#' Read a region of a BAM or CRAM file as JBrowse draws it
#'
#' Runs JBrowse's own alignments engine in V8: the reads it fetches, the
#' pileup rows its display packs them into, their mismatches against the
#' reference, their deletions and splice gaps, and the coverage. Coordinates
#' are 0-based, half-open, on the reference.
#'
#' @param file a BAM or CRAM file, a path or a URL
#' @param chrom,start,end the region, 0-based half-open
#' @param reference the reference FASTA (indexed) or 2bit, path or URL; a CRAM
#'   needs it, and a BAM read without an MD tag is compared against it
#' @param index the BAI, CSI or CRAI, if not beside `file` under its usual name
#' @return a `jbrowse_alignments` list of data frames:
#'   * `reads`: `start`, `end`, `row`, `strand`, `mapq`, `flags`, `name`,
#'     `insert_size`, `pair_orientation`
#'   * `mismatches`: `position`, `base`, `read`, `row`, `frequency` (that base's
#'     share of the depth there), `quality`
#'   * `gaps`: `start`, `end`, `type` (`"deletion"` or `"skip"`), `read`, `row`
#'   * `coverage`: `start`, `end`, `depth`, `forward`, `reverse`
#'
#'   `read` is a row number of `reads`.
#' @export
#' @examples
#' bam <- system.file("extdata", "volvox-sorted.bam", package = "ggjbrowse")
#' a <- jbrowse_alignments(bam, "ctgA", 1000, 3000)
#' head(a$reads)
jbrowse_alignments <- function(file, chrom, start, end, reference = NULL,
                               index = NULL) {
  args <- list(
    adapterConfig = alignments_adapter(file, index),
    region = list(refName = chrom, start = start, end = end))
  if (!is.null(reference)) args$sequenceAdapter <- sequence_adapter(reference)
  t <- call_engine("alignments", args)
  reads <- as_frame(t$reads)
  names(reads)[names(reads) == "insertSize"] <- "insert_size"
  names(reads)[names(reads) == "pairOrientation"] <- "pair_orientation"
  on_read <- function(table) {
    d <- as_frame(table)
    d$read <- d$read + 1
    d$row <- reads$row[d$read]
    d
  }
  structure(
    list(reads = reads, mismatches = on_read(t$mismatches),
      gaps = on_read(t$gaps), coverage = as_frame(t$coverage),
      region = list(chrom = chrom, start = start, end = end)),
    class = "jbrowse_alignments")
}

#' @export
print.jbrowse_alignments <- function(x, ...) {
  cat(sprintf(
    "<jbrowse_alignments> %s:%.0f-%.0f, %d reads on %d rows, %d mismatches, %d gaps\n",
    x$region$chrom, x$region$start, x$region$end, nrow(x$reads),
    length(unique(x$reads$row)), nrow(x$mismatches), nrow(x$gaps)))
  invisible(x)
}

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
#' @param reference_chrom the contig's name in `reference`, where the file
#'   spells it another way (`1` in the BAM, `chr1` in the FASTA)
#' @param index the BAI, CSI or CRAI, if not beside `file` under its usual name
#' @param color what fills each read, as the display's `color` object: a
#'   field name (`"strand"`, `"mapq"`, `"insertSize"`, `"pairOrientation"`,
#'   `"firstOfPairStrand"`, `"tags.XS"`) or a list of the object's slots, such
#'   as `list(field = "insertSize", domain = c(200, 600))` or
#'   `list(field = "strand", range = c("#d33", "#33d"))`; `NULL` fills every
#'   read plain
#' @param filter which reads are kept, as the display's `filter` object: a
#'   list with any of `flagInclude`, `flagExclude` (1540 by default), `readName`,
#'   `tagFilters` (a list of `list(tag = , value = )`), `spliced`,
#'   `properPairs`, `singletons` and `split` (each `"only"` or `"exclude"`)
#' @param base_color the per-base layer painted over the reads, as the
#'   display's `baseColor` object: `"modifications"` for the reads' base
#'   modifications (MM/ML tags), or a list of the object's slots; `NULL`
#'   paints none
#' @param modifications the settings that layer reads, as the display's
#'   `modifications` object: a list with any of `threshold` (the percent
#'   probability under which a call is hidden, 10), `unmodified` (`"hidden"`,
#'   `"calls"` or `"all"`: which unmodified sites are painted too),
#'   `cytosineContext` (`"CG"`, `"CHG"`, `"CHH"` or `"all"`) and
#'   `shownModifications` (the type codes drawn, every one when empty)
#' @param byte_limit the most compressed bytes the region may need, refused off
#'   the index before any read is fetched; the whole table is held in memory
#' @return a `jbrowse_alignments` list of data frames:
#'   * `reads`: `start`, `end`, `row`, `strand`, `mapq`, `flags`, `name`,
#'     `insert_size`, `pair_orientation`, `color` (the fill the display paints,
#'     as `#rrggbb`) and `category` (the bucket the display files the read
#'     under: `fwdStrand`, `mapq`, `unmappedMate`, ...)
#'   * `mismatches`: `position`, `base`, `read`, `row`, `frequency` (how many
#'     reads carry that base there, over the depth), `quality`, `fade`
#'   * `gaps`: `start`, `end`, `type` (`"deletion"` or `"skip"`), `read`, `row`,
#'     `fade`
#'   * `insertions`: `position`, `length`, `sequence`, `read`, `row`, `fade`
#'   * `clips`: `position`, `length`, `type` (`"soft"` or `"hard"`), `read`,
#'     `row`, `fade`
#'   * `modifications`: each base modification call the layer paints:
#'     `position`, `read`, `row`, `type` (the code, `m` for 5mC), `name`,
#'     `probability`, `unmodified` (whether the call is that the base is not
#'     modified), and the `color` and `opacity` the display paints it at;
#'     empty without `base_color`
#'   * `coverage`: `start`, `end`, `depth`, `forward`, `reverse`
#'   * `junctions`: the reads' splice gaps merged, each once: `start`, `end`,
#'     `count` (the reads supporting it), `strand` (1, -1, or 0 where neither
#'     the reads' tags nor the splice motif says) and `motif` (`"GT-AG
#'     (canonical)"`, `"non-canonical"`, or `""` where the region's sequence
#'     does not cover both ends)
#'   * `legend`: `label` and `color`, the key the display shows for `color`,
#'     its title as the frame's `"title"` attribute
#'
#'   `read` is a row number of `reads`. `fade` is the share of the depth a
#'   mark's site holds, zeroed under the display's depth-dependent threshold,
#'   which the pileup fades a sub-pixel mark by.
#' @export
#' @examples
#' bam <- system.file("extdata", "volvox-sorted.bam", package = "ggjbrowse")
#' a <- jbrowse_alignments(bam, "ctgA", 1000, 3000, color = "strand")
#' head(a$reads)
#' a$legend
jbrowse_alignments <- function(file, chrom, start, end, reference = NULL,
                               reference_chrom = chrom, index = NULL,
                               color = NULL, filter = NULL, base_color = NULL,
                               modifications = NULL, byte_limit = 2e8) {
  args <- list(
    adapterConfig = alignments_adapter(file, index),
    region = list(refName = chrom, originalRefName = reference_chrom,
      start = start, end = end),
    byteLimit = byte_limit)
  if (!is.null(reference)) args$sequenceAdapter <- sequence_adapter(reference)
  if (!is.null(color)) args$color <- color_object(color)
  if (!is.null(filter)) args$filter <- filter
  if (!is.null(base_color)) {
    args$baseColor <- if (is.character(base_color)) list(field = base_color)
      else base_color
  }
  if (!is.null(modifications)) {
    if (!is.null(modifications$shownModifications)) {
      modifications$shownModifications <- I(modifications$shownModifications)
    }
    args$modifications <- modifications
  }
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
  legend <- as_frame(t$legend)
  attr(legend, "title") <- t$legendTitle
  structure(
    list(reads = reads, mismatches = on_read(t$mismatches),
      gaps = on_read(t$gaps), insertions = on_read(t$insertions),
      clips = on_read(t$clips), modifications = on_read(t$modifications),
      coverage = as_frame(t$coverage), junctions = as_frame(t$junctions),
      legend = legend,
      region = list(chrom = chrom, start = start, end = end),
      reference = if (is.null(reference)) NULL else
        list(file = reference, chrom = reference_chrom)),
    class = "jbrowse_alignments")
}

#' Read a region of the reference
#'
#' The bases of a region, read through JBrowse's own FASTA or 2bit reader.
#' @param reference the reference FASTA (indexed) or 2bit, path or URL
#' @inheritParams jbrowse_alignments
#' @return one string of `end - start` bases
#' @export
#' @examples
#' fa <- system.file("extdata", "volvox.fa", package = "ggjbrowse")
#' jbrowse_sequence(fa, "ctgA", 1000, 1020)
jbrowse_sequence <- function(reference, chrom, start, end) {
  engine()$call("jbrowse.sequence", as.character(jsonlite::toJSON(list(
    sequenceAdapter = sequence_adapter(reference),
    region = list(refName = chrom, start = start, end = end)),
    auto_unbox = TRUE)), await = TRUE)
}

# The display's color object from a field name or a list of its slots; a
# slot that is a list stays one in JSON, however many values it holds.
color_object <- function(color) {
  if (is.character(color) && length(color) == 1) return(list(field = color))
  for (slot in intersect(c("domain", "range", "labels"), names(color))) {
    color[[slot]] <- I(color[[slot]])
  }
  color
}

#' @export
print.jbrowse_alignments <- function(x, ...) {
  cat(sprintf(
    "<jbrowse_alignments> %s:%.0f-%.0f, %d reads on %d rows, %d mismatches, %d gaps, %d insertions, %d clips, %d modifications, %d junctions\n",
    x$region$chrom, x$region$start, x$region$end, nrow(x$reads),
    length(unique(x$reads$row)), nrow(x$mismatches), nrow(x$gaps),
    nrow(x$insertions), nrow(x$clips), nrow(x$modifications),
    nrow(x$junctions)))
  invisible(x)
}

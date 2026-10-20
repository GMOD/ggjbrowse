#' What a reader returns
#'
#' [jbrowse_alignments()], [jbrowse_features()], [jbrowse_variants()],
#' [jbrowse_wiggle()], [jbrowse_synteny()], [jbrowse_dotplot()] and
#' [jbrowse_multiway()] each return a classed list, `a$reads` a data frame. A
#' fixed set of slots holds what is not a table; every other slot is a table
#' of records or of marks, which each reader's page lists, or a class's
#' drawing state, listed below.
#'
#' `$` matches a slot's name exactly, so `x$region` is `NULL`.
#'
#' @section The slots every result holds:
#' * `regions`: a data frame, a row per region read: `chrom`, `start` and
#'   `end` (0-based, half-open), where the axis lays it (`offset`, the
#'   position its start takes, and `cum_start` and `cum_end`, its ends), and
#'   `missing`, the reason its contig reads empty, `""` where the file has the
#'   contig
#' * `width`: the pixels the regions were laid out across
#' * `display`: the display's config as JBrowse resolves it, every slot
#'   under JBrowse's own name, the schema's defaults filled in
#' * `keys`: the display's color keys, each a data frame of `label` and
#'   `color` titled by its `"title"` attribute, as [guides_jbrowse()]
#'   draws them
#' * `ramps`: the display's color bars, as [guides_jbrowse()] draws them
#' * `y_scale`: the score axis of a [jbrowse_wiggle()] result, or the
#'   coverage band's depth axis of a [jbrowse_alignments()] result, as the
#'   display autoscales it: `domain`, `type` (`"linear"`, `"log"` or
#'   `"symlog"`), `symlog_constant` and `ticks`; `NULL` for the others
#' * `fields`: the file's declarations of the columns it adds, each a data
#'   frame: `info` and `format` (a VCF's header lines), `attributes` (a GFF3's
#'   or BigBed's fields), and `annotations` (the fields of a VCF's `ANN` or
#'   `CSQ`)
#' * `warnings`: what the engine warned of while reading, each also raised as
#'   an R warning
#'
#' A [jbrowse_features()] result also holds `glyph`, `height` and
#' `bp_per_px`, the layout's pixel sizes; a [jbrowse_wiggle()] result
#' `rendering`, `origin`, `row_layout`, `height` and `bp_per_px`; a
#' [jbrowse_alignments()] result `read_pitch`, `glyph`, `arc_span_max` and
#' `sequence`, each region's reference bases where `assembly` has a sequence,
#' which the sequence band draws; and a [jbrowse_variants()] result the
#' `tree` a caller sets.
#'
#' A comparative result lays out several axes: [jbrowse_synteny()]'s two
#' rows, [jbrowse_dotplot()]'s x and y, [jbrowse_multiway()]'s lanes. Its
#' `regions` adds `axis`, naming the axis each region lays out on, and
#' `assembly`, the assembly it is on, and its tables name their regions by
#' axis: `query_region` and `target_region`, `x_region` and `y_region`, or
#' `upper_region` and `lower_region`, each a row of `regions`. A
#' [jbrowse_synteny()] result also holds `offscreen_strip`, `fade_thin`,
#' `has_cigar` and `assemblies`, the axis titles; a [jbrowse_dotplot()]
#' result `height`, `assemblies` and `skipped`; and a [jbrowse_multiway()]
#' result `tick_interval` and `rows_vs_anchor`. `y_scale` is `NULL` and
#' `fields` empty for each. A multi-genome file lists no contigs, so a
#' [jbrowse_multiway()] result's `missing` is always `""`.
#'
#' @section Errors:
#' Each error a reader raises has a class, which a [tryCatch()] handler
#' matches:
#' * `ggjbrowse_too_large`: the regions pass `span_limit` or `byte_limit`,
#'   which the error names, or a whole file or contig passes 5 million rows
#'   in a table
#' * `ggjbrowse_missing_contig`: no region's contig is one the file lists. A
#'   region on such a contig beside others on the file's contigs reads empty
#'   and warns, its reason in `regions$missing`
#' * `ggjbrowse_bad_display`: a `display` slot holds a value the display's
#'   config schema refuses, or a synteny view's `cigarMode` one it does not
#'   take, which the error names
#' * `ggjbrowse_bad_region`: a region is neither a locus nor a contig of the
#'   file, or a data frame of regions reads as 1-based
#' * `ggjbrowse_forked`: the engine was started in a forked worker
#' @name jbrowse_result
NULL

# Each class's parent tables: the column naming each record (`key`), the
# tables naming a record by it (in `child_key`, `key` where unset), and the
# parent's columns a child draws with, joined from the parent at draw time
FAMILIES <- list(
  jbrowse_alignments = list(
    reads = list(key = "read", columns = c("row", "section"),
      children = c("mismatches", "gaps", "insertions", "clips",
        "modifications", "qualities"))),
  jbrowse_features = list(
    features = list(key = "feature",
      columns = c("name", "type", "strand", "start", "end", "section"),
      children = c("subfeatures", "rects", "lines", "chevrons", "arrows",
        "amino_acids", "labels", "samples", "annotations")),
    subfeatures = list(key = "subfeature", columns = character(0),
      children = c("rects", "lines", "chevrons", "arrows", "amino_acids",
        "labels"))),
  jbrowse_variants = list(
    variants = list(key = "variant", columns = c("start", "end", "opacity"),
      children = c("cells", "annotations")),
    rows = list(key = "name", child_key = "row", columns = "sample",
      children = "cells")),
  jbrowse_wiggle = list(
    sources = list(key = "name", child_key = "source", columns = character(0),
      children = c("marks", "bins"))),
  jbrowse_synteny = list(
    alignments = list(key = "alignment", columns = character(0),
      children = "ribbons")),
  jbrowse_dotplot = list(
    alignments = list(key = "alignment", columns = character(0),
      children = "segments")),
  jbrowse_multiway = list(
    alignments = list(key = "alignment", columns = character(0),
      children = "ribbons")))

# The tables an alignments track summarizes its reads in, which the engine
# computed over every read and a filter in R leaves as they were
SUMMARIES <- c("coverage", "allele_coverage", "modification_coverage",
  "interbase_coverage", "indicators", "junctions", "chains", "overlaps",
  "arcs", "arc_ticks")

families_of <- function(x) {
  for (cls in class(x)) if (!is.null(FAMILIES[[cls]])) return(FAMILIES[[cls]])
  list()
}

# The row of `parent` each row of `d` names: the copy in the same region
# where both list regions, as a record read in two regions is listed in each
parent_rows <- function(d, family, parent) {
  key <- d[[family$child_key %||% family$key]]
  own <- parent[[family$key]]
  if (is.numeric(d$region) && is.numeric(parent$region)) {
    if (!is.numeric(key) || !is.numeric(own)) {
      return(match(paste(key, d$region), paste(own, parent$region)))
    }
    n <- max(0, d$region, parent$region) + 1
    match(key * n + d$region, own * n + parent$region)
  } else {
    match(as.character(key), as.character(own))
  }
}

# Each child table of `x`, an unclassed result, with the columns it draws
# with that it lacks, from its parent's row for the same record
joined <- function(x, families) {
  for (family_name in names(families)) {
    family <- families[[family_name]]
    parent <- x[[family_name]]
    if (!is.data.frame(parent) || !length(family$columns)) next
    for (child in intersect(family$children, names(x))) {
      d <- x[[child]]
      need <- intersect(setdiff(family$columns, names(d)), names(parent))
      if (!is.data.frame(d) || !length(need)) next
      at <- parent_rows(d, family, parent)
      for (column in need) d[[column]] <- parent[[column]][at]
      x[[child]] <- d
    }
  }
  x
}

#' @export
`$.jbrowse_result` <- function(x, name) .subset2(x, name, exact = TRUE)

#' @export
`$<-.jbrowse_result` <- function(x, name, value) {
  x[[name]] <- value
  x
}

# Assigning a parent table drops the rows of its children whose record it no
# longer holds, and marks the result `filtered` where it lost records, so
# every table, drawn or analysed, holds the same records
#' @export
`[[<-.jbrowse_result` <- function(x, i, value) {
  cls <- class(x)
  old <- if (is.character(i) && i %in% names(x)) .subset2(x, i)
  x <- unclass(x)
  x[[i]] <- value
  class(x) <- cls
  family <- if (is.character(i)) families_of(x)[[i]]
  if (is.null(family) || !is.data.frame(value)) return(x)
  if (is.data.frame(old) && lost(old, value, family$key)) {
    attr(x, "filtered") <- union(attr(x, "filtered"), i)
  }
  for (child in intersect(family$children, names(x))) {
    d <- .subset2(x, child)
    if (!is.data.frame(d)) next
    key <- d[[family$child_key %||% family$key]]
    keep <- !is.na(parent_rows(d, family, value)) | is.na(key)
    if (!all(keep)) x[[child]] <- d[keep, , drop = FALSE]
  }
  x
}

# Whether `new` lacks a record, or a record's copy in a region, `old` held
lost <- function(old, new, key) {
  id <- function(d) {
    if (is.numeric(d$region)) paste(d[[key]], d$region)
    else as.character(d[[key]])
  }
  !all(id(old) %in% id(new))
}

# A summary table as a geom draws it, warning once where the reads it
# summarizes were filtered after the read
summary_table <- function(x, name) {
  if ("reads" %in% attr(x, "filtered") && nrow(x[[name]])) {
    rlang::warn(paste("the coverage, junctions and other summary tables",
      "count every read the engine read, not the reads filtered in R: read",
      "again with `display = list(filter = ...)` to recompute them"),
      .frequency = "once", .frequency_id = "ggjbrowse_summary_tables")
  }
  x[[name]]
}

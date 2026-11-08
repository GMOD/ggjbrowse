#' A row display's row axis: what a row is, which rows show, and how they read
#'
#' The `rows` setting of a display that draws a row per sample or source, as
#' a track config's `displays` entry writes it, for a reader's `display`:
#' `display = list(rows = jbrowse_rows(field = "haplotype"))`. A plain list
#' works as well; this names each member and writes a named vector as the map
#' JBrowse reads. A row is named by its `name` in the reader's
#' `rows` table, a haplotype row `"<sample> HP1"`, `"HP2"`, where a sample's
#' name also names its haplotypes.
#'
#' @param field what one row is: on [jbrowse_variants()], `"sample"` (the
#'   default) or `"haplotype"`; on a wiggle display, the field whose values
#'   each take a row, such as `"source"`
#' @param sort row names, the rows that lead, in order; the rest follow in
#'   the file's order
#' @param kept row names, the rows shown; a name matching no row draws no
#'   rows, and none listed shows every row
#' @param filter attribute tests the rows shown pass, each as a display's
#'   `filter` entry writes it, such as `list(field = "superpopulation", oneOf
#'   = c("EUR", "AFR"))`, or a list of them, all of which a row passes. It
#'   tests a samples file's columns, before `kept`; a row's name goes in
#'   `kept`
#' @param labels a named vector, the label a row shows by row name; a
#'   sample's label names its haplotypes too, as `"Mother HP1"`
#' @param colors a named vector, a CSS color by row name, over what `color`
#'   gives the row
#' @param color the row attribute the rows are colored by, such as a
#'   samples file column, or a list of `field` and its color `scale`
#' @param step the row height in px, or `"fit"`, the default, the rows
#'   dividing the display's height
#' @param tree a newick tree of the rows, labelled by row name, or a list of
#'   its `newick` and `provenance`
#' @param axis the row axis's guide: a list of any of `tree` and `labels`
#'   (`TRUE` or `FALSE`, whether the sidebar draws each), `width`, and `grid`
#' @return the `rows` list, unset members left out
#' @export
#' @examples
#' vcf <- system.file("extdata", "1000g_dcst2.vcf.gz", package = "ggjbrowse")
#' v <- jbrowse_variants(vcf, "1:155,000,001-155,050,000",
#'   display = list(rows = jbrowse_rows(field = "haplotype", step = 4,
#'     colors = c(HG01879 = "crimson"))))
#' head(v$rows)
jbrowse_rows <- function(field = NULL, sort = NULL, kept = NULL,
                         filter = NULL, labels = NULL, colors = NULL,
                         color = NULL, step = NULL,
                         tree = NULL, axis = NULL) {
  check_row_map(labels, "labels")
  check_row_map(colors, "colors")
  if (!is.null(step) && !identical(step, "fit") &&
    !(is.numeric(step) && length(step) == 1 && step > 0)) {
    stop("`step` is a row height in px above 0, or \"fit\"", call. = FALSE)
  }
  rows <- list(field = field, sort = name_list(sort), kept = name_list(kept),
    filter = if ("field" %in% names(filter)) list(filter) else filter,
    labels = as.list(labels), colors = as.list(colors),
    color = if (is.character(color)) list(field = color) else color,
    step = step, tree = tree, axis = axis)
  rows[!vapply(rows, function(m) is.null(m) || identical(m, list()), NA)]
}

check_row_map <- function(map, slot) {
  if (length(map) && (is.null(names(map)) || any(names(map) == ""))) {
    stop("`", slot, "` names each row it sets, as `c(HG002 = ...)`",
      call. = FALSE)
  }
}

name_list <- function(names) {
  if (length(names)) I(as.character(names))
}

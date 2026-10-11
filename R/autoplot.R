#' The colors JBrowse's alignments display draws in
#'
#' `read` (the plain fill), `forward` and `reverse` (each strand), `skip`,
#' `deletion`, `coverage` and `bases` (a named vector, one per mismatch base),
#' read from the engine's own palette.
#' @export
jbrowse_colors <- local({
  colors <- NULL
  function() {
    if (is.null(colors)) {
      read <- jsonlite::fromJSON(engine()$eval("jbrowse.alignmentColors()"))
      read$bases <- unlist(read$bases)
      colors <<- read
    }
    colors
  }
})

# Each read as its aligned blocks: the read split at its splice skips
read_blocks <- function(x) {
  skips <- x$gaps[x$gaps$type == "skip", , drop = FALSE]
  reads <- x$reads
  reads$read <- seq_len(nrow(reads))
  if (!nrow(skips)) return(reads)
  by_read <- split(skips, skips$read)
  spliced <- reads$read %in% names(by_read)
  blocks <- lapply(which(spliced), function(i) {
    s <- by_read[[as.character(i)]]
    s <- s[order(s$start), ]
    starts <- c(reads$start[i], s$end)
    ends <- c(s$start, reads$end[i])
    keep <- ends > starts
    r <- reads[rep(i, sum(keep)), , drop = FALSE]
    r$start <- starts[keep]
    r$end <- ends[keep]
    r
  })
  do.call(rbind, c(list(reads[!spliced, , drop = FALSE]), blocks))
}

#' Draw alignments as JBrowse's pileup or coverage band
#'
#' @param object a [jbrowse_alignments()] result
#' @param what `"pileup"` (each read on its row, split at its splice gaps, with
#'   its deletions and mismatches, filled as the display's `color` object
#'   paints it and keyed as the display keys it) or `"coverage"`
#' @param ... unused
#' @return a ggplot, to extend with layers, scales and themes or to stack with
#'   patchwork
#' @exportS3Method ggplot2::autoplot
#' @examples
#' bam <- system.file("extdata", "volvox-sorted.bam", package = "ggjbrowse")
#' a <- jbrowse_alignments(bam, "ctgA", 1000, 3000, color = "strand")
#' autoplot(a)
autoplot.jbrowse_alignments <- function(object, what = c("pileup", "coverage"),
                                        ...) {
  what <- match.arg(what)
  colors <- jbrowse_colors()
  xlab <- sprintf("%s (bp)", object$region$chrom)
  xlim <- c(object$region$start, object$region$end)
  if (what == "coverage") {
    return(ggplot2::ggplot(object$coverage) +
      ggplot2::geom_rect(ggplot2::aes(xmin = start, xmax = end, ymin = 0,
        ymax = depth), fill = colors$coverage) +
      ggplot2::coord_cartesian(xlim = xlim, expand = FALSE) +
      ggplot2::labs(x = xlab, y = "depth") +
      ggplot2::theme_minimal())
  }
  blocks <- read_blocks(object)
  legend <- object$legend
  keyed <- legend[!duplicated(legend$color), , drop = FALSE]
  gaps <- object$gaps
  ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes(xmin = start, xmax = end, ymin = row + 0.1,
      ymax = row + 0.9, fill = color), blocks) +
    ggplot2::scale_fill_identity(name = attr(legend, "title"),
      breaks = keyed$color, labels = keyed$label, guide = "legend") +
    ggplot2::geom_segment(ggplot2::aes(x = start, xend = end, y = row + 0.5,
      yend = row + 0.5), gaps[gaps$type == "skip", , drop = FALSE],
      colour = colors$skip, linewidth = 0.3) +
    ggplot2::geom_rect(ggplot2::aes(xmin = start, xmax = end, ymin = row + 0.35,
      ymax = row + 0.65), gaps[gaps$type == "deletion", , drop = FALSE],
      fill = colors$deletion) +
    ggplot2::geom_rect(ggplot2::aes(xmin = position, xmax = position + 1,
      ymin = row + 0.1, ymax = row + 0.9, colour = base), object$mismatches,
      fill = "transparent", linewidth = 0.4, key_glyph = "path") +
    ggplot2::scale_colour_manual(values = colors$bases, name = "Mismatch") +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_cartesian(xlim = xlim, expand = FALSE) +
    ggplot2::labs(x = xlab, y = NULL) +
    ggplot2::theme_minimal()
}

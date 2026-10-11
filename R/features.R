# A feature file's adapter, read off its extension: GFF3, GTF and BED plain or
# bgzipped with a tabix (or CSI) index beside it, or BigBed.
features_adapter <- function(file, index) {
  name <- tolower(sub("\\?.*$", "", file))
  tabix <- function(type, slot) {
    index <- index %||% paste0(file, ".tbi")
    a <- list(type = type, index = list(location = location(index),
      indexType = if (grepl("\\.csi$", index)) "CSI" else "TBI"))
    a[[slot]] <- location(file)
    a
  }
  plain <- function(type, slot) {
    a <- list(type = type)
    a[[slot]] <- location(file)
    a
  }
  if (grepl("\\.(bb|bigbed)$", name)) return(plain("BigBedAdapter", "bigBedLocation"))
  if (grepl("\\.gff3?\\.gz$", name)) return(tabix("Gff3TabixAdapter", "gffGzLocation"))
  if (grepl("\\.gff3?$", name)) return(plain("Gff3Adapter", "gffLocation"))
  if (grepl("\\.gtf\\.gz$", name)) return(tabix("GtfTabixAdapter", "gtfGzLocation"))
  if (grepl("\\.gtf$", name)) return(plain("GtfAdapter", "gtfLocation"))
  if (grepl("\\.bed\\.gz$", name)) return(tabix("BedTabixAdapter", "bedGzLocation"))
  if (grepl("\\.bed$", name)) return(plain("BedAdapter", "bedLocation"))
  stop("cannot tell the format of ", file,
    ": expected .gff3, .gtf or .bed (each optionally .gz with an index) or .bb",
    call. = FALSE)
}

#' Read a region of a feature file as JBrowse draws it
#'
#' Runs JBrowse's own feature display in V8: the features it fetches, the gene
#' glyph it builds from each (transcripts, exons, CDS, UTRs, intron lines,
#' strand chevrons and arrows), the rows its layout packs them into, the
#' labels it shows at this zoom and the colors its `color` object paints.
#' Coordinates are 0-based, half-open, on the reference; vertical positions
#' are the display's pixels, down from the track's top.
#'
#' @param file a GFF3, GTF or BED file (plain, or bgzipped with a tabix or CSI
#'   index beside it) or a BigBed, a path or a URL
#' @param chrom,start,end the region, 0-based half-open
#' @param display the display's settings, as a track config's `displays`
#'   entry writes them for a `LinearBasicDisplay`: a list with any of
#'   `color` (a field name such as `"type"` or `"score"`, or a list of the
#'   color object's slots), `displayMode` (`"normal"`, `"compact"`,
#'   `"superCompact"`, `"collapsed"`), `showLabels` (`"auto"`, `"name"`,
#'   `"nameAndDescription"`, `"description"`, `"none"`), `geneGlyphMode`
#'   (`"auto"`, `"all"`, `"longestCoding"`), `subfeatureLabels` (`"none"`,
#'   `"below"`, `"overlay"`), `featureHeight`, `showOnlyGenes`, `filter` (jexl
#'   expressions) and `displayDirectionalChevrons`
#' @param width the width in pixels the display lays the region out across:
#'   whether a name fits beside its feature, which chevrons and arrows draw and
#'   how the `auto` modes resolve all follow from it
#' @param reference the reference FASTA (indexed) or 2bit, path or URL, for
#'   `color_by_cds`
#' @param reference_chrom the contig's name in `reference`, where the file
#'   spells it another way
#' @param color_by_cds paint each CDS by its reading frame, as the view's
#'   "Color by CDS" does
#' @param index the tabix or CSI index, if not beside `file` under its usual
#'   name
#' @param byte_limit the most compressed bytes the region may need, refused
#'   off the index before any feature is fetched
#' @return a `jbrowse_features` list of data frames:
#'   * `features`: each top-level feature on its row: `id`, `name`, `type`,
#'     `start`, `end`, `strand`, `top` and `bottom` (its row's pixels, labels
#'     included) and `gene`
#'   * `subfeatures`: a gene's transcripts: `id`, `name`, `type`, `start`,
#'     `end`, `top`, `bottom`, `feature`
#'   * `rects`: the boxes the glyph draws (exons, CDS, UTRs, a plain feature):
#'     `start`, `end`, `top`, `height`, `color`, `feature`
#'   * `lines`: the intron and connector lines: `start`, `end`, `y`, `color`,
#'     `feature`
#'   * `chevrons`: the strand chevrons on those lines where the display draws
#'     them: `position`, `y`, `direction`, `color`, `feature`
#'   * `arrows`: the strand arrow off a feature's end where it is wide enough
#'     on screen: `position`, `y`, `direction`, `head_half_height`, `color`,
#'     `feature`
#'   * `labels`: each name, description and transcript label the display
#'     shows: `text`, `kind`, `position` (its left edge, in bp), `y` (its top),
#'     `color`, `feature`
#'   * `legend`: `label` and `color`, the key the display shows for a color
#'     field, its title as the frame's `"title"` attribute
#'
#'   `feature` is a row number of `features`. `glyph` gives the pixel sizes
#'   the display draws text, arrows and chevrons at, and `height` the track's
#'   laid-out height in pixels.
#' @export
#' @examples
#' gff <- system.file("extdata", "volvox.sort.gff3.gz", package = "ggjbrowse")
#' g <- jbrowse_features(gff, "ctgA", 0, 25000)
#' head(g$features)
jbrowse_features <- function(file, chrom, start, end, display = NULL,
                             width = 800, reference = NULL,
                             reference_chrom = chrom, color_by_cds = FALSE,
                             index = NULL, byte_limit = 2e8) {
  args <- list(
    adapterConfig = features_adapter(file, index),
    region = list(refName = chrom, originalRefName = reference_chrom,
      start = start, end = end),
    bpPerPx = (end - start) / width,
    colorByCDS = color_by_cds,
    byteLimit = byte_limit)
  if (!is.null(reference)) args$sequenceAdapter <- sequence_adapter(reference)
  if (!is.null(display)) {
    if (!is.null(display$color)) display$color <- color_object(display$color)
    if (!is.null(display$filter)) display$filter <- I(display$filter)
    args$display <- display
  }
  t <- call_engine("features", args)
  on_feature <- function(table) {
    d <- as_frame(table)
    d$feature <- d$feature + 1
    d
  }
  arrows <- on_feature(t$arrows)
  names(arrows)[names(arrows) == "headHalfHeight"] <- "head_half_height"
  legend <- as_frame(t$legend)
  attr(legend, "title") <- t$legendTitle
  structure(
    list(features = as_frame(t$features),
      subfeatures = on_feature(t$subfeatures), rects = on_feature(t$rects),
      lines = on_feature(t$lines), chevrons = on_feature(t$chevrons),
      arrows = arrows, labels = on_feature(t$labels), legend = legend,
      glyph = t$glyph, height = t$height, bp_per_px = t$bpPerPx,
      region = list(chrom = chrom, start = start, end = end)),
    class = "jbrowse_features")
}

# A line's width in ggplot's linewidth units from millimetres
mm_linewidth <- function(mm) mm * 96 / 25.4 / ggplot2::.pt

#' @export
print.jbrowse_features <- function(x, ...) {
  cat(sprintf(
    "<jbrowse_features> %s:%.0f-%.0f, %d features, %d transcripts, %d boxes, %d labels, %.0f px tall\n",
    x$region$chrom, x$region$start, x$region$end, nrow(x$features),
    nrow(x$subfeatures), nrow(x$rects), nrow(x$labels), x$height))
  invisible(x)
}

#' Draw features as JBrowse's feature track
#'
#' Each feature's glyph on its row: the boxes in the colors the display's
#' `color` object paints, the intron lines with their strand chevrons, the
#' strand arrows, and the labels, at the pixel sizes the display draws them.
#' The y axis is the display's pixels, so the track reads as the browser
#' shows it when the panel is about `height` pixels tall at the figure's
#' scale.
#'
#' @param object a [jbrowse_features()] result
#' @param figure_width_mm the width the figure is printed at, which sets how
#'   large the display's pixels are and so the labels' text size
#' @param ... unused
#' @return a ggplot, to extend or to stack over an alignments track with
#'   patchwork
#' @exportS3Method ggplot2::autoplot
#' @examples
#' gff <- system.file("extdata", "volvox.sort.gff3.gz", package = "ggjbrowse")
#' g <- jbrowse_features(gff, "ctgA", 0, 25000)
#' autoplot(g)
autoplot.jbrowse_features <- function(object, figure_width_mm = 170, ...) {
  bp <- object$bp_per_px
  glyph <- object$glyph
  px_mm <- figure_width_mm * bp / (object$region$end - object$region$start)
  xlim <- c(object$region$start, object$region$end)
  legend <- object$legend
  keyed <- legend[!duplicated(legend$color), , drop = FALSE]
  chev <- object$chevrons
  half_w <- glyph$chevronWidth / 2 * bp
  half_h <- glyph$chevronHeight / 2
  chevron_arms <- rbind(
    data.frame(x = chev$position - half_w * chev$direction, y = chev$y - half_h,
      xend = chev$position + half_w * chev$direction, yend = chev$y,
      color = chev$color),
    data.frame(x = chev$position + half_w * chev$direction, y = chev$y,
      xend = chev$position - half_w * chev$direction, yend = chev$y + half_h,
      color = chev$color))
  a <- object$arrows
  stem_end <- a$position + glyph$arrowLength / 2 * bp * a$direction
  tip <- a$position + glyph$arrowLength * bp * a$direction
  heads <- do.call(rbind, lapply(seq_len(nrow(a)), function(i) data.frame(
    id = i, x = c(stem_end[i], stem_end[i], tip[i]),
    y = a$y[i] + c(-1, 1, 0) * a$head_half_height[i], color = a$color[i])))
  heads <- heads %||% data.frame(id = integer(0), x = numeric(0),
    y = numeric(0), color = character(0))
  labels <- object$labels
  ggplot2::ggplot() +
    ggplot2::geom_segment(ggplot2::aes(x = start, xend = end, y = y, yend = y,
      colour = color), object$lines, linewidth = 0.3) +
    ggplot2::geom_segment(ggplot2::aes(x = x, xend = xend, y = y, yend = yend,
      colour = color), chevron_arms, linewidth = 0.3) +
    ggplot2::geom_rect(ggplot2::aes(xmin = start, xmax = end, ymin = top,
      ymax = top + height, fill = color), object$rects) +
    ggplot2::scale_fill_identity(name = attr(legend, "title"),
      breaks = keyed$color, labels = keyed$label,
      guide = if (nrow(keyed)) "legend" else "none") +
    ggplot2::geom_segment(ggplot2::aes(x = position, xend = stem_end,
      y = y, yend = y, colour = color), cbind(a, stem_end = stem_end),
      linewidth = mm_linewidth(2 * glyph$arrowStemHalfHeight * px_mm)) +
    ggplot2::geom_polygon(ggplot2::aes(x = x, y = y, group = id, fill = color),
      heads) +
    ggplot2::geom_text(ggplot2::aes(x = position, y = y, label = text,
      colour = color), labels, hjust = 0, vjust = 1,
      size = glyph$fontSize * px_mm) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, object$height),
      expand = FALSE) +
    ggplot2::labs(x = sprintf("%s (bp)", object$region$chrom), y = NULL) +
    ggplot2::theme_minimal() +
    no_y
}

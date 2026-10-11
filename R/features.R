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
  if (grepl("\\.vcf\\.gz$", name)) return(tabix("VcfTabixAdapter", "vcfGzLocation"))
  if (grepl("\\.vcf$", name)) return(plain("VcfAdapter", "vcfLocation"))
  if (grepl("\\.bed$", name)) return(plain("BedAdapter", "bedLocation"))
  stop("cannot tell the format of ", file,
    ": expected .gff3, .gtf, .bed or .vcf (each optionally .gz with an index) or .bb",
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
#' @param file a GFF3, GTF, BED or VCF file (plain, or bgzipped with a tabix
#'   or CSI index beside it) or a BigBed, a path or a URL; or a JBrowse track's
#'   `adapter` config, as a list, for an adapter's own settings
#' @param region the region: a locus as JBrowse's location box writes it
#'   (`"ctgA:1,001-9,000"`, 1-based and inclusive), several of them, which
#'   are read in turn and drawn end to end on one axis as a JBrowse view lays
#'   them out, a data frame of `chrom`, `start` and `end`, or a contig name
#'   with `start` and `end`
#' @param start,end with a contig name as `region`, the region's bounds,
#'   0-based half-open
#' @param display the display's settings, as a track config's `displays`
#'   entry writes them for a `LinearBasicDisplay`: a list with any of
#'   `color` (a constant such as `"purple"`, or the color object's slots,
#'   such as `list(field = "type")`), `displayMode` (`"normal"`, `"compact"`,
#'   `"superCompact"`, `"collapsed"`), `showLabels` (`"auto"`, `"name"`,
#'   `"nameAndDescription"`, `"description"`, `"none"`), `geneGlyphMode`
#'   (`"auto"`, `"all"`, `"longestCoding"`), `subfeatureLabels` (`"none"`,
#'   `"below"`, `"overlay"`), `featureHeight`, `showOnlyGenes`, `filter` (jexl
#'   expressions), `displayDirectionalChevrons`, and `heightMode` with
#'   `height`. A figure has no track height to fit, so `heightMode` defaults to
#'   `"grow"`, every feature at full size; `list(heightMode = "fit", height =
#'   100)` gives up the labels, isoforms and row height the browser gives up to
#'   fit a 100 px track
#' @param width the width in pixels the display lays the region out across:
#'   whether a name fits beside its feature, which chevrons and arrows draw and
#'   how the `auto` modes resolve all follow from it. [autoplot()] scales those
#'   pixels to the panel, so a smaller `width` gives larger labels in a figure
#'   of the same size
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
#'     `color`, `feature`, and `subfeature`, the transcript a transcript label
#'     names (a row number of `subfeatures`, `NA` for a feature's own label)
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
jbrowse_features <- function(file, region, start = NULL, end = NULL,
                             display = NULL, width = 800, reference = NULL,
                             reference_chrom = NULL, color_by_cds = FALSE,
                             index = NULL, byte_limit = 2e8) {
  read_regions(region, start, end, function(chrom, start, end, share, i)
    features_region(file, chrom, start, end, display, width * share,
      reference, reference_chrom[i] %||% chrom, color_by_cds, index,
      byte_limit))
}

features_region <- function(file, chrom, start, end, display, width,
                            reference, reference_chrom, color_by_cds, index,
                            byte_limit) {
  args <- list(
    adapterConfig = adapter_of(file, function()
      features_adapter(file, index)),
    region = list(refName = chrom, originalRefName = reference_chrom,
      start = start, end = end),
    bpPerPx = (end - start) / width,
    colorByCDS = color_by_cds,
    byteLimit = byte_limit)
  if (!is.null(reference)) args$sequenceAdapter <- sequence_adapter(reference)
  display <- display %||% list()
  display$heightMode <- display$heightMode %||% "grow"
  display$color <- config_color(display$color)
  if (!is.null(display$filter)) display$filter <- I(display$filter)
  args$display <- display
  t <- call_engine("features", args)
  features <- as_frame(t$features)
  on_feature <- function(table) {
    d <- as_frame(table)
    d$feature <- d$feature + 1
    for (column in setdiff(c("name", "type", "strand", "start", "end"),
      names(d))) {
      d[[column]] <- features[[column]][d$feature]
    }
    if (!is.null(d$color)) d$color <- jbrowse_color(d$color)
    d
  }
  labels <- on_feature(t$labels)
  labels$subfeature <- ifelse(labels$subfeature >= 0, labels$subfeature + 1,
    NA_integer_)
  arrows <- on_feature(t$arrows)
  names(arrows)[names(arrows) == "headHalfHeight"] <- "head_half_height"
  legend <- as_frame(t$legend)
  attr(legend, "title") <- t$legendTitle
  structure(
    list(features = features,
      subfeatures = on_feature(t$subfeatures), rects = on_feature(t$rects),
      lines = on_feature(t$lines), chevrons = on_feature(t$chevrons),
      arrows = arrows, labels = labels, legend = legend,
      keys = list(features = legend),
      glyph = t$glyph, height = t$height, bp_per_px = t$bpPerPx,
      region = list(chrom = chrom, start = start, end = end),
      missing_contig = t$missingContig),
    primary = "features",
    class = c("jbrowse_features", "jbrowse_result"))
}

# Text sized in the display's pixels, scaled to the panel when it draws:
# `px_wide` and `px_tall` are the display's extent in pixels
GeomDisplayText <- ggplot2::ggproto("GeomDisplayText", ggplot2::GeomText,
  default_aes = ggplot2::aes(!!!ggplot2::GeomText$default_aes,
    px_wide = NA, px_tall = NA),
  draw_panel = function(self, data, panel_params, coord) {
    text <- ggplot2::ggproto_parent(ggplot2::GeomText, self)$draw_panel(data,
      panel_params, coord)
    grid::gTree(children = grid::gList(text), px_wide = data$px_wide[1],
      px_tall = data$px_tall[1], cl = "ggjbrowse_display_text")
  })

#' @exportS3Method grid::makeContent
makeContent.ggjbrowse_display_text <- function(x) {
  mm_per_px <- min(
    grid::convertWidth(grid::unit(1, "npc"), "mm", valueOnly = TRUE) /
      x$px_wide,
    grid::convertHeight(grid::unit(1, "npc"), "mm", valueOnly = TRUE) /
      x$px_tall)
  x$children[[1]]$gp$fontsize <- x$children[[1]]$gp$fontsize * mm_per_px
  x
}

#' @export
print.jbrowse_features <- function(x, ...) {
  cat(sprintf(
    "<jbrowse_features> %s:%.0f-%.0f, %d features, %d transcripts, %d boxes, %d labels, %.0f px tall\n",
    x$region$chrom, x$region$start, x$region$end, nrow(x$features),
    nrow(x$subfeatures), nrow(x$rects), nrow(x$labels), x$height))
  invisible(x)
}

# A chevron's two arms, each a segment, at the display's px sizes
chevron_arms <- function(x) {
  chev <- x$chevrons
  half_w <- x$glyph$chevronWidth / 2 * x$bp_per_px
  half_h <- x$glyph$chevronHeight / 2
  rbind(
    transform(chev, x = position - half_w * direction, y0 = y - half_h,
      xend = position + half_w * direction, yend = y),
    transform(chev, x = position + half_w * direction, y0 = y,
      xend = position - half_w * direction, yend = y + half_h))
}

# The glyph's boxes, each at least a display pixel wide, as the display draws
# a 1 bp feature at any zoom
boxes <- function(x) {
  d <- x$rects
  d$end <- pmax(d$end, d$start + x$bp_per_px)
  d
}

# A strand arrow's stem, as a box in data units
arrow_stems <- function(x) {
  a <- x$arrows
  stem_end <- a$position + x$glyph$arrowLength / 2 * x$bp_per_px * a$direction
  half <- x$glyph$arrowStemHalfHeight
  transform(a, xmin = pmin(position, stem_end), xmax = pmax(position, stem_end),
    ymin = y - half, ymax = y + half)
}

# A strand arrow's head, as a triangle's three corners
arrow_heads <- function(x) {
  a <- x$arrows
  length <- x$glyph$arrowLength * x$bp_per_px
  rows <- rep(seq_len(nrow(a)), each = 3)
  h <- a[rows, , drop = FALSE]
  corner <- rep(1:3, nrow(a))
  h$id <- rows
  h$x <- h$position + length * h$direction * ifelse(corner == 3, 1, 0.5)
  h$y <- h$y + c(-1, 1, 0)[corner] * h$head_half_height
  h
}

# The entries of `mapping` a geom draws
for_geom <- function(mapping, geom) {
  if (is.null(mapping)) return(NULL)
  mapping[intersect(names(mapping), geom$aesthetics())]
}

#' A feature track's marks as ggplot layers
#'
#' `geom_jbrowse_glyphs()` draws each feature's glyph as the display does:
#' the intron and connector lines with their strand chevrons, the boxes, and
#' the strand arrows. `geom_jbrowse_labels()` writes the names and
#' descriptions the display shows at this zoom, at its text size scaled to the
#' panel. Each mark table carries its feature's `name`, `type` and `strand`
#' beside the display's `color`, and a mark with no span of its own (a
#' chevron, an arrow, a label) its feature's `start` and `end`, so
#' `aes(fill = type)` colors the boxes by feature type through ggplot's own
#' scales. A `fill` mapping reaches the boxes and arrowheads, and a `colour`
#' mapping every part that has a stroke, the boxes' outlines included.
#' @param mapping aesthetics that override the defaults
#' @param data a [jbrowse_features()] result, or `NULL` for the plot's
#' @param ... passed to [ggplot2::layer()] params
#' @return a list of ggplot layers
#' @export
#' @examples
#' gff <- system.file("extdata", "volvox.sort.gff3.gz", package = "ggjbrowse")
#' g <- jbrowse_features(gff, "ctgA", 0, 25000)
#' library(ggplot2)
#' ggplot(g) + geom_jbrowse_glyphs(aes(fill = type)) + scale_y_reverse()
geom_jbrowse_glyphs <- function(mapping = NULL, data = NULL, ...) {
  list(
    jbrowse_layer(ggplot2::GeomSegment, data, function(x) x$lines,
      for_geom(mapping, ggplot2::GeomSegment),
      ggplot2::aes(x = start, xend = end, y = y, yend = y, colour = color),
      list(linewidth = 0.3, ...), keyed = FALSE),
    jbrowse_layer(ggplot2::GeomSegment, data, chevron_arms,
      for_geom(mapping, ggplot2::GeomSegment),
      ggplot2::aes(x = x, xend = xend, y = y0, yend = yend, colour = color),
      list(linewidth = 0.3, ...), keyed = FALSE),
    jbrowse_layer(ggplot2::GeomRect, data, boxes,
      for_geom(mapping, ggplot2::GeomRect),
      ggplot2::aes(xmin = start, xmax = end, ymin = top, ymax = top + height,
        fill = color), list(...)),
    jbrowse_layer(ggplot2::GeomRect, data, arrow_stems,
      for_geom(mapping, ggplot2::GeomRect),
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = color), list(...), keyed = FALSE),
    jbrowse_layer(ggplot2::GeomPolygon, data, arrow_heads,
      for_geom(mapping, ggplot2::GeomPolygon),
      ggplot2::aes(x = x, y = y, group = id, fill = color), list(...),
      keyed = FALSE))
}

#' @rdname geom_jbrowse_glyphs
#' @export
geom_jbrowse_labels <- function(mapping = NULL, data = NULL, ...) {
  labels <- function(x) {
    d <- x$labels
    n <- nrow(d)
    d$size <- rep(x$glyph$fontSize, n)
    d$px_wide <- rep((x$region$end - x$region$start) / x$bp_per_px, n)
    d$px_tall <- rep(x$height, n)
    d
  }
  jbrowse_layer(GeomDisplayText, data, labels, mapping,
    ggplot2::aes(x = position, y = y, label = text, colour = color,
      size = I(size), px_wide = I(px_wide), px_tall = I(px_tall)),
    list(hjust = 0, vjust = 1, ...), keyed = FALSE)
}

#' Draw features as JBrowse's feature track
#'
#' Each feature's glyph on its row ([geom_jbrowse_glyphs()]) and the labels
#' the display shows ([geom_jbrowse_labels()]), keyed by the display's
#' `color` object. The y axis is the display's pixels. Labels scale with the
#' panel: at draw time they take the display's pixel size across the panel's
#' width or down its height, whichever is smaller, so they keep inside the
#' room the layout left them at any output size.
#'
#' @param object a [jbrowse_features()] result
#' @param ... unused
#' @return a ggplot, to extend or to stack over an alignments track with
#'   patchwork
#' @exportS3Method ggplot2::autoplot
#' @examples
#' gff <- system.file("extdata", "volvox.sort.gff3.gz", package = "ggjbrowse")
#' g <- jbrowse_features(gff, "ctgA", 0, 25000)
#' autoplot(g)
autoplot.jbrowse_features <- function(object, ...) {
  object <- on_axis(object)
  ggplot2::ggplot(object) +
    geom_jbrowse_glyphs() +
    scale_fill_jbrowse(object) +
    geom_jbrowse_labels() +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_cartesian(
      xlim = c(object$region$start, object$region$end),
      ylim = c(0, object$height), expand = FALSE) +
    genomic_x(object) +
    ggplot2::labs(x = axis_label(object), y = NULL) +
    ggplot2::theme_minimal() +
    no_y
}

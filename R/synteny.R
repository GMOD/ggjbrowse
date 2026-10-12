# A synteny file's adapter, read off its extension: PAF (plain or gzipped),
# indexed PIF, MUMmer delta, UCSC chain or MashMap output. The adapter names
# its sides `query` and `target`, as the file's columns do.
synteny_adapter <- function(file, index) {
  name <- tolower(sub("\\?.*$", "", file))
  sides <- list(assemblyNames = c("query", "target"))
  plain <- function(type, slot) {
    a <- c(list(type = type), sides)
    a[[slot]] <- location(file)
    a
  }
  if (grepl("\\.paf(\\.gz)?$", name)) return(plain("PAFAdapter", "pafLocation"))
  if (grepl("\\.delta$", name)) return(plain("DeltaAdapter", "deltaLocation"))
  if (grepl("\\.chain$", name)) return(plain("ChainAdapter", "chainLocation"))
  if (grepl("\\.out$", name)) return(plain("MashMapAdapter", "outLocation"))
  if (grepl("\\.pif\\.gz$", name)) {
    index <- index %||% paste0(file, ".tbi")
    a <- plain("PairwiseIndexedPAFAdapter", "pifGzLocation")
    a$index <- list(location = location(index),
      indexType = if (grepl("\\.csi$", index)) "CSI" else "TBI")
    return(a)
  }
  stop("cannot tell the format of ", file,
    ": expected .paf (optionally .gz), .pif.gz with its index, .delta, .chain or .out",
    call. = FALSE)
}

# One axis of a comparative view for the engine: its regions, the assembly
# they are on where the caller names it, and the pixels it spans
axis_args <- function(regions, assembly, width) {
  a <- list(regions = lapply(seq_len(nrow(regions)), function(i)
    list(refName = regions$chrom[i], start = regions$start[i],
      end = regions$end[i])), width = width)
  if (!is.null(assembly) && !is.na(assembly)) a$assemblyName <- assembly
  a
}

# The engine's camelCase columns as this package's snake_case
snake <- function(d) {
  names(d) <- gsub("([A-Z])", "_\\L\\1", names(d), perl = TRUE)
  d
}

comparative_args <- function(file, index, view) {
  args <- list(adapterConfig = adapter_of(file, function()
    synteny_adapter(file, index)))
  if (!is.null(view)) {
    view$color <- config_color(view$color)
    args$view <- view
  }
  args
}

#' Read a pairwise alignment as JBrowse's linear synteny view draws it
#'
#' Runs JBrowse's own synteny display in V8 over two rows, the query assembly
#' above the target: the alignment blocks the two rows share, the ribbons the
#' band draws between them, with each CIGAR insertion, deletion and skip wide
#' enough to show as its own colored piece, and the color and opacity the
#' view's settings give each. Coordinates are 0-based, half-open, on each
#' assembly.
#'
#' @param file a PAF (plain or gzipped), an indexed PIF (`.pif.gz` with its
#'   tabix index beside it), a MUMmer `.delta`, a UCSC `.chain` or a MashMap
#'   `.out` file, a path or a URL; or a JBrowse synteny track's `adapter`
#'   config, as a list, such as an `MCScanAnchorsAdapter`'s
#' @param query,target the regions each row shows: a locus as JBrowse's
#'   location box writes it (`"ctgA:1-50,001"`, 1-based and inclusive),
#'   several of them, drawn end to end as a view lays them out, or a data
#'   frame of `chrom`, `start` and `end`
#' @param view the view's settings, as a `LinearSyntenyView` snapshot writes
#'   them: a list with any of `color` (a CSS color, or the color object's
#'   slots, such as `list(field = "strand")`, `list(field = "identity")`,
#'   `list(field = "query")` or `list(field = "mapq")`), `opacity` (a number,
#'   or an opacity object with a `field`), `minAlignmentLength`, `cigarMode`
#'   (`"full"` for colored indels, `"matches"` for the matches alone with the
#'   indels left clear, `"off"` for the block alone), `drawLocationMarkers`,
#'   `fadeThinAlignmentsMode` (`"auto"`, `"on"`, `"off"`) and `height`, the
#'   band's pixels
#' @param width the width in pixels both rows draw their regions across,
#'   which decides which CIGAR ops are wide enough to draw and which ribbons
#'   are thin enough to fade
#' @param assemblies for an adapter config naming its assemblies, which of
#'   them each row is on, `c(query, target)`; unset, the adapter's own order.
#'   A pairwise file answers from either side, so `c("volvox", "volvox_del")`
#'   on a PAF of `volvox_del` against `volvox` puts `volvox` on top
#' @param index the tabix or CSI index of a `.pif.gz`, if not beside it
#' @return a `jbrowse_synteny` list:
#'   * `alignments`: each block, in the band's paint order: `id`, `name`,
#'     `query_chrom`, `query_start`, `query_end`, `target_chrom`,
#'     `target_start`, `target_end`, `strand`, `length`, the measurements the
#'     file carries (`identity`, `mean_identity`, `mapq`, and any column the
#'     adapter's `attributeColumns` declares), `color` and `alpha`
#'   * `ribbons`: what the band draws, in paint order, each a trapezoid from
#'     a span on the query row to a span on the target row: `alignment` (a row
#'     number of `alignments`), `kind` (`"match"`, `"insertion"`,
#'     `"deletion"`, `"skip"`, or `"marker"` for a location tick),
#'     `query_region`, `query_chrom`, `query_start`, `query_end`,
#'     `target_region`, `target_chrom`, `target_start`, `target_end`, where
#'     the ribbon joins `query_start` to `target_start` and `query_end` to
#'     `target_end`, so a reverse-strand block crosses; `color`, `alpha`, and
#'     `thin`, whether the band draws it as a one-pixel line
#'   * `offscreen`: the blocks on one `axis` whose mate lies on a contig the
#'     other row is not showing, with their span and their mate's
#'   * `legend`: `label` and `color`, the key the view shows, its title as
#'     the frame's `"title"` attribute
#'
#'   `query` and `target` are the rows' regions as laid out, `height` the
#'   band's pixels, `alpha` the opacity every ribbon draws at, `field` the
#'   color field painted and `ramp` a continuous field's color ramp.
#' @export
#' @examples
#' paf <- system.file("extdata", "volvox_inv_indels.paf", package = "ggjbrowse")
#' s <- jbrowse_synteny(paf, "ctgA:1-49,186", "ctgA:1-50,001")
#' head(s$alignments)
jbrowse_synteny <- function(file, query, target, view = NULL, width = 800,
                            assemblies = NULL, index = NULL) {
  q <- regions_of(query, NULL, NULL)
  t <- regions_of(target, NULL, NULL)
  args <- c(comparative_args(file, index, view),
    list(query = axis_args(q, assemblies[1], width),
      target = axis_args(t, assemblies[2], width)))
  r <- call_engine("synteny", args)
  alignments <- snake(as_frame(r$alignments))
  alignments$color <- jbrowse_color(alignments$color)
  ribbons <- snake(as_frame(r$ribbons))
  ribbons$alignment <- ribbons$alignment + 1
  ribbons$query_region <- ribbons$query_region + 1
  ribbons$target_region <- ribbons$target_region + 1
  ribbons$color <- jbrowse_color(ribbons$color)
  offscreen <- snake(as_frame(r$offscreen))
  offscreen$region <- offscreen$region + 1
  key <- key_frame(unlist(r$legend$label), unlist(r$legend$color),
    r$legendTitle)
  attr(key, "override") <- list(alpha = max(r$alpha, 0.45))
  structure(
    list(alignments = alignments, ribbons = ribbons, offscreen = offscreen,
      legend = key, keys = list(ribbons = key), ramp = r$ramp,
      field = r$field, alpha = r$alpha, height = r$height,
      fade_thin = r$fadeThin, has_cigar = r$hasCigar,
      assemblies = unlist(r$assemblyNames), width = width,
      query = region_layout(q), target = region_layout(t)),
    primary = "ribbons",
    class = c("jbrowse_synteny", "jbrowse_result"))
}

#' @export
print.jbrowse_synteny <- function(x, ...) {
  cat(sprintf(
    "<jbrowse_synteny> %s over %s, %d alignments, %d ribbons, %d off screen\n",
    paste(loci_of_layout(x$query), collapse = ", "),
    paste(loci_of_layout(x$target), collapse = ", "),
    nrow(x$alignments), nrow(x$ribbons), nrow(x$offscreen)))
  invisible(x)
}

# A position on a row as a pixel of the band: its place on the row's regions
# laid end to end, as a share of the row's span, across the view's width
band_px <- function(layout, region, coord, width) {
  shift <- layout$offset[region] - layout$start[region]
  lo <- min(layout$cum_start)
  (coord + shift - lo) / (max(layout$cum_end) - lo) * width
}

# The ribbons and off-screen marks with each foot at its pixel on its row
synteny_on_axis <- function(x) {
  r <- x$ribbons
  if (!is.null(r$.on_axis)) return(x)
  r$query_px_start <- band_px(x$query, r$query_region, r$query_start, x$width)
  r$query_px_end <- band_px(x$query, r$query_region, r$query_end, x$width)
  r$target_px_start <- band_px(x$target, r$target_region, r$target_start,
    x$width)
  r$target_px_end <- band_px(x$target, r$target_region, r$target_end, x$width)
  r$.on_axis <- rep(TRUE, nrow(r))
  x$ribbons <- r
  o <- x$offscreen
  layout <- list(query = x$query, target = x$target)
  o$px_start <- vapply(seq_len(nrow(o)), function(i)
    band_px(layout[[o$axis[i]]], o$region[i], o$start[i], x$width), 0)
  o$px_end <- vapply(seq_len(nrow(o)), function(i)
    band_px(layout[[o$axis[i]]], o$region[i], o$end[i], x$width), 0)
  o$y <- ifelse(o$axis == "query", 1, 0)
  x$offscreen <- o
  x
}

# The points of each ribbon edge from (x0, y0) to (x1, y1): a straight line,
# or the band's cubic Bezier with both control points at mid height
edge_points <- function(x0, x1, y0, y1, curves, n = 16) {
  if (!curves) return(list(x = c(x0, x1), y = c(y0, y1)))
  t <- seq(0, 1, length.out = n)
  s <- t * t * (3 - 2 * t)
  list(x = x0 + (x1 - x0) * s, y = y0 + (y1 - y0) * t)
}

# Each filled ribbon's outline, top to bottom on its first edge and back up
# its second, one `ribbon` group each
ribbon_polygons <- function(r, curves) {
  if (!nrow(r)) {
    return(data.frame(x = numeric(0), y = numeric(0), ribbon = integer(0),
      color = jbrowse_color(character(0)), alpha = numeric(0)))
  }
  pts <- lapply(seq_len(nrow(r)), function(i) {
    down <- edge_points(r$query_px_start[i], r$target_px_start[i], 1, 0,
      curves)
    up <- edge_points(r$target_px_end[i], r$query_px_end[i], 0, 1, curves)
    data.frame(x = c(down$x, up$x), y = c(down$y, up$y), ribbon = i)
  })
  d <- do.call(rbind, pts)
  d$color <- r$color[d$ribbon]
  d$alpha <- r$alpha[d$ribbon]
  d$kind <- r$kind[d$ribbon]
  d$alignment <- r$alignment[d$ribbon]
  d
}

# Each thin ribbon's centerline, from the middle of its top edge to the
# middle of its bottom edge
ribbon_lines <- function(r, curves) {
  if (!nrow(r)) {
    return(data.frame(x = numeric(0), y = numeric(0), ribbon = integer(0),
      color = jbrowse_color(character(0)), alpha = numeric(0)))
  }
  pts <- lapply(seq_len(nrow(r)), function(i) {
    e <- edge_points((r$query_px_start[i] + r$query_px_end[i]) / 2,
      (r$target_px_start[i] + r$target_px_end[i]) / 2, 1, 0, curves)
    data.frame(x = e$x, y = e$y, ribbon = i)
  })
  d <- do.call(rbind, pts)
  d$color <- r$color[d$ribbon]
  d$alpha <- r$alpha[d$ribbon]
  d$kind <- r$kind[d$ribbon]
  d$alignment <- r$alignment[d$ribbon]
  d
}

#' A synteny band's ribbons as ggplot layers
#'
#' `geom_jbrowse_ribbons()` draws each ribbon between the query row at the
#' top of the panel (`y = 1`) and the target row at the bottom (`y = 0`) as
#' the band draws it: a filled trapezoid in its color at its opacity, or a
#' one-pixel line where it is thinner than a pixel, with the band's location
#' ticks where the view draws them. `geom_jbrowse_offscreen_mates()` marks on
#' a row each block whose mate the other row is not showing.
#' @param mapping aesthetics that override the defaults
#' @param data a [jbrowse_synteny()] result, or `NULL` for the plot's
#' @param curves draw each ribbon edge as the view's `drawCurves` does, a
#'   cubic curve vertical at both rows
#' @param ... passed to [ggplot2::layer()] params
#' @return a list of ggplot layers
#' @export
#' @examples
#' paf <- system.file("extdata", "volvox_inv_indels.paf", package = "ggjbrowse")
#' s <- jbrowse_synteny(paf, "ctgA:1-49,186", "ctgA:1-50,001")
#' library(ggplot2)
#' ggplot(s) + geom_jbrowse_ribbons(curves = TRUE) + scale_fill_jbrowse(s)
geom_jbrowse_ribbons <- function(mapping = NULL, data = NULL, curves = FALSE,
                                 ...) {
  filled <- function(x) ribbon_polygons(x$ribbons[!x$ribbons$thin, ,
    drop = FALSE], curves)
  lines <- function(x) ribbon_lines(x$ribbons[x$ribbons$thin, ,
    drop = FALSE], curves)
  list(
    jbrowse_layer(ggplot2::GeomPolygon, data, filled,
      for_geom(mapping, ggplot2::GeomPolygon),
      ggplot2::aes(x = x, y = y, group = ribbon, fill = color,
        alpha = I(alpha)), list(colour = NA, ...)),
    jbrowse_layer(ggplot2::GeomPath, data, lines,
      for_geom(mapping, ggplot2::GeomPath),
      ggplot2::aes(x = x, y = y, group = ribbon, colour = color,
        alpha = I(alpha)), list(linewidth = 0.3, ...), keyed = FALSE))
}

#' @rdname geom_jbrowse_ribbons
#' @export
geom_jbrowse_offscreen_mates <- function(mapping = NULL, data = NULL, ...) {
  jbrowse_layer(ggplot2::GeomSegment, data, function(x) x$offscreen,
    for_geom(mapping, ggplot2::GeomSegment),
    ggplot2::aes(x = px_start, xend = px_end, y = y, yend = y),
    list(colour = "grey30", linewidth = 2, lineend = "butt", ...),
    keyed = FALSE)
}

# The round positions a row's ruler labels, at their pixels on the band
layout_breaks <- function(layout, px) {
  r <- layout
  w <- r$end - r$start
  at <- lapply(seq_len(nrow(r)), function(i) {
    b <- pretty(c(r$start[i], r$end[i]), n = max(2, round(6 * w[i] / sum(w))))
    lo <- if (i == 1) r$start[i] else r$start[i] + 1
    hi <- if (i == nrow(r)) r$end[i] else r$end[i] - 1
    b[b >= lo & b <= hi]
  })
  list(
    breaks = unlist(lapply(seq_len(nrow(r)), function(i)
      px(r, rep(i, length(at[[i]])), at[[i]]))),
    labels = format(unlist(at), big.mark = ",", trim = TRUE,
      scientific = FALSE))
}

loci_of_layout <- function(layout) {
  sprintf("%s:%s-%s", layout$chrom, format(layout$start + 1, big.mark = ",",
    trim = TRUE), format(layout$end, big.mark = ",", trim = TRUE))
}

layout_label <- function(layout, assembly = "") {
  what <- if (nrow(layout) == 1) layout$chrom
    else paste(loci_of_layout(layout), collapse = ", ")
  name <- if (nzchar(assembly)) paste0(assembly, " ") else ""
  sprintf("%s%s (bp%s)", name, what,
    if (nrow(layout) == 1) "" else ", end to end")
}

# The band's two rulers: the target row's positions under the panel and the
# query row's above it, each at its pixel
synteny_x <- function(x) {
  on_row <- function(layout) layout_breaks(layout, function(r, i, at)
    band_px(r, i, at, x$width))
  q <- on_row(x$query)
  t <- on_row(x$target)
  ggplot2::scale_x_continuous(breaks = t$breaks, labels = t$labels,
    name = layout_label(x$target, x$assemblies[2]),
    sec.axis = ggplot2::sec_axis(~., breaks = q$breaks, labels = q$labels,
      name = layout_label(x$query, x$assemblies[1])))
}

synteny_band <- function(x, curves) {
  x <- on_axis(x)
  ggplot2::ggplot(x) +
    geom_jbrowse_ribbons(curves = curves) +
    geom_jbrowse_offscreen_mates() +
    scale_fill_jbrowse(x) +
    ggplot2::coord_cartesian(xlim = c(0, x$width), ylim = c(0, 1),
      expand = FALSE) +
    synteny_x(x) +
    ggplot2::labs(y = NULL) +
    ggplot2::theme_minimal() +
    no_y +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

# The tracks a row shows, each read over the row's regions
row_tracks <- function(tracks, layout, row) {
  if (is.null(tracks)) return(list())
  if (inherits(tracks, "jbrowse_result")) tracks <- list(tracks)
  for (t in tracks) {
    r <- regions_read(t)
    if (!identical(paste(r$start, r$end), paste(layout$start, layout$end))) {
      stop("a ", row, " track draws over the ", row, " row's regions, ",
        paste(loci_of_layout(layout), collapse = ", "), ", but one reads ",
        paste(loci_of(t), collapse = ", "), call. = FALSE)
    }
  }
  tracks
}

#' Draw a pairwise alignment as JBrowse's linear synteny view
#'
#' The band of ribbons between the two rows, the query row's ruler above it
#' and the target row's below, and any tracks over each row: a
#' [jbrowse_features()], [jbrowse_wiggle()], [jbrowse_variants()] or
#' [jbrowse_alignments()] result over that row's regions, drawn as
#' [jbrowse_stack()] draws it, the query's tracks above the band and the
#' target's below.
#'
#' @param object a [jbrowse_synteny()] result
#' @param query,target tracks over each row's regions, a result or a list of
#'   them
#' @param curves draw each ribbon edge as a curve, as the view's `drawCurves`
#'   does
#' @param ... unused
#' @return a ggplot of the band alone, or with tracks a patchwork, whose
#'   height [jbrowse_height()] measures
#' @exportS3Method ggplot2::autoplot
#' @examples
#' e <- function(f) system.file("extdata", f, package = "ggjbrowse")
#' s <- jbrowse_synteny(e("volvox_inv_indels.paf"), "ctgA:1-49,186",
#'   "ctgA:1-50,001", view = list(color = list(field = "strand")))
#' autoplot(s)
#' genes <- jbrowse_features(e("volvox.sort.gff3.gz"), "ctgA:1-50,001",
#'   display = list(showOnlyGenes = TRUE))
#' autoplot(s, target = genes)
autoplot.jbrowse_synteny <- function(object, query = NULL, target = NULL,
                                     curves = FALSE, ...) {
  band <- synteny_band(object, curves)
  above <- row_tracks(query, object$query, "query")
  below <- row_tracks(target, object$target, "target")
  if (!length(above) && !length(below)) return(band)
  bands <- c(lapply(above, track_bands),
    list(list(plots = list(band), heights = object$height)),
    lapply(below, track_bands))
  plots <- do.call(c, lapply(bands, `[[`, "plots"))
  heights <- unlist(lapply(bands, `[[`, "heights"))
  own <- sum(lengths(lapply(bands[seq_along(above)], `[[`, "plots"))) + 1
  plots[-own] <- lapply(plots[-own], function(p) p + no_x)
  plots <- lapply(plots, function(p)
    p + ggplot2::theme(plot.margin = ggplot2::margin(1, 14, 1, 5.5)))
  patchwork::wrap_plots(plots, ncol = 1, heights = heights,
    guides = "collect")
}

#' Read a pairwise alignment as JBrowse's dotplot view draws it
#'
#' Runs JBrowse's own dotplot display in V8: the alignment blocks between the
#' two axes, and the lines the plot draws for them, one per block or, where a
#' block is wide enough on the plot for its CIGAR to show, one per CIGAR
#' step, each in the color the view's settings paint it. Coordinates are
#' 0-based, half-open, on each assembly.
#'
#' @inheritParams jbrowse_synteny
#' @param x,y the regions each axis shows, written as [jbrowse_synteny()]'s
#'   `query` and `target` are: the horizontal axis is the adapter's query
#'   side and the vertical its target, unless `assemblies` says otherwise
#' @param view the view's settings, as a `DotplotView` snapshot writes them:
#'   a list with any of `color`, `opacity`, `minAlignmentLength`,
#'   `minIdentity`, `drawCigar` and `lineWidth`
#' @param width,height the plot's pixels, which decide which blocks are wide
#'   enough for their CIGAR to show
#' @param assemblies which assembly each axis is on, `c(x, y)`
#' @return a `jbrowse_dotplot` list:
#'   * `alignments`: each block: `name`, `x_chrom`, `x_start`, `x_end`,
#'     `y_chrom`, `y_start`, `y_end`, `strand`, `length`, the measurements
#'     the file carries, `color` and `alpha`
#'   * `segments`: the lines the plot draws, in paint order: `alignment` (a
#'     row number of `alignments`), `op` (`"match"`, `"insertion"`,
#'     `"deletion"` or `"skip"`), `x_region`, `x_chrom`, `x_start`, `x_end`,
#'     `y_region`, `y_chrom`, `y_start`, `y_end`, from (`x_start`, `y_start`)
#'     to (`x_end`, `y_end`), `color` and `alpha`
#'   * `legend`: `label` and `color`, the key the view shows
#'
#'   `x` and `y` are the axes' regions as laid out, `line_width` the pixels
#'   each line is drawn at, `alpha` the opacity, and `skipped` counts the
#'   blocks neither axis places.
#' @export
#' @examples
#' paf <- system.file("extdata", "volvox_inv_indels.paf", package = "ggjbrowse")
#' d <- jbrowse_dotplot(paf, "ctgA:1-49,186", "ctgA:1-50,001")
#' head(d$segments)
jbrowse_dotplot <- function(file, x, y, view = NULL, width = 600,
                            height = 600, assemblies = NULL, index = NULL) {
  xr <- regions_of(x, NULL, NULL)
  yr <- regions_of(y, NULL, NULL)
  args <- c(comparative_args(file, index, view),
    list(x = axis_args(xr, assemblies[1], width),
      y = axis_args(yr, assemblies[2], height)))
  r <- call_engine("dotplot", args)
  alignments <- snake(as_frame(r$alignments))
  alignments$color <- jbrowse_color(alignments$color)
  segments <- snake(as_frame(r$segments))
  segments$alignment <- segments$alignment + 1
  segments$x_region <- segments$x_region + 1
  segments$y_region <- segments$y_region + 1
  segments$color <- jbrowse_color(segments$color)
  key <- key_frame(unlist(r$legend$label), unlist(r$legend$color),
    r$legendTitle)
  attr(key, "override") <- list(alpha = max(r$alpha, 0.45))
  structure(
    list(alignments = alignments, segments = segments, legend = key,
      keys = list(segments = key), ramp = r$ramp, field = r$field,
      alpha = r$alpha, line_width = r$lineWidth,
      assemblies = unlist(r$assemblyNames), skipped = r$skipped,
      width = width, height = height,
      x = region_layout(xr), y = region_layout(yr)),
    primary = "segments",
    class = c("jbrowse_dotplot", "jbrowse_result"))
}

#' @export
print.jbrowse_dotplot <- function(x, ...) {
  cat(sprintf(
    "<jbrowse_dotplot> %s against %s, %d alignments, %d segments\n",
    paste(loci_of_layout(x$x), collapse = ", "),
    paste(loci_of_layout(x$y), collapse = ", "),
    nrow(x$alignments), nrow(x$segments)))
  invisible(x)
}

# Each segment's ends on the axes that lay each side's regions end to end
dotplot_on_axis <- function(x) {
  s <- x$segments
  if (!is.null(s$.on_axis)) return(x)
  shift <- function(layout, region) layout$offset[region] - layout$start[region]
  s$x_start <- s$x_start + shift(x$x, s$x_region)
  s$x_end <- s$x_end + shift(x$x, s$x_region)
  s$y_start <- s$y_start + shift(x$y, s$y_region)
  s$y_end <- s$y_end + shift(x$y, s$y_region)
  s$.on_axis <- rep(TRUE, nrow(s))
  x$segments <- s
  x
}

# Pixels at the figure's 150 an inch as ggplot's linewidth, in mm
px_to_mm <- function(px) px / PX_PER_INCH * 25.4

#' A dotplot's lines as a ggplot layer
#'
#' Draws each segment from its start to its end in its color at its
#' opacity, with round ends, as the view does.
#' @param mapping aesthetics that override the defaults
#' @param data a [jbrowse_dotplot()] result, or `NULL` for the plot's
#' @param ... passed to [ggplot2::layer()] params, such as `linewidth`
#' @return a ggplot layer
#' @export
#' @examples
#' paf <- system.file("extdata", "volvox_inv_indels.paf", package = "ggjbrowse")
#' d <- jbrowse_dotplot(paf, "ctgA:1-49,186", "ctgA:1-50,001")
#' library(ggplot2)
#' ggplot(d) + geom_jbrowse_dotplot(aes(colour = op))
geom_jbrowse_dotplot <- function(mapping = NULL, data = NULL, ...) {
  jbrowse_layer(ggplot2::GeomSegment, data, function(x) x$segments,
    for_geom(mapping, ggplot2::GeomSegment),
    ggplot2::aes(x = x_start, xend = x_end, y = y_start, yend = y_end,
      colour = color, alpha = I(alpha)),
    with_defaults(list(...), list(linewidth = px_to_mm(2.5),
      lineend = "round")))
}

# An axis of regions laid end to end: breaks at round positions in each,
# labelled in the region's own coordinates
layout_scale <- function(scale, layout) {
  b <- layout_breaks(layout, function(r, i, at) at - r$start[i] + r$offset[i])
  scale(breaks = b$breaks, labels = b$labels)
}

# The seams between an axis's regions
layout_seams <- function(layout) {
  if (nrow(layout) < 2) return(numeric(0))
  (layout$cum_end[-nrow(layout)] + layout$cum_start[-1]) / 2
}

#' Draw a pairwise alignment as JBrowse's dotplot
#'
#' [geom_jbrowse_dotplot()] at the view's line width, each axis labelled in
#' its regions' positions, with a line at each seam where an axis lays
#' several regions end to end.
#' @param object a [jbrowse_dotplot()] result
#' @param ... unused
#' @return a ggplot
#' @exportS3Method ggplot2::autoplot
#' @examples
#' paf <- system.file("extdata", "volvox_inv_indels.paf", package = "ggjbrowse")
#' autoplot(jbrowse_dotplot(paf, "ctgA:1-49,186", "ctgA:1-50,001",
#'   view = list(color = list(field = "strand"))))
autoplot.jbrowse_dotplot <- function(object, ...) {
  x <- on_axis(object)
  seams <- data.frame(x = layout_seams(x$x))
  rows <- data.frame(y = layout_seams(x$y))
  ggplot2::ggplot(x) +
    ggplot2::geom_vline(data = seams, ggplot2::aes(xintercept = x),
      colour = "grey70") +
    ggplot2::geom_hline(data = rows, ggplot2::aes(yintercept = y),
      colour = "grey70") +
    geom_jbrowse_dotplot(linewidth = px_to_mm(x$line_width)) +
    scale_colour_jbrowse(x) +
    layout_scale(ggplot2::scale_x_continuous, x$x) +
    layout_scale(ggplot2::scale_y_continuous, x$y) +
    ggplot2::coord_cartesian(
      xlim = c(min(x$x$cum_start), max(x$x$cum_end)),
      ylim = c(min(x$y$cum_start), max(x$y$cum_end)), expand = FALSE) +
    ggplot2::labs(x = layout_label(x$x, x$assemblies[1]),
      y = layout_label(x$y, x$assemblies[2])) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey70"))
}

# A synteny file's adapter, read off its extension: PAF (plain or gzipped),
# indexed PIF, MUMmer delta, UCSC chain or MashMap output. The adapter names
# its sides `query` and `target`, as the file's columns do.
synteny_adapter <- function(file, index) {
  a <- guess_adapter(file, index, list(
    list(pattern = "\\.paf(\\.gz)?$", type = "PAFAdapter", slot = "pafLocation"),
    list(pattern = "\\.delta$", type = "DeltaAdapter", slot = "deltaLocation"),
    list(pattern = "\\.chain$", type = "ChainAdapter", slot = "chainLocation"),
    list(pattern = "\\.out$", type = "MashMapAdapter", slot = "outLocation"),
    list(pattern = "\\.pif\\.gz$", type = "PairwiseIndexedPAFAdapter",
      slot = "pifGzLocation", index = "TBI")),
    ".paf (optionally .gz), .pif.gz with its index, .delta, .chain or .out")
  c(a, list(assemblyNames = c("query", "target")))
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

colored_frame <- function(columns) {
  d <- snake(as_frame(columns))
  for (col in intersect(c("color", "halo"), names(d))) {
    d[[col]] <- jbrowse_color(d[[col]])
  }
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

# The names an axis title shows: each assembly name `labels` renames, as the
# named vector `c(query = "peach", target = "grape")` gives it, or an unnamed
# pair in the adapter's query, target order
display_names <- function(names, labels, config) {
  if (is.null(labels)) return(names)
  if (is.null(names(labels))) {
    if (length(labels) != 2) {
      stop("`labels` is a named vector, or a pair in query, target order",
        call. = FALSE)
    }
    names(labels) <- unlist(config$assemblyNames) %||%
      unlist(c(config$queryAssembly, config$targetAssembly))
  }
  hit <- match(names, names(labels))
  ifelse(is.na(hit), names, unname(labels)[hit])
}

# The assembly each axis is on, first and second: `assemblies` where given,
# else the adapter's own names, reversed by `swap`
axis_assemblies <- function(config, assemblies, swap) {
  names <- assemblies %||% unlist(config$assemblyNames) %||%
    unlist(c(config$queryAssembly, config$targetAssembly))
  if (!swap) return(names)
  if (length(names) < 2) {
    stop("swap = TRUE needs two assembly names, the adapter's or `assemblies`",
      call. = FALSE)
  }
  rev(names[1:2])
}

#' The contigs on each side of a pairwise alignment
#'
#' Lists each side's contigs in the order the file first names them, with
#' the length the file states: the regions a whole-genome [jbrowse_dotplot()]
#' or [jbrowse_synteny()] lays end to end where it is given none. A PAF,
#' delta, chain and MashMap file state every contig's length, and an indexed
#' PIF states it on each contig's first row.
#' @inheritParams jbrowse_synteny
#' @return a data frame of `side` (`"query"` or `"target"`), `assembly`,
#'   `chrom` and `length`, `NA` where the file does not say, as in an MCScan
#'   anchors file
#' @export
#' @examples
#' paf <- system.file("extdata", "peach_grape.paf.gz", package = "ggjbrowse")
#' head(jbrowse_contigs(paf))
jbrowse_contigs <- function(file, index = NULL) {
  contigs_of(comparative_args(file, index, NULL)$adapterConfig)
}

contigs_of <- function(config) {
  r <- call_engine("contigs", list(adapterConfig = config))
  side <- function(name) {
    chrom <- as.character(unlist(r[[name]]$refName))
    data.frame(side = rep(name, length(chrom)),
      assembly = rep(r[[name]]$assemblyName, length(chrom)), chrom = chrom,
      length = as.numeric(unlist(r[[name]]$length) %||% numeric(0)))
  }
  rbind(side("query"), side("target"))
}

memo <- function(f) {
  value <- NULL
  function() {
    if (is.null(value)) value <<- f()
    value
  }
}

# The regions an axis shows: each locus given, each contig named alone whole,
# and with none given every contig of the axis's assembly, in natural order
# `side` picks the file's side where the axis names no assembly the file has.
axis_regions <- function(given, contigs, assembly, side) {
  if (is.data.frame(given)) return(regions_of(given, NULL, NULL))
  whole <- !grepl(LOCUS, given)
  if (length(given) && !any(whole)) return(regions_of(given, NULL, NULL))
  all <- contigs()
  named <- length(assembly) && !is.na(assembly) && assembly %in% all$assembly
  own <- all[if (named) all$assembly == assembly else all$side == side, ]
  own <- own[!duplicated(own$chrom), ]
  if (!length(given)) {
    given <- own$chrom[natural_order(own$chrom)]
    whole <- rep(TRUE, length(given))
  }
  where <- if (named) paste("assembly", assembly) else paste(side, "side")
  if (!length(given)) {
    stop("the file names no contig on its ", where, call. = FALSE)
  }
  unknown <- given[whole & !given %in% own$chrom]
  if (length(unknown)) {
    stop(sprintf("\"%s\" is neither a locus such as \"ctgA:1-1,000\" nor a contig of the file's %s, whose contigs are %s",
      unknown[1], where, paste(utils::head(own$chrom, 20), collapse = ", ")),
      call. = FALSE)
  }
  length <- own$length[match(given, own$chrom)]
  if (anyNA(length[whole])) {
    stop("the file does not state the length of ", given[whole][is.na(
      length[whole])][1], ": give the axis its regions, as loci",
      call. = FALSE)
  }
  do.call(rbind, lapply(seq_along(given), function(i)
    if (whole[i]) data.frame(chrom = given[i], start = 0, end = length[i])
    else regions_of(given[i], NULL, NULL)))
}

# Both axes of a comparative view: the assembly each is on and its regions,
# laid end to end as JBrowse's synteny and dotplot workers lay them, with no
# gap
comparative_axes <- function(config, first, second, assemblies, swap) {
  names <- axis_assemblies(config, assemblies, swap)
  contigs <- memo(function() contigs_of(config))
  list(names = names,
    first = axis_regions(first, contigs, names[1], "query"),
    second = axis_regions(second, contigs, names[2], "target"))
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
#'   location box writes it (`"ctgA:1-50,001"`, 1-based and inclusive), a
#'   contig's name alone for the whole of it, several of them, drawn end to
#'   end as a view lays them out, or a data frame of `chrom`, `start` and
#'   `end`. Unset, every contig of that row's assembly, in the order the file
#'   names them ([jbrowse_contigs()]), as JBrowse shows a whole genome
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
#' @param swap put the file's target side on the top row and its query side
#'   on the bottom: the adapter's two assembly names, or `assemblies`,
#'   reversed. A file guessed from its extension names its sides `"query"`
#'   and `"target"`, so `swap = TRUE` is `assemblies = c("target",
#'   "query")`
#' @param labels what each axis title calls its assembly, named by the
#'   assembly: `c(query = "peach", target = "grape")` for a file guessed from
#'   its extension, or the adapter's own assembly names. An unnamed pair is
#'   the adapter's query and target in that order. Unset, a title shows the
#'   assembly name itself, `"query"` and `"target"` for a guessed file
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
#'     other row is not showing, with their span and their mate's, the
#'     `color` and `alpha` the strip along that row's edge paints them, and
#'     `drawn`, whether they clear the strip's length and contig floors
#'   * `offscreen_contigs`: each contig those mates lie on, by `axis`:
#'     `mate_chrom`, `count`, `aligned_bp`, the `label` the strip's hover
#'     prints, `color` and `alpha`
#'   * `offscreen_strip`: what the strips paint, in pixels from the band's
#'     top left: `marks` (`axis`, `x`, `y`, `width`, `height`, `color`,
#'     `alpha`) and `labels` (`axis`, `label`, `x`, `y` at the baseline,
#'     `color`, `alpha`, `halo`), in a `font_px` font
#'   * `legend`: `label` and `color`, the key the view shows, its title as
#'     the frame's `"title"` attribute
#'   * `ramps`: a continuous field's color bar, as [guides_jbrowse()] draws
#'     it
#'
#'   `query` and `target` are the rows' regions as laid out, `height` the
#'   band's pixels, `alpha` the opacity every ribbon draws at and `field` the
#'   color field painted.
#' @export
#' @examples
#' chain <- system.file("extdata", "hg002_chr8.chain", package = "ggjbrowse")
#' s <- jbrowse_synteny(chain, "chr8_PATERNAL:5,250,001-14,250,000", "chr8_MATERNAL:5,250,001-14,250,000")
#' head(s$alignments)
jbrowse_synteny <- function(file, query = NULL, target = NULL, view = NULL,
                            width = 800, assemblies = NULL, swap = FALSE,
                            labels = NULL, index = NULL) {
  args <- comparative_args(file, index, view)
  axes <- comparative_axes(args$adapterConfig, query, target, assemblies, swap)
  q <- axes$first
  t <- axes$second
  args <- c(args, list(query = axis_args(q, axes$names[1], width),
    target = axis_args(t, axes$names[2], width)))
  r <- call_engine("synteny", args)
  alignments <- snake(as_frame(r$alignments))
  alignments$color <- jbrowse_color(alignments$color)
  ribbons <- snake(as_frame(r$ribbons))
  ribbons$alignment <- ribbons$alignment + 1
  ribbons$query_region <- ribbons$query_region + 1
  ribbons$target_region <- ribbons$target_region + 1
  ribbons$color <- jbrowse_color(ribbons$color)
  offscreen <- colored_frame(r$offscreen)
  offscreen$region <- offscreen$region + 1
  key <- key_frame(unlist(r$legend$label), unlist(r$legend$color),
    r$legendTitle)
  attr(key, "override") <- list(alpha = max(r$alpha, 0.45))
  attr(key, "complete") <- TRUE
  structure(
    list(alignments = alignments, ribbons = ribbons, offscreen = offscreen,
      offscreen_contigs = colored_frame(r$offscreenContigs),
      offscreen_strip = list(marks = colored_frame(r$offscreenStrip$marks),
        labels = colored_frame(r$offscreenStrip$labels),
        font_px = r$offscreenStrip$fontPx),
      legend = key, keys = list(ribbons = key), ramps = ramp_frames(r$ramps),
      field = r$field, alpha = r$alpha, height = r$height,
      fade_thin = r$fadeThin, has_cigar = r$hasCigar,
      assemblies = display_names(unlist(r$assemblyNames), labels,
        args$adapterConfig), width = width,
      query = region_layout(q, 0), target = region_layout(t, 0)),
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

# Where each ribbon's feet sit: px along its upper and lower edges, and the
# heights of those edges. The band's query row is at y = 1, its target at 0.
band_feet <- function(r) {
  list(top1 = r$query_px_start, top2 = r$query_px_end,
    bottom1 = r$target_px_start, bottom2 = r$target_px_end,
    y_top = rep(1, nrow(r)), y_bottom = rep(0, nrow(r)))
}

# The pieces `pts` of ribbons `r` as one frame, each piece a `ribbon` group
# in its ribbon's color, alpha, kind and alignment
ribbon_marks <- function(r, pts) {
  if (!nrow(r)) {
    return(data.frame(x = numeric(0), y = numeric(0), ribbon = integer(0),
      color = jbrowse_color(character(0)), alpha = numeric(0)))
  }
  d <- do.call(rbind, pts)
  for (column in c("color", "alpha", "kind", "alignment")) {
    d[[column]] <- r[[column]][d$ribbon]
  }
  d
}

# Each filled ribbon's outline, top to bottom on its first edge and back up
# its second
ribbon_polygons <- function(r, curves, feet = band_feet(r)) {
  ribbon_marks(r, lapply(seq_len(nrow(r)), function(i) {
    down <- edge_points(feet$top1[i], feet$bottom1[i], feet$y_top[i],
      feet$y_bottom[i], curves)
    up <- edge_points(feet$bottom2[i], feet$top2[i], feet$y_bottom[i],
      feet$y_top[i], curves)
    data.frame(x = c(down$x, up$x), y = c(down$y, up$y), ribbon = i)
  }))
}

# Each thin ribbon's centerline, from the middle of its top edge to the
# middle of its bottom edge
ribbon_lines <- function(r, curves, feet = band_feet(r)) {
  ribbon_marks(r, lapply(seq_len(nrow(r)), function(i) {
    e <- edge_points((feet$top1[i] + feet$top2[i]) / 2,
      (feet$bottom1[i] + feet$bottom2[i]) / 2, feet$y_top[i],
      feet$y_bottom[i], curves)
    data.frame(x = e$x, y = e$y, ribbon = i)
  }))
}

#' A synteny band's ribbons as ggplot layers
#'
#' `geom_jbrowse_ribbons()` draws each ribbon between the query row at the
#' top of the panel (`y = 1`) and the target row at the bottom (`y = 0`) as
#' the band draws it: a filled trapezoid in its color at its opacity, or a
#' one-pixel line where it is thinner than a pixel, with the band's location
#' ticks where the view draws them. `geom_jbrowse_offscreen_mates()` draws
#' the strip the band hangs off a row's edge for the blocks whose mate the
#' other row is not showing: the display's own marks, 6 px tall, grey or in
#' the mate contig's ribbon color where the ribbons are colored by that
#' row's contigs, and the contig names it has room for.
#' @param mapping aesthetics that override the defaults
#' @param data a [jbrowse_synteny()] result, `NULL` for the plot's,
#'   or a data frame in place of the table the geom draws from the plot's
#'   result, such as a filtered copy of it
#' @param curves draw each ribbon edge as the view's `drawCurves` does, a
#'   cubic curve vertical at both rows
#' @param ... passed to [ggplot2::layer()] params
#' @return a list of ggplot layers
#' @export
#' @examples
#' chain <- system.file("extdata", "hg002_chr8.chain", package = "ggjbrowse")
#' s <- jbrowse_synteny(chain, "chr8_PATERNAL:5,250,001-14,250,000", "chr8_MATERNAL:5,250,001-14,250,000")
#' library(ggplot2)
#' ggplot(s) + geom_jbrowse_ribbons(curves = TRUE) + scale_fill_jbrowse(s)
geom_jbrowse_ribbons <- function(mapping = NULL, data = NULL, curves = FALSE,
                                 ...) {
  data <- replacing(data, "ribbons", "geom_jbrowse_ribbons")
  filled <- function(x) ribbon_polygons(x$ribbons[!x$ribbons$thin, ,
    drop = FALSE], curves)
  lines <- function(x) ribbon_lines(x$ribbons[x$ribbons$thin, ,
    drop = FALSE], curves)
  list(
    jbrowse_layer(ggplot2::GeomPolygon, data, filled,
      for_geom(mapping, ggplot2::GeomPolygon),
      ggplot2::aes(x = x, y = y, group = ribbon, fill = color,
        alpha = I(alpha)), list(colour = NA, ...), keyed = "complete"),
    jbrowse_layer(ggplot2::GeomPath, data, lines,
      for_geom(mapping, ggplot2::GeomPath),
      ggplot2::aes(x = x, y = y, group = ribbon, colour = color,
        alpha = I(alpha)), list(linewidth = 0.3, ...), keyed = FALSE))
}

#' @rdname geom_jbrowse_ribbons
#' @param labels name each stretch of marks with its contig where the strip
#'   has room, as the band does
#' @export
geom_jbrowse_offscreen_mates <- function(mapping = NULL, data = NULL,
                                         labels = TRUE, ...) {
  data <- replacing(data, NULL, "geom_jbrowse_offscreen_mates")
  text <- function(halo) function(x) strip_labels(x, halo)
  text_layer <- function(halo) {
    jbrowse_layer(GeomDisplayText, data, text(halo), NULL,
      ggplot2::aes(x = x, y = y, label = label, colour = color,
        alpha = I(alpha), size = I(size), px_wide = I(px_wide),
        px_tall = I(px_tall)),
      list(hjust = 0, vjust = 0), keyed = FALSE)
  }
  c(list(jbrowse_layer(ggplot2::GeomRect, data, strip_marks,
      for_geom(mapping, ggplot2::GeomRect),
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = color, alpha = I(alpha)),
      list(colour = NA, ...), keyed = FALSE)),
    if (labels) list(text_layer(TRUE), text_layer(FALSE)))
}

# The strip's marks on the band, which runs from the target row (y = 0) up
# to the query row (y = 1) over `height` pixels
strip_marks <- function(x) {
  m <- x$offscreen_strip$marks
  m$xmin <- m$x
  m$xmax <- m$x + m$width
  m$ymax <- 1 - m$y / x$height
  m$ymin <- 1 - (m$y + m$height) / x$height
  m
}

HALO_PX <- 3

# The strip's labels at their baselines, or their halo: each label again in
# the band's ground, ringed round it half the halo's width out
strip_labels <- function(x, halo) {
  l <- x$offscreen_strip$labels
  n <- nrow(l)
  l$y <- 1 - l$y / x$height
  l$size <- rep(x$offscreen_strip$font_px, n)
  l$px_wide <- rep(x$width, n)
  l$px_tall <- rep(x$height, n)
  if (!halo) return(l)
  angles <- seq(0, 2 * pi, length.out = 9)[-9]
  do.call(rbind, lapply(angles, function(a) {
    h <- l
    h$x <- h$x + HALO_PX / 2 * cos(a)
    h$y <- h$y + HALO_PX / 2 * sin(a) / x$height
    h$color <- h$halo
    h$alpha <- rep(1, n)
    h
  }))
}

OFFSCREEN_KEY_ROWS <- 10

# A key per row whose strip draws marks: a chip of each contig in the color
# the strip paints it, beside the line the strip's hover prints for it, the
# most sequence first
offscreen_guides <- function(x) {
  drawn <- x$offscreen[x$offscreen$drawn, , drop = FALSE]
  if (!nrow(drawn)) return(NULL)
  contigs <- x$offscreen_contigs
  guides <- lapply(intersect(c("query", "target"), drawn$axis), function(ax) {
    on <- contigs[contigs$axis == ax &
      contigs$mate_chrom %in% drawn$mate_chrom[drawn$axis == ax], ]
    on <- on[order(-on$aligned_bp), ]
    labels <- on$label
    fill <- scales::alpha(as.character(on$color), on$alpha)
    rest <- length(labels) - OFFSCREEN_KEY_ROWS
    if (rest > 0) {
      labels <- c(labels[seq_len(OFFSCREEN_KEY_ROWS)], sprintf("+%d more", rest))
      fill <- c(fill[seq_len(OFFSCREEN_KEY_ROWS)], NA)
    }
    chip_guide(labels, fill, sprintf("Off-screen mates, %s row", ax))
  })
  names(guides) <- sprintf("jbrowse_offscreen_%d", seq_along(guides))
  do.call(ggplot2::guides, guides)
}

# Chips and their labels, a row each, as a guide of its own
chip_guide <- function(labels, fill, title) {
  row_in <- 0.18
  chip_in <- 0.14
  n <- length(labels)
  y <- grid::unit(1, "npc") - grid::unit((seq_len(n) - 0.5) * row_in, "in")
  text <- grid::gpar(fontsize = 8, col = "grey30")
  grob <- grid::gTree(children = grid::gList(
    grid::rectGrob(x = 0, y = y, width = grid::unit(chip_in, "in"),
      height = grid::unit(chip_in, "in"), hjust = 0,
      gp = grid::gpar(fill = fill, col = NA)),
    grid::textGrob(labels, x = grid::unit(chip_in + 0.06, "in"), y = y,
      hjust = 0, gp = text)))
  widest <- grid::textGrob(labels[which.max(nchar(labels))], gp = text)
  ggplot2::guide_custom(grob, title = title,
    width = grid::unit(chip_in + 0.06, "in") + grid::grobWidth(widest),
    height = grid::unit(n * row_in, "in"))
}

# Whether an axis names its contigs rather than numbering positions, as
# JBrowse's dotplot labels each region by name: an axis of more than four
# regions, or of several where one spans under a tenth of it
names_contigs <- function(layout) {
  w <- layout$end - layout$start
  nrow(layout) > 4 || (nrow(layout) > 1 && min(w) < sum(w) / 10)
}

# The pixels apart two contig names stand at least, the dotplot's LABEL_PX
LABEL_PX <- 12

# Which labels at `at` px to keep: each its region's, widest region first,
# where it stands LABEL_PX clear of every label kept before it
thin_labels <- function(at, widths) {
  kept <- numeric(0)
  keep <- logical(length(at))
  for (i in order(widths, decreasing = TRUE)) {
    if (all(abs(kept - at[i]) >= LABEL_PX)) {
      kept <- c(kept, at[i])
      keep[i] <- TRUE
    }
  }
  keep
}

# The labels on an axis of `width` px, `px(layout, region, coord)` placing a
# position on it: round positions in each region, labelled in the region's
# own coordinates, or where names_contigs() says, each contig's name at its
# middle
layout_breaks <- function(layout, px, width) {
  r <- layout
  if (names_contigs(r)) {
    n <- nrow(r)
    mid <- px(r, seq_len(n), (r$start + r$end) / 2)
    edges <- px(r, c(1, n), c(r$start[1], r$end[n]))
    keep <- thin_labels((mid - edges[1]) / diff(edges) * width,
      r$end - r$start)
    return(list(breaks = mid[keep], labels = r$chrom[keep]))
  }
  at <- region_breaks(r)
  list(
    breaks = unlist(lapply(seq_len(nrow(r)), function(i)
      px(r, rep(i, length(at[[i]])), at[[i]]))),
    labels = format(unlist(at), big.mark = ",", trim = TRUE,
      scientific = FALSE))
}

layout_label <- function(layout, assembly = "") {
  if (names_contigs(layout)) {
    return(sprintf("%s%d regions, end to end", if (nzchar(assembly))
      paste0(assembly, ", ") else "", nrow(layout)))
  }
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
    band_px(r, i, at, x$width), x$width)
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
    guides_jbrowse(x) +
    offscreen_guides(x) +
    ggplot2::coord_cartesian(xlim = c(0, x$width), ylim = c(0, 1),
      expand = FALSE) +
    synteny_x(x) +
    ggplot2::labs(y = NULL) +
    ggplot2::theme_minimal() +
    no_y +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    upright_names(x$target, "x.bottom") +
    upright_names(x$query, "x.top")
}

# An axis's contig names turned to run along the axis's normal, as the
# dotplot's horizontal axis draws them, so close names do not overlap
upright_names <- function(layout, axis) {
  if (!names_contigs(layout)) return(NULL)
  hjust <- if (axis == "x.top") 0 else 1
  do.call(ggplot2::theme, stats::setNames(list(ggplot2::element_text(
    angle = 90, hjust = hjust, vjust = 0.5)), paste0("axis.text.", axis)))
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
  lapply(tracks, function(t) {
    if (is.null(t$regions)) return(t)
    t$regions <- layout
    t$region$start <- min(layout$cum_start)
    t$region$end <- max(layout$cum_end)
    t
  })
}

#' Draw a pairwise alignment as JBrowse's linear synteny view
#'
#' The band of ribbons between the two rows, the query row's ruler above it
#' and the target row's below, and any tracks over each row: a
#' [jbrowse_features()], [jbrowse_wiggle()], [jbrowse_variants()] or
#' [jbrowse_alignments()] result over that row's regions, drawn as
#' [jbrowse_stack()] draws it, the top row's tracks above the band and the
#' bottom row's below.
#'
#' @param object a [jbrowse_synteny()] result
#' @param above,below tracks over the top row's regions, drawn above the
#'   band, and over the bottom row's, drawn below it: each a result or a
#'   list of them
#' @param curves draw each ribbon edge as a curve, as the view's `drawCurves`
#'   does
#' @param ... must be empty, so a misspelt argument errors
#' @return a ggplot of the band alone, or with tracks a patchwork, whose
#'   height and width [jbrowse_height()] and [jbrowse_width()] measure
#' @exportS3Method ggplot2::autoplot
#' @examples
#' e <- function(f) system.file("extdata", f, package = "ggjbrowse")
#' s <- jbrowse_synteny(e("hg002_chr8.chain"), "chr8_PATERNAL:5,250,001-14,250,000",
#'   "chr8_MATERNAL:5,250,001-14,250,000", view = list(color = list(field = "strand")),
#'   labels = c(query = "paternal", target = "maternal"))
#' autoplot(s)
#' genes <- jbrowse_features(e("hg002_mat_genes.gff3.gz"), "chr8_MATERNAL:5,250,001-14,250,000",
#'   display = list(showOnlyGenes = TRUE))
#' autoplot(s, below = genes)
autoplot.jbrowse_synteny <- function(object, above = NULL, below = NULL,
                                     curves = FALSE, ...) {
  rlang::check_dots_empty()
  band <- synteny_band(object, curves)
  above <- row_tracks(above, object$query, "top")
  below <- row_tracks(below, object$target, "bottom")
  if (!length(above) && !length(below)) return(band)
  bands <- c(lapply(above, track_bands),
    list(list(plots = list(band), heights = object$height)),
    lapply(below, track_bands))
  own <- sum(lengths(lapply(bands[seq_along(above)], `[[`, "plots"))) + 1
  stack_bands(bands, axis_at = own)
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
#'   side and the vertical its target, unless `assemblies` or `swap` says
#'   otherwise. Unset, every contig of that axis's assembly, so
#'   `jbrowse_dotplot(paf)` is the whole genome against the whole genome
#' @param view the view's settings, as a `DotplotView` snapshot writes them:
#'   a list with any of `color`, `opacity`, `minAlignmentLength`,
#'   `minIdentity`, `drawCigar` and `lineWidth`
#' @param width,height the plot's pixels, which decide which blocks are wide
#'   enough for their CIGAR to show
#' @param assemblies which assembly each axis is on, `c(x, y)`
#' @param swap put the file's target side on the horizontal axis and its
#'   query on the vertical, as [jbrowse_synteny()]'s `swap` does the rows
#' @param labels what each axis title calls its assembly, as for
#'   [jbrowse_synteny()]: `c(query = "peach", target = "grape")`
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
#'   * `ramps`: a continuous field's color bar, as [guides_jbrowse()] draws
#'     it
#'
#'   `x` and `y` are the axes' regions as laid out, `line_width` the pixels
#'   each line is drawn at, `alpha` the opacity, and `skipped` counts the
#'   blocks neither axis places.
#' @export
#' @examples
#' chain <- system.file("extdata", "hg002_chr8.chain", package = "ggjbrowse")
#' d <- jbrowse_dotplot(chain, "chr8_PATERNAL:5,250,001-14,250,000", "chr8_MATERNAL:5,250,001-14,250,000")
#' head(d$segments)
jbrowse_dotplot <- function(file, x = NULL, y = NULL, view = NULL,
                            width = 600, height = 600, assemblies = NULL,
                            swap = FALSE, labels = NULL, index = NULL) {
  args <- comparative_args(file, index, view)
  axes <- comparative_axes(args$adapterConfig, x, y, assemblies, swap)
  xr <- axes$first
  yr <- axes$second
  args <- c(args, list(x = axis_args(xr, axes$names[1], width),
    y = axis_args(yr, axes$names[2], height)))
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
  attr(key, "complete") <- TRUE
  structure(
    list(alignments = alignments, segments = segments, legend = key,
      keys = list(segments = key), ramps = ramp_frames(r$ramps),
      field = r$field,
      alpha = r$alpha, line_width = r$lineWidth,
      assemblies = display_names(unlist(r$assemblyNames), labels,
        args$adapterConfig), skipped = r$skipped,
      width = width, height = height,
      x = region_layout(xr, 0), y = region_layout(yr, 0)),
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
#' @param data a [jbrowse_dotplot()] result, `NULL` for the plot's,
#'   or a data frame in place of the table the geom draws from the plot's
#'   result, such as a filtered copy of it
#' @param ... passed to [ggplot2::layer()] params, such as `linewidth`
#' @return a ggplot layer
#' @export
#' @examples
#' paf <- system.file("extdata", "hprc_amylase.paf", package = "ggjbrowse")
#' d <- jbrowse_dotplot(paf, "GRCh38#0#chr1:103,520,894-103,832,637",
#'   "HG00133#1#CM090045.1:103,669,666-103,981,330")
#' library(ggplot2)
#' ggplot(d) + geom_jbrowse_dotplot(aes(colour = op))
geom_jbrowse_dotplot <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "segments", "geom_jbrowse_dotplot")
  jbrowse_layer(ggplot2::GeomSegment, data, function(x) x$segments,
    for_geom(mapping, ggplot2::GeomSegment),
    ggplot2::aes(x = x_start, xend = x_end, y = y_start, yend = y_end,
      colour = color, alpha = I(alpha)),
    with_defaults(list(...), list(linewidth = px_to_mm(2.5),
      lineend = "round")), keyed = "complete")
}

# An axis of regions laid end to end: breaks at round positions in each,
# labelled in the region's own coordinates
layout_scale <- function(scale, layout, width) {
  b <- layout_breaks(layout, function(r, i, at) at - r$start[i] + r$offset[i],
    width)
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
#' @param ... must be empty, so a misspelt argument errors
#' @return a ggplot
#' @exportS3Method ggplot2::autoplot
#' @examples
#' chain <- system.file("extdata", "hg002_chr8.chain", package = "ggjbrowse")
#' autoplot(jbrowse_dotplot(chain, view = list(color = list(field = "strand")),
#'   labels = c(query = "paternal", target = "maternal")))
autoplot.jbrowse_dotplot <- function(object, ...) {
  rlang::check_dots_empty()
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
    guides_jbrowse(x) +
    layout_scale(ggplot2::scale_x_continuous, x$x, x$width) +
    layout_scale(ggplot2::scale_y_continuous, x$y, x$height) +
    ggplot2::coord_cartesian(
      xlim = c(min(x$x$cum_start), max(x$x$cum_end)),
      ylim = c(min(x$y$cum_start), max(x$y$cum_end)), expand = FALSE) +
    ggplot2::labs(x = layout_label(x$x, x$assemblies[1]),
      y = layout_label(x$y, x$assemblies[2])) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey70")) +
    upright_names(x$x, "x.bottom")
}

# Names in the order a reader expects them: chr2 before chr10, each run of
# digits compared as a number
natural_order <- function(x) {
  parts <- regmatches(x, gregexpr("[0-9]+|[^0-9]+", x))
  key <- vapply(parts, function(p) paste(vapply(p, function(s)
    if (grepl("^[0-9]", s)) formatC(as.numeric(s), width = 12, flag = "0",
      format = "d") else s, ""), collapse = ""), "")
  order(key, x)
}

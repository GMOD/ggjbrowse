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

# The engine's alignments, each keyed by `alignment`, the number the tables
# drawn from it name it by
keyed_alignments <- function(columns) {
  d <- snake(as_frame(columns))
  cbind(alignment = seq_len(nrow(d)), d)
}

snake_frame <- function(columns) snake(as_frame(columns))

# The engine's camelCase columns as this package's snake_case
snake <- function(d) {
  names(d) <- gsub("([A-Z])", "_\\L\\1", names(d), perl = TRUE)
  d
}


comparative_args <- function(file, index, display) {
  args <- list(adapterConfig = adapter_of(file, function()
    synteny_adapter(file, index)))
  if (length(display)) {
    display$color <- config_color(display$color)
    args$view <- display
  }
  args
}

# The values a view's enumerated settings take
SYNTENY_CHOICES <- list(cigarMode = c("full", "matches", "off"),
  fadeThinAlignmentsMode = c("auto", "on", "off"))

# Stops where a setting holds a value its view does not take, as the
# linear displays' config schemas refuse one
check_choices <- function(display, choices) {
  for (slot in intersect(names(choices), names(display))) {
    value <- display[[slot]]
    if (length(value) == 1 && value %in% choices[[slot]]) next
    rlang::abort(sprintf("`display$%s` is %s; the view takes %s", slot,
      paste(format(value), collapse = ", "), paste(encodeString(
        choices[[slot]], quote = "\""), collapse = ", ")),
      class = "ggjbrowse_bad_display")
  }
}

# The engine's tables and the warnings it raised reading them, each also
# raised as an R warning
warned_call <- function(table, args) {
  said <- character(0)
  value <- withCallingHandlers(call_engine(table, args), warning = function(w)
    said <<- c(said, conditionMessage(w)))
  list(value = value, warnings = unique(said))
}

# A comparative result: its tables, then the slots every result holds, then
# its own drawing state
comparative_result <- function(cls, primary, tables, regions, width, display,
                               keys, ramps, warnings, state = list()) {
  structure(c(tables, list(regions = regions, width = width,
      display = display, keys = keys, ramps = ramps, y_scale = NULL,
      fields = list(), warnings = warnings), state),
    primary = primary, class = c(cls, "jbrowse_result"))
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

# Each side's contigs in the order the file first names them: `side`
# (`"query"` or `"target"`), `assembly`, `chrom` and `length`, `NA` where the
# file does not say
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

# The regions an axis shows, each with its missing-contig message, `""`
# where the file has its contig: any region a linear reader takes, a data
# frame's rows narrowed to the axis's assembly where it names an `assembly`,
# and with none given every contig of the axis's assembly in natural order.
# `side` picks the file's side where the axis names no assembly the file has.
axis_regions <- function(given, contigs, assembly, side) {
  all <- contigs()
  named <- length(assembly) && !is.na(assembly) && assembly %in% all$assembly
  name <- if (named) assembly else side
  own <- all[if (named) all$assembly == assembly else all$side == side, ]
  own <- own[!duplicated(own$chrom), ]
  where <- if (name %in% c("query", "target")) paste(name, "side")
    else paste("assembly", name)
  bad <- function(...) rlang::abort(paste0(...), class = "ggjbrowse_bad_region")
  if (is.null(given)) {
    if (!nrow(own)) bad("the file names no contig on its ", where)
    given <- own$chrom[natural_order(own$chrom)]
  }
  if (is.data.frame(given) && "assembly" %in% names(given)) {
    mine <- given[as.character(given[["assembly"]]) %in% name, , drop = FALSE]
    if (!nrow(mine)) {
      bad("the region frame names assembly ", paste(unique(
        given[["assembly"]]), collapse = ", "), ", none of them this axis's ",
        name)
    }
    given <- mine
  }
  regions <- regions_of(given, function()
    data.frame(chrom = own$chrom, start = 0, end = own$length))
  if (anyNA(regions$end)) {
    bad("the file does not state the length of ", regions$chrom[
      is.na(regions$end)][1], ": give the axis its regions, as loci")
  }
  regions$missing <- if (!nrow(own)) "" else ifelse(regions$chrom %in%
    own$chrom, "", sprintf("%s is not a contig of the file's %s, whose contigs are %s",
      regions$chrom, where, paste(utils::head(own$chrom, 20), collapse = ", ")))
  regions
}

# Both axes of a comparative view: the assembly each is on and its regions.
# Axes none of whose contigs the file has stop; a region on a contig the file
# lacks beside others warns, its warning in `warnings`. A row with no contig
# of the file still draws the strip of the other row's off-screen mates.
comparative_axes <- function(config, first, second, assemblies, swap) {
  names <- axis_assemblies(config, assemblies, swap)
  contigs <- memo(function() contigs_of(config))
  first <- axis_regions(first, contigs, names[1], "query")
  second <- axis_regions(second, contigs, names[2], "target")
  missing <- check_missing(c(first$missing, second$missing),
    c(first$chrom, second$chrom))
  list(names = names, first = first, second = second,
    warnings = attr(missing, "warning"))
}

# An axis's rows of a result's `regions`: each region named by the `axis` it
# lays out on and the `assembly` it is on, laid end to end from the axis's
# first region as JBrowse's synteny and dotplot workers lay them, with no gap
# unless `gap_frac` gives one
axis_rows <- function(regions, axis, assembly, gap_frac = 0) {
  l <- region_layout(regions[c("chrom", "start", "end")], gap_frac)
  n <- nrow(l)
  data.frame(axis = rep(axis, n),
    assembly = rep(if (length(assembly)) assembly else NA_character_, n), l,
    missing = regions$missing %||% rep("", n))
}

# The rows of `regions` on `axis`
on_axis_of <- function(regions, axis) {
  r <- regions[regions$axis == axis, , drop = FALSE]
  rownames(r) <- NULL
  r
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
#' @param query,target the regions each row shows, drawn end to end as a
#'   view lays them out, in any form a linear reader's `region` takes: loci
#'   as JBrowse's location box writes them (`"ctgA:1-50,001"`, 1-based and
#'   inclusive), contig names alone for whole contigs, a data frame of
#'   `chrom`, `start` and `end` (0-based, half-open), or a GRanges. A frame
#'   with an `assembly` column, as [jbrowse_contigs()] and a result's
#'   `regions` write it, gives the row the regions on its assembly. With
#'   `NULL`, the default, every contig of that row's assembly in natural
#'   order, `chr2` before `chr10`, as JBrowse shows a whole genome. A region
#'   on a contig the file lacks reads empty and warns, and rows with none of
#'   the file's contigs between them stop
#' @param display the view's settings, as a `LinearSyntenyView` snapshot
#'   writes them: a list with any of `color` (a CSS color, or the color object's
#'   slots, such as `list(field = "strand")`, `list(field = "identity")`,
#'   `list(field = "query")` or `list(field = "mapq")`), `opacity` (a number,
#'   or an opacity object with a `field`), `minAlignmentLength`, `cigarMode`
#'   (`"full"` for colored indels, `"matches"` for the matches alone with the
#'   indels left clear, `"off"` for the block alone), `drawLocationMarkers`,
#'   `fadeThinAlignmentsMode` (`"auto"`, `"on"`, `"off"`) and `height`, the
#'   band's pixels
#' @param width the width in pixels both rows draw their regions across,
#'   1,500 as for the linear readers, which decides which CIGAR ops are wide
#'   enough to draw and which ribbons are thin enough to fade
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
#' @return a `jbrowse_synteny` list, beside the slots [jbrowse_result] lists:
#'   * `alignments`: each block, in the band's paint order: `alignment`,
#'     its key, `id`, `name`,
#'     `query_chrom`, `query_start`, `query_end`, `target_chrom`,
#'     `target_start`, `target_end`, `strand`, `length`, the measurements the
#'     file carries (`identity`, `mean_identity`, `mapq`, and any column the
#'     adapter's `attributeColumns` declares), `color` and `alpha`
#'   * `ribbons`: what the band draws, in paint order, each a trapezoid from
#'     a span on the query row to a span on the target row: `alignment` (the
#'     block's key in `alignments`), `kind` (`"match"`, `"insertion"`,
#'     `"deletion"`, `"skip"`, or `"marker"` for a location tick),
#'     `query_region`, `query_chrom`, `query_start`, `query_end`,
#'     `target_region`, `target_chrom`, `target_start`, `target_end`, each
#'     `*_region` the row of `regions` its foot lies in, where
#'     the ribbon joins `query_start` to `target_start` and `query_end` to
#'     `target_end`, so a reverse-strand block crosses; `color`, `alpha`, and
#'     `thin`, whether the band draws it as a one-pixel line
#'   * `offscreen`: the blocks on one `axis` whose mate lies on a contig the
#'     other row is not showing, the `region` they lie in, their span and
#'     their mate's, the
#'     `color` and `alpha` the strip along that row's edge paints them, and
#'     `drawn`, whether they clear the strip's length and contig floors
#'   * `offscreen_contigs`: each contig those mates lie on, by `axis`:
#'     `mate_chrom`, `count`, `aligned_bp`, the `label` the strip's hover
#'     prints, `color` and `alpha`
#'   * `offscreen_strip`: what the strips paint, in pixels from the band's
#'     top left: `marks` (`axis`, `x`, `y`, `width`, `height`, `color`,
#'     `alpha`) and `labels` (`axis`, `label`, `x`, `y` at the baseline,
#'     `color`, `alpha`, `halo`), in a `font_px` font
#'
#'   Assigning `alignments` drops the ribbons of the blocks it leaves out;
#'   the off-screen tables summarize the band's edges as read, and stay.
#'   `keys$ribbons` is the key the view shows, `ramps` a continuous field's
#'   color bar, as [guides_jbrowse()] draws it. `regions` lays out both rows, its `axis` `"query"` for the top row and
#'   `"target"` for the bottom. `display` holds the view's settings as the
#'   engine resolves them, `color` and `opacity` as config objects and
#'   `height` the band's pixels.
#' @export
#' @examples
#' chain <- system.file("extdata", "hg002_chr8.chain", package = "ggjbrowse")
#' s <- jbrowse_synteny(chain,
#'   "chr8_PATERNAL:5,250,001-14,250,000",
#'   "chr8_MATERNAL:5,250,001-14,250,000")
#' head(s$alignments)
jbrowse_synteny <- function(file, query = NULL, target = NULL,
                            display = NULL, width = DEFAULT_WIDTH,
                            assemblies = NULL, swap = FALSE, labels = NULL,
                            index = NULL) {
  check_choices(display, SYNTENY_CHOICES)
  args <- comparative_args(file, index, display)
  axes <- comparative_axes(args$adapterConfig, query, target, assemblies, swap)
  q <- axes$first
  t <- axes$second
  args <- c(args, list(query = axis_args(q, axes$names[1], width),
    target = axis_args(t, axes$names[2], width)))
  read <- warned_call("synteny", args)
  r <- read$value
  regions <- rbind(axis_rows(q, "query", axes$names[1]),
    axis_rows(t, "target", axes$names[2]))
  alignments <- keyed_alignments(r$alignments)
  ribbons <- snake(as_frame(r$ribbons))
  ribbons$alignment <- as.integer(ribbons$alignment) + 1L
  ribbons$query_region <- ribbons$query_region + 1
  ribbons$target_region <- ribbons$target_region + 1 + nrow(q)
  offscreen <- snake_frame(r$offscreen)
  offscreen$region <- offscreen$region + 1 +
    ifelse(offscreen$axis == "target", nrow(q), 0)
  key <- key_frame(unlist(r$legend$label), unlist(r$legend$color),
    r$legendTitle)
  attr(key, "alpha") <- max(r$alpha, 0.45)
  contigs <- snake_frame(r$offscreenContigs)
  comparative_result("jbrowse_synteny", "ribbons",
    list(alignments = alignments, ribbons = ribbons, offscreen = offscreen,
      offscreen_contigs = contigs),
    regions, width, r$config,
    c(list(ribbons = key), offscreen_keys(offscreen, contigs)),
    ramp_frames(r$ramps),
    unique(c(axes$warnings, read$warnings)),
    list(offscreen_strip = list(marks = snake_frame(r$offscreenStrip$marks),
        labels = snake_frame(r$offscreenStrip$labels),
        font_px = r$offscreenStrip$fontPx),
      fade_thin = r$fadeThin, has_cigar = r$hasCigar,
      assemblies = display_names(unlist(r$assemblyNames), labels,
        args$adapterConfig)))
}

#' @export
print.jbrowse_synteny <- function(x, ...) {
  cat(sprintf(
    "<jbrowse_synteny> %s over %s, %d alignments, %d ribbons, %d off screen\n",
    paste(loci_of_layout(on_axis_of(x$regions, "query")), collapse = ", "),
    paste(loci_of_layout(on_axis_of(x$regions, "target")), collapse = ", "),
    nrow(x$alignments), nrow(x$ribbons), nrow(x$offscreen)))
  invisible(x)
}

# A position in a row of `regions` as a pixel of the band: its place on its
# axis's regions laid end to end, as a share of the axis's span, across the
# view's width
band_px <- function(regions, axis, region, coord, width) {
  lo <- tapply(regions$cum_start, regions$axis, min)[axis]
  hi <- tapply(regions$cum_end, regions$axis, max)[axis]
  shift <- regions$offset[region] - regions$start[region]
  unname((coord + shift - lo) / (hi - lo) * width)
}

# The ribbons and off-screen marks with each foot at its pixel on its row
synteny_on_axis <- function(x) {
  r <- x$ribbons
  if (is_on_axis(r)) return(x)
  px <- function(axis, region, coord) band_px(x$regions, axis, region, coord,
    x$width)
  r$query_px_start <- px("query", r$query_region, r$query_start)
  r$query_px_end <- px("query", r$query_region, r$query_end)
  r$target_px_start <- px("target", r$target_region, r$target_start)
  r$target_px_end <- px("target", r$target_region, r$target_end)
  x$ribbons <- structure(r, on_axis = TRUE)
  o <- x$offscreen
  o$px_start <- px(o$axis, o$region, o$start)
  o$px_end <- px(o$axis, o$region, o$end)
  o$y <- ifelse(o$axis == "query", 1, 0)
  x$offscreen <- structure(o, on_axis = TRUE)
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
      color = character(0), alpha = numeric(0)))
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
#' s <- jbrowse_synteny(chain,
#'   "chr8_PATERNAL:5,250,001-14,250,000",
#'   "chr8_MATERNAL:5,250,001-14,250,000")
#' library(ggplot2)
#' ggplot(s) + geom_jbrowse_ribbons(curves = TRUE) + guides_jbrowse(s)
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
      ggplot2::aes(x = x, y = y, group = ribbon, fill = I(color),
        alpha = I(alpha)), list(colour = NA, ...)),
    jbrowse_layer(ggplot2::GeomPath, data, lines,
      for_geom(mapping, ggplot2::GeomPath),
      ggplot2::aes(x = x, y = y, group = ribbon, colour = I(color),
        alpha = I(alpha)), list(linewidth = 0.3, ...)))
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
      ggplot2::aes(x = x, y = y, label = label, colour = I(color),
        alpha = I(alpha), size = I(size), px_wide = I(px_wide),
        px_tall = I(px_tall)),
      list(hjust = 0, vjust = 0))
  }
  c(list(jbrowse_layer(ggplot2::GeomRect, data, strip_marks,
      for_geom(mapping, ggplot2::GeomRect),
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = I(color), alpha = I(alpha)),
      list(colour = NA, ...))),
    if (labels) list(text_layer(TRUE), text_layer(FALSE)))
}

# The strip's marks on the band, which runs from the target row (y = 0) up
# to the query row (y = 1) over `height` pixels
strip_marks <- function(x) {
  m <- x$offscreen_strip$marks
  m$xmin <- m$x
  m$xmax <- m$x + m$width
  px <- x$display$height
  m$ymax <- 1 - m$y / px
  m$ymin <- 1 - (m$y + m$height) / px
  m
}

HALO_PX <- 3

# The strip's labels at their baselines, or their halo: each label again in
# the band's ground, ringed round it half the halo's width out
strip_labels <- function(x, halo) {
  l <- x$offscreen_strip$labels
  n <- nrow(l)
  px <- x$display$height
  l$y <- 1 - l$y / px
  l$size <- rep(x$offscreen_strip$font_px, n)
  l$px_wide <- rep(x$width, n)
  l$px_tall <- rep(px, n)
  if (!halo) return(l)
  angles <- seq(0, 2 * pi, length.out = 9)[-9]
  do.call(rbind, lapply(angles, function(a) {
    h <- l
    h$x <- h$x + HALO_PX / 2 * cos(a)
    h$y <- h$y + HALO_PX / 2 * sin(a) / px
    h$color <- h$halo
    h$alpha <- rep(1, n)
    h
  }))
}

OFFSCREEN_KEY_ROWS <- 10

# A key per row whose strip draws marks: a chip of each contig in the color
# the strip paints it, beside the line the strip's hover prints for it, the
# most sequence first
offscreen_keys <- function(offscreen, contigs) {
  drawn <- offscreen[offscreen$drawn, , drop = FALSE]
  axes <- intersect(c("query", "target"), drawn$axis)
  keys <- lapply(axes, function(ax) {
    on <- contigs[contigs$axis == ax &
      contigs$mate_chrom %in% drawn$mate_chrom[drawn$axis == ax], ]
    on <- on[order(-on$aligned_bp), ]
    labels <- on$label
    fill <- scales::alpha(on$color, on$alpha)
    rest <- length(labels) - OFFSCREEN_KEY_ROWS
    if (rest > 0) {
      labels <- c(labels[seq_len(OFFSCREEN_KEY_ROWS)], sprintf("+%d more", rest))
      fill <- c(fill[seq_len(OFFSCREEN_KEY_ROWS)], NA)
    }
    key_frame(labels, fill, sprintf("Off-screen mates, %s row", ax))
  })
  stats::setNames(keys, sprintf("offscreen_%s", axes))
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
  on_row <- function(axis) layout_breaks(on_axis_of(x$regions, axis),
    function(r, i, at) band_px(r, axis, i, at, x$width), x$width)
  q <- on_row("query")
  t <- on_row("target")
  ggplot2::scale_x_continuous(breaks = t$breaks, labels = t$labels,
    name = layout_label(on_axis_of(x$regions, "target"), x$assemblies[2]),
    sec.axis = ggplot2::sec_axis(~., breaks = q$breaks, labels = q$labels,
      name = layout_label(on_axis_of(x$regions, "query"), x$assemblies[1])))
}

synteny_band <- function(x, curves) {
  x <- on_axis(x)
  ggplot2::ggplot(x) +
    geom_jbrowse_ribbons(curves = curves) +
    geom_jbrowse_offscreen_mates() +
    guides_jbrowse(x) +
    ggplot2::coord_cartesian(xlim = c(0, x$width), ylim = c(0, 1),
      expand = FALSE) +
    synteny_x(x) +
    ggplot2::labs(y = NULL) +
    ggplot2::theme_minimal() +
    no_y +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    upright_names(on_axis_of(x$regions, "target"), "x.bottom") +
    upright_names(on_axis_of(x$regions, "query"), "x.top")
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
    r <- t$regions
    if (!identical(paste(r$start, r$end), paste(layout$start, layout$end))) {
      stop("a ", row, " track draws over the ", row, " row's regions, ",
        paste(loci_of_layout(layout), collapse = ", "), ", but one reads ",
        paste(loci_of(t), collapse = ", "), call. = FALSE)
    }
  }
  lapply(tracks, function(t) {
    axis <- c("offset", "cum_start", "cum_end")
    t$regions[axis] <- layout[axis]
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
#'   "chr8_MATERNAL:5,250,001-14,250,000",
#'   display = list(color = list(field = "strand")),
#'   labels = c(query = "paternal", target = "maternal"))
#' autoplot(s)
#' genes <- jbrowse_features(e("hg002_mat_genes.gff3.gz"), "chr8_MATERNAL:5,250,001-14,250,000",
#'   display = list(filter = list(genes = "only")))
#' autoplot(s, below = genes)
autoplot.jbrowse_synteny <- function(object, above = NULL, below = NULL,
                                     curves = FALSE, ...) {
  rlang::check_dots_empty()
  band <- synteny_band(object, curves)
  above <- row_tracks(above, on_axis_of(object$regions, "query"), "top")
  below <- row_tracks(below, on_axis_of(object$regions, "target"), "bottom")
  if (!length(above) && !length(below)) return(band)
  bands <- c(lapply(above, track_bands),
    list(list(plots = list(band), heights = object$display$height)),
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
#' @param display the view's settings, as a `DotplotView` snapshot writes
#'   them: a list with any of `color`, `opacity`, `minAlignmentLength`,
#'   `minIdentity`, `drawCigar` and `lineWidth`
#' @param width,height the plot's pixels across and up, which decide which
#'   blocks are wide enough for their CIGAR to show: 1,500 across as for the
#'   linear readers, and up as many as across
#' @param assemblies which assembly each axis is on, `c(x, y)`
#' @param swap put the file's target side on the horizontal axis and its
#'   query on the vertical, as [jbrowse_synteny()]'s `swap` does the rows
#' @param labels what each axis title calls its assembly, as for
#'   [jbrowse_synteny()]: `c(query = "peach", target = "grape")`
#' @return a `jbrowse_dotplot` list, beside the slots [jbrowse_result] lists:
#'   * `alignments`: each block: `alignment`, its key, `name`, `x_chrom`, `x_start`, `x_end`,
#'     `y_chrom`, `y_start`, `y_end`, `strand`, `length`, the measurements
#'     the file carries, `color` and `alpha`
#'   * `segments`: the lines the plot draws, in paint order: `alignment` (the
#'     block's key in `alignments`), `op` (`"match"`, `"insertion"`,
#'     `"deletion"` or `"skip"`), `x_region`, `x_chrom`, `x_start`, `x_end`,
#'     `y_region`, `y_chrom`, `y_start`, `y_end`, from (`x_start`, `y_start`)
#'     to (`x_end`, `y_end`), each `*_region` a row of `regions`, `color` and
#'     `alpha`
#'
#'   Assigning `alignments` drops the segments of the blocks it leaves out.
#'   `keys$segments` is the key the view shows, `ramps` a continuous field's
#'   color bar. `regions` lays out both axes, its `axis` `"x"` or `"y"`. `display`
#'   holds the view's settings as the engine resolves them, `color` and
#'   `opacity` as config objects and `lineWidth` a line's pixels. `skipped` counts the blocks
#'   neither axis places.
#' @export
#' @examples
#' chain <- system.file("extdata", "hg002_chr8.chain", package = "ggjbrowse")
#' d <- jbrowse_dotplot(chain,
#'   "chr8_PATERNAL:5,250,001-14,250,000",
#'   "chr8_MATERNAL:5,250,001-14,250,000")
#' head(d$segments)
jbrowse_dotplot <- function(file, x = NULL, y = NULL, display = NULL,
                            width = DEFAULT_WIDTH, height = width,
                            assemblies = NULL, swap = FALSE, labels = NULL,
                            index = NULL) {
  args <- comparative_args(file, index, display)
  axes <- comparative_axes(args$adapterConfig, x, y, assemblies, swap)
  xr <- axes$first
  yr <- axes$second
  args <- c(args, list(x = axis_args(xr, axes$names[1], width),
    y = axis_args(yr, axes$names[2], height)))
  read <- warned_call("dotplot", args)
  r <- read$value
  alignments <- keyed_alignments(r$alignments)
  segments <- snake(as_frame(r$segments))
  segments$alignment <- as.integer(segments$alignment) + 1L
  segments$x_region <- segments$x_region + 1
  segments$y_region <- segments$y_region + 1 + nrow(xr)
  key <- key_frame(unlist(r$legend$label), unlist(r$legend$color),
    r$legendTitle)
  attr(key, "alpha") <- max(r$alpha, 0.45)
  comparative_result("jbrowse_dotplot", "segments",
    list(alignments = alignments, segments = segments),
    rbind(axis_rows(xr, "x", axes$names[1]), axis_rows(yr, "y", axes$names[2])),
    width, r$config, list(segments = key), ramp_frames(r$ramps),
    unique(c(axes$warnings, read$warnings)),
    list(height = height,
      assemblies = display_names(unlist(r$assemblyNames), labels,
        args$adapterConfig), skipped = r$skipped))
}

#' @export
print.jbrowse_dotplot <- function(x, ...) {
  cat(sprintf(
    "<jbrowse_dotplot> %s against %s, %d alignments, %d segments\n",
    paste(loci_of_layout(on_axis_of(x$regions, "x")), collapse = ", "),
    paste(loci_of_layout(on_axis_of(x$regions, "y")), collapse = ", "),
    nrow(x$alignments), nrow(x$segments)))
  invisible(x)
}

# Each segment's ends on the axes that lay each side's regions end to end
dotplot_on_axis <- function(x) {
  s <- x$segments
  if (is_on_axis(s)) return(x)
  r <- x$regions
  shift <- function(region) r$offset[region] - r$start[region]
  s$x_start <- s$x_start + shift(s$x_region)
  s$x_end <- s$x_end + shift(s$x_region)
  s$y_start <- s$y_start + shift(s$y_region)
  s$y_end <- s$y_end + shift(s$y_region)
  x$segments <- structure(s, on_axis = TRUE)
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
      colour = I(color), alpha = I(alpha)),
    with_defaults(list(...), list(linewidth = px_to_mm(2.5),
      lineend = "round")))
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
#' autoplot(jbrowse_dotplot(chain, display = list(color = list(field = "strand")),
#'   labels = c(query = "paternal", target = "maternal")))
autoplot.jbrowse_dotplot <- function(object, ...) {
  rlang::check_dots_empty()
  x <- on_axis(object)
  h <- on_axis_of(x$regions, "x")
  v <- on_axis_of(x$regions, "y")
  seams <- data.frame(x = layout_seams(h))
  rows <- data.frame(y = layout_seams(v))
  ggplot2::ggplot(x) +
    ggplot2::geom_vline(data = seams, ggplot2::aes(xintercept = x),
      colour = "grey70") +
    ggplot2::geom_hline(data = rows, ggplot2::aes(yintercept = y),
      colour = "grey70") +
    geom_jbrowse_dotplot(linewidth = px_to_mm(x$display$lineWidth)) +
    guides_jbrowse(x) +
    layout_scale(ggplot2::scale_x_continuous, h, x$width) +
    layout_scale(ggplot2::scale_y_continuous, v, x$height) +
    ggplot2::coord_cartesian(xlim = c(min(h$cum_start), max(h$cum_end)),
      ylim = c(min(v$cum_start), max(v$cum_end)), expand = FALSE) +
    ggplot2::labs(x = layout_label(h, x$assemblies[1]),
      y = layout_label(v, x$assemblies[2])) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey70")) +
    upright_names(h, "x.bottom")
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

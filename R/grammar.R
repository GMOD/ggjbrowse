# ggjbrowse in grammar-of-graphics terms. Every table the engine returns is
# data, every mark kind is a geom that draws one table, and the display's
# colors are a column of each table, classed so ggplot gives them an identity
# scale. The defaults draw what the display draws; map any other column and
# ggplot's own scales take over. `autoplot()` only composes these pieces.

#' The color the display paints a mark
#'
#' A `#rrggbb` column that ggplot maps through an identity scale whenever an
#' aesthetic names it, so `aes(fill = color)` paints the display's colors with
#' no scale written. [scale_fill_jbrowse()] adds the display's key.
#' @param x color strings
#' @return `x`, classed `jbrowse_color`
#' @export
jbrowse_color <- function(x) {
  structure(tolower(as.character(x)), class = "jbrowse_color")
}

#' @export
`[.jbrowse_color` <- function(x, i) jbrowse_color(unclass(x)[i])

#' @export
as.character.jbrowse_color <- function(x, ...) unclass(x)

#' @export
print.jbrowse_color <- function(x, ...) print(unclass(x), ...)

#' @exportS3Method ggplot2::scale_type
scale_type.jbrowse_color <- function(x) "identity"

# What a geom draws: the jbrowse result its plot or its `data` holds, on its
# axis and read through `pick`, or a data frame the caller handed it, on its
# plot's axis where it is a result's table
table_from <- function(data, pick, finish = function(d, model) d) {
  if (inherits(data, "jbrowse_result")) {
    return(finish(pick(on_axis(data)), data))
  }
  if (is.null(data)) {
    return(function(plot_data) {
      model <- attr(plot_data, "jbrowse")
      if (is.null(model)) {
        stop("a jbrowse geom needs ggplot(<jbrowse result>) or data = <jbrowse result>",
          call. = FALSE)
      }
      finish(pick(on_axis(model)), model)
    })
  }
  if (is.data.frame(data)) {
    return(function(plot_data) {
      model <- attr(plot_data, "jbrowse")
      finish(frame_on_axis(data, model$regions), model)
    })
  }
  data
}

# A table's `color` column back as the display's colors where something on
# the way dropped the class, such as an ifelse() in a mutate(), so ggplot
# still paints it through an identity scale rather than its hue palette
display_colors <- function(d) {
  color <- d$color
  if (is.character(color) && !inherits(color, "jbrowse_color") &&
    all(is.na(color) | grepl("^#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?$", color) |
      color %in% grDevices::colors())) {
    d$color <- jbrowse_color(color)
  }
  d
}

with_defaults <- function(mapping, defaults) {
  for (name in names(mapping)) defaults[[name]] <- mapping[[name]]
  defaults
}

# A layer of a jbrowse table. One whose display colors no key lists stays
# out of the legend until the caller maps an aesthetic of their own; one of
# pileup `rows` draws each at the display's read height.
jbrowse_layer <- function(geom, data, pick, mapping, defaults, params,
                          keyed = TRUE, rows = FALSE) {
  inset <- if (rows) with_inset else function(d, model) d
  finish <- function(d, model) display_colors(inset(d, model))
  ggplot2::layer(geom = geom, stat = "identity", position = "identity",
    data = table_from(data, pick, finish),
    mapping = with_defaults(mapping, defaults),
    inherit.aes = FALSE,
    show.legend = legend_of(with_defaults(mapping, defaults), keyed, mapping),
    params = params)
}

#' @exportS3Method ggplot2::fortify
fortify.jbrowse_result <- function(model, data, ...) {
  structure(on_axis(model)[[attr(model, "primary")]], jbrowse = model)
}

# A key frame as an identity scale's breaks and labels, its line swatches
# drawn solid whatever width the data maps
key_scale <- function(scale, key, title) {
  color <- as.character(key$color)
  breaks <- unique(color)
  labels <- vapply(breaks, function(c) paste(key$label[color == c],
    collapse = " / "), "")
  scale(name = title, breaks = breaks, labels = unname(labels),
    limits = if (isTRUE(attr(key, "complete"))) function(x) union(x, breaks),
    guide = if (length(breaks)) {
      ggplot2::guide_legend(override.aes = c(list(linewidth = 1.5),
        attr(key, "override")))
    } else "none")
}

# The key a result's colors are read through, by name
key_of <- function(object, key) {
  keys <- object$keys
  key <- key %||% names(keys)[1]
  k <- keys[[key]]
  if (is.null(k)) {
    stop("no key ", key, "; this result keys ", paste(names(keys),
      collapse = ", "), call. = FALSE)
  }
  k
}

#' The display's key as a ggplot scale
#'
#' An identity scale over the display's colors whose legend lists the
#' display's key: what each fill or line color means, under the display's
#' title. A result can have several keys (an alignments track keys its
#' reads, its mismatch bases, its read connections and its splice
#' junctions), and `key` picks one.
#' @param object a ggjbrowse result
#' @param key which of its keys: for alignments `"reads"`, `"bases"`,
#'   `"arcs"` or `"junctions"`; `NULL` takes its first
#' @return a ggplot scale
#' @export
scale_fill_jbrowse <- function(object, key = NULL) {
  k <- key_of(object, key)
  key_scale(ggplot2::scale_fill_identity, k, attr(k, "title"))
}

#' @rdname scale_fill_jbrowse
#' @export
scale_colour_jbrowse <- function(object, key = NULL) {
  k <- key_of(object, key)
  key_scale(ggplot2::scale_colour_identity, k, attr(k, "title"))
}

#' @rdname scale_fill_jbrowse
#' @export
scale_color_jbrowse <- scale_colour_jbrowse

# A key frame from labels and colors, under a title
key_frame <- function(label, color, title) {
  structure(data.frame(label = label, color = jbrowse_color(color)),
    title = title)
}

# Whether a layer joins the legends: a layer drawn through a `complete` key
# draws a chip for every entry of it, present in the data or not, as the
# view lists its fixed key; a keyed layer or one the caller maps joins as
# ggplot decides, and one with neither stays out
legend_of <- function(mapping, keyed, own) {
  if (identical(keyed, "complete")) {
    keys <- names(mapping)[vapply(mapping, function(m)
      rlang::is_quosure(m) && identical(rlang::quo_get_expr(m), quote(color)),
      TRUE)]
    if (length(keys)) return(stats::setNames(rep(TRUE, length(keys)), keys))
  }
  if (isTRUE(keyed) || identical(keyed, "complete") || !is.null(own)) NA
  else FALSE
}

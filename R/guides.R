key_frame <- function(label, color, title) {
  structure(data.frame(label = label, color = tolower(color)), title = title)
}

key_shown <- function(key, shown) {
  kept <- key[key$color %in% shown, , drop = FALSE]
  for (a in setdiff(names(attributes(key)), names(attributes(kept)))) {
    attr(kept, a) <- attr(key, a)
  }
  kept
}

ramp_frames <- function(ramps) {
  lapply(ramps, function(r) {
    structure(data.frame(offset = unlist(r$offset),
        color = tolower(unlist(r$color)), opacity = unlist(r$opacity)),
      title = r$title, labels = c(r$minLabel, r$maxLabel))
  })
}

drawn_keys <- function(x) {
  keys <- x$keys
  if (inherits(x, "jbrowse_alignments") && !is.null(keys$junctions) &&
    !isTRUE(attr(keys$junctions, "shared"))) {
    keys$junctions <- key_shown(keys$junctions, shown_junctions(x)$color)
  }
  keys
}

#' The display's keys and color bars as ggplot guides
#'
#' The geoms paint the display's colors as they are, `I(color)`, so they
#' train no scale and leave every scale to your own mappings. The display's
#' keys come back as guides of their own instead: a chip per entry of each
#' key in `object$keys`, under its title, and a bar per color ramp in
#' `object$ramps`, labelled at its ends as the display labels it. They take
#' the plot's theme, so `theme(legend.text = element_text(size = 10))`
#' restyles them, and [jbrowse_stack()] draws a key several tracks share
#' once. Relabel or recolor an entry by editing the key, as in
#' `a$keys$reads$label[1] <- "Forward"`.
#' @param object a ggjbrowse result
#' @param keys the names of the keys to draw, every key by default; for
#'   alignments `"reads"`, `"bases"`, `"junctions"`, `"interbase"` or
#'   `"arcs"`
#' @param ramps draw the display's color bars too
#' @return ggplot guides
#' @export
#' @examples
#' bam <- system.file("extdata", "gapdh_direct_rna.bam",
#'   package = "ggjbrowse")
#' a <- jbrowse_alignments(bam, "chr12:6,534,201-6,538,700",
#'   display = list(color = list(field = "strand")))
#' library(ggplot2)
#' ggplot(a) + geom_jbrowse_reads() + guides_jbrowse(a, "reads") +
#'   scale_y_reverse()
guides_jbrowse <- function(object, keys = names(object$keys), ramps = TRUE) {
  drawn <- Filter(NROW, drawn_keys(object)[intersect(keys,
    names(object$keys))])
  names(drawn) <- sprintf("jbrowse_key_%s", names(drawn))
  bars <- if (ramps) object$ramps else list()
  names(bars) <- sprintf("jbrowse_ramp_%d", seq_along(bars))
  do.call(ggplot2::guides, lapply(c(drawn, bars), guide_jbrowse))
}

guide_jbrowse <- function(key) {
  title <- attr(key, "title")
  ggplot2::new_guide(key = key, title = title, grob = NULL, width = NULL,
    height = NULL, theme = NULL, position = NULL, order = 0,
    hash = rlang::hash(list(title, key)), available_aes = "any",
    super = GuideJbrowse)
}

GuideJbrowse <- ggplot2::ggproto("GuideJbrowse", ggplot2::GuideCustom,
  params = c(ggplot2::GuideCustom$params, list(key = NULL)),
  hashables = rlang::exprs(title, key),
  draw = function(self, theme, position = NULL, direction = NULL,
                  params = self$params) {
    drawn <- if (is.null(params$key$offset)) chips(params$key, theme)
      else ramp_bar(params$key, theme)
    params[names(drawn)] <- drawn
    ggplot2::ggproto_parent(ggplot2::GuideCustom, self)$draw(theme,
      position, direction, params)
  })

legend_text <- function(theme) {
  e <- ggplot2::calc_element("legend.text", theme)
  grid::gpar(fontsize = e$size, col = e$colour, fontfamily = e$family,
    fontface = e$face)
}

in_pt <- function(u) grid::convertUnit(u, "pt", valueOnly = TRUE)

text_pt <- function(labels, gp) {
  max(0, vapply(labels, function(l) grid::convertWidth(grid::grobWidth(
    grid::textGrob(l, gp = gp)), "pt", valueOnly = TRUE), 0))
}

chips <- function(key, theme) {
  gp <- legend_text(theme)
  row <- max(in_pt(ggplot2::calc_element("legend.key.height", theme)),
    1.2 * gp$fontsize)
  chip <- 0.8 * row
  gap <- 5.5
  n <- nrow(key)
  y <- grid::unit(1, "npc") - grid::unit((seq_len(n) - 0.5) * row, "pt")
  fill <- scales::alpha(key$color, attr(key, "alpha") %||% NA)
  list(grob = grid::gTree(children = grid::gList(
      grid::rectGrob(x = grid::unit(row / 2, "pt"), y = y,
        width = grid::unit(chip, "pt"), height = grid::unit(chip, "pt"),
        gp = grid::gpar(fill = fill, col = NA)),
      grid::textGrob(key$label, x = grid::unit(row + gap, "pt"), y = y,
        hjust = 0, gp = gp))),
    width = grid::unit(row + gap + text_pt(key$label, gp), "pt"),
    height = grid::unit(n * row, "pt"))
}

ramp_bar <- function(ramp, theme) {
  gp <- legend_text(theme)
  at <- seq(0, 1, length.out = 128)
  fill <- scales::alpha(
    scales::gradient_n_pal(ramp$color, ramp$offset)(at),
    stats::approx(ramp$offset, ramp$opacity, at)$y)
  bar <- 0.14 * 72.27
  labels <- attr(ramp, "labels")
  list(grob = grid::gTree(children = grid::gList(
      grid::rasterGrob(matrix(fill, nrow = 1), y = grid::unit(1, "npc"),
        height = grid::unit(bar, "pt"), width = grid::unit(1, "npc"),
        vjust = 1, interpolate = TRUE),
      grid::textGrob(labels[1], x = 0, y = 0, hjust = 0, vjust = 0, gp = gp),
      grid::textGrob(labels[2], x = 1, y = 0, hjust = 1, vjust = 0,
        gp = gp))),
    width = grid::unit(max(1.2 * 72.27, text_pt(labels, gp) * 2 + 11), "pt"),
    height = grid::unit(bar + 1.3 * gp$fontsize, "pt"))
}

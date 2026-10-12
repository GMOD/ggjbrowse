# The engine's color bars, each a frame of its stops with its title and the
# end labels the display's bar prints
ramp_frames <- function(ramps) {
  lapply(ramps, function(r) {
    structure(data.frame(offset = unlist(r$offset),
        color = jbrowse_color(unlist(r$color)),
        opacity = unlist(r$opacity)),
      title = r$title, labels = c(r$minLabel, r$maxLabel))
  })
}

#' The display's color bars as ggplot guides
#'
#' A display that paints through a color ramp keys it as a bar labelled at
#' its ends: a density track's scores, reads colored by a linear tag or
#' MAPQ, base quality. The marks carry the colors the ramp gave them, so
#' each bar is a guide of its own beside the identity scale's key. Map a
#' column through a continuous scale of your own to key it the ggplot way.
#' @param object a ggjbrowse result
#' @return ggplot guides, a bar per ramp the display keys
#' @export
#' @examples
#' bw <- system.file("extdata", "gapdh_phylop100way.bw",
#'   package = "ggjbrowse")
#' w <- jbrowse_wiggle(bw, "chr12", 6534200, 6538700,
#'   display = list(mark = "span"))
#' library(ggplot2)
#' ggplot(w) + geom_jbrowse_wiggle() + guides_jbrowse(w)
guides_jbrowse <- function(object) {
  guides <- lapply(object$ramps, ramp_guide)
  names(guides) <- sprintf("jbrowse_ramp_%d", seq_along(guides))
  do.call(ggplot2::guides, guides)
}

# A ramp's bar 1.2 inches wide, its colors read between the stops, and its
# end labels under it as the display prints them
ramp_guide <- function(ramp) {
  at <- seq(0, 1, length.out = 128)
  fill <- scales::alpha(
    scales::gradient_n_pal(as.character(ramp$color), ramp$offset)(at),
    stats::approx(ramp$offset, ramp$opacity, at)$y)
  bar_in <- 0.14
  text <- grid::gpar(fontsize = 8, col = "grey30")
  labels <- attr(ramp, "labels")
  grob <- grid::gTree(children = grid::gList(
    grid::rasterGrob(matrix(fill, nrow = 1), y = grid::unit(1, "npc"),
      height = grid::unit(bar_in, "in"), width = grid::unit(1, "npc"),
      vjust = 1, interpolate = TRUE),
    grid::textGrob(labels[1], x = 0, y = 0, hjust = 0, vjust = 0, gp = text),
    grid::textGrob(labels[2], x = 1, y = 0, hjust = 1, vjust = 0, gp = text)))
  ggplot2::guide_custom(grob, title = attr(ramp, "title"),
    width = grid::unit(1.2, "in"), height = grid::unit(bar_in + 0.16, "in"))
}

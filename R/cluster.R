# A variant matrix's rows in the order it draws them: `rows`' own order, or a
# tree's leaf order where `tree` holds one (anything stats::as.hclust()
# takes), its labels the rows' names; the cells' `row` factor follows
in_row_order <- function(x) {
  if (!is.null(x$tree)) {
    tree <- stats::as.hclust(x$tree)
    leaves <- tree$labels[tree$order]
    if (length(leaves) != nrow(x$rows) || !setequal(leaves, x$rows$name)) {
      stop("the tree's leaves are not the matrix's rows: label it by ",
        "x$rows$name, one leaf per row (", nrow(x$rows), " rows, such as ",
        paste(utils::head(x$rows$name, 2), collapse = ", "), ")",
        call. = FALSE)
    }
    x$rows <- x$rows[match(leaves, x$rows$name), , drop = FALSE]
    rownames(x$rows) <- NULL
  }
  x$cells$row <- factor(as.character(x$cells$row),
    levels = as.character(x$rows$name))
  x
}

# A tree's branches as segments, each leaf at its label's row in `rows`:
# `y` and `yend` in rows, `h` and `hend` the merge heights over the tallest,
# 0 at the leaves
tree_segments <- function(tree, rows) {
  n <- length(tree$order)
  if (n < 2) {
    return(data.frame(y = numeric(0), yend = numeric(0), h = numeric(0),
      hend = numeric(0)))
  }
  leaf <- match(tree$labels, rows)
  height <- tree$height / max(tree$height, .Machine$double.eps)
  y <- numeric(n - 1)
  at <- function(k) if (k < 0) c(leaf[-k], 0) else c(y[k], height[k])
  do.call(rbind, lapply(seq_len(n - 1), function(k) {
    a <- at(tree$merge[k, 1])
    b <- at(tree$merge[k, 2])
    y[k] <<- (a[1] + b[1]) / 2
    data.frame(y = c(a[1], b[1], a[1]), yend = c(a[1], b[1], b[1]),
      h = c(a[2], b[2], height[k]), hend = height[k])
  }))
}

# A tree drawn left of the panel, its leaves `gap` points from the panel's
# edge and its root `span` points further out, each branch at its row
GeomTree <- ggplot2::ggproto("GeomTree", ggplot2::Geom,
  required_aes = c("y", "yend", "h", "hend"),
  default_aes = ggplot2::aes(colour = "grey30", linewidth = 0.3),
  extra_params = c("na.rm", "gap", "span"),
  draw_panel = function(data, panel_params, coord, gap, span) {
    d <- coord$transform(data, panel_params)
    out <- function(h) grid::unit(0, "npc") - grid::unit(gap + h * span, "pt")
    grid::segmentsGrob(x0 = out(d$h), x1 = out(d$hend),
      y0 = grid::unit(d$y, "npc"), y1 = grid::unit(d$yend, "npc"),
      gp = grid::gpar(col = d$colour[1], lwd = d$linewidth[1] * ggplot2::.pt))
  })

#' A variant matrix's row tree as a ggplot layer
#'
#' The tree in a [jbrowse_variants()] result's `tree`, drawn left of the
#' panel beside the rows, each leaf at its row, as the variant display draws
#' its clustering tree in its sidebar. `tree` takes an [stats::hclust()]
#' result, or anything [stats::as.hclust()] takes, labelled by the rows'
#' `name`; setting it also puts the rows in its leaf order. A plot of one's
#' own draws outside the panel only under `coord_cartesian(clip = "off")`.
#' @param mapping aesthetics that override the defaults, `colour` and
#'   `linewidth`
#' @param data a [jbrowse_variants()] result, or `NULL` for the plot's
#' @param gap the points between the panel's left edge and the leaves
#' @param span the points from the leaves to the root
#' @param ... passed to [ggplot2::layer()] params
#' @return a ggplot layer
#' @export
#' @examples
#' vcf <- system.file("extdata", "1000g_dcst2.vcf.gz", package = "ggjbrowse")
#' v <- jbrowse_variants(vcf, "1:155,000,001-155,050,000")
#' v$tree <- hclust(dist(xtabs(dosage ~ row + variant, v$cells)), "average")
#' library(ggplot2)
#' ggplot(v) + geom_jbrowse_genotypes(layout = "columns") +
#'   geom_jbrowse_tree() + coord_cartesian(clip = "off")
geom_jbrowse_tree <- function(mapping = NULL, data = NULL, gap = 2,
                              span = 40, ...) {
  data <- replacing(data, NULL, "geom_jbrowse_tree")
  segments <- function(x) {
    tree <- if (!is.null(x$tree)) stats::as.hclust(x$tree)
    tree_segments(tree, as.character(x$rows$name))
  }
  jbrowse_layer(GeomTree, data, segments, mapping,
    ggplot2::aes(y = y, yend = yend, h = h, hend = hend),
    list(gap = gap, span = span, ...), keyed = FALSE)
}

# The points a tree spans beside a matrix's rows
TREE_PT <- 40

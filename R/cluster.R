#' A variant matrix's genotypes as the display's clustering reads them
#'
#' The matrix the variant display's "Cluster rows by genotype" clusters: a
#' row per row of `x$rows`, in its order, named by its `name`. With a row per
#' sample, each record has a column per ALT allele, scored 0 to 2 as twice
#' that allele's share of the sample's called alleles, so `1/1` and `2/2`
#' stand apart; with a row per haplotype (`unit = "haplotype"`), a column per
#' record, 1 where the haplotype carries an ALT. A no-call is `NA`, as is a
#' haplotype of an unphased call, and a reference call the display skips
#' drawing is 0. `hclust(dist(as.matrix(v)), "average")` is the R script the
#' dialog writes: `dist()` leaves a pair's no-calls out of their distance,
#' where the in-app run fills each with its record's mean.
#' @param x a [jbrowse_variants()] result
#' @param ... unused
#' @return a numeric matrix, its columns named by the record's row number in
#'   `x$variants`, with `.1`, `.2`, ... per ALT of a multi-allelic record
#' @export
#' @examples
#' vcf <- system.file("extdata", "1000g_dcst2.vcf.gz", package = "ggjbrowse")
#' v <- jbrowse_variants(vcf, "1:155,000,001-155,050,000")
#' m <- as.matrix(v)
#' v$tree <- hclust(dist(m), "average")
as.matrix.jbrowse_variants <- function(x, ...) {
  cells <- x$cells[as.character(x$cells$row) %in% x$rows$name, , drop = FALSE]
  i <- match(as.character(cells$row), x$rows$name)
  rows <- as.character(x$rows$name)
  if (identical(x$unit, "haplotype")) {
    m <- matrix(0, length(rows), nrow(x$variants),
      dimnames = list(rows, seq_len(nrow(x$variants))))
    m[cbind(i, cells$variant)] <- cells$dosage
    return(m)
  }
  alts <- lengths(strsplit(x$variants$alt, ",", fixed = TRUE))
  first <- cumsum(c(0, alts))[seq_along(alts)]
  columns <- unlist(lapply(seq_along(alts), function(v)
    if (alts[v] == 1) as.character(v) else paste0(v, ".", seq_len(alts[v]))))
  m <- matrix(0, length(rows), length(columns),
    dimnames = list(rows, columns))
  alleles <- strsplit(cells$genotype, "[/|]")
  called <- lapply(alleles, function(a) suppressWarnings(as.integer(a[a != "."])))
  n <- lengths(called)
  cell <- rep(seq_along(called), n)
  allele <- unlist(called)
  alt <- !is.na(allele) & allele >= 1 & allele <= alts[cells$variant[cell]]
  slots <- max(alts, 0) + 1
  key <- (cell[alt] - 1) * slots + allele[alt]
  hit <- unique(key)
  count <- tabulate(match(key, hit), length(hit))
  hit_cell <- (hit - 1) %/% slots + 1
  hit_allele <- (hit - 1) %% slots + 1
  m[cbind(i[hit_cell], first[cells$variant[hit_cell]] + hit_allele)] <-
    2 * count / n[hit_cell]
  none <- which(n == 0)
  per <- alts[cells$variant[none]]
  m[cbind(rep(i[none], per), rep(first[cells$variant[none]], per) +
    sequence(per))] <- NA
  m
}

# `rows` in a tree's leaf order where `tree` holds one (anything
# stats::as.hclust() takes), its labels the rows' names; else as they are.
# `where` names the rows' table for an error.
by_tree <- function(rows, tree, where) {
  if (is.null(tree)) return(rows)
  tree <- stats::as.hclust(tree)
  leaves <- tree$labels[tree$order]
  if (length(leaves) != nrow(rows) || !setequal(leaves, rows$name)) {
    stop("the tree's leaves are not the rows: label it by ", where,
      "$name, one leaf per row (", nrow(rows), " rows, such as ",
      paste(utils::head(rows$name, 2), collapse = ", "), ")", call. = FALSE)
  }
  rows <- rows[match(leaves, rows$name), , drop = FALSE]
  rownames(rows) <- NULL
  rows
}

# A variant matrix's rows in the order it draws them: `rows`' own order, or
# its tree's leaf order; the cells' `row` factor follows
in_row_order <- function(x) {
  x$rows <- by_tree(x$rows, x$tree, "x$rows")
  x$cells$row <- factor(as.character(x$cells$row),
    levels = as.character(x$rows$name))
  x
}

# A quantitative track's sources in the order it draws them: `sources`' own
# order, or its tree's leaf order, each source's marks on its row, and a
# source left out of `sources` left out of the marks
in_source_order <- function(x) {
  x$sources <- by_tree(x$sources, x$tree, "x$sources")
  for (name in c("marks", "bins")) {
    d <- x[[name]]
    x[[name]] <- d[d$source %in% x$sources$name, , drop = FALSE]
  }
  if (isTRUE(x$row_layout)) {
    x$marks$row <- match(x$marks$source, x$sources$name) - 1
  }
  x
}

# The rows a tree's leaves sit at: their names, and each row's y as its
# plot draws it, from its position in that order
tree_rows <- function(x) {
  if (inherits(x, "jbrowse_wiggle")) {
    n <- nrow(x$sources)
    return(list(name = as.character(x$sources$name),
      at = function(i) n - i + 0.5))
  }
  list(name = as.character(x$rows$name), at = identity)
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

#' A row tree as a ggplot layer
#'
#' The tree in a [jbrowse_variants()] or [jbrowse_wiggle()] result's `tree`,
#' drawn left of the panel beside the rows, each leaf at its row, as the
#' display draws its clustering tree in its sidebar. `tree` takes an
#' [stats::hclust()] result, or anything [stats::as.hclust()] takes, labelled
#' by the rows' `name` (a quantitative track's `sources`' `name`); setting it
#' also puts the rows in its leaf order. A quantitative track's rows are its
#' sources under a row layout. A plot of one's own draws outside the panel
#' only under `coord_cartesian(clip = "off")`.
#' @param mapping aesthetics that override the defaults, `colour` and
#'   `linewidth`
#' @param data a [jbrowse_variants()] or [jbrowse_wiggle()] result, or
#'   `NULL` for the plot's
#' @param gap the points between the panel's left edge and the leaves
#' @param span the points from the leaves to the root
#' @param ... passed to [ggplot2::layer()] params
#' @return a ggplot layer
#' @export
#' @examples
#' vcf <- system.file("extdata", "1000g_dcst2.vcf.gz", package = "ggjbrowse")
#' v <- jbrowse_variants(vcf, "1:155,000,001-155,050,000")
#' v$tree <- hclust(dist(as.matrix(v)), "average")
#' library(ggplot2)
#' ggplot(v) + geom_jbrowse_genotypes(layout = "columns") +
#'   geom_jbrowse_tree() + coord_cartesian(clip = "off")
geom_jbrowse_tree <- function(mapping = NULL, data = NULL, gap = 2,
                              span = 40, ...) {
  data <- replacing(data, NULL, "geom_jbrowse_tree")
  segments <- function(x) {
    tree <- if (!is.null(x$tree)) stats::as.hclust(x$tree)
    rows <- tree_rows(x)
    s <- tree_segments(tree, rows$name)
    s$y <- rows$at(s$y)
    s$yend <- rows$at(s$yend)
    s
  }
  jbrowse_layer(GeomTree, data, segments, mapping,
    ggplot2::aes(y = y, yend = yend, h = h, hend = hend),
    list(gap = gap, span = span, ...), keyed = FALSE)
}

# The points a tree spans beside a matrix's rows
TREE_PT <- 40

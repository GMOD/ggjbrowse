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
#' @return a numeric matrix, its columns named by the record's `variant`
#'   key, with `.1`, `.2`, ... per ALT of a multi-allelic record
#' @export
#' @examples
#' vcf <- system.file("extdata", "1000g_dcst2.vcf.gz", package = "ggjbrowse")
#' v <- jbrowse_variants(vcf, "1:155,000,001-155,050,000")
#' m <- as.matrix(v)
#' v$tree <- hclust(dist(m), "average")
as.matrix.jbrowse_variants <- function(x, ...) {
  variants <- x$variants[!duplicated(x$variants$variant), , drop = FALSE]
  rows <- as.character(x$rows$name)
  cells <- x$cells[as.character(x$cells$row) %in% rows &
    x$cells$variant %in% variants$variant, , drop = FALSE]
  cells <- cells[!duplicated(cells[c("variant", "row")]), , drop = FALSE]
  i <- match(as.character(cells$row), rows)
  v <- match(cells$variant, variants$variant)
  if (identical(x$display$unit, "haplotype")) {
    m <- matrix(0, length(rows), nrow(variants),
      dimnames = list(rows, variants$variant))
    m[cbind(i, v)] <- cells$dosage
    return(m)
  }
  alts <- lengths(strsplit(variants$alt, ",", fixed = TRUE))
  first <- cumsum(c(0, alts))[seq_along(alts)]
  columns <- unlist(lapply(seq_along(alts), function(k)
    if (alts[k] == 1) as.character(variants$variant[k])
    else paste0(variants$variant[k], ".", seq_len(alts[k]))))
  m <- matrix(0, length(rows), length(columns),
    dimnames = list(rows, columns))
  alleles <- strsplit(cells$genotype, "[/|]")
  called <- lapply(alleles, function(a) suppressWarnings(as.integer(a[a != "."])))
  n <- lengths(called)
  cell <- rep(seq_along(called), n)
  allele <- unlist(called)
  alt <- !is.na(allele) & allele >= 1 & allele <= alts[v[cell]]
  slots <- max(alts, 0) + 1
  key <- (cell[alt] - 1) * slots + allele[alt]
  hit <- unique(key)
  count <- tabulate(match(key, hit), length(hit))
  hit_cell <- (hit - 1) %/% slots + 1
  hit_allele <- (hit - 1) %% slots + 1
  m[cbind(i[hit_cell], first[v[hit_cell]] + hit_allele)] <-
    2 * count / n[hit_cell]
  none <- which(n == 0)
  per <- alts[v[none]]
  m[cbind(rep(i[none], per), rep(first[v[none]], per) +
    sequence(per))] <- NA
  m
}

# The trees `tree` holds, each as an hclust: none, the one, or a forest, a
# plain list of them
trees_of <- function(tree) {
  if (is.null(tree)) return(list())
  if (is.list(tree) && !is.object(tree)) lapply(tree, stats::as.hclust)
  else list(stats::as.hclust(tree))
}

# `rows` in a tree's leaf order where `tree` holds one (anything
# stats::as.hclust() takes), its labels the rows' names; else as they are.
# A forest orders the rows its trees hold, tree by tree, and the rest follow.
# `where` names the rows' table for an error.
by_tree <- function(rows, tree, where) {
  trees <- trees_of(tree)
  if (!length(trees)) return(rows)
  leaves <- unlist(lapply(trees, function(t) t$labels[t$order]))
  forest <- is.list(tree) && !is.object(tree)
  if (anyDuplicated(leaves) || !all(leaves %in% rows$name) ||
    (!forest && length(leaves) != nrow(rows))) {
    stop("the tree's leaves are not the rows: label it by ", where,
      "$name, one leaf per row (", nrow(rows), " rows, such as ",
      paste(utils::head(rows$name, 2), collapse = ", "), ")", call. = FALSE)
  }
  at <- match(leaves, rows$name)
  rows <- rows[c(at, setdiff(seq_len(nrow(rows)), at)), , drop = FALSE]
  rownames(rows) <- NULL
  rows
}

# A variant matrix's rows in the order it draws them: `rows`' own order, or
# its tree's leaf order, stacked in `bands`' order where the rows are banded;
# the cells' `row` factor follows
in_row_order <- function(x) {
  x$rows <- by_tree(x$rows, x$tree, "x$rows")
  if (any(x$rows$band != "")) {
    x$rows <- x$rows[order(match(x$rows$band, x$bands$key)), , drop = FALSE]
    rownames(x$rows) <- NULL
  }
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

# Each band's run of rows, top to bottom: its `key`, its `label` from
# `bands`, and its first and last row; NULL while no row is banded
row_bands <- function(rows, bands) {
  band <- rows$band
  if (!length(band) || all(band == "")) return(NULL)
  runs <- rle(band)
  last <- cumsum(runs$lengths)
  label <- bands$label[match(runs$values, bands$key)]
  data.frame(key = runs$values,
    label = ifelse(is.na(label), runs$values, label),
    first = last - runs$lengths + 1, last = last)
}

# The merges of `tree` a banded matrix draws: those inside a clade whose
# leaves are exactly one band's rows, as the display draws a clade per band
band_merges <- function(tree, rows, bands) {
  n <- length(tree$order)
  members <- vector("list", n - 1)
  of <- function(child) if (child < 0) tree$labels[-child] else members[[child]]
  for (k in seq_len(n - 1)) {
    members[[k]] <- c(of(tree$merge[k, 1]), of(tree$merge[k, 2]))
  }
  clades <- lapply(seq_len(nrow(bands)), function(b) {
    names <- rows[bands$first[b]:bands$last[b]]
    Find(function(k) setequal(members[[k]], names), seq_len(n - 1))
  })
  clades <- unlist(clades)
  vapply(members, function(m)
    any(vapply(clades, function(k) all(m %in% members[[k]]), NA)), NA)
}

# A tree's branches as segments, each leaf at its label's row in `rows`:
# `y` and `yend` in rows, `h` and `hend` the merge heights over the tallest,
# 0 at the leaves
tree_segments <- function(tree, rows, top = max(tree$height)) {
  n <- length(tree$order)
  if (n < 2) {
    return(data.frame(y = numeric(0), yend = numeric(0), h = numeric(0),
      hend = numeric(0)))
  }
  leaf <- match(tree$labels, rows)
  height <- tree$height / max(top, .Machine$double.eps)
  merge <- tree$merge
  y <- numeric(n - 1)
  y_of <- function(child) {
    ifelse(child < 0, leaf[pmax(-child, 1)], y[pmax(child, 1)])
  }
  h_of <- function(child) ifelse(child < 0, 0, height[pmax(child, 1)])
  for (k in seq_len(n - 1)) y[k] <- mean(y_of(merge[k, ]))
  a <- y_of(merge[, 1])
  b <- y_of(merge[, 2])
  data.frame(y = c(rbind(a, b, a)), yend = c(rbind(a, b, b)),
    h = c(rbind(h_of(merge[, 1]), h_of(merge[, 2]), height)),
    hend = rep(height, each = 3))
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
#' also puts the rows in its leaf order. A list of trees is a forest, which
#' orders the rows tree by tree. Under a variant display's `facet`, only a
#' clade whose leaves are exactly one band's rows draws, as in the display.
#' A quantitative track's rows are its sources under a row layout. A plot of
#' one's own draws outside the panel only under
#' `coord_cartesian(clip = "off")`.
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
    trees <- trees_of(x$tree)
    rows <- tree_rows(x)
    top <- max(0, unlist(lapply(trees, `[[`, "height")))
    bands <- if (inherits(x, "jbrowse_variants")) row_bands(x$rows, x$bands)
    s <- do.call(rbind, c(list(tree_segments(NULL, rows$name)),
      lapply(trees, function(tree) {
        s <- tree_segments(tree, rows$name, top)
        if (is.null(bands) || !nrow(s)) return(s)
        s[rep(band_merges(tree, rows$name, bands), each = 3), , drop = FALSE]
      })))
    s$y <- rows$at(s$y)
    s$yend <- rows$at(s$yend)
    s
  }
  jbrowse_layer(GeomTree, data, segments, mapping,
    ggplot2::aes(y = y, yend = yend, h = h, hend = hend),
    list(gap = gap, span = span, ...))
}

# The points a tree spans beside a matrix's rows
TREE_PT <- 40

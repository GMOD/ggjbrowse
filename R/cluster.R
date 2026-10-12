#' Order a variant matrix's rows by their genotypes
#'
#' Clusters the rows of a [jbrowse_variants()] result as the variant display's
#' "Cluster rows by genotype" does: each row's dosage of each ALT allele at
#' each record (0, 1 or 2 on a diploid scale, a no-call taking the record's
#' mean), Euclidean distances, and average linkage, the R script JBrowse's
#' clustering dialog writes. [autoplot()] draws the tree beside the rows.
#' @param x a [jbrowse_variants()] result
#' @param method the linkage, as [stats::hclust()] takes it; `"average"` is
#'   the display's
#' @return `x` with its rows in the tree's leaf order and `tree`, the
#'   [stats::hclust()] result, its labels the row names
#' @export
#' @examples
#' vcf <- system.file("extdata", "dog10k_cyp1a2_snvs.vcf.gz",
#'   package = "ggjbrowse")
#' v <- jbrowse_cluster(jbrowse_variants(vcf, "chr30:38,258,001-38,265,000"))
#' head(v$rows$name)
jbrowse_cluster <- function(x, method = "average") {
  if (!inherits(x, "jbrowse_variants")) {
    stop("jbrowse_cluster orders a jbrowse_variants() result's rows",
      call. = FALSE)
  }
  if (nrow(x$rows) < 2) return(x)
  x$tree <- stats::hclust(stats::dist(genotype_matrix(x)), method = method)
  in_row_order(x)
}

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

# Each row's dosage of each ALT at each record, a column per ALT: 2 times the
# share of the row's called alleles that are that ALT, a haplotype row
# counting its own allele. A row with no cell at a record is homozygous
# reference there, and a no-call takes its column's mean.
genotype_matrix <- function(x) {
  alts <- lengths(strsplit(x$variants$alt, ","))
  first <- cumsum(c(0, alts))[seq_along(alts)]
  m <- matrix(0, nrow(x$rows), sum(alts), dimnames = list(x$rows$name, NULL))
  cells <- x$cells
  row <- match(as.character(cells$row), x$rows$name)
  alleles <- strsplit(cells$genotype, "[/|]")
  if (any(x$rows$name != x$rows$sample)) {
    alleles <- Map(function(a, h) a[h + 1], alleles, x$rows$haplotype[row])
  }
  k <- alts[cells$variant]
  each <- rep(seq_len(nrow(cells)), k)
  alt <- sequence(k)
  m[cbind(row[each], first[cells$variant[each]] + alt)] <- mapply(
    function(a, j) {
      called <- a[!is.na(a) & a != "."]
      if (length(called)) 2 * sum(called == j) / length(called) else NA
    }, alleles[each], as.character(alt))
  missing <- which(is.na(m), arr.ind = TRUE)
  means <- colMeans(m, na.rm = TRUE)
  m[missing] <- ifelse(is.nan(means[missing[, 2]]), 0, means[missing[, 2]])
  m
}

# A tree's branches as segments: `y` and `yend` in rows, `h` and `hend` the
# merge heights over the tallest, 0 at the leaves
tree_segments <- function(tree) {
  n <- length(tree$order)
  leaf <- integer(n)
  leaf[tree$order] <- seq_len(n)
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

# A tree drawn left of the panel, its leaves `inner_pt` from the panel's
# edge and its root `width_pt` further out, each branch at its row
GeomTree <- ggplot2::ggproto("GeomTree", ggplot2::Geom,
  required_aes = c("y", "yend", "h", "hend"),
  default_aes = ggplot2::aes(colour = "grey30", linewidth = 0.3),
  extra_params = c("na.rm", "inner_pt", "width_pt"),
  draw_panel = function(data, panel_params, coord, inner_pt, width_pt) {
    d <- coord$transform(data, panel_params)
    out <- function(h) grid::unit(0, "npc") -
      grid::unit(inner_pt + h * width_pt, "pt")
    grid::segmentsGrob(x0 = out(d$h), x1 = out(d$hend),
      y0 = grid::unit(d$y, "npc"), y1 = grid::unit(d$yend, "npc"),
      gp = grid::gpar(col = d$colour[1], lwd = d$linewidth[1] * ggplot2::.pt))
  })

tree_layer <- function(tree, inner_pt, width_pt) {
  ggplot2::layer(geom = GeomTree, stat = "identity", position = "identity",
    data = tree_segments(stats::as.hclust(tree)),
    mapping = ggplot2::aes(y = y, yend = yend, h = h, hend = hend),
    inherit.aes = FALSE, show.legend = FALSE,
    params = list(inner_pt = inner_pt, width_pt = width_pt))
}

# The width a tree takes beside a matrix's rows
TREE_PT <- 40

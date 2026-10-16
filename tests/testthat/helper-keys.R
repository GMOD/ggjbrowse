# The row of `parent` each row of `child` names by `key`, in the same region
# where both list regions
of <- function(child, parent, key, child_key = key) {
  same <- !is.null(child$region) && !is.null(parent$region)
  id <- function(d, k) if (same) paste(d[[k]], d$region) else d[[k]]
  parent[match(id(child, child_key), id(parent, key)), , drop = FALSE]
}

# The keys and color bars a plot's jbrowse guides draw, by title
drawn_guides <- function(p) {
  keys <- lapply(ggplot2::ggplot_build(p)$plot$guides$params, `[[`, "key")
  stats::setNames(keys, vapply(keys, attr, "", "title"))
}

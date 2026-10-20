# Quantitative tracks

`jbrowse_wiggle()` reads a BigWig through JBrowse's quantitative display: the
bins the file serves at that zoom, the marks the display draws from them, and
the score axis it autoscales to. `display` takes its settings (`mark`,
`y`, `extent`, `color`, `scales`, ...), and `width` picks the zoom
level that answers ([bars, a line and density](gallery.md#conservation)).
`w$y_scale` is the score axis, which `scale_y_jbrowse(w)` draws.

- A `MultiWiggleAdapter` config draws a row per source in `w$sources`' order.
  Reorder or filter that table to rearrange them; filtering drops the left-out
  sources' bins and marks ([twelve cell lines](gallery.md#dnase_rows)).
- `as.matrix(w)` is the matrix the display's "Cluster rows" clusters: a row
  per source, a column per pixel.
- `w$tree <- hclust(dist(as.matrix(w)), "average")` puts the rows in the
  tree's order, and `autoplot()` draws the tree beside them
  ([clustered](gallery.md#dnase_clustered)).

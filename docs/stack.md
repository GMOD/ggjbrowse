# A view of several tracks

The code below reads `e()` from the [quick start](../README.md#quick-start).

`jbrowse_stack()` stacks tracks as a JBrowse view does: each track's plot, or
an alignments track's bands, top to bottom in the order given, on one x axis,
with the keys collected.

- `jbrowse_stack()` lays several loci end to end, as a view lays out several
  regions, and labels each tick in its own region's positions
  ([two regions](gallery.md#two_regions)).
- A transcript's exons, each with a few bases of flank, lay out side by side
  with the introns collapsed, as JBrowse's **Collapse introns** does, and each
  splice junction arcs from one exon straight to the next
  ([GAPDH's exons](gallery.md#collapsed_introns)).
- Every track draws at one width: `jbrowse_stack(..., width = 1200)`, or the
  widest `width` any track was read at. A feature or quantitative track drawn
  narrower than it was laid out for warns.
- Each band is an ordinary ggplot, so patchwork's `&` restyles them all, as
  in `jbrowse_stack(genes, signal, reads) & theme_classic()`. A ggplot of
  your own stacks as a band on the same axis.
- A track's name titles it, as the view labels a track:
  `jbrowse_stack(Genes = genes, reads)`, or `setNames(reads, contexts)` for
  a list of them.
- An alignments track read with `display = list(coverage = list(show = FALSE))`
  stacks no coverage band.
- A variant matrix gives each row the display's `rowHeight` pixels, or
  without one fits its rows into the display's `height`, 200 pixels by
  default, as the display does.
- A key several tracks share appears once, listing every entry any of them
  holds. A tumor with hard clips over a normal without them gets one
  "Insertions and clips" key, "Hard clip" included.
- Files that spell a contig differently stack on one axis when each reader
  takes the same `assembly`: the readers then name the contig in every table
  as the region does.

Each band takes the pixels the display gives it. Reads draw at the display's
`pileup$readHeight` (7, or 3 and 1 for its compact and super-compact settings),
and `jbrowse_height()` gives the figure height that keeps those proportions,
so a compact pileup saves as a shorter figure. `jbrowse_width()` gives the
width at the same scale, so a track read at a wide view's `width` draws its
labels at the display's size:

```r
loci <- c("chr12:6,534,201-6,535,700", "chr12:6,537,201-6,538,700")
genes <- jbrowse_features(e("gapdh_gencode.gff3.gz"), loci)
reads <- jbrowse_alignments(e("gapdh_direct_rna.bam"), loci,
  display = list(pileup = list(readHeight = 3)))
fig <- jbrowse_stack(genes, reads)
ggsave("stack.png", fig, width = jbrowse_width(fig),
  height = jbrowse_height(fig))
```

## From a JBrowse view

JBrowse's linear genome view is gaining an **Export R script** menu item, not
yet in a release. It writes one ggjbrowse reader per open track, then
`jbrowse_stack()` and `ggsave()` at `jbrowse_width(fig)` by
`jbrowse_height(fig)`.

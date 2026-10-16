# The grammar

The examples read `a` and `e()` from the
[quick start](../README.md#quick-start).

`autoplot()` only composes pieces you can use yourself:

- **Data:** `ggplot(a)` holds a result, and each table is a data frame
  (`a$reads`, `a$mismatches`, `a$coverage`, ...).
- **Geoms:** one per table: `geom_jbrowse_reads()`,
  `geom_jbrowse_mismatches()`, `geom_jbrowse_coverage()`,
  `geom_jbrowse_arcs()`, `geom_jbrowse_glyphs()`, `geom_jbrowse_wiggle()`,
  `geom_jbrowse_genotypes()` and the rest. Each draws as the display does
  and takes any column of its table in `aes()`.
- **Colors:** each table's `color` column is the display's color, which
  each geom paints as it is, `I(color)`, so it trains no scale. Map another
  column and ggplot's own scales take it, beside the other geoms' display
  colors.
- **Keys:** `guides_jbrowse(a)` adds the display's keys and color bars as
  guides that take the plot's theme. Edit `a$keys` to relabel an entry.
- **Scales:** `scale_y_jbrowse(w)` adds a quantitative display's score axis,
  and `scale_x_jbrowse(a)` the axis of several regions laid end to end,
  which `ggplot(a)` adds unasked, moving the result's tables in your own
  layers onto it. A frame of your own with `chrom` and `start` and `end`,
  or `position`, moves to the region it falls in.
- **Decisions:** a geom draws what its table says, such as
  `a$junctions$shown`, so editing the column draws otherwise.
- **Filters:** assign a filtered table into the result, as in
  `a$reads <- subset(a$reads, strand == 1)`, and its marks follow, since
  each names its read by a key. A geom takes a filtered mark table as
  `data`, but not a filtered `reads`.
- **Facets:** the `section` column is the display's facet, so
  `facet_grid(section ~ .)` splits a plot the way `facet =` splits the track,
  on an alignments or a feature track.

```r
library(ggplot2)

ggplot(a) +
  geom_jbrowse_reads(aes(fill = mapq)) +
  geom_jbrowse_mismatches() +
  scale_fill_viridis_c() +
  scale_y_reverse()
```

The tables work without the geoms as well:

```r
cram <- jbrowse_alignments(e("gapdh_direct_rna.cram"),
  "chr12:6,534,201-6,535,200", assembly = e("gapdh_hg38.fa.gz"))
ggplot(cram$mismatches, aes(position, frequency, colour = I(color))) +
  geom_point() +
  guides_jbrowse(cram, "bases", ramps = FALSE)
```

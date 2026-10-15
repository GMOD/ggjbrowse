# Synteny

`jbrowse_synteny()`, `jbrowse_dotplot()` and `jbrowse_multiway()` read
alignments (PAF, indexed PIF, MUMmer delta, UCSC chain, MashMap) through
JBrowse's synteny displays.

- A file guessed from its extension names its sides `query` and `target`. A
  PAF names its sequences by assembly when they read `sample#haplotype#contig`
  ([the amylase locus](gallery.md#amylase_synteny)), and a JBrowse
  adapter config names assemblies, which `assemblies = c(query, target)` puts
  on rows.
- `swap = TRUE` puts the file's target on the top row.
- `labels = c(query = "peach", target = "grape")` writes those names in the
  axis titles, whichever side each is on.
- `guides_jbrowse()` keys a continuous field such as
  `list(field = "identity")` as the display's `0%` to `100%` bar. A threshold
  field such as `mapq` keys its bins as chips.

## Two rows

`jbrowse_synteny()` returns the blocks two rows of regions share and the
ribbons the linear synteny display draws between them, in its colors and
opacity, with each CIGAR insertion, deletion and skip wider than a pixel
([HG002's chr8 inversion](gallery.md#inversion_synteny),
[a CIGAR zoom](gallery.md#amylase_zoom)).

- `view` takes the `LinearSyntenyView`'s settings (`color`, `opacity`,
  `cigarMode`, `minAlignmentLength`, `drawLocationMarkers`, ...).
- `autoplot(s, above = genes, below = genes2)` stacks tracks over either row
  ([genes on both rows](gallery.md#synteny_genes)).
- `geom_jbrowse_ribbons()` draws `s$ribbons`, each foot's span on its row,
  with `thin` set where the band draws a one-pixel line. `s$alignments` holds
  the blocks.
- `geom_jbrowse_offscreen_mates()` draws `s$offscreen_strip`, the marks and
  contig names the band hangs off a row's edge for blocks whose mate the other
  row does not show (`s$offscreen`). `autoplot()` keys each contig from
  `s$offscreen_contigs`.

## Dotplots

`jbrowse_dotplot()` reads the same files through JBrowse's dotplot display,
which draws a line per block, or per CIGAR step where a block is wide enough
on the plot, and `view` takes the `DotplotView`'s settings
([HG002's chr8](gallery.md#inversion_dotplot)). With no regions given,
each axis shows every contig of its assembly in natural order, so a
whole-genome dotplot is one call
([peach against grape](gallery.md#whole_genome)).
`geom_jbrowse_dotplot()` draws its lines.

## Several genomes

`jbrowse_multiway()` reads an all-vs-all PAF, an indexed PIF or a
multi-genome adapter config through JBrowse's multi-way synteny display. It
returns an anchor lane over the locus given, a lane per other genome in the
alignments' own coordinates, and the ribbons between adjacent lanes
([TP53 in four primates](gallery.md#multiway)).

- `lanes` picks and orders the genomes, or pins a lane to a contig.
- `display` takes the display's settings (`ribbonColor`, `drawCurves`,
  `bridgeSkippedLanes`, ...).
- `autoplot(m, above = genes)` stacks tracks over the anchor lane.

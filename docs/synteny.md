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
- `display` takes the view's or display's settings as JBrowse writes them,
  and the result's `display` holds them as the engine resolves them, every
  default filled in: `s$display$encoding$opacity$value` is the level the ribbons
  draw at.
- Each reader lays its regions out 1,500 px wide unless given `width`, as
  the linear readers do.

## Regions on each axis

A result holds every axis's regions in one `regions` frame, its `axis`
naming the row, axis or lane a region lays out on and `assembly` the
assembly it is on.

- Each row of a synteny view (`query`, `target`) and each axis of a dotplot
  (`x`, `y`) takes the forms a linear reader's `region` takes: loci, bare
  contigs, a frame of `chrom`, `start` and `end`, a GRanges, or `NULL` for
  every contig of its assembly.
- A frame with an `assembly` column gives each axis its own rows, so
  `jbrowse_synteny(paf, jbrowse_contigs(paf), jbrowse_contigs(paf))` reads
  both whole genomes, as does `jbrowse_synteny(paf)`.
- A table names a region per axis, `query_region` and `target_region` or
  `x_region` and `y_region`, each a row of `regions`.
- A region on a contig the file lacks reads empty, warns and names why in
  `regions$missing`. Rows with none of the file's contigs between them stop
  with `ggjbrowse_missing_contig`.

## Two rows

`jbrowse_synteny()` returns the blocks two rows of regions share and the
ribbons the linear synteny display draws between them, in its colors and
opacity, with each CIGAR insertion, deletion and skip wider than a pixel
([HG002's chr8 inversion](gallery.md#inversion_synteny),
[a CIGAR zoom](gallery.md#amylase_zoom)).

- `display` takes the `LinearSyntenyView`'s settings (`color`, `opacity`,
  `cigarMode`, `minAlignmentLength`, `drawLocationMarkers`, ...).
- `autoplot(s, above = genes, below = genes2)` stacks tracks over either row
  ([genes on both rows](gallery.md#synteny_genes)).
- `geom_jbrowse_ribbons()` draws `s$ribbons`, each foot's span on its row,
  with `thin` set where the band draws a one-pixel line. `s$alignments` holds
  the blocks, each keyed by `alignment`, and
  `s$alignments <- subset(s$alignments, identity > 0.99)` drops the ribbons
  of the blocks it leaves out.
- `geom_jbrowse_offscreen_mates()` draws `s$offscreen_strip`, the marks and
  contig names the band hangs off a row's edge for blocks whose mate the other
  row does not show (`s$offscreen`). `s$keys$offscreen_query` and
  `s$keys$offscreen_target` key each contig from `s$offscreen_contigs`.

## Dotplots

`jbrowse_dotplot()` reads the same files through JBrowse's dotplot display,
which draws a line per block, or per CIGAR step where a block is wide enough
on the plot, and `display` takes the `DotplotView`'s settings
([HG002's chr8](gallery.md#inversion_dotplot)). With no regions given,
each axis shows every contig of its assembly in natural order, so a
whole-genome dotplot is one call
([peach against grape](gallery.md#whole_genome)).
`geom_jbrowse_dotplot()` draws its lines, and assigning `d$alignments` drops
the segments of the blocks it leaves out. `height` sets the plot's pixels
up, as many as `width` unless given.

## Several genomes

`jbrowse_multiway()` reads an all-vs-all PAF, an indexed PIF or a
multi-genome adapter config through JBrowse's multi-way synteny display. It
returns an anchor lane over the locus given, a lane per other genome in the
alignments' own coordinates, and the ribbons between adjacent lanes
([TP53 in four primates](gallery.md#multiway)).

- `lanes` picks and orders the genomes, or pins a lane to a contig.
- `display` takes the display's settings (`ribbons` with its `encoding` and
  `mark`, `hideUnlabelled`, ...).
- `autoplot(m, above = genes)` stacks tracks over the anchor lane.
- `region` takes the anchor's loci, a frame or a GRanges, but not a bare
  contig: a multi-genome file states no contig's length. `regions` holds the
  anchor's regions and each other lane's window, its `axis` the lane's
  assembly.
- `jbrowse_stack()` refuses a comparative result, since it lays out several
  axes: `autoplot()` with `above` and `below` stacks tracks on its rows.

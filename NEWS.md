# ggjbrowse (development version)

- The engine decides what R copied: the read-connection band's scale
  (`arc_scale`, so an arc's apex keeps clear of the band's edge by half its
  widest stroke and the read cloud's axis reaches 1,000 bp at least), a
  deletion label's and a mismatch letter's opacity, a comparative key's chip
  alpha, the genotype rows' background (drawn under a picked cell color too)
  and height, and the anchor region a multi-way block lies in.
- A feature listed for a label that reaches in from the margin keeps its
  file's fields, and a VCF record there reads, where it errored.
- A record two regions list has its `samples` and `annotations` once.

- The read-connection band hangs from its top under `readConnections`'
  default placement, `"down"`, the read cloud included, and rises from its
  bottom under `"up"`, as the display draws it; the cloud drew upside down.
- `geom_jbrowse_letters()` fades each mismatch letter as its box fades under
  the base-quality encoding, and writes none under 5% opacity; it takes
  `width`, as the box geoms do.
- A multi-way key's chips are opaque, as the display draws them.
- `jbrowse_dotplot()`'s `skipped` is `count`, `total`, `x_chroms` and
  `y_chroms`, each contig list a character vector.

- A locus past its contig's end stops at the end the file or assembly
  states, as the location box cuts it, and one starting past the end is a
  `ggjbrowse_bad_region` error.
- The ruler ticks base N at its left edge, 0-based N - 1, as the view does,
  where the tick sat a base to the right.
- A deletion's length label fades and centres on the part of it in view,
  so one running off the panel keeps its label and a sliver of one in a
  second region loses it.
- `jbrowse_wiggle()` keys the engine's categorical colors, such as a
  threshold scale's bands, in `keys$scores`, and the one-pixel floor of a
  mark is a pixel of the width the track draws at.
- `as.matrix()` of a wiggle result is `NA` in a pixel no bin covers, as the
  display's clustering script writes it, snaps a column count a few ulps
  short of a whole number, and of a result with no bins, or a variants
  result with no cells, returns the matrix where it errored. A record with
  no ALT takes one column, as the display's matrix gives it.
- The genotype rows' grey background draws only while reference calls draw
  no cell, and survives a key two stacked tracks share.
- A VCF record two regions list has its `samples` listed once.
- `geom_jbrowse_ribbons()` and its multi-way kin map any ribbon column, and
  `kind` and `alignment` when no ribbon is thin.
- `jbrowse_multiway()`'s `lanes` keep the row `labels` in `display`, and its
  new `genomes` takes the lanes' assemblies.
- A `MultiWiggleAdapter`'s `bigWigs` may name local files.
- `jbrowse_rows()` loses `groups`, which no reader's display accepts.

- `jbrowse_variants()` on a one-sample VCF lists its records as the display
  draws them, with no genotype row: every record the filters keep, as on a
  sites-only VCF, so a calls VCF with no GT lists its records where it
  listed none. `rows` and `cells` stay empty unless `genotypes$show` is
  `TRUE`, and while no genotype row draws, `minorAlleleFrequency` and
  `missingness` filters test nothing.

- `jbrowse_features()` warns when the layout's rows cannot hold every
  feature in view and some are left out of its tables, with the count; a
  dense track read whole under `heightMode = "grow"` dropped them silently.
- A table the engine always returns that comes back missing is an error
  naming it, where it read as an empty table.

- A multi-way read frames its anchor locus in its contig as the display
  does: it fetches the view's blocks around the window, so `blocks` and
  `ribbons` can run past 0 and `width` for the plot to clip, an alignment's
  `id` and `group` name the fetched block, and the legend and lane order
  cover what was fetched. `jbrowse_multiway(genomes = )` takes the lanes'
  assemblies, which give the anchor's contig its length and stop each lane's
  line at its contig's ends.
- A synteny read drops the ribbons the band culls past its overdraw, and
  `offscreen`, `offscreenContigs` and `offscreenStrip` mark a ribbon whose
  mate lies past the facing row's view, as the band's strip does.

- A locus reads as JBrowse frames it: a window into its whole contig, with
  the view's margin fetched around it, so rows lay out as the display lays
  them out. The tables list what the window shows, plus what a listed row
  draws from just outside it (a read a connector joins, a feature whose
  label reaches in, the bins a line runs in from), and every index between
  tables names a listed row.
- Synteny and dotplot ribbons run to the window's edges, clipped there as
  the band clips them, where they stopped at the requested loci; a ribbon's
  px may fall outside `[0, width]`. The contig comes from the assembly, or
  from the lengths the PAF or chain file states.
- `jbrowse_wiggle()`'s `marks$clip_side` flags a bar its axis clips (1 at
  the top, -1 the bottom) and `clip_strip` gives the strip JBrowse paints
  there; `geom_jbrowse_wiggle()` draws it.
- A deletion that starts before the window fades as it does inside it,
  where it was drawn opaque.
- Alignments rows, connectors and dip reserve come from the display's own
  layout code; a faceted, collapsed pileup puts every read on one row.
- A malformed region (no end, an end before its start) is an error, and a
  refusal over `byte_limit` names the region asked for.

- `jbrowse_variants()`'s `variants$color` is what the records band paints:
  the record's structural-variant class color, or goldenrod (`"#daa520"`)
  for a record with none, where it was the full-dosage genotype color.
  `cells$color` keeps the genotype colors `keys$genotypes` lists.
- A synteny or dotplot result with no color field keys nothing, as the view
  shows no key there; `keys$ribbons` is empty.
- With `coverage = list(show = FALSE)`, `jbrowse_alignments()`'s coverage
  scale is `[0, 0]` with no ticks, as the hidden band draws no axis.
- Rows given as a tree (`rows$tree`) come in the tree's leaf order, as the
  display draws them, where they came in the file's order.
- A multi-way result labels its lanes by `rows$labels` and keeps a chosen
  lane with no alignments in the window.

- Under `unit = "chain"` the `chains` table gives way to `connectors`, every
  connector the pileup draws, and `geom_jbrowse_connectors()` draws each as
  the display does: straight, a curve dipping under its rows or bowing over
  them, or a loop with its arrowhead, placed by the engine at the figure's
  pixels per base, across regions too. `geom_jbrowse_chains()` is its
  deprecated name.
- `jbrowse_alignments()`'s `reads` table drops `category` for `level`, the
  read's value as the color's `scale$domain` names it (`"1"`, `"-1"`,
  `"RL"`, `"short"`, a MAPQ, a tag's value, or `""`). The `arcs` table
  drops `category` for `level`, the `insertSizeAndOrientation` level whose
  color it paints, and `class`, a split read's junction class
  (`"deletion"`, `"duplication"`, `"inversion"`).
- `insertSizeAndOrientation` is the one paired-end read color, ordinal over
  its levels: the `insertSize` and `pairOrientation` color fields are gone,
  and `display = list(normalInsertSize = list(min = 200, max = 600))` pins the
  insert sizes it colors normal.
- `jbrowse_variants()` reads its `display` as a `LinearVariantDisplay`'s,
  the one variant display, and its result's `display$type` says so.
- A variant display's row settings live under `rows`: `unit` is
  `rows$field` (`"sample"` or `"haplotype"`), `rowColor` is `rows$color` for
  a samples file column and `rows$colors` for a color by row name, and
  `rowHeight` is `rows$step`. A haplotype row is named `"<sample> HP1"`,
  `"HP2"`, from 1, and a sample's label or color reaches its haplotypes.
- `jbrowse_rows()` writes a display's `rows`, a named vector of `labels` or
  `colors` as the map the display reads.
- A named vector anywhere in `display` is written as a JSON object, where it
  had lost its names, and a `TRUE` or `FALSE` in a slot named like an array
  slot, such as `rows$axis$labels`, stays a boolean.
- `geom_jbrowse_ribbons()` and `geom_jbrowse_lane_ribbons()` draw a thin
  ribbon one pixel wide, as the band does, where they drew it nearly two;
  `linewidth` widens it.
- A multi-way lane's ruler leaves out a tick label that would run past the
  view's edge, where it drew it cut in half, and a lane names the other
  contigs placing a comparable share (`lanes$also_on`), as the display's
  header does.
- `autoplot()` of a synteny or multi-way result titles each track a named
  list passes as `above` or `below`, as `jbrowse_stack()` titles a named
  track, where it dropped the names.
- The engine is `@jbrowse/tables` 5.0.0-beta.20, bundled in
  `inst/js/jbrowse-tables.js`.

# ggjbrowse 0.11.0

- A CRAM read again in the same session under another assembly decodes its
  bases against that assembly's sequence. A first read against a FASTA that
  spells the contig differently, with no alias table, had left every
  substituted base as N for the rest of the session.
- A CRAM read against a sequence that has no name for its contig warns that
  its bases decode as N, and names the sequence's contigs.
- The engine is `@jbrowse/tables` 5.0.0-beta.18 at jbrowse-components
  d5a0299d00, bundled in `inst/js/jbrowse-tables.js`.

# ggjbrowse 0.10.0

The first tagged release, for reference on GitHub; it is not on CRAN. Install
it with `pak::pak("GMOD/ggjbrowse@v0.10.0")`.

- Readers for BAM and CRAM (`jbrowse_alignments()`), GFF3, GTF, BED and BigBed
  (`jbrowse_features()`), BigWig (`jbrowse_wiggle()`), VCF
  (`jbrowse_variants()`), and PAF, PIF, delta, chain and MashMap
  (`jbrowse_synteny()`), each returning the tables JBrowse draws from.
- A geom per table, `autoplot()` for a whole track, and `jbrowse_stack()` for
  several tracks over one region.
- Display settings in JBrowse 2 v5's spellings, the wiggle display's
  `mark = list(extent = )` included.
- The engine is `@jbrowse/tables` 5.0.0-beta.18, bundled in
  `inst/js/jbrowse-tables.js`.

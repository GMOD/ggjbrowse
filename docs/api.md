# API

- Each reader returns a classed list of data frames, one per kind of record
  or mark the JBrowse display draws, beside a fixed set of other slots.
  `?jbrowse_alignments` and its siblings list every column; this page maps
  the pieces and states the contract.
- Coordinates are 0-based and half-open, as in the engine. A locus string such
  as `"chr12:6,534,201-6,538,700"` is 1-based and inclusive, as JBrowse's
  location box writes it.

## Readers

| Function | Reads | Returns a `jbrowse_*` list of |
|---|---|---|
| `jbrowse_alignments()` | BAM, CRAM | `reads`, `mismatches`, `gaps`, `insertions`, `clips`, `coverage`, `allele_coverage`, `junctions`, `arcs`, ... |
| `jbrowse_features()` | GFF3, GTF, BED, BigBed, VCF | `features`, `subfeatures`, `rects`, `lines`, `chevrons`, `arrows`, `amino_acids`, `labels`, `sections`; for a VCF, `samples` and `annotations` |
| `jbrowse_wiggle()` | BigWig | `bins`, `marks`, `sources`, `rules` |
| `jbrowse_variants()` | VCF | `variants`, `cells`, `rows`, `bands`, `annotations` |
| `jbrowse_synteny()` | PAF, PIF, delta, chain, MashMap | `alignments`, `ribbons`, `offscreen*` |
| `jbrowse_dotplot()` | the same | `alignments`, `segments` |
| `jbrowse_multiway()` | an all-vs-all PAF | `lanes`, `segments`, `blocks`, `ribbons`, `alignments` |

- `file` takes a path, a URL, or a JBrowse track's adapter config, including
  a `FromConfigAdapter` that lists its features inline.
- `region` takes:
  - `NULL`, the default: the whole file, every contig it lists but those its
    index shows empty, each to the length the file or `assembly` states
  - loci, `"chr1:1,001-2,000"`, several drawn end to end
  - a contig name alone, `"chr1"`, for the whole contig
  - a data frame of `chrom`, `start` and `end`, 0-based half-open, such as
    `jbrowse_contigs()` returns
  - a GRanges, which GenomicRanges counts from 1
- A frame with `seqnames`, or with a `width` of `end - start + 1`, reads as
  1-based and is refused.
- `jbrowse_synteny()` takes a `region` per row, `query` and `target`, and
  `jbrowse_dotplot()` one per axis, `x` and `y`, each in any of these forms.
  `NULL` gives an axis every contig of its assembly. A frame with an
  `assembly` column, as `jbrowse_contigs()` and a result's `regions` write
  it, gives each axis the rows on its assembly.
- `jbrowse_multiway()`'s `region` takes loci, a frame or a GRanges for its
  anchor lane. A multi-genome file states no contig's length, so a bare
  contig or `NULL` is refused.
- A whole file or contig is refused past 5 million rows in any one table,
  and `jbrowse_alignments()` past `span_limit`, asking for a region.
- `display` takes a track config's `displays` entry, or for
  `jbrowse_synteny()` and `jbrowse_dotplot()` the view's settings as its
  snapshot writes them.
- `width` is the pixels the reader lays the regions out across, 1,500 unless
  given: the labels and zoom of a feature or quantitative track, the fade of
  a variant matrix, and the CIGAR steps a synteny band or dotplot draws
  follow from it. `jbrowse_dotplot()`'s `height` is its plot's pixels up,
  its `width` unless given.
- `assembly` takes a `jbrowse_assembly()`, an assembly as a JBrowse config's
  `assemblies` entry writes it, or a FASTA or 2bit path alone. Its sequence
  is the reference every reader reads against. A region may name a contig by
  any alias in the assembly's alias table. The engine reads each file under
  the file's own alias, and the reader names the contig in its tables as the
  region does.

## The result

Every reader's result holds these slots, which `?jbrowse_result` documents.
`$` matches a slot's name exactly, so `g$region` is `NULL`, not
`g$regions`:

| Slot | Holds |
|---|---|
| `regions` | a data frame, a row per region: `chrom`, `start`, `end` (0-based, half-open), the axis layout `offset`, `cum_start` and `cum_end`, and `missing`, why the region's contig reads empty, `""` where the file has it |
| `width` | the pixels the regions were laid out across |
| `display` | the display's or view's resolved config, under JBrowse's own camelCase names |
| `keys` | the display's color keys, each a frame of `label` and `color` titled by its `"title"` attribute |
| `ramps` | the display's color bars |
| `y_scale` | the score axis (wiggle) or depth axis (alignments): `domain`, `type`, `symlog_constant`, `ticks`; `NULL` for the others |
| `fields` | the file's declarations of its columns: `info`, `format`, `attributes`, `annotations` |
| `warnings` | what the engine warned of, each also raised as an R warning |

Every other slot is a table, but for a few that hold a class's drawing
state:

- features: `glyph`, `height`, `bp_per_px`
- wiggle: `rendering`, `origin`, `row_layout`, `height`, `bp_per_px`
- alignments: `read_pitch`, `glyph`, `arc_span_max`, and `sequence`, each
  region's reference bases where `assembly` has a sequence
- variants: `tree`, where a caller sets it
- synteny: `offscreen_strip`, `fade_thin`, `has_cigar`, `assemblies`
- dotplot: `height`, `assemblies`, `skipped`
- multiway: `tick_interval`, `rows_vs_anchor`

A comparative result lays out several axes in its one `regions` frame:

- `axis` names the axis a row lays out on: `"query"` (the top row) or
  `"target"` for synteny, `"x"` or `"y"` for a dotplot, and the lane's
  assembly for multiway. `assembly` names the assembly the region is on.
- `display` holds the view's settings, or for multiway the display's
  config, as the engine resolves them, `color` and `opacity` as config
  objects.

## Records and their keys

Every table has a `region` column, the 1-based row of `regions` each row was
read in, and a table with coordinates of its own has `chrom`. A record read
in two regions is listed in each, its copies sharing its key.

A comparative table names a region per axis: a ribbon's `query_region` and
`target_region`, a segment's `x_region` and `y_region`, a multi-way
ribbon's `upper_region` and `lower_region`, each a row of `regions`.
`alignments` names no region, since a block can cross several.

| Parent | Key | Children |
|---|---|---|
| `reads` | `read` | `mismatches`, `gaps`, `insertions`, `clips`, `modifications`, `qualities` |
| `features` | `feature` | `subfeatures`, `rects`, `lines`, `chevrons`, `arrows`, `amino_acids`, `labels`, `samples`, `annotations` |
| `subfeatures` | `subfeature` | `rects`, `lines`, `chevrons`, `arrows`, `amino_acids`, `labels` |
| `variants` | `variant` | `cells`, `annotations` |
| `rows` | `name` | `cells`, by `row` |
| `sources` | `name` | `bins` and `marks`, by `source` |
| `alignments` (synteny, multiway) | `alignment` | `ribbons` |
| `alignments` (dotplot) | `alignment` | `segments` |

- A child row holds its parent's key, and its `region` identifies the
  parent's copy in that region.
- The geoms join a parent's columns onto its children when they draw: a
  read's `row` and `section`, a feature's `name`, `type`, `strand` and
  `section` (and a chevron's, arrow's or label's `start` and `end`), and a
  cell's record `start`, `end` and `opacity` and row `sample`. Reassigning
  `a$reads$row` moves a read's marks.
- Assigning a parent, `a$reads <- subset(a$reads, mapq >= 20)`, drops the
  children's rows whose record it leaves out, so drawing, `as.matrix()` and
  `saveRDS()` all see the same records. A dplyr tibble assigns the same way.
- A geom takes a data frame in place of a child or summary table, but not a
  parent's: `geom_jbrowse_reads(data = ...)` errors, pointing to assignment.
- The summary tables an alignments track computes over every read
  (`coverage`, `allele_coverage`, `modification_coverage`,
  `interbase_coverage`, `indicators`, `junctions`, `chains`, `overlaps`,
  `arcs`, `arc_ticks`) don't recompute after `reads` is filtered. Drawing one
  then warns once: read again with `display = list(filter = ...)`.

## Errors

| Class | When |
|---|---|
| `ggjbrowse_too_large` | the regions pass `span_limit` or `byte_limit`, or a whole file or contig passes 5 million rows in a table |
| `ggjbrowse_missing_contig` | no region's contig is one the file lists; a region on a missing contig beside others warns and reads empty |
| `ggjbrowse_bad_display` | a `display` slot holds a value the config schema refuses, or a synteny view's `cigarMode` or `fadeThinAlignmentsMode` one it does not take, which the error names |
| `ggjbrowse_bad_region` | a region is neither a locus nor a contig of the file, or a frame reads as 1-based |
| `ggjbrowse_parent_data` | a geom is handed a parent table as `data` |
| `ggjbrowse_forked` | the engine was started in a forked worker |

## Helpers

- `jbrowse_assembly(sequence, aliases)`: an assembly, its reference sequence
  and its contigs' alias table, for any reader's `assembly`
- `jbrowse_contigs(x)`: an assembly's, a file's, or each side of a pairwise
  file's contigs as `assembly`, `chrom`, `start` and `end`, a region frame
  any reader takes; a linear file adds `has_data`
- `jbrowse_sequence(assembly, region)`: each region's reference bases as a
  string
- `jbrowse_reset()`: restarts the engine and clears its caches

## Drawing

- `autoplot(result)` composes the geoms into the browser's track.
  `autoplot(a, "pileup")` draws one band of an alignments track.
- `jbrowse_stack(...)` stacks results, and your own ggplots, on one x axis,
  each named one titled. A comparative result, which lays out several axes,
  stops the stack; `autoplot(s, above = genes)` stacks tracks over a
  synteny or multi-way result's rows. A key several tracks share is drawn
  once, listing every entry any of them holds. `jbrowse_height()` and
  `jbrowse_width()` measure the figure.
- A data frame of your own with `chrom` and `start` and `end`, or
  `position`, in a `ggplot(<result>)` layer or a plot in a stack, lands in
  the region it falls in, a row over two regions in each.
- `autoplot()` and `jbrowse_stack()` take `width` to draw at another size.
  A stack otherwise draws at the widest `width` any of its tracks was read
  at. Drawing a feature or quantitative track narrower than it was laid out
  for warns, since its labels and zoom were chosen for more pixels.
- One geom per table, all named `geom_jbrowse_*()`: `reads`, `mismatches`,
  `deletions`, `insertions`, `clips`, `coverage`, `alleles`, `junctions`,
  `glyphs`, `amino_acids`, `labels`, `wiggle`, `genotypes`, `ribbons`,
  `dotplot`, `lanes` and the rest. A geom's `width` defaults to the
  result's.
- Each geom paints its table's `color` as it is, `I(color)`, so the
  display's colors train no scale and a mapping of your own takes its scale
  beside them: `geom_jbrowse_reads(aes(fill = mapq)) +
  scale_fill_viridis_c()`.
- `guides_jbrowse(x, keys, ramps)` draws `x$keys` and `x$ramps` as guides
  that take the plot's theme. Edit a key to relabel or recolor an entry.
- A geom draws the decisions its table holds: `junctions$shown`,
  `allele_coverage$shown`, each mark's `fade` and `variants$opacity`.
  Editing the column draws otherwise, with no read:
  `a$junctions$shown <- a$junctions$count >= 5`.
- `scale_x_jbrowse()` and `scale_y_jbrowse()` add the regions' axis and the
  `y_scale`.
- A saved result draws without its files: the sequence band reads
  `sequence`, and the genotype fade `variants$opacity`.

## Notes

- `reads$id` is the record's own id, the same in every call, where mates and
  supplementary segments share a `name`.
- Every reader reads several regions in one engine call, which lays them out
  as one view, with one pileup layout, one set of variant rows and one depth
  axis.
- `coverage` and the band's other tables stop at the region's ends, though
  the reads reach past them. A band the `display` hides (`showCoverage =
  FALSE`, ...) reads as empty tables.
- Under `display$facet`, an alignments or feature result lists its
  `sections` in stacking order, and `section` names the one a row sits in.
- `fade` is the share of the depth a mark's site holds, zero under the
  display's threshold. R applies it at the figure's pixels per base.
- `jbrowse_variants()` returns its `rows` in the display's arrangement, by
  the `rows` config (`domain`, `labels`, `kept`) and the `facet` bands.
  `tree` takes an `hclust`, or a list of them, one per band.
- `jbrowse_variants()` gives each record its `filter` and a column per INFO
  key, named as the key, or `info_<key>` where the key names a record column
  such as `type`. The header's `Number` decides the shape: a flag is logical,
  `Number=1` a plain column, and `A`, `R`, `G` or `.` a list column of
  vectors. `fields$info` names each key's column.
- `jbrowse_variants()` keeps a record only where a sample calls an allele, as
  the genotype matrix does, so a calls VCF with no GT keeps none.
  `jbrowse_features()` reads every record, as the plain variant track draws
  them: `features` gains `ref`, `alt`, `filter` and the INFO columns, and
  `samples` holds a row per record and sample with a column per FORMAT key,
  typed as INFO is.
- A VCF's SnpEff `ANN` or VEP `CSQ` reads as `annotations`, a row per record
  and consequence, with a column per field the key's header line names, all
  character. `variant` (or `feature`) names the record, and the raw string
  leaves the INFO columns.
- `jbrowse_features()` gives each feature and transcript a column per field
  the file states, keyed as the file writes it: a GFF3's `ID`, `Name`,
  `Parent` and the rest beside the table's own `id` and `name`, and a BigBed's
  autoSql columns, typed as it declares them. A key naming a table column
  takes `attr_<key>`. A GFF3 has no header, so its columns are the keys the
  regions' records use, `NA` where a record lacks one.
  `fields$attributes` lists each key's type and column.
- `variants$opacity` is the fade the display paints a record's cells at
  where records share a pixel at the result's `width`, through the engine's
  `variantDensityAlpha`. `v$variants$opacity <- 1` turns it
  off.
- `byte_limit` refuses regions whose compressed size, read off the index,
  passes it. The default is 200 MB. `jbrowse_alignments()`'s `span_limit`
  refuses regions spanning more than 10 Mb in all, since the coverage band
  holds every base.
- A setting the engine drops without an error, such as a column sort with no
  `pos`, raises an R warning naming it, and lands in `warnings`.
- A CRAM needs an `assembly` with a sequence, because it stores reads
  against the FASTA it was encoded with.
- The bisulfite layer needs `assembly`, because it calls each cytosine by
  comparing the read's base to the reference's.
- `jbrowse_features()` with `assembly` translates each CDS once `width`
  gives a base a pixel or more, striping `rects` by codon and listing each
  residue in `amino_acids`. A residue's `text` is the letter the display
  writes at 8 pixels a base or more, numbered once the codon fits the number.

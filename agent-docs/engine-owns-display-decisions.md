# Handoff: the engine owns the display's decisions

- **Status:** done (2026-10-15). Since 2026-10-19 the engine's
  `alignmentMarkAlpha` sizes and fades the pileup's marks, so R copies two
  small rules, listed below, and still draws some insertion details
  differently.
  jbrowse-components 42060cf6cf returns the colors, keys, junction flags and
  read pitch; ggjbrowse's step-2 branch reads them and deletes its copies.
- **Payoff:** R stops copying rules JBrowse writes once, so a change there
  reaches ggjbrowse through the bundle, and JBrowse's parity tests cover it.

## The rule

- The engine decides anything JBrowse decides from the track's config and
  palette: a mark's color, a constant of the display's geometry.
- R decides anything that depends on the figure: pixels per base, inches,
  label sizes, which ggplot layer draws what.
- A geom draws the decision a column holds (`shown`, `fade`, `opacity`)
  and takes no argument to re-decide it. A caller overrides one by editing
  the column, which needs no read.
- A value that needs both comes from an engine call that takes the
  figure's side as an argument. `variantDensityAlpha` takes bp per px.
  `alignmentMarkAlpha` takes px per bp and each mark's `fade`, `length` and
  `quality` columns, and R's `mark_alpha()` calls it once per geom build, so
  a caller who edits `fade` changes the draw.

## What the file's fields are: nobody's decision

> A table is the region's records in tidy form: every per-record field the
> file states is a column, every per-base field the display decodes is its
> own long table, and the display's decisions (row, color, section, fade)
> ride along as extra columns. JBrowse chooses what it adds, never which of
> the file's fields R sees, and R never re-parses a string the engine
> already parsed.

What follows from it, settled 2026-10-15:

- **No `fields =` argument.** No reader projects the file's fields. CIGAR,
  `seq`, `qual`, MD, MM and ML stay out, since they are per-base: the engine
  already decodes them into `gaps`, `mismatches`, `modifications` and
  `qualities`.
- **BAM and CRAM tags are the one exception, for now.** They arrive once
  binary transport makes them free, or when someone asks for one.
- **A GFF has no header**, so its attribute columns are the keys the
  region's records use. Binding regions fills a missing column with `NA`, as
  INFO does.
- **A whole file is the same reader with no region**, under `byte_limit` and
  a cap of 5 million rows in any one table (genotype cells above all). Past
  the cap it refuses, asking for a region.
- **A name collision takes the flat prefix INFO uses:** `info_type`,
  `attr_type`.
- **VEP's `CSQ` and SnpEff's `ANN` parse into `annotations`**, a row per
  record and consequence, its columns the header's `Format:` field list. The
  raw strings no longer reach R.

How each piece reaches R:

- `features` and `subfeatures` carry a column per field the file states,
  keyed as written: a GFF3's `ID` and `Name` sit beside the engine's own `id`
  and `name`. `fields$attributes` declares each one's type.
- `reads$id` is the record's own id, the same in every call, where `read` is
  the read's key within one call, shared by its copies in several regions.
- A gene's marks name their transcript by its `subfeature` key, so an exon
  joins to its transcript.
- The engine sends numbers as bytes and repeated text as codes
  (`from_wire()` in `R/engine.R`), so a column's shape is the engine's, never
  jsonlite's guess.

## What the engine now decides

- **Colors:** every alignments table's `color`, base colors resolved
  through `effectiveBaseColors()`, so the muted bases under a modifications
  or bisulfite layer come along.
- **Keys:** the base key (empty while muted), the junction key and the
  interbase key, in the browser's own wording.
- **Junctions:** `shown` under `filter.junctions`, `arcHeight` from
  `sashimiArcHeightFraction()` and `strokeWidth` from
  `sashimiStrokeWidth()`, in px, which R draws at 150 px an inch.
- **Alleles:** `shown` under `coverage.alleles.filter`.
- **Read pitch:** `pileup.mark.height` and its spacing, which R's read inset
  takes.
- **Defaults:** each track's resolved config, as the result's `display`,
  and for the synteny, dotplot and multi-way views their resolved settings.

R keeps `base_colors()`, the engine's palette, for the sequence band's
reference bases.

## Rules R still copies

- **A one-pixel floor.** `geom_jbrowse_deletions()` widens a deletion to a
  pixel, as gap.slang does, and `position_rects()` draws a one-base mark a
  pixel wide at least.
- **`label_fade()`:** `labelFadeOpacity`, the ramp that fades a deletion's
  length label, at the thresholds the engine returns in `glyph$label`.
- **A key chip's alpha.** `jbrowse_synteny()` and `jbrowse_dotplot()` draw
  their key's chips at the ribbons' alpha floored at 0.45, synteny-core's
  `LEGEND_CHIP_ALPHA_FLOOR`, and `jbrowse_multiway()` at 0.6. The engine
  could return the chip alpha beside `legend`.

## Where R still draws differently from JBrowse

- **Insertion details.** R draws no serifs (on an insertion under 10 bp at
  3 px per bp or more), and `geom_jbrowse_letters()` writes every insertion's
  length beside it, where JBrowse writes `(N)` beside a small one and `N`
  inside a boxed one.

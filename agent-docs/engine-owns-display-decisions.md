# Handoff: the engine owns the display's decisions

- **Status:** done (2026-10-15), with three engine fixes found since listed
  under "Waiting on the engine". Since 2026-10-19 the engine's
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

- **The genotype rows' background.** `variants_result()` keeps, as the
  genotype key's `background` attribute, the color of the entry the engine
  labels "Homozygous reference" (or "Reference" for haplotype rows), the
  display's `REFERENCE_COLOR`, since a cell color field's entries come
  before it; `reference_rows()` fills each row with it. The engine could return the row
  background beside `legend`.

## Where R still draws differently from JBrowse

- **Insertion details.** R draws no serifs (on an insertion under 10 bp at
  3 px per bp or more), and `geom_jbrowse_letters()` writes every insertion's
  length beside it, where JBrowse writes `(N)` beside a small one and `N`
  inside a boxed one.
- **A multi-row wiggle's range.** JBrowse writes the rows' shared domain,
  `[0, 200]`, at the top right; R's rows show no scale. This is R's to draw,
  from `y_scale$domain`.
- **A transcript past the view's edge.** JBrowse marks it with `≫` at the
  edge (the `amino_acids` gallery figure); R does not. Not yet traced to
  where JBrowse decides it.

## Waiting on the engine

Found 2026-11-02 comparing the gallery halves. Each is JBrowse's decision,
so R does not copy it: the fix is in jbrowse-components, then an engine bump
here. On 2026-11-03 local jbrowse-components `main` was 195 commits past the
pinned v5.0.0-beta.20, none of them fixing these, and a bump from it takes
those too.

- **Synteny ribbons stop short of the view's edges.** `inversion_synteny`'s
  top-right corner and `amylase_zoom`'s corners are white where JBrowse
  fills them. In chr8, alignment 2 runs to 25,080,663, but its ribbon ends
  at 13,998,216, where the bottom row's region ends.
  - Cause: the live display's displayed region is the whole contig, and
    `syntenyFetchRegions()` buffers the visible window inside it. The
    tables engine (`plugins/linear-comparative-view/src/tables.ts`) makes
    the requested loci the displayed regions and passes them as the fetch
    and target windows, so `buildSyntenyGeometry` trims each ribbon where
    either row's window ends.
  - `place()` clamps each foot to its region on its own, so a wider window
    alone would skew a partly off-screen ribbon, where it should clip.
    Fixing it changes how the tables lay out an axis: each ribbon's region
    and its px.
- **A dotplot with no color field keys "Default / alignment".** The views'
  `colorScales` (synteny-core `TrackColorsMixin`) show no key unless
  `hasLegendKey`, which needs a track, strand, measurement or column field.
  `comparativeLegend()` returns chips for every field, and
  `packages/tables/src/synteny.test.ts` asserts the default
  `['match', 'insertion', 'deletion']`, so the tables' key is deliberate
  there. Settle which rule the tables follow before changing it.
- **Clipped wiggle bars carry no flag.** JBrowse paints a `CLIP_STRIP_PX`
  strip of `CLIP_STRIP_COLOR` (render-core `clipStrip`) across a row's edge
  where a bar passes its domain; R paints none. In `dnase_rows`, K562 peaks
  at 427 on a 0–260 domain. The engine could flag each clipped mark in
  `marks` and return the strip's color and px.
- **A deletion's fade depends on where the view starts.** In
  `hg008t_cuzd1_deletion.bam` the 28 reads' 1,798 bp deletion at
  chr10:122,835,344 fade 0.451 read over `122,835,000-122,836,200`, 1 over
  `122,836,201-122,837,500` and 0.847 over `122,835,345-122,837,500`, and
  each region's copy of it differs in a multi-region read.
  - Cause: `deletionStartFrequencies()` anchors a deletion's frequency at its
    start, but `computePositionFrequencies()` takes the depth through
    `depthNear()`, which clamps the position to `coverageStartPos`, while
    the count comes from the `rare` map of starts outside the window. A
    deletion starting before the window is counted against the depth inside
    itself, where the deleted reads are missing, so it saturates.
    (`plugins/alignments/src/shared/computeFrequenciesAndThresholds.ts`.)
  - The browser runs the same pipeline per region, so a block that starts
    inside a deletion likely draws it opaque too.

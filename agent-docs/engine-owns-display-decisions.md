# Handoff: the engine owns the display's decisions

- **Status:** proposed, unstarted (2026-10-12). Today's figures match the
  browser.
- **When:** with the next change to jbrowse-components'
  `plugins/alignments/src/tables.ts`, so it rides that change's release.
- **Payoff:** R stops copying rules JBrowse writes once, so a change there
  reaches ggjbrowse through the bundle, and JBrowse's parity tests cover it.

## The rule

- The engine decides anything JBrowse decides from the track's config and
  palette: a mark's color, a constant of the display's geometry.
- R decides anything that depends on the figure: pixels per base, inches,
  label sizes, which ggplot layer draws what.
- A value that needs both comes back from the engine as its config-side
  input, and R finishes it. The sub-pixel fade works this way: the engine
  returns each mark's `fade` share, and R's `mark_alpha()` applies it at
  the figure's pixels per base.

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
  and `name`. `attribute_fields` declares each one's type.
- `reads$id` is the record's own id, the same in every call, where `read` is
  a row number within one call.
- A gene's marks name their transcript in `subfeature`, so an exon joins to
  its transcript.
- The engine sends numbers as bytes and repeated text as codes
  (`from_wire()` in `R/engine.R`), so a column's shape is the engine's, never
  jsonlite's guess.

## What R copies today

`alignments_region()` in `R/alignments.R` colors these tables itself, from
`jbrowse_colors()` (`R/autoplot.R`), the palette that
`packages/tables/src/colors.ts` exports:

| Table | R's rule |
|---|---|
| `mismatches`, `allele_coverage` | `colors$bases[base]` |
| `gaps` | skip or deletion color by `type` |
| `insertions`, `clips` | insertion, softclip or hardclip color |
| `interbase_coverage`, `indicators` | color by `type` |
| `coverage` | the coverage grey |
| `junctions` | forward, reverse or read color by `strand` |

The engine already returns `color` for `reads`, `modifications`,
`qualities`, `modificationCoverage`, `arcs` and `arcTicks`.

R also copies these rules:

- **The muted bases.** `paints_modifications()` and `MUTED_BASE <- "#888"`
  in `R/alignments.R` grey every mismatched base while the per-base layer
  paints modifications. JBrowse writes that rule once, in
  `effectiveBaseColors()` (`plugins/alignments/src/features/mismatch/baseColors.ts`),
  gated on `isModificationScheme()`. R gates it on a list of field names,
  and reads `#888` where JBrowse reads `palette.mutedSnpBase`.
- **The keys.** `alignments_region()` builds the Base key and the junction
  strand key from the palette.
- **Sashimi arc heights.** `junction_arcs()` in `R/autoplot.R` copies
  `sashimiArcHeightFraction()` from
  `plugins/alignments/src/features/sashimi/bandFeed.ts`: 0.3 to 0.95 of the
  band, on a log scale from 50 bp to 100 kb. The two agree today.
- **The read inset.** `read_inset()` in `R/geoms-alignments.R` copies
  `featureSpacingForHeight()` from
  `LinearAlignmentsDisplay/menus/compactnessPresets.ts`: a 1 px gap under a
  read taller than 3 px.

## Where R draws differently from JBrowse

R's drawing of these two marks departs from the display's. Each wants a
decision, and a test once decided:

- **Junction stroke width.** JBrowse's `sashimiStrokeWidth()` draws
  `log(count + 1)` px, clamped to 1–16. R's `junction_paths()` draws
  `log(count + 1) / 2` in ggplot linewidth units (about 3 px each at the
  stack's 150 px an inch), clamped to 0.4–4.
- **Long insertions.** JBrowse's insertion shader never fades an insertion
  of 10 bp or more, and widens it to `length × px per bp / 3`, 1 to 5 px,
  or to a box for its length label once it spans 15 px at a feature height
  of 5 or more. R's `geom_jbrowse_insertions()` fades every insertion by
  px per bp squared and draws all of them as one 0.4 linewidth tick.

## The change

1. **jbrowse-components**, `plugins/alignments/src/tables.ts`:
   - Return a `color` array for each table in the first table above,
     resolving base colors through `effectiveBaseColors()` so the mute
     rule comes along.
   - Return the base key (empty while muted) and the junction strand key
     as the read key is returned.
   - Return each junction's arc height fraction from
     `sashimiArcHeightFraction()` and its stroke width from
     `sashimiStrokeWidth()`.
   - Return the read pitch (`featureHeight` and its spacing) beside
     `settings`.
   - Test each in `packages/tables`; adding columns breaks no host.
2. **The bundle:** `scripts/build-js.sh` from that checkout, which stamps
   it "unreleased", or `scripts/take-release.sh` once `@jbrowse/tables`
   publishes it.
3. **ggjbrowse:**
   - Read the colors and keys in `alignments_region()`, and delete
     `MUTED_BASE` and `paints_modifications()`.
   - Keep `jbrowse_colors()` only for `sequence_plot()`'s reference bases,
     or have the engine return those too.
   - Take arc heights in `junction_arcs()`, stroke widths in
     `junction_paths()` and the inset in `read_inset()` from the result.
   - Settle the long-insertion divergence, here or in a follow-up.

## Order and checks

- Land step 1 first. The current R code ignores the new columns, so
  jbrowse-components' downstream canary, which runs
  `scripts/check-against-jbrowse.sh`, stays green.
- Step 3 needs the new bundle. Land it together with step 2, never ahead.
- Expect the README figures to come out pixel-identical, except where a
  divergence above gets fixed. Redraw them with
  `Rscript -e 'devtools::load_all(); source("scripts/figures.R")'` and
  compare against main before committing.

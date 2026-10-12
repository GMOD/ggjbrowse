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

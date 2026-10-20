# Alignments

`jbrowse_alignments(file, region, display, width, assembly)` reads a BAM or
CRAM region through JBrowse's alignments display.

- **Files:** a path or a URL, which `curl` reads in byte ranges. CRAM needs
  `assembly`: a FASTA with its `.fai`, a bgzipped one with its `.fai` and
  `.gzi`, a 2bit, or a `jbrowse_assembly()` naming one.
- **Bands:** the track stacks the display's splice-junction arcs, its
  coverage band with each mismatch position's alleles, and its pileup, each
  read's mismatches, indels and clips faded where the browser fades them.
  `autoplot(a, "pileup")` draws one band alone, as do `"coverage"`,
  `"junctions"` and `"arcs"`.

## Display settings

`display` is the display's config, written as a JBrowse track config writes
it, so `jbrowse_alignments()` colors, filters, sorts and splits reads as
JBrowse does:

- `color = list(field = "mapq")`, or
  `list(field = "insertSize", domain = c(200, 600))`. A bare string such as
  `color = "purple"` is one constant color.
- `filter = list(flagExclude = 1796, tags = list(list(tag = "HP", value = "1")))`
  keeps reads by their flags and tags; its `junctions` and `connections`
  pick the splice junctions and read connections drawn:
  `filter = list(junctions = list(minReads = 5, nonCanonical = "exclude"))`.
- `sort = list(field = "base", pos = 122836300)` sorts by the base at a
  position ([a somatic deletion](gallery.md#deletion)).
- `pileup = list(readHeight = 3)` packs reads compact
  ([compact reads](gallery.md#compact_reads)).
- `baseColor = list(field = "modifications")` paints MM/ML modification calls
  ([CpG methylation](gallery.md#methylation)).
- `baseColor = list(field = "bisulfite")` paints the methylation bisulfite or
  EM-seq reads show against `assembly`, in the context
  `modifications = list(cytosineContext = "CHH")` names
  ([CHH methylation](gallery.md#bisulfite)).
- `readConnections = list(show = TRUE)` adds the read-pair band, an arc per
  mate pair and split read colored by insert size and pair orientation.
  `mark = "cloud"` plots each abnormal pair at its insert size instead.
- `unit = "chain"` puts a read's mate or split segments on one row, joined by
  a hairline ([COLO829 split reads](gallery.md#long_reads_sv)).
- `facet = list(field = "split")` splits the reads into the display's stacked
  sections, each packed separately
  ([an inversion](gallery.md#inversion)).
- `coverage`, `pileup`, `sashimi` and `readConnections` each take the
  band's own `show` and height, as a track config writes them:
  `coverage = list(show = FALSE)`,
  `sashimi = list(placement = "down")`,
  `readConnections = list(show = TRUE, showOrdinaryPairs = FALSE)`.

Each read comes back on the row the display lays it on, with the fill the
display paints it, and `a$keys$reads` is the display's key. A mismatch,
gap, insertion or clip names its read by `read`, so
`a$reads <- subset(a$reads, mapq >= 20)` drops the filtered reads' marks
too. The coverage band and the junctions keep counting every read; read
again with `display = list(filter = ...)` to recompute them.

The engine decides what each band shows, and `a$display` holds the
resolved config it decided from:

- Every mark comes back with the `color` the display paints it, and
  `a$keys` with the browser's keys and titles, such as "Splice junctions".
- `junctions$shown` marks the arcs the sashimi band draws, by
  `filter$junctions`. Each arc rises to its `arcHeight` share of the band and
  strokes `strokeWidth` px, which R draws at 150 px an inch.
- `allele_coverage$shown` marks the alleles the coverage band stacks, by
  `coverage$minAlleleFraction`.
- `autoplot()` and the geoms draw what these columns say, so editing them
  draws otherwise: `a$junctions$shown <- a$junctions$count >= 5`, or
  `a$mismatches$fade <- 1` to stop a sub-pixel mismatch fading.

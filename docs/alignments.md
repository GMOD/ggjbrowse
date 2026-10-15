# Alignments

`jbrowse_alignments(file, region, display, assembly)` reads a BAM or CRAM
region through JBrowse's alignments display.

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
- `filter = list(flagExclude = 1796, tagFilters = list(list(tag = "HP", value = "1")))`
- `sort = list(type = "basePair", pos = 122836300)` sorts by the base at a
  position ([a somatic deletion](gallery.md#deletion)).
- `featureHeight = 3` packs reads compact
  ([compact reads](gallery.md#compact_reads)).
- `baseColor = list(field = "modifications")` paints MM/ML modification calls
  ([CpG methylation](gallery.md#methylation)).
- `baseColor = list(field = "bisulfite")` paints the methylation bisulfite or
  EM-seq reads show against `assembly`, in the context
  `modifications = list(cytosineContext = "CHH")` names
  ([CHH methylation](gallery.md#bisulfite)).
- `readConnections = "arc"` adds the read-pair band, an arc per mate pair and
  split read colored by insert size and pair orientation. `"cloud"`
  plots each abnormal pair at its insert size instead.
- `unit = "chain"` puts a read's mate or split segments on one row, joined by
  a hairline ([COLO829 split reads](gallery.md#long_reads_sv)).
- `facet = list(field = "split")` splits the reads into the display's stacked
  sections, each packed separately
  ([an inversion](gallery.md#inversion)).
- Any other slot, such as `showOrdinaryPairs = FALSE`.

Each read comes back on the row the display lays it on, with the fill the
display paints it, and `a$legend` is the display's key.

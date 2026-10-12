# API

- Each reader returns a classed list of data frames, one per kind of mark the
  JBrowse display draws. `?jbrowse_alignments` and its siblings list every
  column; this page maps the pieces.
- Coordinates are 0-based and half-open, as in the engine. A locus string such
  as `"ctgA:1,001-9,000"` is 1-based and inclusive, as JBrowse's location box
  writes it.

## Readers

| Function | Reads | Returns a `jbrowse_*` list of |
|---|---|---|
| `jbrowse_alignments()` | BAM, CRAM | `reads`, `mismatches`, `gaps`, `insertions`, `clips`, `coverage`, `allele_coverage`, `junctions`, `arcs`, ... |
| `jbrowse_features()` | GFF3, GTF, BED, BigBed, VCF | `features`, `subfeatures`, `rects`, `lines`, `chevrons`, `arrows`, `labels` |
| `jbrowse_wiggle()` | BigWig | `bins`, `marks`, `sources`, `rules` |
| `jbrowse_variants()` | VCF | `variants`, `cells`, `rows` |
| `jbrowse_synteny()` | PAF, PIF, delta, chain, MashMap | `alignments`, `ribbons`, `offscreen*` |
| `jbrowse_dotplot()` | the same | `alignments`, `segments` |
| `jbrowse_multiway()` | an all-vs-all PAF | `lanes`, `segments`, `blocks`, `ribbons`, `alignments` |

- Every result also has `legend` and `ramps`, the keys the display shows.
  `guides_jbrowse()` draws them.
- `file` takes a path, a URL, or a JBrowse track's adapter config.
- `region` takes a locus string, a contig name with `start` and `end`, or a
  data frame of `chrom`, `start` and `end` for several regions.
- `display` takes a track config's `displays` entry. `jbrowse_synteny()` and
  `jbrowse_dotplot()` take the view's snapshot as `view`.

## Helpers

- `jbrowse_contigs(file)`: each side's contigs and lengths in a pairwise
  alignment file
- `jbrowse_sequence(reference, chrom, start, end)`: the reference's bases as
  one string
- `jbrowse_reset()`: restarts the engine and clears its caches

## Drawing

- `autoplot(result)` composes the geoms into the browser's track.
  `autoplot(a, "pileup")` draws one band of an alignments track.
- `jbrowse_stack(...)` stacks results, and your own ggplots, on one x axis.
  `jbrowse_height()` measures the figure.
- One geom per table, all named `geom_jbrowse_*()`: `reads`, `mismatches`,
  `deletions`, `insertions`, `clips`, `coverage`, `alleles`, `junctions`,
  `glyphs`, `labels`, `wiggle`, `genotypes`, `ribbons`, `dotplot`, `lanes` and
  the rest.
- `jbrowse_color()` marks a column as the colors the display paints, so
  `aes(fill = color)` needs no scale. `scale_fill_jbrowse()`,
  `scale_color_jbrowse()` and `scale_x_jbrowse()` add the display's keys and
  axes.

## Notes

- `reads$read` is the id every other alignments table names a read by.
  Filtering `reads` keeps the id.
- `fade` is the share of the depth a mark's site holds, zero under the
  display's threshold. R applies it at the figure's pixels per base.
- `byte_limit` refuses a region whose compressed size, read off the index,
  passes it. The default is 200 MB.
- A CRAM needs `reference`, because it stores reads against the FASTA it was
  encoded with.

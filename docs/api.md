# API

- Each reader returns a classed list of data frames, one per kind of mark the
  JBrowse display draws. `?jbrowse_alignments` and its siblings list every
  column; this page maps the pieces.
- Coordinates are 0-based and half-open, as in the engine. A locus string such
  as `"chr12:6,534,201-6,538,700"` is 1-based and inclusive, as JBrowse's
  location box writes it.

## Readers

| Function | Reads | Returns a `jbrowse_*` list of |
|---|---|---|
| `jbrowse_alignments()` | BAM, CRAM | `reads`, `mismatches`, `gaps`, `insertions`, `clips`, `coverage`, `allele_coverage`, `junctions`, `arcs`, ... |
| `jbrowse_features()` | GFF3, GTF, BED, BigBed, VCF | `features`, `subfeatures`, `rects`, `lines`, `chevrons`, `arrows`, `amino_acids`, `labels`, `sections`; for a VCF, `samples`, `info_fields`, `format_fields` |
| `jbrowse_wiggle()` | BigWig | `bins`, `marks`, `sources`, `rules` |
| `jbrowse_variants()` | VCF | `variants`, `cells`, `rows`, `bands`, `info_fields` |
| `jbrowse_synteny()` | PAF, PIF, delta, chain, MashMap | `alignments`, `ribbons`, `offscreen*` |
| `jbrowse_dotplot()` | the same | `alignments`, `segments` |
| `jbrowse_multiway()` | an all-vs-all PAF | `lanes`, `segments`, `blocks`, `ribbons`, `alignments` |

- Every result also has `legend` and `ramps`, the keys the display shows.
  `guides_jbrowse()` draws them.
- `file` takes a path, a URL, or a JBrowse track's adapter config, including
  a `FromConfigAdapter` that lists its features inline.
- `region` takes a locus string, a contig name with `start` and `end`, or a
  data frame of `chrom`, `start` and `end` for several regions.
- `display` takes a track config's `displays` entry. `jbrowse_synteny()` and
  `jbrowse_dotplot()` take the view's snapshot as `view`.
- `assembly` takes a `jbrowse_assembly()`, or an assembly as a JBrowse
  config's `assemblies` entry writes it. A region may name a contig by any
  alias in the assembly's alias table. The engine reads each file under the
  file's own alias, and the reader names the contig in its tables as the
  region does. The assembly's sequence is the `reference` where none is
  given.

## Helpers

- `jbrowse_assembly(sequence, aliases)`: an assembly, its reference sequence
  and its contigs' alias table, for any reader's `assembly`
- `jbrowse_regions(assembly)`: the assembly's contigs as regions, to read a
  file whole
- `jbrowse_contigs(file)`: each side's contigs and lengths in a pairwise
  alignment file
- `jbrowse_sequence(reference, chrom, start, end)`: the reference's bases as
  one string
- `jbrowse_reset()`: restarts the engine and clears its caches

## Drawing

- `autoplot(result)` composes the geoms into the browser's track.
  `autoplot(a, "pileup")` draws one band of an alignments track.
- `jbrowse_stack(...)` stacks results, and your own ggplots, on one x axis,
  each named one titled. A key several tracks share is drawn once, listing
  every entry any of them holds.
  `jbrowse_height()` and `jbrowse_width()` measure the figure.
- One geom per table, all named `geom_jbrowse_*()`: `reads`, `mismatches`,
  `deletions`, `insertions`, `clips`, `coverage`, `alleles`, `junctions`,
  `glyphs`, `amino_acids`, `labels`, `wiggle`, `genotypes`, `ribbons`,
  `dotplot`, `lanes` and the rest.
- `jbrowse_color()` marks a column as the colors the display paints, so
  `aes(fill = color)` needs no scale. `scale_fill_jbrowse()`,
  `scale_color_jbrowse()` and `scale_x_jbrowse()` add the display's keys and
  axes.

## Notes

- `reads$read` is the id every other alignments table names a read by.
  Filtering `reads` keeps the id.
- Under `display$facet`, an alignments or feature result lists its
  `sections` in stacking order, and each table's `section` column names the
  one a row sits in.
- `fade` is the share of the depth a mark's site holds, zero under the
  display's threshold. R applies it at the figure's pixels per base.
- `jbrowse_variants()` returns its `rows` in the display's arrangement, by
  the `rows` config (`domain`, `labels`, `kept`) and the `facet` bands.
  `tree` takes an `hclust`, or a list of them, one per band.
- `jbrowse_variants()` gives each record its `filter` and a column per INFO
  key, named as the key, or `info_<key>` where the key names a record column
  such as `type`. The header's `Number` decides the shape: a flag is logical,
  `Number=1` a plain column, and `A`, `R`, `G` or `.` a list column of
  vectors. `info_fields` names each key's column.
- `jbrowse_variants()` keeps a record only where a sample calls an allele, as
  the genotype matrix does, so a calls VCF with no GT keeps none.
  `jbrowse_features()` reads every record, as the plain variant track draws
  them: `features` gains `ref`, `alt`, `filter` and the INFO columns, and
  `samples` holds a row per record and sample with a column per FORMAT key,
  typed as INFO is.
- `geom_jbrowse_genotypes()` fades records that share a pixel at the figure's
  `width`, through the engine's `variantDensityAlpha`. `fade = FALSE` turns
  it off.
- `byte_limit` refuses a region whose compressed size, read off the index,
  passes it. The default is 200 MB.
- A CRAM needs `reference`, or an `assembly` with a sequence, because it
  stores reads against the FASTA it was encoded with.
- The bisulfite layer needs `reference`, because it calls each cytosine by
  comparing the read's base to the reference's.
- `jbrowse_features()` with `reference` translates each CDS once `width`
  gives a base a pixel or more, striping `rects` by codon and listing each
  residue in `amino_acids`. A residue's `text` is the letter the display
  writes at 8 pixels a base or more, numbered once the codon fits the number.

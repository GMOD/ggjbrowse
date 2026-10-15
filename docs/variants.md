# Variant matrices

`jbrowse_variants()` reads a multi-sample VCF through JBrowse's variant
display: the records its filters keep, each sample's genotype at each, and
the fill the display paints it, its hue shaded by the alt dosage
([by superpopulation](gallery.md#variants_by_population)).

- `display` takes the display's settings: `unit = "haplotype"` for a row per
  haplotype, `color`, `filter`, `minorAlleleFrequencyFilter`, ...
- `samples` is a table of the samples' attributes, which `rowColor` marks the
  rows by.
- `v$variants` gives each record's `filter`, its FILTER column, and a column
  per INFO key (`AC`, `AN`, ...), typed by the header. A key with a value per
  alt, such as `AC`, is a list column. `v$fields$info` lists each key's
  header line.
- The matrix keeps a record only where a sample calls an allele. A calls VCF
  with no GT, such as a somatic SV caller's, reads through
  `jbrowse_features()` instead, as JBrowse's plain variant track draws it.
  Its `features` hold each record's `filter` and INFO columns, and its
  `samples` each sample's FORMAT values (`AD`, `DP`, ...).
- `rows` arranges the rows as a saved session does. `domain` lists the rows
  that lead, `labels` renames rows, and `kept` names the rows shown.
- `facet = "superpopulation"` stacks the rows in a band per value, named
  beside them ([by band](gallery.md#variants_banded)).
- Reordering or filtering `v$rows` rearranges the matrix's rows, each
  staying in its band. Assigning a filtered `v$rows` or `v$variants` drops
  the cells of the rows or records it leaves out, so `as.matrix(v)` follows.
- `as.matrix(v)` is the genotype matrix the display's "Cluster rows by
  genotype" clusters, a column per ALT allele with `NA` for a no-call.
- A tree in `v$tree` puts the rows in its leaf order, and `autoplot()` draws
  it beside them. The display's dialog writes
  `v$tree <- hclust(dist(as.matrix(v)), "average")`, and any `hclust`, or
  anything `as.hclust()` takes, labelled by row name, will do
  ([clustered](gallery.md#variants_clustered)).
- Under bands, a list of trees, one per band, clusters each band apart, as
  the display does.
- Zoomed out, records sharing a pixel fade by how many share it, as the
  display fades them, so a row shades by its share of alt records. The fade
  is `v$variants$opacity`, at the `width` the reader read at.

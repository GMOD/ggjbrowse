# ggjbrowse (development version)

- `jbrowse_alignments()`'s `reads` table drops `category` for `level`, the
  read's value as the color's `scale$domain` names it (`"1"`, `"-1"`,
  `"RL"`, `"short"`, a MAPQ, a tag's value, or `""`). The `arcs` table
  drops `category` for `level`, the `insertSizeAndOrientation` level whose
  color it paints, and `class`, a split read's junction class
  (`"deletion"`, `"duplication"`, `"inversion"`).
- `insertSizeAndOrientation` is the one paired-end read color, ordinal over
  its levels: the `insertSize` and `pairOrientation` color fields are gone,
  and `display = list(normalInsertSize = list(min = 200, max = 600))` pins the
  insert sizes it colors normal.
- `jbrowse_variants()` reads its `display` as a `LinearVariantDisplay`'s,
  the one variant display, and its result's `display$type` says so.
- The engine is `@jbrowse/tables` 5.0.0-beta.19 at jbrowse-components
  f83577e48e, bundled in `inst/js/jbrowse-tables.js`.

# ggjbrowse 0.11.0

- A CRAM read again in the same session under another assembly decodes its
  bases against that assembly's sequence. A first read against a FASTA that
  spells the contig differently, with no alias table, had left every
  substituted base as N for the rest of the session.
- A CRAM read against a sequence that has no name for its contig warns that
  its bases decode as N, and names the sequence's contigs.
- The engine is `@jbrowse/tables` 5.0.0-beta.18 at jbrowse-components
  d5a0299d00, bundled in `inst/js/jbrowse-tables.js`.

# ggjbrowse 0.10.0

The first tagged release, for reference on GitHub; it is not on CRAN. Install
it with `pak::pak("GMOD/ggjbrowse@v0.10.0")`.

- Readers for BAM and CRAM (`jbrowse_alignments()`), GFF3, GTF, BED and BigBed
  (`jbrowse_features()`), BigWig (`jbrowse_wiggle()`), VCF
  (`jbrowse_variants()`), and PAF, PIF, delta, chain and MashMap
  (`jbrowse_synteny()`), each returning the tables JBrowse draws from.
- A geom per table, `autoplot()` for a whole track, and `jbrowse_stack()` for
  several tracks over one region.
- Display settings in JBrowse 2 v5's spellings, the wiggle display's
  `mark = list(extent = )` included.
- The engine is `@jbrowse/tables` 5.0.0-beta.18, bundled in
  `inst/js/jbrowse-tables.js`.

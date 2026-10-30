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

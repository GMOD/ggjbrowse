# ggjbrowse (development version)

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

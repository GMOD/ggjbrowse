# Feature tracks

`jbrowse_features()` reads a GFF3, GTF, BED or BigBed region through
JBrowse's feature display and returns what it draws:

- each feature on its row, its gene glyph's boxes, intron lines, strand
  chevrons and arrows in their colors, and its transcripts
- the labels the display shows at that zoom, scaled to keep to the room the
  layout left them at any output size

`display` takes the display's settings (`color`, `displayMode`,
`showLabels`, `geneGlyphMode`, `filter`, ...), and `width` is the
pixel width the display lays the region out across, which decides which
names fit. `jbrowse_stack(genes, reads)` draws genes over reads
([GENCODE genes](gallery.md#genes_and_reads)).

- Under `facet = "strand"`, or `list(field = "type", domain = ...,
  hidden = ...)`, the display packs a section per value, and `autoplot()`
  draws each section in its own panel.

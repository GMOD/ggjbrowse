#' What a linear reader returns
#'
#' [jbrowse_alignments()], [jbrowse_features()], [jbrowse_variants()] and
#' [jbrowse_wiggle()] each return a classed list, `a$reads` a data frame. A
#' fixed set of slots holds what is not a table; every other slot is a table
#' of records or of marks, which each reader's page lists.
#'
#' @section The slots every result holds:
#' * `regions`: a data frame, a row per region read: `chrom`, `start` and
#'   `end` (0-based, half-open), where the axis lays it (`offset`, the
#'   position its start takes, and `cum_start` and `cum_end`, its ends), and
#'   `missing`, the reason its contig reads empty, `""` where the file has the
#'   contig
#' * `width`: the pixels the regions were laid out across
#' * `display`: the display's config as JBrowse resolves it, every slot
#'   under JBrowse's own name, the schema's defaults filled in
#' * `keys`: the display's color keys, each a data frame of `label` and
#'   `color` titled by its `"title"` attribute, as [scale_fill_jbrowse()]
#'   draws them
#' * `ramps`: the display's color bars, as [guides_jbrowse()] draws them
#' * `y_scale`: the score axis of a [jbrowse_wiggle()] result, or the
#'   coverage band's depth axis of a [jbrowse_alignments()] result, as the
#'   display autoscales it: `domain`, `type` (`"linear"`, `"log"` or
#'   `"symlog"`), `symlog_constant` and `ticks`; `NULL` for the others
#' * `fields`: the file's declarations of the columns it adds, each a data
#'   frame: `info` and `format` (a VCF's header lines), `attributes` (a GFF3's
#'   or BigBed's fields), and `annotations` (the fields of a VCF's `ANN` or
#'   `CSQ`)
#' * `warnings`: what the engine warned of while reading, each also raised as
#'   an R warning
#'
#' A [jbrowse_features()] result also holds `glyph`, `height` and
#' `bp_per_px`, the layout's pixel sizes; a [jbrowse_wiggle()] result
#' `rendering`, `origin`, `row_layout`, `height` and `bp_per_px`; a
#' [jbrowse_alignments()] result `read_pitch`, `glyph`, `arc_span_max` and
#' `reference`; and a [jbrowse_variants()] result the `tree` a caller sets.
#' @name jbrowse_result
NULL

#!/bin/sh
# inst/extdata/hprc_amylase.paf: minimap2 alignments of four HPRC haplotypes
# and GRCh38 over the amylase locus, each aligned to the next, as
# jbrowse-components builds them for its amylase tutorial.
# Usage: scripts/extdata-hprc.sh [path to jbrowse-components]; needs only cp.
set -eu
jbrowse=${1:-$HOME/src/jbrowse-components}
out=$(cd "$(dirname "$0")" && pwd)/../inst/extdata
cp "$jbrowse/test_data/amylase/amylase_adjacent.paf" "$out/hprc_amylase.paf"

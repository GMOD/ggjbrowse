#!/bin/sh
# inst/extdata/volvox-simple-inv-bow.bam: the four long reads of
# jbrowse-components' test_data/volvox/volvox-simple-inv.bam whose chains join
# ctgA:15,000-19,950 and ctgA:30,010-33,000, two of them bowing over a hidden
# segment. Usage: JBROWSE_COMPONENTS=... scripts/extdata-volvox-bow.sh; needs
# samtools.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
jbc=${JBROWSE_COMPONENTS:-$here/../../jbrowse-components}
out=$here/../inst/extdata/volvox-simple-inv-bow.bam
names=$(mktemp)
cat > "$names" <<NAMES
3be100ba-dc79-4424-5b37-51e58fcec283
6b760d21-9e1a-5179-9d95-d556f12b4d78
6b7af653-b61a-9922-60e5-b41464bdebf3
4f718ecd-53b5-0569-a21b-cd160c9b4a10
NAMES
samtools view -b --no-PG -N "$names" -o "$out" "$jbc/test_data/volvox/volvox-simple-inv.bam"
samtools index "$out"
rm "$names"

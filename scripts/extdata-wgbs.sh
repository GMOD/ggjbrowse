#!/bin/sh
# inst/extdata's bisulfite clip, TAIR10 NC_003070.9:4,398,001-4,412,000, the
# window of JBrowse's bisulfite tutorial:
#   arabidopsis_wgbs.bam    Col-0 WGBS reads (DRR029742, bwameth) over the
#                           window, read sequences kept for the C-to-T calls
#   arabidopsis_tair10.fa.gz NC_003070.9 at its true coordinates, Ns but for
#                           the window and 1 kb either side
# Usage: scripts/extdata-wgbs.sh; needs samtools, bgzip and twoBitToFa.
set -eu
out=$(pwd)/inst/extdata
cram=https://jbrowse.org/demos/bisulfite/arabidopsis_wgbs_bisulfite.cram
twobit=https://hgdownload.soe.ucsc.edu/hubs/GCF/000/001/735/GCF_000001735.4/GCF_000001735.4.2bit
contig=NC_003070.9
start=4398001
end=4412000
tmp=$(mktemp -d)
cd "$tmp"

{
  echo ">$contig"
  {
    head -c $((start - 1001)) /dev/zero | tr '\0' N
    twoBitToFa -seq=$contig -start=$((start - 1001)) -end=$((end + 1000)) \
      "$twobit" /dev/stdout | tail -n +2 | tr -d '\n'
  } | fold -w 60
} | bgzip -i -I "$out/arabidopsis_tair10.fa.gz.gzi" \
  > "$out/arabidopsis_tair10.fa.gz"
samtools faidx "$out/arabidopsis_tair10.fa.gz"

# the CRAM's slices reach past the window, so decoding needs the whole contig
twoBitToFa -seq=$contig "$twobit" contig.fa
samtools view -b -T contig.fa "$cram" \
  "$contig:$start-$end" > reads.bam
samtools calmd -b reads.bam "$out/arabidopsis_tair10.fa.gz" \
  > "$out/arabidopsis_wgbs.bam" 2> /dev/null
samtools index "$out/arabidopsis_wgbs.bam"
cd "$out"
rm -r "$tmp"

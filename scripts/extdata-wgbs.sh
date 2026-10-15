#!/bin/sh
# inst/extdata's bisulfite clip, TAIR10 NC_003070.9:4,398,001-4,412,000, the
# window of JBrowse's bisulfite tutorial:
#   arabidopsis_wgbs.bam     Col-0 WGBS reads (DRR029742, bwameth) over the
#                            window, read sequences kept for the C-to-T calls
#   arabidopsis_tair10.fa.gz NC_003070.9 at its true coordinates, Ns but for
#                            the window and 1 kb either side
#   arabidopsis_refseq.bb    the TAIR10 hub's NCBI RefSeq genes over the window
#   arabidopsis_rmsk.bb      the hub's RepeatMasker repeats over the window
#   arabidopsis_methyldackel_{CpG,CHG,CHH}.bw
#                            MethylDackel's percent methylated per cytosine
#                            over the window, one file per context
# Usage: scripts/extdata-wgbs.sh; needs samtools, bgzip, twoBitToFa,
# bigBedInfo, bigBedToBed, bedToBigBed, bigWigToBedGraph and bedGraphToBigWig.
set -eu
out=$(pwd)/inst/extdata
cram=https://jbrowse.org/demos/bisulfite/arabidopsis_wgbs_bisulfite.cram
hub=https://hgdownload.soe.ucsc.edu/hubs/GCF/000/001/735/GCF_000001735.4
twobit=$hub/GCF_000001735.4.2bit
bigwigs=https://jbrowse.org/demos/bisulfite/arabidopsis_wgbs_
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

bbi=$hub/bbi/GCF_000001735.4_TAIR10.1.
curl -fsS -o refseq.bb "${bbi}ncbiRefSeq.bb"
bigBedInfo -chroms refseq.bb |
  awk -v c=$contig '$1 == c { print $1 "\t" $3 }' > chrom.sizes

# the features overlapping the window, each whole, under the file's own autoSql
clip_bigbed() {
  [ -f "$3.bb" ] || curl -fsS -o "$3.bb" "$1"
  bigBedInfo -as "$3.bb" |
    awk '/^table/ { on = 1 } on { print } on && /^[ \t]*\)[ \t]*$/ { exit }' \
      > "$3.as"
  bigBedToBed -chrom=$contig -start=$((start - 1)) -end=$end "$3.bb" "$3.bed"
  bedToBigBed -type="$4" -as="$3.as" -tab "$3.bed" chrom.sizes "$out/$2"
}
clip_bigbed "${bbi}ncbiRefSeq.bb" arabidopsis_refseq.bb refseq bed12+9
clip_bigbed "${bbi}rmsk.bb" arabidopsis_rmsk.bb rmsk bed9+5

for context in CpG CHG CHH; do
  bigWigToBedGraph -chrom=$contig -start=$((start - 1)) -end=$end \
    "$bigwigs$context.bw" stdout | sort -k1,1 -k2,2n > "$context.bedGraph"
  bedGraphToBigWig "$context.bedGraph" chrom.sizes \
    "$out/arabidopsis_methyldackel_$context.bw"
done
cd "$out"
rm -r "$tmp"

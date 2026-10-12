#!/bin/sh
# inst/extdata's GAPDH files, hg38 chr12:6,500,001-6,570,000:
#   gapdh_direct_rna.bam  NA12878 nanopore direct RNA over GAPDH, about 1 read
#                         in 1,250, MD tags added so no reference is needed
#   gapdh_hg38.fa.gz      chr12 at its true coordinates, Ns but for the window
#   gapdh_gencode.gff3.gz GENCODE v47's features in the window
# Usage: scripts/extdata-gapdh.sh; needs samtools, bgzip and tabix.
set -eu
out=$(pwd)/inst/extdata
bam=https://jbrowse.org/genomes/GRCh38/NA12878/NA12878-DirectRNA.pass.dedup.NoU.fastq.hg38.minimap2.sorted.bam
fasta=https://jbrowse.org/genomes/GRCh38/fasta/hg38.prefix.fa.gz
gff=https://jbrowse.org/genomes/GRCh38/gencode/gencode.v47.chr_patch_hapl_scaff.annotation.sorted.gff3.gz
start=6500001
end=6570000
tmp=$(mktemp -d)
# samtools and tabix save a remote file's index in the working directory
cd "$tmp"

{
  echo ">chr12"
  {
    head -c $((start - 1)) /dev/zero | tr '\0' N
    samtools faidx "$fasta" "12:$start-$end" | tail -n +2 | tr -d '\n'
  } | fold -w 60
} | bgzip -i -I "$out/gapdh_hg38.fa.gz.gzi" > "$out/gapdh_hg38.fa.gz"
samtools faidx "$out/gapdh_hg38.fa.gz"

samtools view -b -F 256 -s 7.0008 "$bam" chr12:6534512-6538374 > "$tmp/reads.bam"
samtools calmd -b "$tmp/reads.bam" "$out/gapdh_hg38.fa.gz" \
  > "$out/gapdh_direct_rna.bam" 2> /dev/null
samtools index "$out/gapdh_direct_rna.bam"

# GENCODE writes gene_name and transcript_name where a display labels by Name
tabix "$gff" "chr12:$start-$end" | sort -k1,1 -k4,4n | awk -F'\t' -v OFS='\t' '
  $9 !~ /(^|;)Name=/ && match($9, /(transcript|gene)_name=[^;]*/) {
    if ($3 != "gene" && match($9, /transcript_name=[^;]*/)) {
      $9 = "Name=" substr($9, RSTART + 16, RLENGTH - 16) ";" $9
    } else if (match($9, /gene_name=[^;]*/)) {
      $9 = "Name=" substr($9, RSTART + 10, RLENGTH - 10) ";" $9
    }
  } { print }' | bgzip > "$out/gapdh_gencode.gff3.gz"
tabix -f -p gff "$out/gapdh_gencode.gff3.gz"
cd "$out"
rm -r "$tmp"

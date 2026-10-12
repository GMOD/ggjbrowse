#!/bin/sh
# inst/extdata's human structural variant clips. Each keeps the reads over its
# view, cut to the view plus 250 bp either side by scripts/trim-window.awk and
# given MD tags, so no reference is needed:
#   colo829_tumor_chr3.bam      COLO829 melanoma nanopore reads over the two chr3
#                               breakpoints of the der(3) chain,
#                               chr3:25,357,600-25,361,000
#   hg00151_chr1_inversion.bam  HG00151 nanopore reads over a 1.2 kb chr1
#                               inversion, chr1:197,786,900-197,789,700
#   hg008t_cuzd1_deletion.bam   HG008-T PacBio HiFi reads over a 1.8 kb somatic
#                               deletion in CUZD1, chr10:122,831,700-122,840,800
# Usage: scripts/extdata-sv.sh; needs samtools and awk.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
out=$here/../inst/extdata
fasta=https://jbrowse.org/genomes/GRCh38/fasta/hg38.prefix.fa.gz
tmp=$(mktemp -d)
cd "$tmp"

clip() {
  name=$1; src=$2; contig=$3; start=$4; end=$5; sample=$6
  samtools faidx "$fasta" "${contig#chr}" | sed "1s/.*/>$contig/" > "$contig.fa"
  samtools faidx "$contig.fa"
  samtools view -h -T "$contig.fa" $sample "$src" "$contig:$start-$end" |
    awk -v ws=$((start - 250)) -v we=$((end + 250)) -f "$here/trim-window.awk" |
    samtools sort -O bam - |
    samtools calmd -b - "$contig.fa" > "$out/$name.bam" 2> /dev/null
  samtools index "$out/$name.bam"
}

ont=https://ont-open-data.s3.amazonaws.com/colo829_2024.03
clip colo829_tumor_chr3 "$ont/wf_somatic_variation/sup/COLO829_tumor.ht.cram" \
  chr3 25357600 25361000 "-s 1.5"
clip hg00151_chr1_inversion \
  https://jbrowse.org/demos/ont/HG00151-ONT-hg38.chr1_inversion.bam \
  chr1 197786900 197789700 "-s 1.7"
clip hg008t_cuzd1_deletion \
  https://jbrowse.org/demos/cgiab/HG008-T_chr10_CUZD1_deletion.bam \
  chr10 122831700 122840800 ""

cd "$out"
rm -r "$tmp"

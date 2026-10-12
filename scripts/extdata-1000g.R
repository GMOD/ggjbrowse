#!/usr/bin/env Rscript
# Writes the samples file scripts/extdata-1000g.sh clips the VCF to: 20
# people from each 1000 Genomes superpopulation, taken from its populations
# in turn. Usage: scripts/extdata-1000g.R <population tsv> <vcf sample list>
args <- commandArgs(TRUE)
populations <- utils::read.delim(args[1], header = FALSE,
  col.names = c("name", "population", "url"))
superpopulations <- list(
  African = c("YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"),
  European = c("CEU", "TSI", "FIN", "GBR", "IBS"),
  `East Asian` = c("CHB", "JPT", "CHS", "CDX", "KHV"),
  `South Asian` = c("GIH", "PJL", "BEB", "STU", "ITU"),
  `Admixed American` = c("MXL", "PUR", "CLM", "PEL"))
people <- populations[populations$name %in% readLines(args[2]), ]
picked <- do.call(rbind, lapply(names(superpopulations), function(s) {
  d <- people[people$population %in% superpopulations[[s]], ]
  d <- d[order(d$population, d$name), ]
  turn <- stats::ave(seq_len(nrow(d)), d$population, FUN = seq_along)
  d <- utils::head(d[order(turn, d$population), ], 20)
  data.frame(name = d$name, superpopulation = s, population = d$population)
}))
utils::write.table(picked, stdout(), sep = "\t", quote = FALSE,
  row.names = FALSE)


library(SeqArray)
library(LDshrink)
library(readr)
library(dplyr)
library(RSSp)
library(purrr)

inputf <- snakemake@input[["traitf"]]
geno_gdsf <- snakemake@input[["geno_gdsf"]]
geno_ogdsf <- snakemake@output[["geno_gdsf"]]

haplo_gdsf <- snakemake@input[["haplo_gdsf"]]
haplo_ogdsf <- snakemake@output[["haplo_gdsf"]]






input_df <- read_delim(inputf,delim="\t")


geno_gds <- seqOpen(geno_gdsf)
output_df <- subset_gds(geno_gds,input_df,region_id=T)

p <- calc_p(geno_gds)

stopifnot(nrow(output_df)==p)

seqExport(geno_gds,geno_ogdsf)
seqClose(geno_gds)


geno_gds <- seqOpen(haplo_gdsf)
output_df <- subset_gds(geno_gds,input_df,region_id=T)

p <- calc_p(geno_gds)

stopifnot(nrow(output_df)==p)

seqExport(geno_gds,haplo_ogdsf)
seqClose(geno_gds)


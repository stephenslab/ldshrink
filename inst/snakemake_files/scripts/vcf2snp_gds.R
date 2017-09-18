library(SNPRelate)

library(tidyverse)

inpf <- snakemake@input[["vcff"]]
mapf <- snakemake@input[["mapf"]]
outf <- snakemake@output[["gdsf"]]
cores <- as.integer(snakemake@threads)
stop()
map_data <- map_df(mapf,read_delim,delim=" ")

snpVCF2GDS(inpf, outf, parallel = cores,storage.option="LZ4_RA.fast")


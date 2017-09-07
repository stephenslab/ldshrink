library(SeqArray)
library(tidyverse)
library(LDshrink)


inpf <- snakemake@input[["vcff"]]
toutf <- snakemake@output[["temp_gds"]]


cores <- as.integer(snakemake@threads)

seqVCF2GDS(vcf.fn = inpf, out.fn = toutf,parallel = cores,
           storage.option = "LZ4_RA.fast")                                                


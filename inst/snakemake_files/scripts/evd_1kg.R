library(SeqArray)
library(tidyverse)
library(LDshrink)


gdsf <- snakemake@input[["gdsf"]]
outf <- snakemake@output[["evdf"]]
region_id <- as.integer(as.numeric(snakemake@params[["region_id"]]))


 save.image()
stopifnot(!is.null(region_id), !is.null(gdsf), !is.null(outf))

stopifnot(file.exists(gdsf), !file.exists(outf))

chunkwise_LDshrink(gds_file = gdsf, region_id = region_id, outfile = outf)


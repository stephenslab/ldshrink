library(LDshrink)
library(tidyverse)
gdsf <- snakemake@input[["gdsf"]]
chrom <-as.integer(as.numeric(snakemake@params[["chrom"]]))
out_dir <-snakemake@params[["outdir"]]

chunkwise_LDshrink_ldsc(gdsf,chrom,out_dir)


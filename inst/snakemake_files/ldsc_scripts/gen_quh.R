library(RSSp)
library(readr)
library(dplyr)
library(purrr)

gwasf <- snakemake@input[["gwasf"]]
evdf <- snakemake@params[["evdf"]]

stopifnot(file.exists(evdf))
gwas_df <- read_delim(gwasf,delim="\t") %>% mutate(chr=as.character(chr))


write_delim(
    gen_quh_chunk(gwas_df,evdf),
    snakemake@output[["output_gwasf"]],
    delim="\t",col_names=T)



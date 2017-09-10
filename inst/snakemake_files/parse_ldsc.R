
library(dplyr)
library(readr)
library(RSSp)

inf <- snakemake@input[["logf"]]
tpf <- snakemake@input[["tparamf"]]
fgeneid <- snakemake@params[["fgeneid"]]
outf <- snakemake@output[["logf"]]

#save.image()
res_df <- parse_ldsc_h2log(inf) %>%
    mutate(fgeneid = as.character(fgeneid))

res_df <- read_delim(tpf,delim="\t") %>%
    mutate(fgeneid = as.character(fgeneid)) %>%
    inner_join(res_df)
write_delim(res_df, path = outf, delim = "\t")

library(readr)
library(dplyr)
library(purrr)



write_delim(
    map_df(snakemake@input[["quhf"]],read_delim,delim="\t"),
    snakemake@output[["quhf"]],delim="\t",col_names=T)
            


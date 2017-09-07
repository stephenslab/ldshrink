library(SeqArray)
library(tidyverse)
library(LDshrink)



outf <- snakemake@output[["gdsf"]]
toutf <- snakemake@input[["temp_gds"]]
mapf <- snakemake@input[["mapf"]]
panelf <- snakemake@input[["panelf"]]
breakf <- snakemake@input[["breakf"]]

cores <- as.integer(snakemake@threads)


map_df <- readRDS(mapf)
import_panel_data(temp_gds=toutf,
                  map_df=map_df,
                  output_file = outf,
                  ld_break_file=breakf,
                  overwrite=F,
                  parallel=cores,
                  panel_ind_file=panelf,
                  subset_MAF=0.05)
                                                


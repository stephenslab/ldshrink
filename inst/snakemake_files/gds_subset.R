library(SeqArray)
library(tidyverse)
library(LDshrink)



outf <- snakemake@output[["gdsf"]]
mapdf <- readRDS(snakemake@input[["mapf"]])
toutf <- snakemake@input[["temp_gds"]]
breakf <- snakemake@input[["breakf"]]



save.image()
import_panel_data(temp_gds=toutf,
                  map_df=mapdf,
                  output_file = outf,
                  ld_break_file=breakf,
                  overwrite=T)
                                                

